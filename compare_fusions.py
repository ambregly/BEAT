#!/usr/bin/env python3
"""
Comparaison des fusions detectees par les kmers avec les fusions de reference BEAT AML.

On construit, pour chaque source, un ensemble de "cles de fusion" definies par :
    (SampleID, left_gene, left_chr, right_gene, right_chr, fusion_index)

Puis on calcule :
    - Vrais positifs (TP)  : cles presentes dans BEAT AML ET dans kmer
    - Faux positifs (FP)   : cles presentes uniquement dans kmer
    - Faux negatifs (FN)   : cles presentes uniquement dans BEAT AML
    - Precision = TP / (TP + FP)
    - Recall    = TP / (TP + FN)

Particularites du fichier kmer :
    - colonne 1 : "left_gene_left_chr_pos_pos_right_gene_right_chr_pos_pos_idx|idx|idx"
                  -> une ligne peut porter PLUSIEURS fusion_index (separes par '|'),
                     chacun donne une cle distincte.
    - colonne 2 : SampleID avec un suffixe 'R' (ex: BA2409R) qui correspond a
                  BA2409 dans BEAT AML -> on retire le 'R' final.

Usage :
    python compare_fusions.py BEAT_AML.csv kmer.tsv
    python compare_fusions.py BEAT_AML.csv kmer.tsv --outdir resultats
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


# colonne 1 du fichier kmer :
#   left_gene _ left_chr _ pos _ pos _ right_gene _ right_chr _ pos _ pos _ idx|idx|idx
# Les noms de genes peuvent contenir des '_', donc on s'ancre sur les champs
# "chr..." et sur les positions numeriques plutot que de faire un simple split('_').
KMER_COL1_RE = re.compile(
    r"^(?P<left_gene>.+?)_"
    r"(?P<left_chr>chr[^_]+)_"
    r"(?P<left_pos1>\d+)_"
    r"(?P<left_pos2>\d+)_"
    r"(?P<right_gene>.+?)_"
    r"(?P<right_chr>chr[^_]+)_"
    r"(?P<right_pos1>\d+)_"
    r"(?P<right_pos2>\d+)_"
    r"(?P<indices>[\d|]+)$"
)


def norm(value):
    """Normalisation commune d'un champ (gene, chr, index) : trim + majuscules."""
    return str(value).strip().upper()


def normalize_sample_id(sample_id):
    """SampleID kmer (BA2409R) -> SampleID BEAT AML (BA2409) : on retire le 'R' final."""
    sid = str(sample_id).strip()
    if sid.endswith("R"):
        sid = sid[:-1]
    return norm(sid)


def make_key(sample_id, left_gene, left_chr, right_gene, right_chr, fusion_index):
    """Cle de fusion comparable entre les deux sources."""
    return (
        norm(sample_id),
        norm(left_gene),
        norm(left_chr),
        norm(right_gene),
        norm(right_chr),
        norm(fusion_index),
    )


def load_beat_aml(path):
    """Lit BEAT_AML.csv et renvoie un dict {cle: [lignes brutes]}."""
    keys = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        required = {
            "SampleID", "left_gene", "left_chr",
            "right_gene", "right_chr", "fusion_index",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"[BEAT AML] colonnes manquantes : {sorted(missing)}")

        for row in reader:
            key = make_key(
                row["SampleID"],
                row["left_gene"], row["left_chr"],
                row["right_gene"], row["right_chr"],
                row["fusion_index"],
            )
            keys[key].append(row)
    return keys


def load_kmer(path):
    """Lit le fichier kmer (tsv) et renvoie un dict {cle: [lignes brutes]}.

    Chaque ligne peut generer plusieurs cles (un fusion_index par index liste).
    """
    keys = defaultdict(list)
    skipped = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for lineno, fields in enumerate(reader, start=1):
            if not fields or not fields[0].strip():
                continue
            col1 = fields[0].strip()
            sample_raw = fields[1].strip() if len(fields) > 1 else ""

            m = KMER_COL1_RE.match(col1)
            if not m:
                skipped.append((lineno, col1))
                continue

            sample_id = normalize_sample_id(sample_raw)
            for idx in m.group("indices").split("|"):
                idx = idx.strip()
                if not idx:
                    continue
                key = make_key(
                    sample_id,
                    m.group("left_gene"), m.group("left_chr"),
                    m.group("right_gene"), m.group("right_chr"),
                    idx,
                )
                keys[key].append({"line": lineno, "col1": col1, "sample": sample_raw})

    if skipped:
        sys.stderr.write(
            f"[kmer] {len(skipped)} ligne(s) non parsable(s) ignoree(s) "
            f"(ex. ligne {skipped[0][0]}: {skipped[0][1]!r})\n"
        )
    return keys


def write_keys(path, keys):
    """Ecrit un ensemble de cles dans un fichier tsv trie."""
    header = ["SampleID", "left_gene", "left_chr",
              "right_gene", "right_chr", "fusion_index"]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        for key in sorted(keys):
            writer.writerow(key)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("beat_aml", help="fichier BEAT_AML.csv")
    parser.add_argument("kmer", help="fichier kmer (tsv)")
    parser.add_argument("--outdir", default=None,
                        help="repertoire ou ecrire les listes TP/FP/FN (optionnel)")
    args = parser.parse_args()

    beat = load_beat_aml(args.beat_aml)
    kmer = load_kmer(args.kmer)

    beat_keys = set(beat)
    kmer_keys = set(kmer)

    tp = beat_keys & kmer_keys      # dans les deux
    fp = kmer_keys - beat_keys      # uniquement kmer
    fn = beat_keys - kmer_keys      # uniquement BEAT AML

    n_tp, n_fp, n_fn = len(tp), len(fp), len(fn)
    precision = n_tp / (n_tp + n_fp) if (n_tp + n_fp) else 0.0
    recall = n_tp / (n_tp + n_fn) if (n_tp + n_fn) else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) else 0.0)

    print("=== Comparaison fusions kmer vs BEAT AML ===")
    print(f"Cles BEAT AML : {len(beat_keys)}")
    print(f"Cles kmer     : {len(kmer_keys)}")
    print()
    print(f"Vrais positifs (TP) : {n_tp}")
    print(f"Faux positifs  (FP) : {n_fp}")
    print(f"Faux negatifs  (FN) : {n_fn}")
    print()
    print(f"Precision = TP/(TP+FP) = {precision:.4f}")
    print(f"Recall    = TP/(TP+FN) = {recall:.4f}")
    print(f"F1-score               = {f1:.4f}")

    if args.outdir:
        import os
        os.makedirs(args.outdir, exist_ok=True)
        write_keys(os.path.join(args.outdir, "vrais_positifs.tsv"), tp)
        write_keys(os.path.join(args.outdir, "faux_positifs.tsv"), fp)
        write_keys(os.path.join(args.outdir, "faux_negatifs.tsv"), fn)
        print(f"\nListes ecrites dans : {args.outdir}/")


if __name__ == "__main__":
    main()
