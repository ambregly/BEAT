#!/usr/bin/env python3
"""
Comparaison des fusions detectees par les kmers avec les fusions de reference BEAT AML.

Calcule vrais positifs (TP), faux positifs (FP), faux negatifs (FN), puis
precision = TP/(TP+FP) et recall = TP/(TP+FN).

Une fusion est identifiee par :
    (SampleID, left_gene, left_chr, right_gene, right_chr, fusion_index)

Particularites du fichier kmer :
    - colonne 1 : "left_gene_left_chr_pos_pos_right_gene_right_chr_pos_pos_idx|idx|idx"
                  -> une ligne porte la liste des fusion_index BEAT AML associes a
                     cette signature (separes par '|').
    - colonne 2 : SampleID avec un suffixe 'R' (ex: BA2409R) qui correspond a
                  BA2409 dans BEAT AML -> on retire le 'R' final.

Trois criteres de comparaison (--match) :
    - row (DEFAUT) : l'unite comptee est la LIGNE kmer. Ses index cherchent leur
                     paire (echantillon + fusion_index) dans BEAT AML. Des qu'UN
                     index matche -> 1 vrai positif (on ignore les autres index de
                     la ligne). Si AUCUN ne matche -> 1 seul faux positif.
                     FN = enregistrements BEAT AML jamais retrouves.
    - index        : chaque couple (echantillon, fusion_index) est compte
                     individuellement (un index non matche = un FP).
    - genepair     : match par (echantillon + paire de genes), fusion_index ignore.

Usage :
    python compare_fusions.py BEAT_AML.csv kmer.tsv
    python compare_fusions.py BEAT_AML.csv kmer.tsv --match index --outdir resultats
"""

import argparse
import csv
import os
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
    """Normalisation commune d'un champ (gene, index) : trim + majuscules."""
    return str(value).strip().upper()


def norm_chr(value):
    """Normalisation d'un chromosome : trim, majuscules, sans prefixe 'chr'.

    BEAT AML stocke le chromosome en nombre nu (ex: '6', '18'), tandis que le
    fichier kmer le prefixe par 'chr' (ex: 'chr6'). On harmonise les deux.
    """
    chrom = str(value).strip().upper()
    if chrom.startswith("CHR"):
        chrom = chrom[3:]
    return chrom


def normalize_sample_id(sample_id):
    """SampleID kmer (BA2409R) -> SampleID BEAT AML (BA2409) : on retire le 'R' final."""
    sid = str(sample_id).strip()
    if sid.endswith("R"):
        sid = sid[:-1]
    return norm(sid)


def make_key(sample_id, left_gene, left_chr, right_gene, right_chr,
             fusion_index, use_index=True):
    """Cle de fusion comparable entre les deux sources.

    use_index=True  -> cle = echantillon + paire de genes + fusion_index
    use_index=False -> cle = echantillon + paire de genes (index ignore)
    """
    key = (
        norm(sample_id),
        norm(left_gene),
        norm_chr(left_chr),
        norm(right_gene),
        norm_chr(right_chr),
    )
    if use_index:
        key = key + (norm(fusion_index),)
    return key


def load_beat_aml(path, use_index=True):
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
                use_index=use_index,
            )
            keys[key].append(row)
    return keys


def parse_kmer_rows(path):
    """Lit le fichier kmer (tsv) et renvoie la liste des lignes parsees.

    Chaque element est un dict : line, sample_id (normalise), sample_raw,
    left_gene, left_chr, right_gene, right_chr, indices (liste de str), col1.
    """
    rows = []
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

            indices = [i.strip() for i in m.group("indices").split("|") if i.strip()]
            rows.append({
                "line": lineno,
                "sample_id": normalize_sample_id(sample_raw),
                "sample_raw": sample_raw,
                "left_gene": m.group("left_gene"),
                "left_chr": m.group("left_chr"),
                "right_gene": m.group("right_gene"),
                "right_chr": m.group("right_chr"),
                "indices": indices,
                "col1": col1,
            })

    if skipped:
        sys.stderr.write(
            f"[kmer] {len(skipped)} ligne(s) non parsable(s) ignoree(s) "
            f"(ex. ligne {skipped[0][0]}: {skipped[0][1]!r})\n"
        )
    return rows


def row_candidate_keys(row, use_index=True):
    """Cles candidates generees par une ligne kmer (une par fusion_index)."""
    return [
        make_key(row["sample_id"], row["left_gene"], row["left_chr"],
                 row["right_gene"], row["right_chr"], idx, use_index=use_index)
        for idx in (row["indices"] or [""])
    ]


def compare_row(beat_keys, kmer_rows):
    """Mode 'row' : l'unite comptee est la ligne kmer.

    - TP : ligne kmer dont AU MOINS un (echantillon + index) est dans BEAT AML.
    - FP : ligne kmer dont AUCUN candidat n'est dans BEAT AML (comptee une fois).
    - FN : enregistrements BEAT AML jamais retrouves par aucune ligne kmer.
    """
    matched_beat = set()
    tp_rows, fp_rows = [], []
    for row in kmer_rows:
        hits = [c for c in row_candidate_keys(row, use_index=True) if c in beat_keys]
        if hits:
            tp_rows.append((row, hits))
            matched_beat.update(hits)
        else:
            fp_rows.append(row)

    fn_keys = beat_keys - matched_beat
    return {
        "tp": len(tp_rows),
        "fp": len(fp_rows),
        "fn": len(fn_keys),
        "matched_beat": len(matched_beat),
        "tp_rows": tp_rows,
        "fp_rows": fp_rows,
        "fn_keys": fn_keys,
    }


def compare_set(beat_keys, kmer_rows, use_index):
    """Modes 'index' / 'genepair' : comparaison ensembliste des cles."""
    kmer_keys = set()
    for row in kmer_rows:
        kmer_keys.update(row_candidate_keys(row, use_index=use_index))

    tp = beat_keys & kmer_keys
    fp = kmer_keys - beat_keys
    fn = beat_keys - kmer_keys
    return {
        "tp": len(tp), "fp": len(fp), "fn": len(fn),
        "tp_keys": tp, "fp_keys": fp, "fn_keys": fn,
    }


def write_keys(path, keys, use_index=True):
    """Ecrit un ensemble de cles dans un fichier tsv trie."""
    header = ["SampleID", "left_gene", "left_chr", "right_gene", "right_chr"]
    if use_index:
        header.append("fusion_index")
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        for key in sorted(keys):
            writer.writerow(key)


def write_rows(path, rows, with_hit=False):
    """Ecrit des lignes kmer (mode 'row') dans un tsv."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        header = ["kmer_line", "SampleID", "left_gene", "left_chr",
                  "right_gene", "right_chr", "indices"]
        if with_hit:
            header.append("index_matche")
        writer.writerow(header)
        for item in rows:
            row, hits = item if with_hit else (item, None)
            base = [row["line"], row["sample_id"], row["left_gene"],
                    norm_chr(row["left_chr"]), row["right_gene"],
                    norm_chr(row["right_chr"]), "|".join(row["indices"])]
            if with_hit:
                base.append("|".join(sorted({h[-1] for h in hits})))
            writer.writerow(base)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("beat_aml", help="fichier BEAT_AML.csv")
    parser.add_argument("kmer", help="fichier kmer (tsv)")
    parser.add_argument("--outdir", default=None,
                        help="repertoire ou ecrire les listes TP/FP/FN (optionnel)")
    parser.add_argument("--match", choices=["row", "index", "genepair"],
                        default="row",
                        help="critere de comparaison (defaut: row). Voir l'aide en "
                             "tete de fichier.")
    args = parser.parse_args()

    use_index = (args.match != "genepair")
    beat = load_beat_aml(args.beat_aml, use_index=use_index)
    beat_keys = set(beat)
    kmer_rows = parse_kmer_rows(args.kmer)

    if args.match == "row":
        res = compare_row(beat_keys, kmer_rows)
        denom_recall_total = len(beat_keys)
    else:
        res = compare_set(beat_keys, kmer_rows, use_index=use_index)
        denom_recall_total = len(beat_keys)

    n_tp, n_fp, n_fn = res["tp"], res["fp"], res["fn"]
    precision = n_tp / (n_tp + n_fp) if (n_tp + n_fp) else 0.0
    # recall base sur les enregistrements BEAT AML retrouves
    found_beat = res.get("matched_beat", n_tp)
    recall = found_beat / denom_recall_total if denom_recall_total else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) else 0.0)

    print("=== Comparaison fusions kmer vs BEAT AML ===")
    print(f"Critere de match : {args.match}")
    print(f"Lignes kmer       : {len(kmer_rows)}")
    print(f"Fusions BEAT AML  : {len(beat_keys)}")
    print()
    print(f"Vrais positifs (TP) : {n_tp}")
    print(f"Faux positifs  (FP) : {n_fp}")
    print(f"Faux negatifs  (FN) : {n_fn}")
    if args.match == "row":
        print(f"(enregistrements BEAT AML retrouves : {found_beat}/{denom_recall_total})")
    print()
    print(f"Precision = TP/(TP+FP) = {precision:.4f}")
    print(f"Recall    = TP/(TP+FN) = {recall:.4f}")
    print(f"F1-score               = {f1:.4f}")

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
        if args.match == "row":
            write_rows(os.path.join(args.outdir, "vrais_positifs.tsv"),
                       res["tp_rows"], with_hit=True)
            write_rows(os.path.join(args.outdir, "faux_positifs.tsv"),
                       res["fp_rows"], with_hit=False)
            write_keys(os.path.join(args.outdir, "faux_negatifs.tsv"),
                       res["fn_keys"], use_index=True)
        else:
            write_keys(os.path.join(args.outdir, "vrais_positifs.tsv"),
                       res["tp_keys"], use_index)
            write_keys(os.path.join(args.outdir, "faux_positifs.tsv"),
                       res["fp_keys"], use_index)
            write_keys(os.path.join(args.outdir, "faux_negatifs.tsv"),
                       res["fn_keys"], use_index)
        print(f"\nListes ecrites dans : {args.outdir}/")


if __name__ == "__main__":
    main()
