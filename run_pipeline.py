#!/usr/bin/env python3
"""
Orchestrateur : enchaine compare_fusions.py puis visualize_fusions.py pour un
(ou plusieurs) jeu(x) de donnees, afin de comparer facilement des parametres.

Pour chaque jeu, cree l'arborescence :
    <outdir>/<label>/resultats/   (TP/FP/FN)
    <outdir>/<label>/figures/     (tableaux CSV/Excel/PDF, Venn, scatter, hist F1,
                                   liste FN enrichie)
Enchaine compare_fusions.py -> visualize_fusions.py -> extra_plots.py.

Deux modes :

1) Un seul jeu (en ligne de commande) :
    python run_pipeline.py --beat BEAT_AML2.csv --kmer kmer2.tsv \
        --label run1 --match row --rows-per-page 30

2) Plusieurs jeux via un manifeste TSV (--manifest). Colonnes attendues
   (entete obligatoire ; colonnes optionnelles si absentes -> valeurs par defaut) :
       label  beat  kmer  match  rows_per_page  vizome_count  kmer_normal
   exemple de ligne :
       run1   BEAT_AML2.csv  kmer2.tsv  row  30  junction_read_count  kmer_normal.tsv
   (kmer_normal facultatif : laisser vide pour ne pas filtrer)

    python run_pipeline.py --manifest jeux.tsv --outdir analyses
"""

import argparse
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COMPARE = os.path.join(HERE, "compare_fusions.py")
VISUALIZE = os.path.join(HERE, "visualize_fusions.py")
EXTRA = os.path.join(HERE, "extra_plots.py")


def run(cmd):
    print("  $ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def process(label, beat, kmer, outdir, match, rows_per_page, vizome_count,
            kmer_normal=None):
    base = os.path.join(outdir, label)
    res_dir = os.path.join(base, "resultats")
    fig_dir = os.path.join(base, "figures")
    os.makedirs(res_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    norm_info = f" kmer_normal={kmer_normal}" if kmer_normal else ""
    print(f"\n=== [{label}] beat={beat} kmer={kmer} match={match}{norm_info} ===")

    extra_normal = ["--kmer-normal", kmer_normal] if kmer_normal else []

    # 1. comparaison -> TP/FP/FN
    run([sys.executable, COMPARE, beat, kmer, "--outdir", res_dir,
         "--match", match] + extra_normal)

    # 2. tableaux + Venn (depuis le dossier resultats)
    run([sys.executable, VISUALIZE, res_dir, "--outdir", fig_dir,
         "--rows-per-page", str(rows_per_page)])

    # 3. graphes supplementaires (depuis les fichiers sources, besoin des index)
    if match in ("row", "index"):
        run([sys.executable, EXTRA, beat, kmer, "--outdir", fig_dir,
             "--match", match, "--vizome-count", vizome_count] + extra_normal)
    else:
        print("  (extra_plots ignore : --match genepair n'a pas d'index pour "
              "relier les comptages)")

    print(f"=== [{label}] termine -> {base}/ ===")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", help="fichier TSV listant plusieurs jeux (voir aide)")
    p.add_argument("--beat", help="fichier BEAT_AML.csv (mode jeu unique)")
    p.add_argument("--kmer", help="fichier kmer.tsv (mode jeu unique)")
    p.add_argument("--label", default="run", help="nom du jeu (mode jeu unique)")
    p.add_argument("--match", default="row", choices=["row", "index", "genepair"])
    p.add_argument("--rows-per-page", type=int, default=40)
    p.add_argument("--vizome-count", default="junction_read_count",
                   choices=["junction_read_count", "spanning_frag_count",
                            "junction_plus_spanning", "ffpm"],
                   help="comptage vizome pour le nuage de points (defaut "
                        "junction_read_count)")
    p.add_argument("--kmer-normal", default=None,
                   help="fichier kmer normal : retire ces paires de genes de BEAT "
                        "et kmer (liste noire d'artefacts)")
    p.add_argument("--outdir", default="analyses", help="repertoire racine de sortie")
    args = p.parse_args()

    if args.manifest:
        with open(args.manifest, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            jobs = list(reader)
        if not jobs:
            sys.exit("manifeste vide")
        for j in jobs:
            kn = (j.get("kmer_normal") or "").strip() or args.kmer_normal
            process(
                label=j["label"].strip(),
                beat=j["beat"].strip(),
                kmer=j["kmer"].strip(),
                outdir=args.outdir,
                match=j.get("match", "row").strip() or "row",
                rows_per_page=int(j.get("rows_per_page") or args.rows_per_page),
                vizome_count=(j.get("vizome_count") or args.vizome_count).strip(),
                kmer_normal=kn or None,
            )
    else:
        if not (args.beat and args.kmer):
            sys.exit("Mode jeu unique : --beat et --kmer sont requis "
                     "(ou utilisez --manifest).")
        process(args.label, args.beat, args.kmer, args.outdir,
                args.match, args.rows_per_page, args.vizome_count,
                kmer_normal=args.kmer_normal)

    print(f"\nTout est ecrit sous : {args.outdir}/")


if __name__ == "__main__":
    main()
