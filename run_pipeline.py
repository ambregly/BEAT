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
        --label run_common --match row --common-samples --rows-per-page 30

2) Plusieurs jeux via un manifeste TSV (--manifest). Colonnes attendues
   (entete obligatoire ; colonnes optionnelles si absentes -> valeurs par defaut) :
       label  beat  kmer  match  common_samples  rows_per_page  vizome_count
   exemple de ligne :
       run1   BEAT_AML2.csv  kmer2.tsv  row  yes  30  junction_read_count
   (common_samples : yes/true/1 pour activer)

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


def process(label, beat, kmer, outdir, match, common_samples, rows_per_page,
            vizome_count):
    base = os.path.join(outdir, label)
    res_dir = os.path.join(base, "resultats")
    fig_dir = os.path.join(base, "figures")
    os.makedirs(res_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    print(f"\n=== [{label}] beat={beat} kmer={kmer} "
          f"match={match} common_samples={common_samples} ===")

    # 1. comparaison -> TP/FP/FN
    cmp_cmd = [sys.executable, COMPARE, beat, kmer, "--outdir", res_dir,
               "--match", match]
    if common_samples:
        cmp_cmd.append("--common-samples")
    run(cmp_cmd)

    # 2. tableaux + Venn (depuis le dossier resultats)
    run([sys.executable, VISUALIZE, res_dir, "--outdir", fig_dir,
         "--rows-per-page", str(rows_per_page)])

    # 3. graphes supplementaires (depuis les fichiers sources, besoin des index)
    if match in ("row", "index"):
        extra_cmd = [sys.executable, EXTRA, beat, kmer, "--outdir", fig_dir,
                     "--match", match, "--vizome-count", vizome_count]
        if common_samples:
            extra_cmd.append("--common-samples")
        run(extra_cmd)
    else:
        print("  (extra_plots ignore : --match genepair n'a pas d'index pour "
              "relier les comptages)")

    print(f"=== [{label}] termine -> {base}/ ===")


def truthy(v):
    return str(v).strip().lower() in ("yes", "true", "1", "oui", "y")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", help="fichier TSV listant plusieurs jeux (voir aide)")
    p.add_argument("--beat", help="fichier BEAT_AML.csv (mode jeu unique)")
    p.add_argument("--kmer", help="fichier kmer.tsv (mode jeu unique)")
    p.add_argument("--label", default="run", help="nom du jeu (mode jeu unique)")
    p.add_argument("--match", default="row", choices=["row", "index", "genepair"])
    p.add_argument("--common-samples", action="store_true",
                   help="restreindre aux echantillons communs BEAT/kmer")
    p.add_argument("--rows-per-page", type=int, default=40)
    p.add_argument("--vizome-count", default="junction_read_count",
                   choices=["junction_read_count", "spanning_frag_count",
                            "junction_plus_spanning", "ffpm"],
                   help="comptage vizome pour le nuage de points (defaut "
                        "junction_read_count)")
    p.add_argument("--outdir", default="analyses", help="repertoire racine de sortie")
    args = p.parse_args()

    if args.manifest:
        with open(args.manifest, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            jobs = list(reader)
        if not jobs:
            sys.exit("manifeste vide")
        for j in jobs:
            process(
                label=j["label"].strip(),
                beat=j["beat"].strip(),
                kmer=j["kmer"].strip(),
                outdir=args.outdir,
                match=j.get("match", "row").strip() or "row",
                common_samples=truthy(j.get("common_samples", "no")),
                rows_per_page=int(j.get("rows_per_page") or args.rows_per_page),
                vizome_count=(j.get("vizome_count") or args.vizome_count).strip(),
            )
    else:
        if not (args.beat and args.kmer):
            sys.exit("Mode jeu unique : --beat et --kmer sont requis "
                     "(ou utilisez --manifest).")
        process(args.label, args.beat, args.kmer, args.outdir,
                args.match, args.common_samples, args.rows_per_page,
                args.vizome_count)

    print(f"\nTout est ecrit sous : {args.outdir}/")


if __name__ == "__main__":
    main()
