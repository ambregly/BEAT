#!/usr/bin/env python3
"""
Graphiques supplementaires sur la comparaison fusions kmer vs vizome (BEAT AML).

A partir des fichiers SOURCES (pour acceder aux comptages), produit :
  1. faux_negatifs_vizome.csv : fusions presentes uniquement dans vizome (FN),
     enrichies du comptage vizome.
  2. scatter_kmer_vs_vizome.png : sur l'ensemble des VRAIS POSITIFS, nuage de
     points (1 point = 1 paire fusion/echantillon) du comptage kmer (X) vs
     comptage vizome (Y), avec droite de regression et correlations.
     Donnees brutes dans tp_comptages.csv.
  3. hist_f1_par_fusion.png : barplot de repartition des scores F1 par fusion.

Le matching (TP/FP/FN) est recalcule exactement comme compare_fusions.py
(mode 'row' par defaut), avec l'option --common-samples.

Usage :
  python extra_plots.py BEAT_AML2.csv kmer2.tsv --outdir figures --common-samples
  python extra_plots.py BEAT_AML2.csv kmer2.tsv --outdir figures \
      --vizome-count junction_read_count
"""

import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import compare_fusions as cf

VIZOME_COUNT_CHOICES = ["junction_read_count", "spanning_frag_count",
                        "junction_plus_spanning", "ffpm"]


def vizome_count(beat_row, mode):
    """Extrait le comptage vizome d'une ligne BEAT AML selon le mode choisi."""
    def num(col):
        try:
            return float(beat_row.get(col, "") or 0)
        except ValueError:
            return 0.0
    if mode == "junction_plus_spanning":
        return num("junction_read_count") + num("spanning_frag_count")
    return num(mode)


def fusion_label(lg, rg):
    return f"{cf.norm(lg)}_{cf.norm(rg)}"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("beat_aml")
    p.add_argument("kmer")
    p.add_argument("--outdir", default="figures")
    p.add_argument("--match", choices=["row", "index"], default="row",
                   help="critere de matching (defaut row). 'genepair' non gere ici "
                        "car on a besoin des index pour relier les comptages.")
    p.add_argument("--common-samples", action="store_true")
    p.add_argument("--vizome-count", choices=VIZOME_COUNT_CHOICES,
                   default="junction_read_count",
                   help="colonne BEAT AML utilisee comme comptage vizome (defaut "
                        "junction_read_count)")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ---- chargement + matching (comme compare_fusions, use_index=True) -------
    beat = cf.load_beat_aml(args.beat_aml, use_index=True)   # {key: [rows]}
    beat_keys = set(beat)
    kmer_rows = cf.parse_kmer_rows(args.kmer)

    if args.common_samples:
        bs = {k[0] for k in beat_keys}
        ks = {r["sample_id"] for r in kmer_rows}
        common = bs & ks
        beat_keys = {k for k in beat_keys if k[0] in common}
        kmer_rows = [r for r in kmer_rows if r["sample_id"] in common]

    res = cf.compare_row(beat_keys, kmer_rows)

    # ---- 1. liste des faux negatifs (uniquement vizome), enrichie -----------
    fn_path = os.path.join(args.outdir, "faux_negatifs_vizome.csv")
    with open(fn_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["SampleID", "fusion", "left_gene", "left_chr",
                    "right_gene", "right_chr", "fusion_index", args.vizome_count])
        for key in sorted(res["fn_keys"]):
            sample, lg, lc, rg, rc, idx = key
            brow = beat[key][0]
            w.writerow([sample, fusion_label(lg, rg), lg, lc, rg, rc, idx,
                        vizome_count(brow, args.vizome_count)])
    n_fn = len(res["fn_keys"])

    # ---- 2. comptages sur les vrais positifs --------------------------------
    # un point = une paire (echantillon, fusion) ; on somme les comptages des TP
    pair_kmer = defaultdict(float)
    pair_vizome = defaultdict(float)
    for row, hits in res["tp_rows"]:
        hit = hits[0]                      # 'on oublie les autres' -> 1er match
        brow = beat[hit][0]
        pair = (row["sample_id"], fusion_label(row["left_gene"], row["right_gene"]))
        if row["count"] is not None:
            pair_kmer[pair] += row["count"]
        pair_vizome[pair] += vizome_count(brow, args.vizome_count)

    pairs = sorted(set(pair_kmer) | set(pair_vizome))
    tp_csv = os.path.join(args.outdir, "tp_comptages.csv")
    with open(tp_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["SampleID", "fusion", "comptage_kmer", "comptage_vizome"])
        for pr in pairs:
            w.writerow([pr[0], pr[1], pair_kmer.get(pr, 0), pair_vizome.get(pr, 0)])

    x = np.array([pair_kmer.get(pr, 0) for pr in pairs], dtype=float)
    y = np.array([pair_vizome.get(pr, 0) for pr in pairs], dtype=float)

    scatter_path = os.path.join(args.outdir, "scatter_kmer_vs_vizome.png")
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.scatter(x, y, s=14, alpha=0.5, color="#3b7dd8", edgecolor="none")
    if len(x) >= 2 and np.ptp(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 100)
        ax.plot(xs, slope * xs + intercept, color="#d83b3b", lw=2,
                label=f"y = {slope:.2f}x + {intercept:.2f}")
        r, p_r = stats.pearsonr(x, y)
        rho, p_s = stats.spearmanr(x, y)
        ax.legend(loc="upper left")
        ax.text(0.02, 0.88,
                f"Pearson r = {r:.3f}\nSpearman rho = {rho:.3f}\nn = {len(x)} paires",
                transform=ax.transAxes, va="top",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax.set_xlabel("Comptage kmer (3e colonne)")
    ax.set_ylabel(f"Comptage vizome ({args.vizome_count})")
    ax.set_title("Vrais positifs : comptage kmer vs vizome\n(1 point = 1 paire fusion/echantillon)")
    fig.tight_layout()
    fig.savefig(scatter_path, dpi=150)
    plt.close(fig)

    # ---- 3. histogramme des scores F1 par fusion ----------------------------
    tp_by_fus = defaultdict(int); fp_by_fus = defaultdict(int); fn_by_fus = defaultdict(int)
    for row, _ in res["tp_rows"]:
        tp_by_fus[fusion_label(row["left_gene"], row["right_gene"])] += 1
    for row in res["fp_rows"]:
        fp_by_fus[fusion_label(row["left_gene"], row["right_gene"])] += 1
    for key in res["fn_keys"]:
        fn_by_fus[fusion_label(key[1], key[3])] += 1

    all_fus = set(tp_by_fus) | set(fp_by_fus) | set(fn_by_fus)
    f1_values = []
    for f in all_fus:
        TP, FP, FN = tp_by_fus[f], fp_by_fus[f], fn_by_fus[f]
        prec = TP / (TP + FP) if (TP + FP) else 0.0
        rec = TP / (TP + FN) if (TP + FN) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        f1_values.append(f1)
    f1_values = np.array(f1_values)

    hist_path = os.path.join(args.outdir, "hist_f1_par_fusion.png")
    bins = np.linspace(0, 1, 11)  # tranches de 0.1
    counts, edges = np.histogram(f1_values, bins=bins)
    fig, ax = plt.subplots(figsize=(8, 5))
    centers = (edges[:-1] + edges[1:]) / 2
    ax.bar(centers, counts, width=0.09, color="#5aa469", edgecolor="white")
    ax.set_xticks(np.round(bins, 1))
    ax.set_xlabel("Score F1 (par fusion)")
    ax.set_ylabel("Nombre de fusions")
    ax.set_title(f"Repartition des scores F1 par fusion (n = {len(f1_values)})")
    for c, ce in zip(counts, centers):
        if c > 0:
            ax.text(ce, c, str(c), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(hist_path, dpi=150)
    plt.close(fig)

    # ---- resume -------------------------------------------------------------
    print(f"Vrais positifs (lignes kmer) : {res['tp']}")
    print(f"Paires fusion/echantillon (points du scatter) : {len(pairs)}")
    print(f"Faux negatifs (uniquement vizome) : {n_fn}")
    print(f"Fusions (pour F1) : {len(f1_values)} | F1 median = {np.median(f1_values):.3f}")
    print("Fichiers ecrits :")
    print(f"  - {fn_path}")
    print(f"  - {tp_csv}")
    print(f"  - {scatter_path}")
    print(f"  - {hist_path}")


if __name__ == "__main__":
    main()
