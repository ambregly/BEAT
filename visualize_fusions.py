#!/usr/bin/env python3
"""
Visualisations de la comparaison fusions kmer vs BEAT AML.

Produit deux sorties :
  1. Un tableau de confusion PAR FUSION (paire de genes) : TP, FP, FN, TN,
     precision, recall. Ecrit en CSV (toutes les fusions) et en image PNG
     (les --top fusions les plus actives).
  2. Un diagramme de Venn : fusions detectees par kmer, fusions de BEAT AML,
     et leur intersection.

Granularite : niveau (paire de genes x echantillon).
  - TP : echantillons ou la fusion est dans BEAT AML ET detectee par kmer
  - FP : echantillons ou kmer detecte la fusion mais absente de BEAT AML
  - FN : echantillons ou la fusion est dans BEAT AML mais non detectee par kmer
  - TN : echantillons ou la fusion n'est ni dans BEAT AML ni dans kmer
         (TN = nombre total d'echantillons - TP - FP - FN)

Usage :
    python visualize_fusions.py BEAT_AML.csv kmer.tsv --outdir figures
    python visualize_fusions.py BEAT_AML.csv kmer.tsv --outdir figures --top 30
"""

import argparse
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")  # backend sans affichage (sauvegarde fichier)
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2

import compare_fusions as cf


def build_genepair_sets(beat_path, kmer_path):
    """Construit, par paire de genes, les ensembles d'echantillons BEAT et kmer.

    Renvoie :
        beat_samples : dict {gene_pair_key: set(samples)}
        kmer_samples : dict {gene_pair_key: set(samples)}
        labels       : dict {gene_pair_key: "LEFT_RIGHT"}
        all_samples  : set de tous les echantillons observes (BEAT U kmer)
    """
    beat_samples = defaultdict(set)
    kmer_samples = defaultdict(set)
    labels = {}
    all_samples = set()

    # BEAT AML : cle paire de genes = (lg, lc, rg, rc) ; sample en colonne SampleID
    beat = cf.load_beat_aml(beat_path, use_index=False)
    for key in beat:                       # key = (sample, lg, lc, rg, rc)
        sample, lg, lc, rg, rc = key
        gp = (lg, lc, rg, rc)
        beat_samples[gp].add(sample)
        labels.setdefault(gp, f"{lg}_{rg}")
        all_samples.add(sample)

    # kmer : une ligne = une detection d'une paire de genes dans un echantillon
    for row in cf.parse_kmer_rows(kmer_path):
        gp = (cf.norm(row["left_gene"]), cf.norm_chr(row["left_chr"]),
              cf.norm(row["right_gene"]), cf.norm_chr(row["right_chr"]))
        kmer_samples[gp].add(row["sample_id"])
        labels.setdefault(gp, f"{cf.norm(row['left_gene'])}_{cf.norm(row['right_gene'])}")
        all_samples.add(row["sample_id"])

    return beat_samples, kmer_samples, labels, all_samples


def confusion_per_fusion(beat_samples, kmer_samples, labels, all_samples):
    """Tableau (DataFrame) de confusion par paire de genes."""
    n_total = len(all_samples)
    records = []
    for gp in set(beat_samples) | set(kmer_samples):
        b = beat_samples.get(gp, set())
        k = kmer_samples.get(gp, set())
        tp = len(b & k)
        fp = len(k - b)
        fn = len(b - k)
        tn = n_total - (tp + fp + fn)
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        records.append({
            "fusion": labels[gp],
            "left_chr": gp[1], "right_chr": gp[3],
            "TP": tp, "FP": fp, "FN": fn, "TN": tn,
            "precision": round(precision, 4),
            "recall": round(recall, 4),
        })
    df = pd.DataFrame.from_records(records)
    # tri par activite (fusions les plus presentes en premier)
    df["_activite"] = df["TP"] + df["FP"] + df["FN"]
    df = df.sort_values(["_activite", "fusion"], ascending=[False, True])
    df = df.drop(columns="_activite").reset_index(drop=True)
    return df


def plot_table(df, path, top):
    """Rend les 'top' premieres fusions sous forme de tableau image."""
    cols = ["fusion", "left_chr", "right_chr", "TP", "FP", "FN", "TN",
            "precision", "recall"]
    sub = df[cols].head(top)

    fig_h = 1.0 + 0.32 * len(sub)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.axis("off")
    ax.set_title(f"Confusion par fusion (top {len(sub)} sur {len(df)})",
                 fontsize=13, pad=12)

    table = ax.table(cellText=sub.values, colLabels=cols,
                     cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.3)

    # en-tete en gras + couleur
    for j in range(len(cols)):
        cell = table[0, j]
        cell.set_facecolor("#40466e")
        cell.set_text_props(color="white", fontweight="bold")

    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_venn(beat_path, kmer_path, path, use_index):
    """Diagramme de Venn des fusions kmer vs BEAT AML."""
    beat = set(cf.load_beat_aml(beat_path, use_index=use_index))
    kmer = set()
    for row in cf.parse_kmer_rows(kmer_path):
        kmer.update(cf.row_candidate_keys(row, use_index=use_index))

    only_kmer = len(kmer - beat)
    only_beat = len(beat - kmer)
    common = len(kmer & beat)

    fig, ax = plt.subplots(figsize=(7, 6))
    v = venn2(subsets=(only_kmer, only_beat, common),
              set_labels=("Fusions kmer", "Fusions BEAT AML"), ax=ax)
    for region, color in (("10", "#66c2a5"), ("01", "#fc8d62"), ("11", "#8da0cb")):
        if v.get_patch_by_id(region):
            v.get_patch_by_id(region).set_color(color)
            v.get_patch_by_id(region).set_alpha(0.7)
    niveau = "echantillon + fusion_index" if use_index else "echantillon + paire de genes"
    ax.set_title(f"Fusions kmer vs BEAT AML\n(cle : {niveau})", fontsize=12)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return only_kmer, common, only_beat


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("beat_aml", help="fichier BEAT_AML.csv")
    parser.add_argument("kmer", help="fichier kmer (tsv)")
    parser.add_argument("--outdir", default="figures",
                        help="repertoire de sortie (defaut: figures)")
    parser.add_argument("--top", type=int, default=25,
                        help="nombre de fusions affichees dans l'image du tableau")
    parser.add_argument("--venn-match", choices=["index", "genepair"],
                        default="genepair",
                        help="cle du diagramme de Venn (defaut: genepair)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Tableau de confusion par fusion
    beat_s, kmer_s, labels, all_samples = build_genepair_sets(args.beat_aml, args.kmer)
    df = confusion_per_fusion(beat_s, kmer_s, labels, all_samples)

    csv_path = os.path.join(args.outdir, "confusion_par_fusion.csv")
    df.to_csv(csv_path, index=False)

    table_png = os.path.join(args.outdir, "tableau_confusion.png")
    plot_table(df, table_png, args.top)

    # 2. Diagramme de Venn
    venn_png = os.path.join(args.outdir, "venn_kmer_beataml.png")
    only_kmer, common, only_beat = plot_venn(
        args.beat_aml, args.kmer, venn_png, use_index=(args.venn_match == "index"))

    print(f"Echantillons consideres (BEAT U kmer) : {len(all_samples)}")
    print(f"Fusions (paires de genes) distinctes  : {len(df)}")
    print()
    print("Fichiers ecrits :")
    print(f"  - {csv_path}        (tableau complet par fusion)")
    print(f"  - {table_png}       (image top {args.top})")
    print(f"  - {venn_png}        (Venn : kmer-only={only_kmer}, "
          f"commun={common}, BEAT-only={only_beat})")


if __name__ == "__main__":
    main()
