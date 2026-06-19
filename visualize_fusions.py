#!/usr/bin/env python3
"""
Visualisations de la comparaison fusions kmer vs vizome (BEAT AML), construites
DIRECTEMENT a partir des fichiers produits par compare_fusions.py :
    <resultats>/vrais_positifs.tsv   (TP)
    <resultats>/faux_positifs.tsv    (FP)
    <resultats>/faux_negatifs.tsv    (FN)
Les chiffres sont donc strictement coherents avec ces fichiers.

Produit (dans --outdir) :
  1. confusion_par_fusion.csv / .pdf       : une ligne PAR FUSION (paire de genes)
  2. confusion_par_echantillon.csv / .pdf  : une ligne PAR ECHANTILLON
     -> TP, FP, FN, TN, precision, recall ; TOUTES les lignes, paginees.
  3. confusion.xlsx : classeur Excel regroupant les deux tableaux (un onglet chacun)
  4. venn_kmer_vizome.png : Venn (kmer = TP+FP, vizome = TP+FN, intersection = TP)

Definitions (granularite = ligne de fichier, comme compare_fusions.py) :
  Par fusion      : TP/FP/FN = nb de lignes de la fusion dans chaque fichier ;
                    TN = nb total d'echantillons - nb d'echantillons de la fusion.
  Par echantillon : TP/FP/FN = nb de lignes de l'echantillon dans chaque fichier ;
                    TN = nb total de fusions - nb de fusions de l'echantillon.

Usage :
    python visualize_fusions.py <dossier_resultats> --outdir figures
    python visualize_fusions.py resultats --outdir figures --rows-per-page 40
"""

import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn2

COMMON_COLS = ["SampleID", "left_gene", "left_chr", "right_gene", "right_chr"]


def read_tsv(path):
    if not os.path.exists(path):
        raise SystemExit(f"fichier introuvable : {path}")
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        miss = set(COMMON_COLS) - set(reader.fieldnames or [])
        if miss:
            raise SystemExit(f"[{path}] colonnes manquantes : {sorted(miss)}")
        rows = []
        for r in reader:
            r["fusion"] = f'{r["left_gene"]}_{r["right_gene"]}'
            rows.append(r)
        return rows


def build_confusion(tp, fp, fn, group_col, other_col, universe):
    """Tableau de confusion regroupe par 'group_col'. TN base sur 'universe'."""
    def counts(rows):
        d = defaultdict(int)
        for r in rows:
            d[r[group_col]] += 1
        return d

    ctp, cfp, cfn = counts(tp), counts(fp), counts(fn)

    distinct_other = defaultdict(set)
    for rows in (tp, fp, fn):
        for r in rows:
            distinct_other[r[group_col]].add(r[other_col])

    keys = sorted(set(ctp) | set(cfp) | set(cfn))
    records = []
    for k in keys:
        TP, FP, FN = ctp.get(k, 0), cfp.get(k, 0), cfn.get(k, 0)
        TN = universe - len(distinct_other[k])
        precision = TP / (TP + FP) if (TP + FP) else 0.0
        recall = TP / (TP + FN) if (TP + FN) else 0.0
        records.append([k, TP, FP, FN, TN, round(precision, 4), round(recall, 4)])
    # tri par activite decroissante
    records.sort(key=lambda x: (-(x[1] + x[2] + x[3]), x[0]))
    header = [group_col, "TP", "FP", "FN", "TN", "precision", "recall"]
    return header, records


def write_csv(path, header, records):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(records)


def write_xlsx(path, sheets):
    """Ecrit un classeur Excel. 'sheets' = liste de (nom, header, records)."""
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment

    wb = Workbook()
    wb.remove(wb.active)
    head_fill = PatternFill("solid", fgColor="40466E")
    head_font = Font(color="FFFFFF", bold=True)
    for name, header, records in sheets:
        ws = wb.create_sheet(title=name[:31])  # Excel limite a 31 caracteres
        ws.append(header)
        for cell in ws[1]:
            cell.fill = head_fill
            cell.font = head_font
            cell.alignment = Alignment(horizontal="center")
        for rec in records:
            ws.append(rec)
        ws.freeze_panes = "A2"  # fige la ligne d'entete
        # largeur de colonnes approximative
        for j, col in enumerate(header, start=1):
            width = max(len(str(col)),
                        *(len(str(r[j - 1])) for r in records)) + 2 if records else len(col) + 2
            ws.column_dimensions[ws.cell(row=1, column=j).column_letter].width = min(width, 40)
    wb.save(path)


def write_table_pdf(path, header, records, titre, rows_per_page):
    n = len(records)
    n_pages = max(1, (n + rows_per_page - 1) // rows_per_page)
    with PdfPages(path) as pdf:
        for pg in range(n_pages):
            i0, i1 = pg * rows_per_page, min((pg + 1) * rows_per_page, n)
            sub = records[i0:i1]
            fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 portrait
            ax.axis("off")
            ax.set_title(f"{titre} - page {pg + 1}/{n_pages} "
                         f"(lignes {i0 + 1}-{i1} sur {n})",
                         fontsize=12, fontweight="bold")
            table = ax.table(cellText=sub, colLabels=header,
                             cellLoc="center", loc="upper center")
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.25)
            for j in range(len(header)):
                c = table[0, j]
                c.set_facecolor("#40466e")
                c.set_text_props(color="white", fontweight="bold")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    return n_pages


def plot_venn(path, n_tp, n_fp, n_fn):
    fig, ax = plt.subplots(figsize=(7, 6))
    v = venn2(subsets=(n_fp, n_fn, n_tp),
              set_labels=("Fusions kmer", "Fusions vizome"), ax=ax)
    for region, color in (("10", "#66c2a5"), ("01", "#fc8d62"), ("11", "#8da0cb")):
        if v.get_patch_by_id(region):
            v.get_patch_by_id(region).set_color(color)
            v.get_patch_by_id(region).set_alpha(0.7)
    ax.set_title("Fusions kmer vs vizome", fontsize=13, fontweight="bold")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("resultats", help="dossier contenant vrais/faux_positifs/negatifs.tsv")
    parser.add_argument("--outdir", default="figures", help="repertoire de sortie")
    parser.add_argument("--rows-per-page", type=int, default=40,
                        help="lignes par page dans les PDF (defaut 40)")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    tp = read_tsv(os.path.join(args.resultats, "vrais_positifs.tsv"))
    fp = read_tsv(os.path.join(args.resultats, "faux_positifs.tsv"))
    fn = read_tsv(os.path.join(args.resultats, "faux_negatifs.tsv"))
    n_tp, n_fp, n_fn = len(tp), len(fp), len(fn)

    n_samples = len({r["SampleID"] for rows in (tp, fp, fn) for r in rows})
    n_fusions = len({r["fusion"] for rows in (tp, fp, fn) for r in rows})

    # 1. par fusion
    h_f, rec_f = build_confusion(tp, fp, fn, "fusion", "SampleID", n_samples)
    write_csv(os.path.join(args.outdir, "confusion_par_fusion.csv"), h_f, rec_f)
    npg_f = write_table_pdf(os.path.join(args.outdir, "tableau_confusion_par_fusion.pdf"),
                            h_f, rec_f, "Confusion par fusion", args.rows_per_page)

    # 2. par echantillon
    h_s, rec_s = build_confusion(tp, fp, fn, "SampleID", "fusion", n_fusions)
    write_csv(os.path.join(args.outdir, "confusion_par_echantillon.csv"), h_s, rec_s)
    npg_s = write_table_pdf(os.path.join(args.outdir, "tableau_confusion_par_echantillon.pdf"),
                            h_s, rec_s, "Confusion par echantillon", args.rows_per_page)

    # classeur Excel regroupant les deux tableaux (un onglet chacun)
    xlsx_path = os.path.join(args.outdir, "confusion.xlsx")
    write_xlsx(xlsx_path, [("par_fusion", h_f, rec_f),
                           ("par_echantillon", h_s, rec_s)])

    # 3. venn
    venn_path = os.path.join(args.outdir, "venn_kmer_vizome.png")
    plot_venn(venn_path, n_tp, n_fp, n_fn)

    print(f"TP={n_tp}  FP={n_fp}  FN={n_fn}")
    print(f"Precision = {n_tp/(n_tp+n_fp):.4f} | Recall = {n_tp/(n_tp+n_fn):.4f}")
    print(f"Echantillons : {n_samples} | Fusions : {n_fusions}")
    print("Fichiers ecrits :")
    print(f"  - {args.outdir}/confusion_par_fusion.csv")
    print(f"  - {args.outdir}/confusion_par_echantillon.csv")
    print(f"  - {xlsx_path}  (onglets : par_fusion, par_echantillon)")
    print(f"  - {args.outdir}/tableau_confusion_par_fusion.pdf ({npg_f} pages)")
    print(f"  - {args.outdir}/tableau_confusion_par_echantillon.pdf ({npg_s} pages)")
    print(f"  - {venn_path}")


if __name__ == "__main__":
    main()
