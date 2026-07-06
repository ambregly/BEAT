#!/usr/bin/env python3
"""
Visualisations de la comparaison fusions kmer vs vizome (BEAT AML), construites
DIRECTEMENT a partir des fichiers produits par compare_fusions.py :
    <resultats>/vrais_positifs.tsv   (TP)
    <resultats>/faux_positifs.tsv    (FP)
    <resultats>/faux_negatifs.tsv    (FN)
Les chiffres sont donc strictement coherents avec ces fichiers.

Produit (dans --outdir) :
  1. confusion_par_fusion.csv / .pdf       : une ligne PAR JONCTION UNIQUE
     (fusion + chromosomes + positions) -> chaque breakpoint distinct est compte
     separement ; la somme TP+FP+FN+TN d'une ligne vaut le nombre d'echantillons.
  2. confusion_par_echantillon.csv / .pdf  : une ligne PAR ECHANTILLON
     -> TP, FP, FN, TN, precision, recall, F1 ; TOUTES les lignes, paginees.
        Une ligne TOTAL (en tete) somme TP/FP/FN/TN et recalcule prec/recall/F1.
     (les Venn, scatter et histogramme F1 restent au niveau paire de genes)
  3. confusion.xlsx : classeur Excel regroupant les deux tableaux (un onglet chacun)
  4. venn_kmer_vizome.png : Venn des DETECTIONS (kmer = TP+FP, vizome = TP+FN,
     intersection = TP).
  5. venn_fusions_uniques.png : Venn des NOMS de fusion uniques (geneA_geneB
     compte une seule fois : kmer seul / vizome seul / les deux).

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
            # nom de fusion = paire de genes (majuscules) -> pour les REPRESENTATIONS
            r["fusion"] = f'{r["left_gene"].upper()}_{r["right_gene"].upper()}'
            # jonction = fusion + chromosomes + POSITIONS -> pour la TABLE de confusion
            # (chaque breakpoint distinct est une fusion unique)
            r["left_pos"] = r.get("left_pos", "")
            r["right_pos"] = r.get("right_pos", "")
            r["jonction"] = (r["fusion"], r["left_chr"], r["left_pos"],
                             r["right_chr"], r["right_pos"])
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
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
        records.append([k, TP, FP, FN, TN, round(precision, 4), round(recall, 4),
                        round(f1, 4)])
    # tri par activite decroissante
    records.sort(key=lambda x: (-(x[1] + x[2] + x[3]), x[0]))

    # ligne TOTAL : somme des TP/FP/FN/TN puis precision/recall/F1 recalcules
    sTP = sum(r[1] for r in records)
    sFP = sum(r[2] for r in records)
    sFN = sum(r[3] for r in records)
    sTN = sum(r[4] for r in records)
    tprec = sTP / (sTP + sFP) if (sTP + sFP) else 0.0
    trec = sTP / (sTP + sFN) if (sTP + sFN) else 0.0
    tf1 = 2 * tprec * trec / (tprec + trec) if (tprec + trec) else 0.0
    total_row = ["TOTAL", sTP, sFP, sFN, sTN,
                 round(tprec, 4), round(trec, 4), round(tf1, 4)]
    records = [total_row] + records   # TOTAL en premiere ligne

    header = [group_col, "TP", "FP", "FN", "TN", "precision", "recall", "F1"]
    return header, records


def build_confusion_junction(tp, fp, fn, n_samples):
    """Table de confusion regroupee par JONCTION UNIQUE (fusion + positions).

    Chaque breakpoint distinct est une fusion a part entiere ; la somme
    TP+FP+FN+TN de chaque ligne vaut n_samples.
    """
    def counts(rows):
        d = defaultdict(int)
        for r in rows:
            d[r["jonction"]] += 1
        return d

    ctp, cfp, cfn = counts(tp), counts(fp), counts(fn)
    distinct_samples = defaultdict(set)
    for rows in (tp, fp, fn):
        for r in rows:
            distinct_samples[r["jonction"]].add(r["SampleID"])

    keys = set(ctp) | set(cfp) | set(cfn)
    records = []
    for k in keys:
        fusion, lchr, lpos, rchr, rpos = k
        TP, FP, FN = ctp.get(k, 0), cfp.get(k, 0), cfn.get(k, 0)
        TN = n_samples - len(distinct_samples[k])
        precision = TP / (TP + FP) if (TP + FP) else 0.0
        recall = TP / (TP + FN) if (TP + FN) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
        records.append([fusion, lchr, lpos, rchr, rpos, TP, FP, FN, TN,
                        round(precision, 4), round(recall, 4), round(f1, 4)])
    # tri par activite decroissante puis nom de fusion
    records.sort(key=lambda x: (-(x[5] + x[6] + x[7]), x[0], x[2], x[4]))

    sTP = sum(r[5] for r in records); sFP = sum(r[6] for r in records)
    sFN = sum(r[7] for r in records); sTN = sum(r[8] for r in records)
    tprec = sTP / (sTP + sFP) if (sTP + sFP) else 0.0
    trec = sTP / (sTP + sFN) if (sTP + sFN) else 0.0
    tf1 = 2 * tprec * trec / (tprec + trec) if (tprec + trec) else 0.0
    total = ["TOTAL", "", "", "", "", sTP, sFP, sFN, sTN,
             round(tprec, 4), round(trec, 4), round(tf1, 4)]
    records = [total] + records

    header = ["fusion", "left_chr", "left_pos", "right_chr", "right_pos",
              "TP", "FP", "FN", "TN", "precision", "recall", "F1"]
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
    ax.set_title("Fusions kmer vs vizome\n(comptage des detections : TP/FP/FN)",
                 fontsize=12, fontweight="bold")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_venn_noms_uniques(path, tp, fp, fn):
    """Venn sur les NOMS de fusion uniques (geneA_geneB compte une seule fois).

    Logique paire-de-genes (inchangee) : une fusion est 'dans kmer' si elle
    apparait en TP ou FP, 'dans vizome' si elle apparait en TP ou FN.
    """
    kmer_names = {r["fusion"] for r in tp} | {r["fusion"] for r in fp}
    vizome_names = {r["fusion"] for r in tp} | {r["fusion"] for r in fn}
    only_kmer = len(kmer_names - vizome_names)
    only_vizome = len(vizome_names - kmer_names)
    both = len(kmer_names & vizome_names)

    fig, ax = plt.subplots(figsize=(7, 6))
    v = venn2(subsets=(only_kmer, only_vizome, both),
              set_labels=("Fusions kmer", "Fusions vizome"), ax=ax)
    for region, color in (("10", "#66c2a5"), ("01", "#fc8d62"), ("11", "#8da0cb")):
        if v.get_patch_by_id(region):
            v.get_patch_by_id(region).set_color(color)
            v.get_patch_by_id(region).set_alpha(0.7)
    ax.set_title("Noms de fusion uniques (geneA_geneB)\n"
                 "1 fusion = 1 fois (sans dedoublement)",
                 fontsize=12, fontweight="bold")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return only_kmer, both, only_vizome


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
    n_junctions = len({r["jonction"] for rows in (tp, fp, fn) for r in rows})

    # 1. par fusion = par JONCTION UNIQUE (fusion + positions) ; somme/ligne = n_samples
    h_f, rec_f = build_confusion_junction(tp, fp, fn, n_samples)
    write_csv(os.path.join(args.outdir, "confusion_par_fusion.csv"), h_f, rec_f)
    npg_f = write_table_pdf(os.path.join(args.outdir, "tableau_confusion_par_fusion.pdf"),
                            h_f, rec_f, "Confusion par fusion (jonction unique)",
                            args.rows_per_page)

    # 2. par echantillon (autre dimension = jonction ; somme/ligne = n_junctions)
    h_s, rec_s = build_confusion(tp, fp, fn, "SampleID", "jonction", n_junctions)
    write_csv(os.path.join(args.outdir, "confusion_par_echantillon.csv"), h_s, rec_s)
    npg_s = write_table_pdf(os.path.join(args.outdir, "tableau_confusion_par_echantillon.pdf"),
                            h_s, rec_s, "Confusion par echantillon", args.rows_per_page)

    # classeur Excel regroupant les deux tableaux (un onglet chacun)
    xlsx_path = os.path.join(args.outdir, "confusion.xlsx")
    write_xlsx(xlsx_path, [("par_fusion", h_f, rec_f),
                           ("par_echantillon", h_s, rec_s)])

    # 3. venn des detections (TP/FP/FN)
    venn_path = os.path.join(args.outdir, "venn_kmer_vizome.png")
    plot_venn(venn_path, n_tp, n_fp, n_fn)

    # 4. venn des noms de fusion uniques (geneA_geneB compte une seule fois)
    venn_noms = os.path.join(args.outdir, "venn_fusions_uniques.png")
    ok, both, ov = plot_venn_noms_uniques(venn_noms, tp, fp, fn)

    print(f"TP={n_tp}  FP={n_fp}  FN={n_fn}")
    print(f"Precision = {n_tp/(n_tp+n_fp):.4f} | Recall = {n_tp/(n_tp+n_fn):.4f}")
    print(f"Echantillons : {n_samples} | Paires de genes : {n_fusions} | "
          f"Jonctions uniques : {n_junctions}")
    print(f"Noms de fusion uniques : kmer seul={ok}, communs={both}, vizome seul={ov}")
    print("Fichiers ecrits :")
    print(f"  - {args.outdir}/confusion_par_fusion.csv")
    print(f"  - {args.outdir}/confusion_par_echantillon.csv")
    print(f"  - {xlsx_path}  (onglets : par_fusion, par_echantillon)")
    print(f"  - {args.outdir}/tableau_confusion_par_fusion.pdf ({npg_f} pages)")
    print(f"  - {args.outdir}/tableau_confusion_par_echantillon.pdf ({npg_s} pages)")
    print(f"  - {venn_path}")
    print(f"  - {venn_noms}")


if __name__ == "__main__":
    main()
