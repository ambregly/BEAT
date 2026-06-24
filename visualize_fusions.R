#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Visualisations de la comparaison fusions kmer vs vizome (BEAT AML), construites
# DIRECTEMENT a partir des fichiers produits par compare_fusions.py :
#     <resultats>/vrais_positifs.tsv   (TP)
#     <resultats>/faux_positifs.tsv    (FP)
#     <resultats>/faux_negatifs.tsv    (FN)
# Les chiffres sont donc strictement coherents avec ces fichiers.
#
# Produit :
#   1. confusion_par_fusion.csv / .pdf       : une ligne PAR FUSION (paire de genes)
#   2. confusion_par_echantillon.csv / .pdf  : une ligne PAR ECHANTILLON
#      -> TP, FP, FN, TN, precision, recall, F1 ; TOUTES les lignes, paginees.
#         Une ligne TOTAL (en tete) somme TP/FP/FN/TN et recalcule prec/recall/F1.
#   3. confusion.xlsx : classeur Excel regroupant les deux tableaux (un onglet chacun)
#   4. venn_kmer_vizome.png : Venn des detections (kmer=TP+FP, vizome=TP+FN, inter=TP)
#   5. venn_fusions_uniques.png : Venn des NOMS de fusion uniques (1 fusion = 1 fois)
#
# Definitions (granularite = ligne de fichier, comme compare_fusions.py) :
#   Par fusion      : TP/FP/FN = nb de lignes de la fusion dans chaque fichier ;
#                     TN = nb total d'echantillons - nb d'echantillons de la fusion.
#   Par echantillon : TP/FP/FN = nb de lignes de l'echantillon dans chaque fichier ;
#                     TN = nb total de fusions - nb de fusions de l'echantillon.
#
# Usage :
#   Rscript visualize_fusions.R <dossier_resultats> [dossier_sortie] [lignes_par_page]
#   Rscript visualize_fusions.R resultats figures 40
# ---------------------------------------------------------------------------

## ---- dependances (installees automatiquement si absentes) -----------------
need <- c("VennDiagram", "gridExtra", "grid", "openxlsx")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressMessages({ library(VennDiagram); library(gridExtra); library(grid); library(openxlsx) })

## ---- arguments ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript visualize_fusions.R <dossier_resultats> [sortie] [lignes_par_page]")
}
res_dir       <- args[1]
out_dir       <- ifelse(length(args) >= 2, args[2], "figures")
rows_per_page <- ifelse(length(args) >= 3, as.integer(args[3]), 40L)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- lecture des fichiers -------------------------------------------------
read_tsv <- function(name) {
  path <- file.path(res_dir, name)
  if (!file.exists(path)) stop(paste("fichier introuvable :", path))
  df <- read.delim(path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, colClasses = "character")
  need_cols <- c("SampleID", "left_gene", "left_chr", "right_gene", "right_chr")
  miss <- setdiff(need_cols, colnames(df))
  if (length(miss) > 0)
    stop(paste0("[", name, "] colonnes manquantes : ", paste(miss, collapse = ", ")))
  # noms de genes en majuscules pour fusionner les variantes de casse
  df$fusion <- paste(toupper(df$left_gene), toupper(df$right_gene), sep = "_")
  df
}

tp <- read_tsv("vrais_positifs.tsv")
fp <- read_tsv("faux_positifs.tsv")
fn <- read_tsv("faux_negatifs.tsv")
n_tp <- nrow(tp); n_fp <- nrow(fp); n_fn <- nrow(fn)

n_total_samples <- length(unique(c(tp$SampleID, fp$SampleID, fn$SampleID)))
n_total_fusions <- length(unique(c(tp$fusion,   fp$fusion,   fn$fusion)))

## ---- construction d'un tableau de confusion -------------------------------
# group_col   : colonne de regroupement (lignes du tableau)
# other_col   : autre dimension (sert au calcul du TN)
# universe    : taille de l'univers de l'autre dimension (pour TN)
build_confusion <- function(group_col, other_col, universe) {
  ctp <- table(tp[[group_col]]); cfp <- table(fp[[group_col]]); cfn <- table(fn[[group_col]])
  g_all <- c(tp[[group_col]], fp[[group_col]], fn[[group_col]])
  o_all <- c(tp[[other_col]], fp[[other_col]], fn[[other_col]])
  distinct_other <- tapply(o_all, g_all, function(x) length(unique(x)))

  keys <- sort(unique(c(names(ctp), names(cfp), names(cfn))))
  g <- function(tab, k) ifelse(k %in% names(tab), as.integer(tab[k]), 0L)

  TP <- vapply(keys, function(k) g(ctp, k), integer(1))
  FP <- vapply(keys, function(k) g(cfp, k), integer(1))
  FN <- vapply(keys, function(k) g(cfn, k), integer(1))
  TN <- universe - as.integer(distinct_other[keys])
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall    <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  f1        <- ifelse((precision + recall) > 0,
                      2 * precision * recall / (precision + recall), 0)

  df <- data.frame(key = keys, TP = TP, FP = FP, FN = FN, TN = TN,
                   precision = round(precision, 4), recall = round(recall, 4),
                   F1 = round(f1, 4),
                   stringsAsFactors = FALSE)
  names(df)[1] <- group_col
  df <- df[order(-(df$TP + df$FP + df$FN), df[[group_col]]), ]
  rownames(df) <- NULL

  # ligne TOTAL : somme des TP/FP/FN/TN puis precision/recall/F1 recalcules
  sTP <- sum(df$TP); sFP <- sum(df$FP); sFN <- sum(df$FN); sTN <- sum(df$TN)
  tprec <- ifelse((sTP + sFP) > 0, sTP / (sTP + sFP), 0)
  trec  <- ifelse((sTP + sFN) > 0, sTP / (sTP + sFN), 0)
  tf1   <- ifelse((tprec + trec) > 0, 2 * tprec * trec / (tprec + trec), 0)
  total <- data.frame(key = "TOTAL", TP = sTP, FP = sFP, FN = sFN, TN = sTN,
                      precision = round(tprec, 4), recall = round(trec, 4),
                      F1 = round(tf1, 4), stringsAsFactors = FALSE)
  names(total)[1] <- group_col
  df <- rbind(total, df)            # TOTAL en premiere ligne
  rownames(df) <- NULL
  df
}

## ---- ecriture d'un tableau en PDF multi-pages -----------------------------
write_table_pdf <- function(df, pdf_path, titre) {
  n   <- nrow(df)
  npg <- max(1L, ceiling(n / rows_per_page))
  tt  <- ttheme_default(base_size = 9,
            colhead = list(fg_params = list(col = "white"),
                           bg_params = list(fill = "#40466e")))
  pdf(pdf_path, width = 8.27, height = 11.69)  # A4 portrait
  for (pg in seq_len(npg)) {
    i0 <- (pg - 1L) * rows_per_page + 1L
    i1 <- min(pg * rows_per_page, n)
    grid.newpage()
    grid.draw(tableGrob(df[i0:i1, ], rows = NULL, theme = tt))
    grid.draw(textGrob(sprintf("%s - page %d/%d (lignes %d-%d sur %d)",
                               titre, pg, npg, i0, i1, n),
                       gp = gpar(fontsize = 12, fontface = "bold"), y = 0.985))
  }
  dev.off()
  npg
}

## ---- 1. tableau PAR FUSION ------------------------------------------------
conf_fus <- build_confusion("fusion", "SampleID", n_total_samples)
write.csv(conf_fus, file.path(out_dir, "confusion_par_fusion.csv"), row.names = FALSE)
npg_fus <- write_table_pdf(conf_fus, file.path(out_dir, "tableau_confusion_par_fusion.pdf"),
                           "Confusion par fusion")

## ---- 2. tableau PAR ECHANTILLON -------------------------------------------
conf_smp <- build_confusion("SampleID", "fusion", n_total_fusions)
write.csv(conf_smp, file.path(out_dir, "confusion_par_echantillon.csv"), row.names = FALSE)
npg_smp <- write_table_pdf(conf_smp, file.path(out_dir, "tableau_confusion_par_echantillon.pdf"),
                           "Confusion par echantillon")

## ---- classeur Excel regroupant les deux tableaux --------------------------
write.xlsx(list(par_fusion = conf_fus, par_echantillon = conf_smp),
           file = file.path(out_dir, "confusion.xlsx"),
           headerStyle = createStyle(fgFill = "#40466E", fontColour = "#FFFFFF",
                                     textDecoration = "bold"))

## ---- 3. diagramme de Venn (detections : TP/FP/FN) -------------------------
venn_png <- file.path(out_dir, "venn_kmer_vizome.png")
png(venn_png, width = 1400, height = 1200, res = 200)
grid.newpage()
draw.pairwise.venn(
  area1 = n_tp + n_fp, area2 = n_tp + n_fn, cross.area = n_tp,
  category = c("Fusions kmer", "Fusions vizome"),
  fill = c("#66c2a5", "#fc8d62"), alpha = c(0.7, 0.7),
  lty = "blank", cex = 1.4, cat.cex = 1.2, cat.pos = c(-30, 30), ind = TRUE)
grid.draw(textGrob("Fusions kmer vs vizome (detections : TP/FP/FN)",
                   y = 0.95, gp = gpar(fontsize = 12, fontface = "bold")))
dev.off()

## ---- 4. Venn des noms de fusion uniques (geneA_geneB compte 1 fois) --------
cf_nofotal <- conf_fus[conf_fus$fusion != "TOTAL", ]
in_kmer   <- (cf_nofotal$TP + cf_nofotal$FP) > 0
in_vizome <- (cf_nofotal$TP + cf_nofotal$FN) > 0
nu_both   <- sum(in_kmer & in_vizome)
nu_kmer   <- sum(in_kmer & !in_vizome)
nu_vizome <- sum(in_vizome & !in_kmer)
venn_noms <- file.path(out_dir, "venn_fusions_uniques.png")
png(venn_noms, width = 1400, height = 1200, res = 200)
grid.newpage()
draw.pairwise.venn(
  area1 = nu_kmer + nu_both, area2 = nu_vizome + nu_both, cross.area = nu_both,
  category = c("Fusions kmer", "Fusions vizome"),
  fill = c("#66c2a5", "#fc8d62"), alpha = c(0.7, 0.7),
  lty = "blank", cex = 1.4, cat.cex = 1.2, cat.pos = c(-30, 30), ind = TRUE)
grid.draw(textGrob("Noms de fusion uniques (1 fusion = 1 fois)",
                   y = 0.95, gp = gpar(fontsize = 12, fontface = "bold")))
dev.off()

## ---- resume console -------------------------------------------------------
cat(sprintf("TP=%d  FP=%d  FN=%d\n", n_tp, n_fp, n_fn))
cat(sprintf("Precision = %.4f | Recall = %.4f\n",
            n_tp / (n_tp + n_fp), n_tp / (n_tp + n_fn)))
cat(sprintf("Echantillons : %d | Fusions : %d\n", n_total_samples, n_total_fusions))
cat(sprintf("Noms de fusion uniques : kmer seul=%d, communs=%d, vizome seul=%d\n",
            nu_kmer, nu_both, nu_vizome))
cat("Fichiers ecrits :\n")
cat(sprintf("  - %s/confusion_par_fusion.csv\n", out_dir))
cat(sprintf("  - %s/confusion_par_echantillon.csv\n", out_dir))
cat(sprintf("  - %s/confusion.xlsx (onglets : par_fusion, par_echantillon)\n", out_dir))
cat(sprintf("  - %s/tableau_confusion_par_fusion.pdf (%d pages)\n", out_dir, npg_fus))
cat(sprintf("  - %s/tableau_confusion_par_echantillon.pdf (%d pages)\n", out_dir, npg_smp))
cat(sprintf("  - %s\n", venn_png))
cat(sprintf("  - %s\n", venn_noms))
