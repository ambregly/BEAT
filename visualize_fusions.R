#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Visualisations de la comparaison fusions kmer vs BEAT AML, construites
# DIRECTEMENT a partir des fichiers produits par compare_fusions.py :
#     <resultats>/vrais_positifs.tsv   (TP)
#     <resultats>/faux_positifs.tsv    (FP)
#     <resultats>/faux_negatifs.tsv    (FN)
#
# Les chiffres sont donc strictement coherents avec ces fichiers.
#
# Produit :
#   1. confusion_par_fusion.csv : TP, FP, FN, TN, precision, recall par fusion
#      (fusion = left_gene_right_gene).
#   2. tableau_confusion.pdf    : TOUTES les fusions, paginees (N par page).
#   3. venn_kmer_vizome.png     : diagramme de Venn kmer vs vizome
#      (kmer = TP+FP, vizome = TP+FN, intersection = TP).
#
# Granularite : ligne du fichier (= meme unite de comptage que compare_fusions.py).
#   TN par fusion = nb total d'echantillons - nb d'echantillons ou la fusion
#   apparait (dans TP, FP ou FN).
#
# Usage :
#   Rscript visualize_fusions.R <dossier_resultats> [dossier_sortie] [lignes_par_page]
#   Rscript visualize_fusions.R resultats figures 40
# ---------------------------------------------------------------------------

## ---- dependances (installees automatiquement si absentes) -----------------
need <- c("VennDiagram", "gridExtra", "grid")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressMessages({
  library(VennDiagram)
  library(gridExtra)
  library(grid)
})

## ---- arguments ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript visualize_fusions.R <dossier_resultats> [sortie] [lignes_par_page]")
}
res_dir        <- args[1]
out_dir        <- ifelse(length(args) >= 2, args[2], "figures")
rows_per_page  <- ifelse(length(args) >= 3, as.integer(args[3]), 40L)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- lecture des fichiers -------------------------------------------------
read_tsv <- function(name) {
  path <- file.path(res_dir, name)
  if (!file.exists(path)) stop(paste("fichier introuvable :", path))
  df <- read.delim(path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, colClasses = "character")
  # colonnes communes aux 3 fichiers, quel que soit le mode (row/index/genepair)
  need_cols <- c("SampleID", "left_gene", "left_chr", "right_gene", "right_chr")
  miss <- setdiff(need_cols, colnames(df))
  if (length(miss) > 0)
    stop(paste0("[", name, "] colonnes manquantes : ", paste(miss, collapse = ", ")))
  df$fusion <- paste(df$left_gene, df$right_gene, sep = "_")
  df
}

tp <- read_tsv("vrais_positifs.tsv")
fp <- read_tsv("faux_positifs.tsv")
fn <- read_tsv("faux_negatifs.tsv")

n_tp <- nrow(tp); n_fp <- nrow(fp); n_fn <- nrow(fn)

## ---- tableau de confusion par fusion --------------------------------------
all_samples <- unique(c(tp$SampleID, fp$SampleID, fn$SampleID))
n_total <- length(all_samples)

count_by_fusion <- function(df) table(df$fusion)
ctp <- count_by_fusion(tp); cfp <- count_by_fusion(fp); cfn <- count_by_fusion(fn)
all_fusions <- sort(unique(c(names(ctp), names(cfp), names(cfn))))

get <- function(tab, key) ifelse(key %in% names(tab), as.integer(tab[key]), 0L)

# echantillons distincts par fusion (toutes categories confondues) -> pour TN
samp_by_fusion <- tapply(
  c(tp$SampleID, fp$SampleID, fn$SampleID),
  c(tp$fusion,   fp$fusion,   fn$fusion),
  function(s) length(unique(s))
)

rows <- lapply(all_fusions, function(f) {
  TP <- get(ctp, f); FP <- get(cfp, f); FN <- get(cfn, f)
  TN <- n_total - as.integer(samp_by_fusion[f])
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  data.frame(fusion = f, TP = TP, FP = FP, FN = FN, TN = TN,
             precision = round(precision, 4), recall = round(recall, 4),
             stringsAsFactors = FALSE)
})
conf <- do.call(rbind, rows)
conf <- conf[order(-(conf$TP + conf$FP + conf$FN), conf$fusion), ]
rownames(conf) <- NULL

csv_path <- file.path(out_dir, "confusion_par_fusion.csv")
write.csv(conf, csv_path, row.names = FALSE)

## ---- tableau : TOUTES les fusions en PDF multi-pages -----------------------
n_fus    <- nrow(conf)
n_pages  <- max(1L, ceiling(n_fus / rows_per_page))
table_pdf <- file.path(out_dir, "tableau_confusion.pdf")
tt <- ttheme_default(base_size = 9,
                     colhead = list(fg_params = list(col = "white"),
                                    bg_params = list(fill = "#40466e")))
pdf(table_pdf, width = 8.27, height = 11.69)  # format A4 portrait
for (pg in seq_len(n_pages)) {
  i0 <- (pg - 1L) * rows_per_page + 1L
  i1 <- min(pg * rows_per_page, n_fus)
  sub <- conf[i0:i1, ]
  grid.newpage()
  title <- textGrob(sprintf("Confusion par fusion - page %d/%d (fusions %d-%d sur %d)",
                            pg, n_pages, i0, i1, n_fus),
                    gp = gpar(fontsize = 12, fontface = "bold"), y = 0.985)
  tg <- tableGrob(sub, rows = NULL, theme = tt)
  grid.draw(tg)
  grid.draw(title)
}
dev.off()

## ---- diagramme de Venn ----------------------------------------------------
# kmer = TP + FP ; vizome = TP + FN ; intersection = TP  (= les fichiers)
venn_png <- file.path(out_dir, "venn_kmer_vizome.png")
png(venn_png, width = 1400, height = 1200, res = 200)
grid.newpage()
vp <- draw.pairwise.venn(
  area1 = n_tp + n_fp,          # fusions kmer
  area2 = n_tp + n_fn,          # fusions vizome
  cross.area = n_tp,            # communes
  category = c("Fusions kmer", "Fusions vizome"),
  fill = c("#66c2a5", "#fc8d62"), alpha = c(0.7, 0.7),
  lty = "blank", cex = 1.4, cat.cex = 1.2,
  cat.pos = c(-30, 30), ind = TRUE)
grid.draw(textGrob("Fusions kmer vs vizome",
                   y = 0.95, gp = gpar(fontsize = 13, fontface = "bold")))
dev.off()

## ---- resume console -------------------------------------------------------
cat(sprintf("TP=%d  FP=%d  FN=%d\n", n_tp, n_fp, n_fn))
cat(sprintf("Precision = %.4f | Recall = %.4f\n",
            n_tp / (n_tp + n_fp), n_tp / (n_tp + n_fn)))
cat(sprintf("Echantillons (union des 3 fichiers) : %d\n", n_total))
cat(sprintf("Fusions distinctes : %d (tableau sur %d page(s))\n", nrow(conf), n_pages))
cat("Fichiers ecrits :\n")
cat(sprintf("  - %s\n  - %s\n  - %s\n", csv_path, table_pdf, venn_png))
