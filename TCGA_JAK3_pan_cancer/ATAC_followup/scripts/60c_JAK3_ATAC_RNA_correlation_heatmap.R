# ============================================================
# 60c — JAK3 matched ATAC-RNA correlation heatmap
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(
  library(ggplot2)
)

infile <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results/JAK3_ATAC_RNA_within_cancer_correlations.csv"

figdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/figures"

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

x <- read.csv(
  infile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cancer_order <- c(
  "CHOL","KIRC","KIRP","THCA","HNSC",
  "STAD","LUAD","COAD","LIHC"
)

enhancer_order <- c(
  "JAK3_m1",
  "JAK3_p1",
  "JAK3_p2"
)

# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------

x$Cell_label <- ""

# Primary tests
idx_primary <- x$Analysis_status == "Primary"

x$Cell_label[idx_primary] <- sprintf(
  "ρ %.2f\nN=%d",
  x$Spearman_rho[idx_primary],
  x$Matched_N[idx_primary]
)

# Add * only for GLOBAL FDR significant tests
sig <- idx_primary &
  !is.na(x$Global_FDR) &
  x$Global_FDR < 0.05

x$Cell_label[sig] <- paste0(
  x$Cell_label[sig],
  "\n*"
)

# Exploratory low-N
idx_exp <- x$Analysis_status == "Exploratory_low_N"

x$Cell_label[idx_exp] <- sprintf(
  "ρ %.2f\nN=%d\n†",
  x$Spearman_rho[idx_exp],
  x$Matched_N[idx_exp]
)

# Not tested
idx_not <- x$Analysis_status == "Not_tested_N_lt_5"

x$Cell_label[idx_not] <- sprintf(
  "N=%d\nNA",
  x$Matched_N[idx_not]
)

# Keep fill blank for untested cases
x$Plot_rho <- x$Spearman_rho
x$Plot_rho[idx_not] <- NA_real_

x$Cancer <- factor(
  x$Cancer,
  levels = rev(cancer_order)
)

x$Enhancer_ID <- factor(
  x$Enhancer_ID,
  levels = enhancer_order,
  labels = c(
    "JAK3_m1\n-2.6 kb",
    "JAK3_p1\n+241.6 kb",
    "JAK3_p2\n+0.36 kb"
  )
)

# ------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------

p <- ggplot(
  x,
  aes(
    x = Enhancer_ID,
    y = Cancer,
    fill = Plot_rho
  )
) +
  geom_tile(
    linewidth = 0.7
  ) +
  geom_text(
    aes(label = Cell_label),
    size = 3.5,
    lineheight = 0.9
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90",
    name = "Spearman\nrho"
  ) +
  labs(
    title =
      "Patient-matched JAK3 chromatin accessibility and RNA expression",
    subtitle =
      "Within-cancer correlations for official TCGA Data S7 JAK3-linked elements",
    x = NULL,
    y = NULL,
    caption = paste0(
      "* global BH-FDR < 0.05 across primary tests (N >= 10). ",
      "† exploratory low-N analysis (5-9 matched patients). ",
      "NA = fewer than 5 matched patients."
    )
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.y = element_text(
      size = 10
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 8.5
    )
  )

ggsave(
  file.path(
    figdir,
    "Figure_10_JAK3_ATAC_RNA_correlation_heatmap.png"
  ),
  p,
  width = 8,
  height = 7,
  dpi = 320
)

ggsave(
  file.path(
    figdir,
    "Figure_10_JAK3_ATAC_RNA_correlation_heatmap.pdf"
  ),
  p,
  width = 8,
  height = 7
)

cat("\n========================================\n")
cat("FIGURE 10 COMPLETE\n")
cat("========================================\n\n")

cat("Global-FDR significant primary associations:\n\n")

print(
  x[
    sig,
    c(
      "Cancer",
      "Enhancer_ID",
      "Matched_N",
      "Spearman_rho",
      "Global_FDR"
    )
  ],
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    figdir,
    "Figure_10_JAK3_ATAC_RNA_correlation_heatmap.png"
  ),
  "\n"
)
