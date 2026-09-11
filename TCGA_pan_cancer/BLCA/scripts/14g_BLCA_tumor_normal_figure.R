# ============================================================
# 14g — TCGA-BLCA tumor vs normal figure
# JAK3 + EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

required <- c("ggplot2", "patchwork")

missing <- required[
  !sapply(required, requireNamespace, quietly = TRUE)
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
  )
}

library(ggplot2)
library(patchwork)

expr_file <-
  "results/14_TCGA_BLCA/14a_BLCA_JAK3_EPHB2_expression.csv"

stats_file <-
  "results/14_TCGA_BLCA/14b_BLCA_tumor_normal_summary.csv"

prev_file <-
  "results/14_TCGA_BLCA/14b_BLCA_prevalence_summary.csv"

outfile <-
  "TCGA_pan_cancer/BLCA/figures/01_tumor_normal_JAK3_EPHB2.png"

x <- read.csv(
  expr_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stats <- read.csv(
  stats_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

prev <- read.csv(
  prev_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

x <- x[
  x$extraction_status == "OK" &
  x$sample_type %in%
    c("Primary Tumor", "Solid Tissue Normal"),
]

x$Group <- ifelse(
  x$sample_type == "Primary Tumor",
  "Tumor",
  "Normal"
)

x$Group <- factor(
  x$Group,
  levels = c("Normal", "Tumor")
)

x$log2TPM1 <- log2(x$TPM + 1)

# ------------------------------------------------------------
# Plot function
# ------------------------------------------------------------

make_plot <- function(gene) {

  z <- x[
    x$gene_name == gene,
  ]

  s <- stats[
    stats$Gene == gene,
  ]

  p <- prev[
    prev$Gene == gene,
  ]

  if (gene == "JAK3") {

    interpretation <-
      "No significant tumor-normal difference"

  } else {

    interpretation <-
      "Significant tumor-associated upregulation"
  }

  fdr_text <- if (
    s$Wilcoxon_FDR < 0.001
  ) {
    format(
      s$Wilcoxon_FDR,
      scientific = TRUE,
      digits = 3
    )
  } else {
    sprintf("%.3f", s$Wilcoxon_FDR)
  }

  prevalence_text <- sprintf(
    "Aberrant prevalence: %d/%d tumors (%.1f%%)",
    p$Aberrant_tumor_n,
    p$Tumor_patient_n,
    p$Aberrant_tumor_percent
  )

  subtitle <- sprintf(
    "Median TPM: Tumor %.2f | Normal %.2f",
    s$Tumor_median_TPM,
    s$Normal_median_TPM
  )

  caption <- paste0(
    "Tumor n = ", s$Tumor_n,
    " | Normal n = ", s$Normal_n,
    "\nTumor-normal FDR = ", fdr_text,
    "\n", prevalence_text,
    "\n", interpretation
  )

  ggplot(
    z,
    aes(
      x = Group,
      y = log2TPM1,
      fill = Group
    )
  ) +

    geom_violin(
      trim = FALSE,
      scale = "width",
      alpha = 0.65,
      linewidth = 0.4
    ) +

    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.9,
      linewidth = 0.45
    ) +

    scale_fill_manual(
      values = c(
        "Normal" = "#0072B2",
        "Tumor" = "#D55E00"
      )
    ) +

    labs(
      title = gene,
      subtitle = subtitle,
      x = NULL,
      y = expression(log[2] * "(TPM + 1)"),
      caption = caption
    ) +

    theme_classic(
      base_size = 13
    ) +

    theme(
      legend.position = "none",

      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      ),

      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5
      ),

      axis.text.x = element_text(
        size = 12,
        face = "bold"
      ),

      plot.caption = element_text(
        size = 9.5,
        hjust = 0,
        lineheight = 1.15,
        margin = margin(t = 12)
      )
    )
}

p1 <- make_plot("JAK3")
p2 <- make_plot("EPHB2")

final_plot <-
  p1 + p2 +
  plot_annotation(
    title =
      "TCGA-BLCA: JAK3 and EPHB2 Tumor–Normal Expression",

    subtitle =
      "Bladder urothelial carcinoma | TCGA STAR-count RNA-seq",

    caption =
      paste0(
        "Expression unit: TPM. ",
        "Tumor-normal analysis uses all available samples. ",
        "Directional prevalence is patient-level and based on ",
        "the normal expression distribution."
      ),

    theme = theme(
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      ),

      plot.subtitle = element_text(
        size = 12,
        hjust = 0.5
      ),

      plot.caption = element_text(
        size = 9,
        hjust = 0
      )
    )
  )

ggsave(
  filename = outfile,
  plot = final_plot,
  width = 12,
  height = 7.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("\n========================================\n")
cat("BLCA TUMOR-NORMAL FIGURE COMPLETE\n")
cat("========================================\n")
cat("Saved:\n ", outfile, "\n")
