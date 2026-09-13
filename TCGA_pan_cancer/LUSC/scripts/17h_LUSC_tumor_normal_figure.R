# ============================================================
# 17h — TCGA-LUSC tumor-normal expression figure
# JAK3 + EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(patchwork)

expr_file <-
  "results/17_TCGA_LUSC/17a_LUSC_JAK3_EPHB2_expression.csv"

stats_file <-
  "results/17_TCGA_LUSC/17b_LUSC_tumor_normal_summary.csv"

prev_file <-
  "results/17_TCGA_LUSC/17b_LUSC_prevalence_summary.csv"

outfile <-
  "TCGA_pan_cancer/LUSC/figures/01_tumor_normal_JAK3_EPHB2.png"

x <- read.csv(expr_file, stringsAsFactors = FALSE)
stats <- read.csv(stats_file, stringsAsFactors = FALSE)
prev <- read.csv(prev_file, stringsAsFactors = FALSE)

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

make_plot <- function(gene) {

  z <- x[x$gene_name == gene, ]
  s <- stats[stats$Gene == gene, ]
  p <- prev[prev$Gene == gene, ]

  fdr <- s$Wilcoxon_FDR

  fdr_text <- if (
    !is.na(fdr) && fdr < 0.001
  ) {
    format(
      fdr,
      scientific = TRUE,
      digits = 3
    )
  } else {
    sprintf("%.3f", fdr)
  }

  interpretation <- if (
    !is.na(fdr) &&
    fdr < 0.05 &&
    s$Tumor_median_TPM > s$Normal_median_TPM
  ) {

    "Significant tumor-associated upregulation"

  } else if (
    !is.na(fdr) &&
    fdr < 0.05 &&
    s$Tumor_median_TPM < s$Normal_median_TPM
  ) {

    "Significant tumor-associated downregulation"

  } else {

    "No significant tumor-normal difference"
  }

  prevalence_note <- if (
    fdr < 0.05
  ) {

    sprintf(
      "High-expression prevalence: %d/%d tumors (%.1f%%)",
      p$Aberrant_tumor_n,
      p$Tumor_patient_n,
      p$Aberrant_tumor_percent
    )

  } else {

    sprintf(
      "Above normal P95: %d/%d tumors (%.1f%%; descriptive)",
      p$Aberrant_tumor_n,
      p$Tumor_patient_n,
      p$Aberrant_tumor_percent
    )
  }

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
      alpha = 0.65
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.9
    ) +
    scale_fill_manual(
      values = c(
        "Normal" = "#0072B2",
        "Tumor" = "#D55E00"
      )
    ) +
    labs(
      title = gene,
      subtitle = sprintf(
        "Median TPM: Normal %.2f | Tumor %.2f",
        s$Normal_median_TPM,
        s$Tumor_median_TPM
      ),
      x = NULL,
      y = expression(log[2] * "(TPM + 1)"),
      caption = paste0(
        "Tumor n = ", s$Tumor_n,
        " | Normal n = ", s$Normal_n,
        "\nTumor-normal FDR = ", fdr_text,
        "\n", prevalence_note,
        "\n", interpretation
      )
    ) +
    theme_classic(base_size = 13) +
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
        face = "bold"
      ),
      plot.caption = element_text(
        size = 9.5,
        hjust = 0,
        lineheight = 1.15,
        margin = margin(t = 10)
      )
    )
}

p1 <- make_plot("JAK3")
p2 <- make_plot("EPHB2")

final_plot <-
  p1 + p2 +
  plot_annotation(
    title =
      "TCGA-LUSC: JAK3 and EPHB2 Tumor–Normal Expression",
    subtitle =
      "Lung squamous cell carcinoma | TCGA STAR-count RNA-seq",
    caption =
      paste0(
        "Expression unit: TPM. ",
        "Directional prevalence uses the normal-derived ",
        "5th/95th percentile threshold."
      ),
    theme = theme(
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )
  )

ggsave(
  outfile,
  final_plot,
  width = 12,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

cat("\nLUSC figure saved:\n")
cat(outfile, "\n")
