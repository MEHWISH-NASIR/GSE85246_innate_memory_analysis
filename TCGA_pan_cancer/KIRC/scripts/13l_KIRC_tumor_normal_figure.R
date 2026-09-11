# ============================================================
# 13l — TCGA-KIRC tumor vs normal figure
# JAK3 + EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

required <- c("ggplot2", "patchwork")

missing <- required[
  !sapply(
    required,
    requireNamespace,
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  install.packages(
    missing,
    repos = "https://cloud.r-project.org"
  )
}

library(ggplot2)
library(patchwork)

infile <-
  "results/13_TCGA_priority_genes/13g_KIRC_JAK3_EPHB2_expression.csv"

outfile <-
  "TCGA_pan_cancer/KIRC/figures/01_tumor_normal_JAK3_EPHB2.png"

x <- read.csv(
  infile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

x <- x[
  x$extraction_status == "OK" &
  x$sample_type %in%
    c(
      "Primary Tumor",
      "Solid Tissue Normal"
    ),
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
# Reusable plot function
# ------------------------------------------------------------

make_plot <- function(
  gene,
  fdr_text,
  prevalence_text,
  interpretation
) {

  z <- x[x$Gene == gene, ]

  tumor_tpm <-
    z$TPM[z$Group == "Tumor"]

  normal_tpm <-
    z$TPM[z$Group == "Normal"]

  tumor_med <-
    median(tumor_tpm, na.rm = TRUE)

  normal_med <-
    median(normal_tpm, na.rm = TRUE)

  n_tumor <-
    sum(z$Group == "Tumor")

  n_normal <-
    sum(z$Group == "Normal")

  subtitle <- sprintf(
    "Median TPM: Tumor %.2f | Normal %.2f",
    tumor_med,
    normal_med
  )

  caption <- paste0(
    "Tumor n = ", n_tumor,
    " | Normal n = ", n_normal,
    "\n",
    "Tumor vs normal FDR: ", fdr_text,
    "\n",
    prevalence_text,
    "\n",
    interpretation
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
      alpha = 0.85,
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
      y = expression(
        log[2] * "(TPM + 1)"
      ),
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

      axis.title.y = element_text(
        size = 12
      ),

      plot.caption = element_text(
        size = 9.5,
        hjust = 0,
        lineheight = 1.15,
        margin = margin(t = 12)
      ),

      plot.margin = margin(
        10, 15, 10, 15
      )
    )
}

# ------------------------------------------------------------
# JAK3
# ------------------------------------------------------------

p1 <- make_plot(
  gene = "JAK3",

  fdr_text = "< 1e-15",

  prevalence_text =
    paste0(
      "Directional prevalence: ",
      "401/533 tumors (75.2%) above ",
      "the normal 95th percentile"
    ),

  interpretation =
    "Strong tumor-associated upregulation."
)

# ------------------------------------------------------------
# EPHB2
# ------------------------------------------------------------

p2 <- make_plot(
  gene = "EPHB2",

  fdr_text = "3.13e-15",

  prevalence_text =
    paste0(
      "Directional prevalence: ",
      "187/533 tumors (35.1%) below ",
      "the normal 5th percentile"
    ),

  interpretation =
    "Significant tumor-associated downregulation."
)

# ------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------

final_plot <-
  p1 + p2 +
  plot_annotation(
    title =
      "TCGA-KIRC: JAK3 and EPHB2 Tumor–Normal Expression",

    subtitle =
      paste0(
        "Kidney renal clear cell carcinoma | ",
        "TCGA STAR-count RNA-seq"
      ),

    caption =
      paste0(
        "Expression distributions are based on TPM. ",
        "Prevalence estimates are patient-level and ",
        "defined relative to the normal expression distribution."
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
cat("KIRC TUMOR-NORMAL FIGURE COMPLETE\n")
cat("========================================\n")
cat("Saved:\n ", outfile, "\n")
