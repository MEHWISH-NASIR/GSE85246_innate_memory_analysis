
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setDTthreads(1)

outdir <- "macrophage_transcriptomics_validation/E-MTAB-9214/results"
figdir <- "macrophage_transcriptomics_validation/E-MTAB-9214/figures"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

x <- data.table(
  sample = c(
    "36_D7HDM","52_D7HDM","123_D7HDM",
    "36_ctr","52_ctr","123_ctr"
  ),
  group = c(
    "HDM","HDM","HDM",
    "Control","Control","Control"
  ),
  JAK3 = c(
    107.9161109,
    85.10733582,
    104.023369,
    33.02359916,
    24.43957257,
    28.04502478
  )
)

x[, group := factor(group, levels = c("Control","HDM"))]

fwrite(
  x,
  file.path(outdir, "31a_EMTAB9214_JAK3_normalized_values.csv")
)

summary <- x[
  ,
  .(
    n = .N,
    mean = mean(JAK3),
    median = median(JAK3),
    sd = sd(JAK3)
  ),
  by = group
]

fwrite(
  summary,
  file.path(outdir, "31a_EMTAB9214_JAK3_group_summary.csv")
)

formal <- data.table(
  dataset = "E-MTAB-9214",
  gene = "JAK3",
  contrast_original = "Control_vs_D7HDM",
  original_log2FC = -1.700294815,
  standardized_contrast = "D7HDM_vs_Control",
  standardized_log2FC = 1.700294815,
  P_value = 9.96023e-10,
  FDR = 6.4398e-07
)

fwrite(
  formal,
  file.path(outdir, "31a_EMTAB9214_JAK3_formal_result.csv")
)

cat("\n===== JAK3 FORMAL RESULT =====\n")
print(formal)

cat("\n===== GROUP SUMMARY =====\n")
print(summary)

p <- ggplot(
  x,
  aes(
    x = group,
    y = JAK3
  )
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.05,
    size = 3
  ) +
  labs(
    title = "E-MTAB-9214: persistent JAK3 expression after HDM training",
    subtitle = "Resting human macrophages after the training/rest period",
    x = NULL,
    y = "Deposited normalized JAK3 expression"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(figdir, "31a_EMTAB9214_JAK3_HDM_vs_control.png"),
  p,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figdir, "31a_EMTAB9214_JAK3_HDM_vs_control.pdf"),
  p,
  width = 6,
  height = 5
)

cat("\n31a COMPLETE\n")
