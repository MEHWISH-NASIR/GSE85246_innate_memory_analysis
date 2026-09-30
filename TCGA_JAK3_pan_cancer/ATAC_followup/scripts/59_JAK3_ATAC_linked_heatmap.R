# ============================================================
# 59 — JAK3-linked ATAC regulatory-element heatmap
#
# Combines:
# - RNA-based JAK3 priority score
# - median accessibility of 3 official Data S7 JAK3 links
#
# Heatmap color is scaled separately within each feature.
# Numbers printed in cells are the original values.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(
  library(ggplot2)
)

atac_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results/JAK3_linked_ATAC_summary_9cancers.csv"

ranking_file <-
  "TCGA_JAK3_pan_cancer/results/release/JAK3_tumor_priority_ranking.csv"

outdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results"

figdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/figures"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

atac <- read.csv(
  atac_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rank <- read.csv(
  ranking_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

priority <- rank[
  rank$Priority_rank <= 9,
  c(
    "Cancer",
    "Priority_rank",
    "Priority_score"
  )
]

# ------------------------------------------------------------
# Create wide ATAC summary
# ------------------------------------------------------------

atac_small <- atac[
  ,
  c(
    "Cancer",
    "Enhancer_ID",
    "Median_accessibility"
  )
]

atac_wide <- reshape(
  atac_small,
  idvar = "Cancer",
  timevar = "Enhancer_ID",
  direction = "wide"
)

names(atac_wide) <- sub(
  "^Median_accessibility\\.",
  "",
  names(atac_wide)
)

dat <- merge(
  priority,
  atac_wide,
  by = "Cancer",
  all.x = TRUE,
  sort = FALSE
)

dat <- dat[
  order(dat$Priority_rank),
]

# ------------------------------------------------------------
# Feature metadata
# ------------------------------------------------------------

feature_info <- unique(
  atac[
    ,
    c(
      "Enhancer_ID",
      "Enhancer_Distance",
      "Enhancer_Correlation",
      "Enhancer_FDR"
    )
  ]
)

write.csv(
  feature_info,
  file.path(
    outdir,
    "JAK3_ATAC_linked_element_metadata.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Long format
# ------------------------------------------------------------

features <- c(
  "Priority_score",
  "JAK3_m1",
  "JAK3_p1",
  "JAK3_p2"
)

labels <- c(
  "RNA priority",
  "JAK3_m1\n-2.6 kb",
  "JAK3_p1\n+241.6 kb",
  "JAK3_p2\n+0.36 kb"
)

long <- do.call(
  rbind,
  lapply(
    seq_along(features),
    function(i) {

      data.frame(
        Cancer = dat$Cancer,
        Priority_rank = dat$Priority_rank,
        Feature = labels[i],
        Raw_value = dat[[features[i]]],
        stringsAsFactors = FALSE
      )
    }
  )
)

# ------------------------------------------------------------
# Scale each feature independently from 0 to 1
#
# This is only for heatmap color.
# Cell labels retain original values.
# ------------------------------------------------------------

long$Scaled_value <- NA_real_

for (f in unique(long$Feature)) {

  idx <- long$Feature == f

  x <- long$Raw_value[idx]

  rng <- range(
    x,
    na.rm = TRUE
  )

  if (
    all(is.finite(rng)) &&
    diff(rng) > 0
  ) {

    long$Scaled_value[idx] <-
      (x - rng[1]) /
      (rng[2] - rng[1])

  } else {

    long$Scaled_value[idx] <- 0.5
  }
}

# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------

long$Cell_label <- ifelse(
  long$Feature == "RNA priority",
  sprintf("%.1f", long$Raw_value),
  sprintf("%.2f", long$Raw_value)
)

cancer_levels <- dat$Cancer

long$Cancer <- factor(
  long$Cancer,
  levels = rev(cancer_levels)
)

long$Feature <- factor(
  long$Feature,
  levels = labels
)

# ------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------

p <- ggplot(
  long,
  aes(
    x = Feature,
    y = Cancer,
    fill = Scaled_value
  )
) +
  geom_tile(
    linewidth = 0.65
  ) +
  geom_text(
    aes(label = Cell_label),
    size = 3.8
  ) +
  scale_fill_gradient(
    limits = c(0, 1),
    name = "Relative\nsignal"
  ) +
  labs(
    title =
      "JAK3 RNA prioritization and linked chromatin accessibility",
    subtitle =
      "Official TCGA Data S7 JAK3-linked ATAC elements across prioritized tumor types",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Numbers are original RNA priority scores or median normalized ATAC accessibility. ",
      "Heatmap colors are scaled independently within each column. ",
      "m1/p1/p2 denote official pan-cancer JAK3-linked regulatory elements."
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

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

ggsave(
  file.path(
    figdir,
    "Figure_9_JAK3_ATAC_linked_elements_heatmap.png"
  ),
  p,
  width = 9,
  height = 7,
  dpi = 320
)

ggsave(
  file.path(
    figdir,
    "Figure_9_JAK3_ATAC_linked_elements_heatmap.pdf"
  ),
  p,
  width = 9,
  height = 7
)

# ------------------------------------------------------------
# Save combined table
# ------------------------------------------------------------

write.csv(
  dat,
  file.path(
    outdir,
    "JAK3_RNA_priority_ATAC_linked_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

cat("\n========================================\n")
cat("JAK3 RNA + ATAC INTEGRATED SUMMARY\n")
cat("========================================\n\n")

print(
  dat,
  row.names = FALSE
)

cat("\nOfficial linked-element metadata:\n\n")

print(
  feature_info,
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    figdir,
    "Figure_9_JAK3_ATAC_linked_elements_heatmap.png"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "JAK3_RNA_priority_ATAC_linked_summary.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("59 COMPLETE\n")
cat("========================================\n")
