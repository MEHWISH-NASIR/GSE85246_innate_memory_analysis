
rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
})

outdir <- "figures/18_JAK3_EPHB2_reconciliation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

h27 <- read.csv(
  "results/18_JAK3_EPHB2_reconciliation/H3K27ac/18b_JAK3_EPHB2_H3K27ac_replicate_values.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

h4 <- read.csv(
  "results/18_JAK3_EPHB2_reconciliation/H3K4me1/18b_JAK3_EPHB2_H3K4me1_replicate_values.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

h27$Mark <- "H3K27ac"
h4$Mark <- "H3K4me1"

x <- rbind(h27, h4)

# Primary scaling only
x <- x[
  x$scaling_method == "PeerMedian_primary",
  ,
  drop = FALSE
]

x$Day <- ifelse(
  grepl("_d1_", x$sample),
  "Day 1",
  "Day 6"
)

x$Condition <- ifelse(
  grepl("^LPS", x$sample),
  "LPS",
  "RPMI"
)

x$Replicate <- ifelse(
  grepl("rep1$", x$sample),
  "rep1",
  "rep2"
)

x$Day <- factor(
  x$Day,
  levels = c("Day 1", "Day 6")
)

x$Condition <- factor(
  x$Condition,
  levels = c("RPMI", "LPS")
)

p <- ggplot(
  x,
  aes(
    x = Day,
    y = log2_signal,
    group = Condition,
    shape = Condition
  )
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  stat_summary(
    aes(group = Condition),
    fun = mean,
    geom = "line",
    position = position_dodge(width = 0.35),
    linewidth = 0.7
  ) +
  facet_grid(
    SYMBOL ~ Mark,
    scales = "free_y"
  ) +
  labs(
    title = "Individual replicate chromatin signals for JAK3 and EPHB2",
    subtitle = "Primary peer-median scaling; points represent individual biological replicates",
    x = NULL,
    y = "log2(signal + 1)",
    shape = "Condition"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    outdir,
    "Figure_17_JAK3_EPHB2_individual_replicates.png"
  ),
  p,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    outdir,
    "Figure_17_JAK3_EPHB2_individual_replicates.pdf"
  ),
  p,
  width = 9,
  height = 6
)

cat("\n18e COMPLETE\n")

