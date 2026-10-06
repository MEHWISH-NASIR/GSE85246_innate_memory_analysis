
library(data.table)
library(ggplot2)

setDTthreads(1)

root <- "TCGA_multikinase_analysis"

indir <- file.path(
  root,
  "results/final_integration"
)

figdir <- file.path(
  root,
  "figures/final_integration"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("STEP 14: PUBLICATION-READY FIGURES\n")
cat("========================================\n")

# ============================================================
# 1. Load data
# ============================================================

compact <- fread(
  file.path(
    indir,
    "13_supervisor_final_multikinase_comparison.csv"
  )
)

domains <- fread(
  file.path(
    indir,
    "13_final_domain_support_matrix.csv"
  )
)

atac_sig <- fread(
  file.path(
    indir,
    "13_significant_ATAC_RNA_associations.csv"
  )
)

cat("Compact table rows:", nrow(compact), "\n")
cat("Domain table rows:", nrow(domains), "\n")
cat("Significant ATAC rows:", nrow(atac_sig), "\n")

# ============================================================
# 2. Shared formatting
# ============================================================

tier_order <- c(
  "REFERENCE",
  "VERY_STRONG",
  "STRONG",
  "SUPPORTED",
  "EXPLORATORY"
)

compact[, Final_evidence_tier :=
  factor(
    Final_evidence_tier,
    levels = tier_order
  )
]

compact <- compact[
  order(
    Final_evidence_tier,
    -Final_multiomic_index,
    gene
  )
]

gene_order <- compact$gene

publication_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey50"),
    strip.text = element_text(face = "bold")
  )

save_both <- function(plot_obj, filename, width, height) {
  ggsave(
    file.path(figdir, paste0(filename, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600
  )
  ggsave(
    file.path(figdir, paste0(filename, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600
  )
}

# ============================================================
# 3. Figure 1: Six-domain evidence heatmap
# ============================================================

domain_cols <- c(
  "Primary_RNA_support",
  "Early_stage_support",
  "Progression_support",
  "Survival_support",
  "Immune_adjusted_support",
  "ATAC_evidence_support"
)

domain_labels <- c(
  Primary_RNA_support = "Primary RNA",
  Early_stage_support = "Early-stage RNA",
  Progression_support = "Stage/grade progression",
  Survival_support = "Survival",
  Immune_adjusted_support = "Immune-adjusted progression",
  ATAC_evidence_support = "ATAC-RNA coupling"
)

dom_plot <- merge(
  domains[, c("gene", domain_cols), with = FALSE],
  compact[, .(gene, Final_multiomic_index, Final_evidence_tier)],
  by = "gene",
  all.x = TRUE
)

dom_long <- melt(
  dom_plot,
  id.vars = c("gene", "Final_multiomic_index", "Final_evidence_tier"),
  variable.name = "Domain",
  value.name = "Support"
)

dom_long[, Domain :=
  factor(
    domain_labels[as.character(Domain)],
    levels = unname(domain_labels)
  )
]

dom_long[, gene :=
  factor(
    gene,
    levels = rev(gene_order)
  )
]

fig1 <- ggplot(
  dom_long,
  aes(x = Domain, y = gene, fill = Support)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", Support)),
    size = 3.1
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08306B",
    limits = c(0, 1),
    name = "Support"
  ) +
  labs(
    title = "Figure 1. Six-domain multi-omic support matrix",
    subtitle = "Gene-wise support across RNA and ATAC evidence domains",
    x = "",
    y = ""
  ) +
  publication_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

save_both(
  fig1,
  "Figure_1_six_domain_support_heatmap",
  width = 10,
  height = 6.5
)

# ============================================================
# 4. Figure 2: Final multikinase ranking
# ============================================================

plot2 <- copy(compact)
plot2[, gene :=
  factor(
    gene,
    levels = rev(gene_order)
  )
]

fig2 <- ggplot(
  plot2,
  aes(
    x = gene,
    y = Final_multiomic_index,
    fill = Final_evidence_tier
  )
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  coord_flip() +
  geom_text(
    aes(
      label = paste0(
        "Domains=", JAK3_domains_final
      )
    ),
    hjust = -0.1,
    size = 3.3
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Figure 2. Final multikinase ranking",
    subtitle = "Bars show final multi-omic index; labels show number of domains at least as strong as JAK3",
    x = "",
    y = "Final multi-omic index",
    fill = "Evidence tier"
  ) +
  publication_theme +
  theme(
    legend.position = "bottom"
  )

save_both(
  fig2,
  "Figure_2_final_multikinase_ranking",
  width = 9.5,
  height = 6.5
)

# ============================================================
# 5. Figure 3: JAK3 benchmark profile
# ============================================================

selected_genes <- c(
  "JAK3",
  "HCK",
  "EPHB2",
  "LYN",
  "IRAK2"
)

profile_long <- dom_long[
  gene %in% selected_genes
]

profile_long[, gene :=
  factor(
    as.character(gene),
    levels = selected_genes
  )
]

fig3 <- ggplot(
  profile_long,
  aes(
    x = Domain,
    y = Support,
    color = gene,
    group = gene
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    title = "Figure 3. JAK3 benchmark domain profile",
    subtitle = "Comparison of the main follow-up candidates against the JAK3 reference profile",
    x = "",
    y = "Domain support",
    color = "Gene"
  ) +
  publication_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

save_both(
  fig3,
  "Figure_3_JAK3_benchmark_profile",
  width = 11,
  height = 6.5
)

# ============================================================
# 6. Figure 4: Significant ATAC-RNA heatmap
# ============================================================

atac_heat <- atac_sig[
  ,
  .(
    Best_rho =
      Spearman_rho[
        which.max(abs(Spearman_rho))
      ],
    Significant_enhancers = .N
  ),
  by = .(gene, Cancer)
]

cancer_order <- atac_heat[
  ,
  .(
    Total_sig = sum(Significant_enhancers)
  ),
  by = Cancer
][
  order(-Total_sig, Cancer)
]$Cancer

atac_heat[, Cancer :=
  factor(
    Cancer,
    levels = cancer_order
  )
]

atac_heat[, gene :=
  factor(
    gene,
    levels = rev(gene_order)
  )
]

fig4 <- ggplot(
  atac_heat,
  aes(
    x = Cancer,
    y = gene,
    fill = Best_rho
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = Significant_enhancers),
    size = 3.0
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Best\nrho"
  ) +
  labs(
    title = "Figure 4. Significant ATAC-RNA associations across cancers",
    subtitle = "Cell values show number of significant enhancer associations; color shows strongest correlation per gene-cancer pair",
    x = "Cancer type",
    y = ""
  ) +
  publication_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

save_both(
  fig4,
  "Figure_4_significant_ATAC_RNA_heatmap",
  width = 12,
  height = 6.8
)

# ============================================================
# 7. Figure 5: RNA vs ATAC landscape
# ============================================================

jak3_rna <- compact[
  gene == "JAK3",
  RNA_evidence_index
]

jak3_atac <- compact[
  gene == "JAK3",
  ATAC_evidence_support
]

fig5 <- ggplot(
  compact,
  aes(
    x = RNA_evidence_index,
    y = ATAC_evidence_support,
    color = Final_evidence_tier,
    size = JAK3_domains_final
  )
) +
  geom_vline(
    xintercept = jak3_rna,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = jak3_atac,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_point(alpha = 0.9) +
  geom_text(
    aes(label = gene),
    vjust = -0.9,
    size = 3.5,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "Figure 5. RNA versus ATAC evidence landscape",
    subtitle = "Dashed lines mark the JAK3 benchmark for RNA and ATAC evidence",
    x = "RNA evidence index",
    y = "ATAC evidence support",
    color = "Evidence tier",
    size = "JAK3-level\ndomains"
  ) +
  publication_theme +
  theme(
    legend.position = "right"
  )

save_both(
  fig5,
  "Figure_5_RNA_vs_ATAC_landscape",
  width = 8.5,
  height = 6.8
)

cat("\nFigures written to:\n")
cat(figdir, "\n")

cat("\nGenerated files:\n")
print(list.files(figdir, full.names = FALSE))

cat("\n14 COMPLETE\n")

