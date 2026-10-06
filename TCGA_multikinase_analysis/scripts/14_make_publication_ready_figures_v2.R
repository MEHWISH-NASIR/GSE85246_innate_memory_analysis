
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
  "figures/final_integration_v2"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("STEP 14 V2: PUBLICATION-READY FIGURES\n")
cat("========================================\n")

# ============================================================
# 1. Load data
# ============================================================

compact <- fread(
  file.path(
    indir,
    "13_supervisor_final_multikinase_comparison.csv"
  ),
  nThread = 1
)

domains <- fread(
  file.path(
    indir,
    "13_final_domain_support_matrix.csv"
  ),
  nThread = 1
)

atac_sig <- fread(
  file.path(
    indir,
    "13_significant_ATAC_RNA_associations.csv"
  ),
  nThread = 1
)

cat("Compact rows:", nrow(compact), "\n")
cat("Domain rows:", nrow(domains), "\n")
cat("Significant ATAC rows:", nrow(atac_sig), "\n")

# ============================================================
# 2. Discovery classes
# ============================================================

compact[, Discovery_class :=
  fifelse(
    gene == "JAK3",
    "Reference",
    fifelse(
      grepl("Washout", Discovery_support, ignore.case = TRUE),
      "Washout-emergent",
      "Strict persistent"
    )
  )
]

compact[, Discovery_class :=
  factor(
    Discovery_class,
    levels = c(
      "Reference",
      "Strict persistent",
      "Washout-emergent"
    )
  )
]

# ============================================================
# 3. Gene ordering
# ============================================================

gene_order <- compact[
  order(
    -Final_multiomic_index,
    gene
  ),
  gene
]

# ============================================================
# 4. Shared theme
# ============================================================

publication_theme <- theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(12, 16, 12, 12)
  )

save_both <- function(plot_obj, filename, width, height) {

  ggsave(
    filename = file.path(
      figdir,
      paste0(filename, ".png")
    ),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = file.path(
      figdir,
      paste0(filename, ".pdf")
    ),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
}

# ============================================================
# 5. Figure 1
# Six-domain support heatmap
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
  Progression_support = "Stage/grade",
  Survival_support = "Survival",
  Immune_adjusted_support = "Immune-adjusted",
  ATAC_evidence_support = "ATAC-RNA"
)

dom <- merge(
  domains[, c("gene", domain_cols), with = FALSE],
  compact[, .(gene, Discovery_class)],
  by = "gene"
)

dom_long <- melt(
  dom,
  id.vars = c(
    "gene",
    "Discovery_class"
  ),
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
  aes(
    x = Domain,
    y = gene,
    fill = Support
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", Support)
    ),
    size = 3.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#154360",
    limits = c(0, 1),
    name = "Support"
  ) +
  labs(
    title = "Six-domain multi-omic support across candidate genes",
    subtitle = "Normalized support across RNA, clinical and ATAC-RNA evidence domains",
    x = NULL,
    y = NULL
  ) +
  publication_theme +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

save_both(
  fig1,
  "Figure_1_six_domain_support_heatmap_v2",
  10.5,
  6.5
)

# ============================================================
# 6. Figure 2
# Final multiomic index
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
    fill = Discovery_class
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.25
  ) +
  coord_flip() +
  geom_text(
    aes(
      label = sprintf("%.3f", Final_multiomic_index)
    ),
    hjust = -0.15,
    size = 3.4
  ) +
  scale_fill_manual(
    values = c(
      "Reference" = "#D95F5F",
      "Strict persistent" = "#3A86A8",
      "Washout-emergent" = "#8C6BB1"
    ),
    name = "Discovery class"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.15)
    )
  ) +
  labs(
    title = "Integrated multi-omic evidence across candidate genes",
    subtitle = "Final index combines five RNA/clinical domains with ATAC-RNA regulatory coupling",
    x = NULL,
    y = "Final multi-omic evidence index"
  ) +
  publication_theme +
  theme(
    legend.position = "bottom"
  )

save_both(
  fig2,
  "Figure_2_final_multiomic_index_v2",
  9.5,
  6.5
)

# ============================================================
# 7. Figure 3
# Benchmark profile
# ============================================================

selected_genes <- c(
  "JAK3",
  "HCK",
  "EPHB2",
  "LYN",
  "IRAK2"
)

profile <- dom_long[
  as.character(gene) %in% selected_genes
]

profile[, gene :=
  factor(
    as.character(gene),
    levels = selected_genes
  )
]

fig3 <- ggplot(
  profile,
  aes(
    x = Domain,
    y = Support,
    group = gene,
    color = gene
  )
) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    size = 2.8
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    title = "Multi-domain evidence profiles relative to JAK3",
    subtitle = "Comparison of the principal follow-up candidates across six evidence domains",
    x = NULL,
    y = "Normalized domain support",
    color = "Gene"
  ) +
  publication_theme +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "bottom"
  )

save_both(
  fig3,
  "Figure_3_JAK3_candidate_profiles_v2",
  11,
  6.5
)

# ============================================================
# 8. Figure 4
# Cancer-specific ATAC-RNA coupling
# ============================================================

atac_heat <- atac_sig[
  ,
  .(
    Best_rho =
      Spearman_rho[
        which.max(abs(Spearman_rho))
      ],
    Significant_elements = .N
  ),
  by = .(
    gene,
    Cancer
  )
]

cancer_order <- atac_heat[
  ,
  .(
    Total_sig = sum(Significant_elements)
  ),
  by = Cancer
][
  order(
    -Total_sig,
    Cancer
  ),
  Cancer
]

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
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = Significant_elements
    ),
    size = 3.2
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  labs(
    title = "Cancer-specific ATAC-RNA coupling of candidate regulatory elements",
    subtitle = "Color shows strongest significant correlation; numbers show significant linked elements. Blank cells indicate no FDR-significant association.",
    x = "Cancer type",
    y = NULL
  ) +
  publication_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_both(
  fig4,
  "Figure_4_cancer_specific_ATAC_RNA_coupling_v2",
  12.5,
  7
)

# ============================================================
# 9. Figure 5
# RNA versus ATAC landscape
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
    y = ATAC_evidence_support
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
  geom_point(
    aes(
      fill = Discovery_class,
      size = JAK3_domains_final
    ),
    shape = 21,
    color = "black",
    stroke = 0.4
  ) +
  geom_text(
    aes(
      label = gene
    ),
    vjust = -1.0,
    size = 3.6,
    check_overlap = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Reference" = "#D95F5F",
      "Strict persistent" = "#3A86A8",
      "Washout-emergent" = "#8C6BB1"
    ),
    name = "Discovery class"
  ) +
  scale_size_continuous(
    range = c(3, 7),
    breaks = 1:6,
    name = "Domains ≥ JAK3"
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.15, 0.15)
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.12, 0.15)
    )
  ) +
  labs(
    title = "RNA versus ATAC evidence landscape",
    subtitle = "Dashed lines indicate the JAK3 benchmark in each evidence axis",
    x = "RNA evidence index",
    y = "ATAC evidence support"
  ) +
  publication_theme +
  theme(
    legend.position = "right"
  )

save_both(
  fig5,
  "Figure_5_RNA_vs_ATAC_landscape_v2",
  9.5,
  7.2
)

# ============================================================
# 10. Finish
# ============================================================

cat("\nFigures written to:\n")
cat(figdir, "\n")

cat("\nGenerated files:\n")
print(
  list.files(
    figdir,
    full.names = FALSE
  )
)

cat("\n14 V2 COMPLETE\n")

