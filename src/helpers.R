# helpers.R - Shared helper functions for RNA-seq differential expression analysis

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(stringr)
    library(forcats)
    library(ggplot2)
    library(ggrepel)
})

FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 1

cond_cols <- c("healthy" = "grey40", 
               "ad_non_lesion" = "#8DA0CB", "ad_lesion" = "#FC8D62",
               "pso_non_lesion" = "#66C2A5", "pso_lesion" = "#E78AC3"
)
disease_cols <- c("CTRL" = "grey30", "AD" = "#E69F00", "PSO" = "#56B4E9")
lesion_cols <- c("healthy" = "grey80", "non-lesional" = "skyblue3", "lesional" = "firebrick3")

check_sample_alignment <- function(counts, meta) {
    count_samples <- colnames(counts)
    meta_samples <- meta$sample

    aligned <- all(count_samples == meta_samples)
    print("Sample alignment: ", aligned)

    if (!aligned) {
        common <- intersect(count_samples, meta_samples)
        print("Reordering counts to match metadata order...")
        counts <- counts[, meta_samples[meta_samples %in% count_samples]]
    }

    invisible(counts)
}

mark_deg <- function(res) {
    res %>%
        mutate(
            DEG = (adj.P.Val <= FDR_CUTOFF & abs(logFC) >= LFC_CUTOFF),
            direction = case_when(DEG & logFC > 0 ~ "Up", DEG & logFC < 0 ~ "Down", TRUE ~ "NS")
        )
}

run_go_enrich <- function(genes, universe) {
  ego <- clusterProfiler::enrichGO(
    gene = genes, universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 0.999, qvalueCutoff = 0.999, readable = TRUE
  )

  ego@result %>%
    as_tibble() %>%
    filter(p.adjust <= FDR_CUTOFF) %>%
    mutate(
      neglog10_padj = -log10(p.adjust),
      Description = stringr::str_trunc(Description, 70)
    )
}
fisher_enrichment <- function(deg_genes, gene_set, universe) {
    gene_set_u <- intersect(gene_set, universe)
    deg_u <- intersect(deg_genes, universe)
    in_set <- universe %in% gene_set_u
    is_deg <- universe %in% deg_u
    ft <- fisher.test(matrix(c(
        sum(is_deg & in_set), sum(is_deg & !in_set),
        sum(!is_deg & in_set), sum(!is_deg & !in_set)
    ), nrow = 2))
    tibble(p_value = ft$p.value, odds_ratio = ft$estimate)
}

plot_go_bar <- function(go_res, title) {
    top_go <- go_res %>%
        arrange(p.adjust) %>%
        slice_head(n = 15) %>%
        mutate(Description = fct_reorder(Description, neglog10_padj))
    ggplot(top_go, aes(x = neglog10_padj, y = Description)) +
        geom_col() +
        theme_classic(base_size = 12) +
        labs(x = expression(-log[10] ~ adjusted ~ italic(p)), y = NULL, title = title)
}

plot_volcano <- function(res, title) {
    top_labels <- res %>%
        filter(DEG) %>%
        arrange(adj.P.Val) %>%
        slice_head(n = 10)
    ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
        geom_point(size = 0.8, alpha = 0.6) +
        scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70")) +
        geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
        geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey40") +
        ggrepel::geom_text_repel(data = top_labels, aes(label = gene), size = 2.5, max.overlaps = 15) +
        theme_classic(base_size = 12) +
        labs(x = expression(log[2] ~ FC), 
             y = expression(-log[10] ~ adjusted ~ italic(p)),
            title = title, color = "Regulation") +
        theme(legend.position = "top")
}

plot_ma <- function(res, title) {
    ggplot(res, aes(x = AveExpr, y = logFC, color = direction)) +
        geom_point(size = 0.6, alpha = 0.5) +
        scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70")) +
        geom_hline(
            yintercept = c(-LFC_CUTOFF, 0, LFC_CUTOFF),
            linetype = c("dashed", "solid", "dashed"), color = c("grey40", "black", "grey40")
        ) +
        theme_classic(base_size = 12) +
        labs(
            x = "Average Expression (log2 CPM)", y = expression(log[2] ~ FC),
            title = title, color = "Regulation"
        ) +
        theme(legend.position = "top")
}

get_deg_genes <- function(res, n = NULL) {
    out <- res %>%
        filter(DEG) %>%
        arrange(adj.P.Val)
    if (!is.null(n)) out <- slice_head(out, n = n)
    pull(out, gene)
}

plot_pca <- function(pca_df, pc_x, pc_y, title) {
    ggplot(pca_df, aes(x = .data[[pc_x]], y = .data[[pc_y]], color = condition)) +
        geom_point(size = 3, alpha = 0.9) +
        scale_color_manual(values = cond_cols) +
        theme_classic(base_size = 13) +
        labs(title = title, color = "Condition")
}

plot_fc_concordance <- function(res_x, res_y, x_label, y_label, title_prefix) {
    joined <- res_x %>%
        select(gene, logFC_x = logFC, DEG_x = DEG) %>%
        inner_join(res_y %>% select(gene, logFC_y = logFC, DEG_y = DEG), by = "gene") %>%
        mutate(class = case_when(
            DEG_x & DEG_y ~ "DE in both",
            DEG_x & !DEG_y ~ paste0("DE in ", gsub(" vs.*", "", x_label), " only"),
            !DEG_x & DEG_y ~ paste0("DE in ", gsub(" vs.*", "", y_label), " only"),
            TRUE ~ "Not DE"
        ))
    rho <- cor(joined$logFC_x, joined$logFC_y, method = "spearman")
    ggplot(joined, aes(x = logFC_x, y = logFC_y, color = class)) +
        geom_point(size = 0.8, alpha = 0.7) +
        scale_color_manual(values = c(
            "DE in both" = "red", "Not DE" = "grey80",
            setNames("orange", paste0("DE in ", gsub(" vs.*", "", x_label), " only")),
            setNames("blue", paste0("DE in ", gsub(" vs.*", "", y_label), " only"))
        )) +
        theme_classic(base_size = 13) +
        labs(
            x = paste0("log2 FC (", x_label, ")"), y = paste0("log2 FC (", y_label, ")"),
            title = paste0(title_prefix, ": Spearman r = ", round(rho, 2))
        )
}

plot_deg_counts <- function(res_list, labels, title) {
    counts_df <- tibble(
        contrast = labels,
        up = sapply(res_list, function(r) sum(r$direction == "Up")),
        down = sapply(res_list, function(r) sum(r$direction == "Down"))
    ) %>%
        tidyr::pivot_longer(cols = c(up, down), names_to = "regulation", values_to = "n") %>%
        mutate(regulation = recode(regulation, "up" = "Upregulated", "down" = "Downregulated"))

    ggplot(counts_df, aes(x = contrast, y = n, fill = regulation)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = c("Upregulated" = "#D73027", "Downregulated" = "#4575B4")) +
        theme_classic(base_size = 13) +
        labs(x = NULL, y = "Number of DEGs (FDR <= 0.05, |log2FC| >= 1)", fill = "Regulation", title = title)
}

plot_gene_boxplot <- function(logexpr, genes, meta, title) {
    expr_long <- t(logexpr[genes, meta$sample, drop = FALSE]) %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
        left_join(meta, by = "sample") %>%
        mutate(gene = factor(gene, levels = genes))

    ggplot(expr_long, aes(x = condition, y = expr, color = condition)) +
        geom_boxplot(outlier.size = 0.7, alpha = 0.8) +
        scale_color_manual(values = cond_cols) +
        facet_wrap(~gene, scales = "free_y", ncol = 4) +
        theme_classic(base_size = 11) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
        labs(x = NULL, y = "log2 expression (voom)", title = title)
}

plot_expr_heatmap <- function(logexpr, genes, meta, title) {
    mat <- logexpr[genes, meta$sample, drop = FALSE]
    mat_scaled <- t(scale(t(mat)))

    ha <- ComplexHeatmap::HeatmapAnnotation(
        Condition = meta$condition,
        col = list(Condition = cond_cols)
    )

    ComplexHeatmap::Heatmap(mat_scaled,
        name = "Z-score",
        top_annotation = ha,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 7),
        column_title = title,
        cluster_columns = TRUE,
        cluster_rows = TRUE
    )
}
