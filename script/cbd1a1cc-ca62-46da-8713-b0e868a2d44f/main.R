library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)

# 读取参数
params <- jsonlite::fromJSON("params.json", simplifyVector = FALSE)
file_path <- params$input_file$content
contrast_str <- params$contrast
pvalue_cutoff <- params$pvalue_cutoff
padj_cutoff <- params$padj_cutoff
lfc_cutoff <- params$lfc_cutoff
output_prefix <- params$output_name
plot_width <- params$plot_width
plot_height <- params$plot_height

# 读取计数矩阵和分组信息
# 假设 input_file 中 count_matrix 指向 TSV 计数表（行=基因，列=样本），sample_group 指向分组向量（TSV 或 RDS）
count_mat_path <- params$input_file$count_matrix$content
sample_group_path <- params$input_file$sample_group$content

counts_df <- readr::read_tsv(count_mat_path, show_col_types = FALSE)
sample_group_df <- readr::read_tsv(sample_group_path, show_col_types = FALSE)

# 构建 DESeqDataSet
rownames(counts_df) <- counts_df[[1]]
counts_df <- counts_df[-1]
colnames(counts_df) <- colnames(sample_group_df)

# 确保列名一致
stopifnot(all(colnames(counts_df) %in% colnames(sample_group_df)))

# 创建 sample table
sample_table <- data.frame(
  row.names = colnames(counts_df),
  condition = sample_group_df[match(colnames(counts_df), colnames(sample_group_df)), 1],
  stringsAsFactors = FALSE
)

# 转为 matrix（整数）
count_matrix <- as.matrix(counts_df)

# 创建 DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_table,
  design = ~ condition
)

# 运行 DESeq
message("Running DESeq pipeline...")
dds <- DESeq(dds)

# 提取结果
res <- results(dds, contrast = strsplit(contrast_str, "_vs_")[[1]])
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# 输出差异结果
res_tsv_path <- file.path("output", paste0(output_prefix, "_results.tsv"))
readr::write_tsv(res_df, res_tsv_path)
message(sprintf("DESeq2 results saved to: %s", res_tsv_path))

# 绘制 PCA
message("Plotting PCA...")
rld <- rlog(dds, blind = FALSE)
pca_data <- plotPCA(rld, intgroup = "condition") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path("output", paste0(output_prefix, "_pca.png")), plot = pca_data, width = plot_width, height = plot_height, dpi = 300)

# MA plot
message("Plotting MA plot...")
ma_plot <- plotMA(res, alpha = pvalue_cutoff, ylim = c(-5, 5)) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path("output", paste0(output_prefix, "_ma.png")), plot = ma_plot, width = plot_width, height = plot_height, dpi = 300)

# Volcano plot
message("Plotting volcano plot...")
volcano <- EnhancedVolcano(res_df,
  lab = res_df$gene_id,
  x = "log2FoldChange",
  y = "padj",
  xlim = c(-max(abs(res_df$log2FoldChange)), max(abs(res_df$log2FoldChange))),
  title = paste0("Volcano Plot (", contrast_str, ")"),
  FCcutoff = lfc_cutoff,
  pCutoff = padj_cutoff,
  drawConnectors = TRUE,
  widthConnector = 0.5
)
ggsave(filename = file.path("output", paste0(output_prefix, "_volcano.png")), plot = volcano, width = plot_width, height = plot_height, dpi = 300)

# Heatmap of top 50 DEGs
message("Plotting heatmap...")
top_genes <- head(rownames(res), 50)
mat <- assay(rld)[top_genes, ]
mat_z <- t(apply(mat, 1, scale))

pheatmap(mat_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_col = sample_table,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = paste0("Top 50 DEGs Heatmap (", contrast_str, ")")
)

pdf_heatmap_path <- file.path("output", paste0(output_prefix, "_heatmap.pdf"))
dev.copy2pdf(file = pdf_heatmap_path)
message(sprintf("Heatmap PDF saved to: %s", pdf_heatmap_path))

# Summary table
summary_df <- data.frame(
  TotalGenes = nrow(count_matrix),
  DEGs_Total = sum(!is.na(res_df$padj) & res_df$padj <= padj_cutoff),
  Upregulated = sum(!is.na(res_df$padj) & res_df$padj <= padj_cutoff & res_df$log2FoldChange > lfc_cutoff),
  Downregulated = sum(!is.na(res_df$padj) & res_df$padj <= padj_cutoff & res_df$log2FoldChange < -lfc_cutoff)
)
summary_tsv_path <- file.path("output", paste0(output_prefix, "_summary.tsv"))
readr::write_tsv(summary_df, summary_tsv_path)
message(sprintf("Summary saved to: %s", summary_tsv_path))

message("DESeq2 analysis completed successfully.")
