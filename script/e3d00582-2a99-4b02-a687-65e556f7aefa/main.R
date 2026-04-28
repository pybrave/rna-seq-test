library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

# 读取参数wwwww
params <- jsonlite::fromJSON("params.json", simplifyVector = FALSE)

count_matrix_path <- params$count_matrix$content
sample_annotation_path <- params$sample_annotation$content

# 读取输入数据
count_matrix <- readr::read_tsv(count_matrix_path, show_col_types = FALSE)
sample_annotation <- readr::read_tsv(sample_annotation_path, show_col_types = FALSE)

# 构建DESeqDataSet
# 假设 count_matrix 第一列为基因名，其余列为样本
rownames(count_matrix) <- count_matrix[[1]]
count_matrix <- count_matrix[,-1]

# 确保列名与 sample_annotation 中的 sample_id 匹配
colnames(count_matrix) <- gsub("\\.\\d+$", "", colnames(count_matrix))
sample_annotation$sample_id <- gsub("\\.\\d+$", "", sample_annotation$sample_id)

# 按 sample_id 排序以对齐
count_matrix <- count_matrix[, match(sample_annotation$sample_id, colnames(count_matrix))]

# 创建 DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_matrix),
  colData = sample_annotation,
  design = ~ condition
)

# 运行 DESeq
dds <- DESeq(dds)

# 提取 results
res <- results(dds,
               contrast = c("condition", params$test_condition, params$ref_condition),
               alpha = params$pvalue_cutoff
)

# 转换为 data.frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]  # 移除 NA
res_df <- res_df[order(res_df$padj), ]   # 按 padj 排序

# 输出差异结果表格
output_prefix <- params$output_name
results_tsv <- file.path("output", paste0(output_prefix, "_diff_genes.tsv"))
readr::write_tsv(res_df, results_tsv)
message(sprintf("Differential expression results saved to: %s", results_tsv))

# 绘制 MA plot
ma_plot <- plotMA(res, alpha = params$pvalue_cutoff, 
                   ylim = c(-5, 5), 
                   main = "MA Plot",
                   alpha = 0.7)
ma_png <- file.path("output", paste0(output_prefix, "_MA_plot.png"))
ggsave(ma_png, plot = ma_plot, width = params$plot_width, height = params$plot_height, dpi = 300)

# 绘制火山图
volcano_plot <- EnhancedVolcano(res_df,
                                lab = res_df$gene,
                                x = "log2FoldChange",
                                y = "padj",
                                xlim = c(-5, 5),
                                title = "Volcano Plot",
                                pCutoff = params$padj_cutoff,
                                FCcutoff = params$lfc_cutoff,
                                pointSize = 3.0,
                                labSize = 3.0)
volcano_png <- file.path("output", paste0(output_prefix, "_volcano_plot.png"))
ggsave(volcano_png, plot = volcano_plot, width = params$plot_width, height = params$plot_height, dpi = 300)

# 绘制热图（前50个显著差异基因）
res_sig <- subset(res_df, padj < params$padj_cutoff & abs(log2FoldChange) > params$lfc_cutoff)
if (nrow(res_sig) > 0) {
  top_genes <- head(rownames(res_sig), 50)
  if (length(top_genes) == 0) top_genes <- rownames(res_sig)[1]
  
  # 获取原始 counts 并标准化
  rld <- rlog(dds, blind = FALSE)
  mat <- assay(rld)[top_genes, ]
  
  # 样本聚类热图
  pheatmap(mat,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           annotation_col = sample_annotation[, c("condition")],
           fontsize_row = 8,
           fontsize_col = 10)
  
  heatmap_png <- file.path("output", paste0(output_prefix, "_heatmap.png"))
  dev.copy(png, heatmap_png, width = params$plot_width * 100, height = params$plot_height * 100)
  dev.off()
}

# 输出 normalized counts
norm_counts <- counts(dds, normalized = TRUE)
norm_tsv <- file.path("output", paste0(output_prefix, "_normalized_counts.tsv"))
readr::write_tsv(as.data.frame(norm_counts), norm_tsv)
message(sprintf("Normalized counts saved to: %s", norm_tsv))

# 输出 size factors
sf_df <- data.frame(sample = colnames(dds), size_factor = sizeFactors(dds))
sf_tsv <- file.path("output", paste0(output_prefix, "_size_factors.tsv"))
readr::write_tsv(sf_df, sf_tsv)
message(sprintf("Size factors saved to: %s", sf_tsv))

message("DESeq2 analysis completed successfully.")