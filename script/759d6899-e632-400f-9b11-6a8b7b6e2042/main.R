## BWA比对脚本 (R环境封装调用，实际执行使用系统shell)

# 读取参数
params <- jsonlite::fromJSON("params.json", simplifyVector = FALSE)

# 提取参数
input_fastq1 <- params$input_fastq1
input_fastq2 <- params$input_fastq2
reference_genome <- params$reference_genome
output_sam <- params$output_sam
output_prefix <- params$output_name
threads <- params$threads
mem_limit_gb <- params$mem_limit_gb
aligner_options <- params$aligner_options

# 验证必要参数
if (is.null(input_fastq1) || is.null(input_fastq2) || is.null(reference_genome)) {
  stop("Missing required input: input_fastq1, input_fastq2, or reference_genome")
}

# 构建BWA命令
bwa_cmd <- sprintf(
  "bwa mem -t %d -M %s %s %s %s > %s",
  threads,
  aligner_options,
  shQuote(reference_genome),
  shQuote(input_fastq1),
  shQuote(input_fastq2),
  shQuote(output_sam)
)

message("Running BWA command:", bwa_cmd)

# 执行比对（需确保系统已安装bwa）
system(bwa_cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)

# 检查输出文件是否存在
if (!file.exists(output_sam)) {
  stop(sprintf("BWA alignment failed: output SAM file '%s' not generated.", output_sam))
}

# 可选：转换为 BAM 并排序（若 samtools 可用）
if (requireNamespace("system2", quietly = TRUE)) {
  bam_file <- gsub("\\.sam$", "\\.bam", output_sam)
  sorted_bam <- gsub("\\.bam$", "_sorted.bam", bam_file)
  
  # samtools view
  view_cmd <- sprintf("samtools view -Sb %s > %s", shQuote(output_sam), shQuote(bam_file))
  message("Converting SAM to BAM: ", view_cmd)
  system(view_cmd)
  
  # samtools sort
  sort_cmd <- sprintf("samtools sort -@ %d -m %dg %s -o %s", 
                      threads, mem_limit_gb, shQuote(bam_file), shQuote(sorted_bam))
  message("Sorting BAM: ", sort_cmd)
  system(sort_cmd)
  
  # samtools index
  index_cmd <- sprintf("samtools index %s", shQuote(sorted_bam))
  message("Indexing BAM: ", index_cmd)
  system(index_cmd)
  
  if (file.exists(sorted_bam)) {
    message(sprintf("Sorted & indexed BAM saved to: %s", sorted_bam))
  }
}

message(sprintf("BWA alignment completed. Output saved to: %s", output_sam))
