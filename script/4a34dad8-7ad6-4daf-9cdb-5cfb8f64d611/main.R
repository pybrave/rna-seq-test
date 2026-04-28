library(tidyverse)
library(DESeq2)
params <- jsonlite::fromJSON("params.json")

df <- read_tsv(params$count$content)  

feature_var <- params$count$feature_var
control <- params$count$control
treatment <- params$count$treatment
exp <- df |>
  select(c(feature_var, control, treatment)) |>
  column_to_rownames(feature_var)

metadata <- data.frame(
  sample=c(control,treatment),
  group=c(rep("control",length(control)),rep("treatment",length(treatment)))
) |> column_to_rownames("sample")


exp <- exp[rowSums(exp) > 0, ]
all(rownames(metadata) %in% colnames(exp))
if(!identical(rownames(metadata) ,colnames(exp))){
    stop("样本名称不匹配，请检查metadata和表达矩阵的列名是否一致")
}

dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = metadata,
                              design = ~ group)
dds <- DESeq(dds)
res <- results(dds)
as.data.frame(res) |>
  rownames_to_column("feature") |>
  left_join(rownames_to_column(exp,"feature"),by="feature") |>
  write_tsv(str_glue("output/all.deg.tsv"))


