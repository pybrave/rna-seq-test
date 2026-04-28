library(tidyverse)
library(jsonlite)

params <- jsonlite::fromJSON("params.json")

gene_id_mapping <- read_tsv(params$gene_id_mapping)

kegg_enrichment <- read_tsv(params$kegg_enrichment)

# species <- 

df_sig_comp_vec <- setNames(gene_id_mapping$logfc, gene_id_mapping$kegg_id) |> as.list()

# enrich_res_df_json <- kegg_enrichment |>
#   mutate(organism ="mmu" ) |>
#   mutate(pathwayId=str_replace(ID,"mmu",""))

enrich_res_df_json <- kegg_enrichment |>
  mutate(pathwayId= ID )

list(compound = df_sig_comp_vec,list=enrich_res_df_json) |> toJSON(auto_unbox = TRUE) |>
  write("output/kegg_map.vis")
