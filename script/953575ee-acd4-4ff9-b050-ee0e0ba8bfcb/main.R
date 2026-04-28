library(tidyverse)
library(jsonlite)
library(clusterProfiler)
library(KEGGREST)

# 2: In readChar(con, 5L, useBytes = TRUE) :
#   cannot open compressed file '//.local/share/clusterProfiler/kegg_category.rda', probable reason 'No such file or directory'
# Output saved to: output

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

to_number <- function(x, default, lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(x %||% default))
  if (is.na(value)) value <- default
  value <- max(lower, min(upper, value))
  value
}

normalize_filter_by <- function(x) {
  value <- tolower(as.character(x %||% "p"))
  if (!value %in% c("p", "q")) "p" else value
}

normalize_species <- function(x) {
  allowed <- c("hsa", "mmu", "rno", "dre", "ath", "sce", "cel", "dme")
  value <- tolower(as.character(x %||% "hsa"))
  if (!value %in% allowed) {
    warning(sprintf("Unsupported species '%s', fallback to '%s'", value, "hsa"))
    value <- "hsa"
  }
  value
}

normalize_gene_id_type <- function(x) {
  allowed <- c("symbol", "entrez", "kegg", "uniprot")
  value <- tolower(as.character(x %||% "symbol"))
  if (!value %in% allowed) {
    warning(sprintf("Unsupported gene_id_type '%s', fallback to '%s'", value, "symbol"))
    value <- "symbol"
  }
  value
}
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

get_orgdb_package <- function(species) {
  species_pkg <- c(
    hsa = "org.Hs.eg.db",
    mmu = "org.Mm.eg.db",
    rno = "org.Rn.eg.db",
    dre = "org.Dr.eg.db",
    ath = "org.At.tair.db",
    sce = "org.Sc.sgd.db",
    cel = "org.Ce.eg.db",
    dme = "org.Dm.eg.db"
  )
  species_pkg[[species]] %||% NA_character_
}

convert_symbol_to_entrez <- function(symbol_ids, species) {
  input_ids <- unique(as.character(symbol_ids))
  input_ids <- input_ids[!is.na(input_ids) & input_ids != ""]

  org_pkg <- get_orgdb_package(species)
  if (is.na(org_pkg)) {
    return(list(
      ids = character(),
      table = tibble(feature = input_ids, entrezid = NA_character_),
      status = "no_orgdb_for_species",
      package = "none",
      input_count = length(input_ids),
      mapped_count = 0,
      unmapped_count = length(input_ids)
    ))
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace(org_pkg, quietly = TRUE)) {
    # return(list(
    #   ids = character(),
    #   table = tibble(feature = input_ids, entrezid = NA_character_),
    #   status = "orgdb_missing",
    #   package = org_pkg,
    #   input_count = length(input_ids),
    #   mapped_count = 0,
    #   unmapped_count = length(input_ids)
    # ))
    install_if_missing(org_pkg)
    
  }

  orgdb <- getExportedValue(org_pkg, org_pkg)
  mapped <- suppressWarnings(
    AnnotationDbi::mapIds(
      x = orgdb,
      keys = input_ids,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
  )

  map_table <- tibble(
    feature = names(mapped),
    entrezid = as.character(unname(mapped))
  )

  mapped_ids <- map_table |>
    dplyr::filter(!is.na(entrezid) & entrezid != "") |>
    dplyr::pull(entrezid) |>
    unique()

  mapped_count <- map_table |>
    dplyr::filter(!is.na(entrezid) & entrezid != "") |>
    nrow()

  list(
    ids = mapped_ids,
    table = map_table,
    status = "ok",
    package = org_pkg,
    input_count = nrow(map_table),
    mapped_count = mapped_count,
    unmapped_count = nrow(map_table) - mapped_count
  )
}

parse_ratio <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  parts <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NA_real_)
  num <- suppressWarnings(as.numeric(parts[1]))
  den <- suppressWarnings(as.numeric(parts[2]))
  if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
  num / den
}

run_kegg_enrichment <- function(gene_ids, species, candidates = c("ncbi-geneid", "kegg", "uniprot")) {
  for (key_type in candidates) {
    res <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = gene_ids,
        organism = species,
        keyType = key_type,
        pvalueCutoff = 1,
        qvalueCutoff = 1
      ),
      error = function(e) NULL
    )
    org_pkg <- get_orgdb_package(species)
    res <- setReadable(res, OrgDb = org_pkg, keyType="ENTREZID")
    

    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      return(list(result = res, key_type = key_type))
    }
  }
  list(result = NULL, key_type = "none")
}

params <- jsonlite::fromJSON("params.json", simplifyVector = FALSE)
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

deg_path <- as.character(params$deg %||% "")
if (deg_path == "" || !file.exists(deg_path)) {
  stop(sprintf("Input deg file not found: %s", deg_path))
}

deg <- readr::read_tsv(deg_path, show_col_types = FALSE)

required_cols <- c("feature", "log2FoldChange", "pvalue", "padj")
missing_cols <- setdiff(required_cols, colnames(deg))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

filter_by <- normalize_filter_by(params$filter_by)
species <- normalize_species(params$species)
gene_id_type <- normalize_gene_id_type(params$gene_id_type)
top_n <- as.integer(to_number(params$top_n, default = 20, lower = 1, upper = 100))

sig_column <- if (filter_by == "q") "padj" else "pvalue"

cutoff <- to_number(params$cutoff, default = 0.05, lower = 0, upper = 1)
filter_deg <- deg |>
  dplyr::filter(!is.na(.data[[sig_column]]) & .data[[sig_column]] < cutoff) |>
  dplyr::select("feature", "log2FoldChange", "pvalue", "padj")

feature_list <- unique(as.character(filter_deg$feature))
feature_list <- feature_list[!is.na(feature_list) & feature_list != ""]

filtered_deg_file <- file.path(output_dir, "filtered_deg.tsv")
enrich_file <- file.path(output_dir, "kegg_enrichment.tsv")
plot_file <- file.path(output_dir, "kegg_enrichment_barplot.pdf")
id_map_file <- file.path(output_dir, "gene_id_mapping.tsv")

readr::write_tsv(filter_deg, filtered_deg_file)

enrich_df <- tibble()
used_key_type <- "none"
plot_status <- "not generated"
conversion_status <- "not_needed"
conversion_package <- "none"
conversion_input_count <- length(feature_list)
conversion_mapped_count <- length(feature_list)
conversion_unmapped_count <- 0

id_map_table <- tibble(feature = feature_list, entrezid = NA_character_)
# https://rest.kegg.jp/link/mmu/pathway
# https://rest.kegg.jp/list/pathway/mmu

# https://rest.kegg.jp/conv/ncbi-geneid/mmu
# R convert entrez to kegg: mmu:100000029



enrich_candidates <- c("ncbi-geneid", "kegg", "uniprot")
ids_for_enrich <- feature_list

if (gene_id_type == "symbol") {
  conversion <- convert_symbol_to_entrez(feature_list, species)
  id_map_table <- conversion$table
  conversion_status <- conversion$status
  conversion_package <- conversion$package
  conversion_input_count <- conversion$input_count
  conversion_mapped_count <- conversion$mapped_count
  conversion_unmapped_count <- conversion$unmapped_count
  ids_for_enrich <- conversion$ids
  enrich_candidates <- c("ncbi-geneid")
} else if (gene_id_type == "entrez") {
  enrich_candidates <- c("ncbi-geneid")
} else if (gene_id_type == "kegg") {
  enrich_candidates <- c("kegg")
} else if (gene_id_type == "uniprot") {
  enrich_candidates <- c("uniprot")
}

kegg_id <- keggConv(species, "ncbi-geneid")
kegg_df <- tibble::tibble(
  ncbi_geneid = sub("ncbi-geneid:", "", names(kegg_id)),
  kegg_id = unname(kegg_id)
)
# add kegg id to mapping table if possible
if (nrow(id_map_table) > 0 && "entrezid" %in% colnames(id_map_table)) {
  filter_deg_fc <- filter_deg |>
    select(feature,logfc = log2FoldChange   )
  id_map_table <- id_map_table |>
    dplyr::left_join(kegg_df, by = c("entrezid" = "ncbi_geneid")) |>
    dplyr::left_join(filter_deg_fc, by="feature")
  
  id_map_table  <- id_map_table |> na.omit()
}

if (nrow(enrich_df) == 0 && length(ids_for_enrich) > 0) {
  enrich_run <- run_kegg_enrichment(ids_for_enrich, species, candidates = enrich_candidates)
  saveRDS(enrich_run$result, file = "output/kegg_enrich.rds")
  # barplot_p <- barplot(enrich_run$result, showCategory=20)
  # ggsave(str_glue("output/barplot.png"),barplot_p,  width = 10, height = 7,units = "in")
  # ggsave(str_glue("output/barplot.download.pdf"),barplot_p,  width = 10, height = 7,units = "in")
  # 
  
  used_key_type <- enrich_run$key_type

  if (!is.null(enrich_run$result)) {
    enrich_df <- as_tibble(as.data.frame(enrich_run$result)) |>
      dplyr::mutate(
        GeneRatioNumeric = vapply(GeneRatio, parse_ratio, numeric(1)),
        Log10Padjust = -log10(p.adjust)
      ) |>
      dplyr::arrange(p.adjust)
  }
}




readr::write_tsv(id_map_table, id_map_file)

readr::write_tsv(enrich_df, enrich_file)

# if (nrow(enrich_df) > 0) {
#   plot_df <- enrich_df |>
#     dplyr::slice_head(n = top_n) |>
#     dplyr::mutate(Description = forcats::fct_reorder(Description, Log10Padjust))
# 
#   p <- ggplot(plot_df, aes(x = Description, y = Log10Padjust, fill = GeneRatioNumeric)) +
#     geom_col(width = 0.75) +
#     coord_flip() +
#     scale_fill_gradient(low = "#9ecae1", high = "#08519c", na.value = "#bdbdbd") +
#     labs(
#       title = sprintf("KEGG enrichment (%s)", species),
#       x = "Pathway",
#       y = "-log10(adjusted p-value)",
#       fill = "Gene ratio"
#     ) +
#     theme_bw(base_size = 12) +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       axis.text.y = element_text(size = 9),
#       legend.position = "right"
#     )
# 
#   ggsave(plot_file, p, width = 10, height = 7)
#   plot_status <- "generated"
# }

info_lines <- c(
  "# Analysis Output",
  "",
  "## Run Info",
  sprintf("- run_time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("- params_path: %s", "params.json"),
  sprintf("- output_path: %s", output_dir),
  "",
  "## Params",
  sprintf("- deg: %s", deg_path),
  sprintf("- filter_by: %s", filter_by),
  sprintf("- cutoff: %.6g", cutoff),
  sprintf("- species: %s", species),
  sprintf("- gene_id_type: %s", gene_id_type),
  sprintf("- top_n: %d", top_n),
  "",
  "## Stats",
  sprintf("- input_row_count: %d", nrow(deg)),
  sprintf("- significant_feature_count: %d", nrow(filter_deg)),
  sprintf("- unique_feature_count: %d", length(feature_list)),
  sprintf("- id_conversion_status: %s", conversion_status),
  sprintf("- id_conversion_package: %s", conversion_package),
  sprintf("- id_conversion_input_count: %d", conversion_input_count),
  sprintf("- id_conversion_mapped_count: %d", conversion_mapped_count),
  sprintf("- id_conversion_unmapped_count: %d", conversion_unmapped_count),
  sprintf("- enrichment_terms_count: %d", nrow(enrich_df)),
  sprintf("- enrichment_key_type: %s", used_key_type),
  sprintf("- plot_status: %s", plot_status),
  "",
  "## Outputs",
  sprintf("- filtered_deg_file: %s", filtered_deg_file),
  sprintf("- id_mapping_file: %s", id_map_file),
  sprintf("- enrichment_file: %s", enrich_file),
  sprintf("- enrichment_plot: %s", plot_file)
)

readr::write_lines(info_lines, file.path(output_dir, "output.md"))
message(sprintf("Output saved to: %s", output_dir))











