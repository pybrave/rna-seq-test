library(tidyverse)
library(clusterProfiler)
library(ggtree)


params <- jsonlite::fromJSON("params.json")

enrich <- readRDS(params$rds)

barplot_p <- barplot(enrich, showCategory=20)
ggsave(str_glue("output/barplot.png"),barplot_p,  width = 10, height = 7,units = "in")
ggsave(str_glue("output/barplot.download.pdf"),barplot_p,  width = 10, height = 7,units = "in")

dotplot_p <- dotplot(enrich, showCategory=20)
ggsave(str_glue("output/dotplot.png"),dotplot_p,  width = 10, height = 7,units = "in")
ggsave(str_glue("output/dotplot.download.pdf"),dotplot_p,  width = 10, height = 7,units = "in")

# edox2 <- pairwise_termsim(enrich)
# p1 <- treeplot(edox2, cladelab_offset=8, tiplab_offset=.3, fontsize_cladelab =5) + 
#   hexpand(.2)
# p2 <- treeplot(edox2, cluster_method = "average", 
#                cladelab_offset=14, tiplab_offset=.3, fontsize_cladelab =5) + 
#   hexpand(.3)
# aplot::plot_list(p1, p2, tag_levels='A', ncol=2)