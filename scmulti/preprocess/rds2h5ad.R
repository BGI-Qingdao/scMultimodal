library(Seurat)
library(SeuratDisk)

args <- commandArgs()
fn <- args[6]
base_name <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(fn))
rdata <- readRDS(fn)
SaveH5Seurat(rdata, filename = paste0(base_name, ".h5Seurat"))
Convert(paste0(base_name, ".h5Seurat"), dest = "h5ad")
