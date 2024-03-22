# Title     : TODO
# Objective : TODO
# Created by: Oreo
# Created on: 2024-03-04

library(Seurat)
library(Signac)
library(Matrix)
library(argparse)

# Create an argument parser object
parser <- argparse::ArgumentParser()

# Add the input directory argument
parser$add_argument("--rna", help = "scRNA rds file")
parser$add_argument("--atac", required = TRUE, help = "Path to the atac rds file")
parser$add_argument("-m", "--metadata_fn", help = "Path to ATAC metadata.tsv file")
parser$add_argument("-f", "--fragments", help = "Path to ATAC fragments.tsv.gz file")
parser$add_argument("-g", "--gff", help = "Path to GFF/GTF gene annotation file")
parser$add_argument("--gene_name_col", help = "gene name columns in to GFF/GTF gene annotation file")
parser$add_argument("--name", help = "Project name")
parser$add_argument("--celltype", default='celltype', help = "cell type label in meta.data")
parser$add_argument("-o", "--output", help = "output directory")
# Parse the command-line arguments
args <- parser$parse_args()

# Load scATAC data
print('Loading scATAC data...')
data.atac <- readRDS(args$atac)
print(data.atac)
# Load scRNA data
print('Loading scRNA data...')
scrna_fn <- args$rna
data.rna <- readRDS(scrna_fn)
print(str(data.rna))
print(colnames(data.rna@meta.data))
print(head(data.rna@meta.data[[args$celltype]]))

# Identifying anchors between scRNA-seq and scATAC-seq datasets
# transfer label
# quantify gene activity
print('Calculating Gene Activity...')
gene.activities <- GeneActivity(data.atac, features = VariableFeatures(data.rna))
print('Gene Activity Done')

# add gene activities as a new assay
data.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(data.atac) <- "ACTIVITY"
print('Normalizing and Scaling data...')
data.atac <- NormalizeData(data.atac)
data.atac <- ScaleData(data.atac, features = rownames(data.atac))

# Identify anchors
print('Finding Transfer Anchors...')
transfer.anchors <- FindTransferAnchors(reference = data.rna,
                                        query = data.atac,
                                        features = VariableFeatures(object = data.rna),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca")
print('-----Finished finding transfer anchors-----')

# Annotate scATAC-seq cells via label transfer
print('Transfering data...')
celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = data.rna@meta.data[[args$celltype]],
                                     weight.reduction = data.atac[["integrated_lsi"]],
                                     dims = 2:30)
print('-----Finished transferring data-----')
data.atac <- AddMetaData(data.atac, metadata = celltype.predictions)

saveRDS(data.atac, file.path(args$output, paste0(args$name, '.merged.atac.rds')))
print('-----Saved results to disk-----')


plot_annotation <- function(merged.atac, data.rna, celltype_label){
  p1 <- DimPlot(merged.atac, reduction = 'integrated_lsi', group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
  p2 <- DimPlot(data.rna, group.by = celltype_label, label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
  p1 + p2
  ggsave(filename = file.path(args$output, "atac_predicted_id.png"))
}

plot_umap <- function(data.rna, merged.atac, celltype_label){
  p1 <- DimPlot(data.rna, group.by = celltype_label, label = TRUE) + NoLegend() + ggtitle("RNA")
  p2 <- DimPlot(merged.atac, reduction = 'integrated_lsi', group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
  p1 + p2
  plot <- (p1 + p2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
  ggsave(filename = file.path(args$output, "atacseq_b4_integration.png"))#, height = 7, width = 12, plot = plot, quality = 50)
}

plot_annotation(data.atac, data.rna, args$celltype)
plot_umap(data.rna, data.atac, args$celltype)
