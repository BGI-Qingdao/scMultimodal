# Title     : label transfer
# Objective : Transfer cell type labels from scRNA data to scATAC data
# Created by: Yao Li
# Created on: 2024-03-04
# Identifying anchors between scRNA-seq and scATAC-seq datasets to Transfer Label

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
parser$add_argument("--list", required = FALSE, help = "library id list text file")
# Parse the command-line arguments
args <- parser$parse_args()


# -------- 0. Load in data
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

# -------- 1. Quantify gene activity
print('Calculating Gene Activity...')
GeneActivity2 <- function(sample_ids, data.atac){
  feature_names <- list()
  for (sidx in 1:length(sample_ids)){
    sid <- sample_ids[[sidx]]
    subdata <- subset(data.atac, subset=sample==sid)
    gene.activities <- GeneActivity(subdata, features = VariableFeatures(data.rna))  # only consider HVGs according to scRNA data
    print(paste0(sid, ': ',length(rownames(gene.activities))))
    feature_names[[length(feature_names)+1]] <- rownames(gene.activities)  # append list into a list
  }
  print(length(feature_names))
  shared_feature_names <- Reduce(intersect,feature_names)
  print(paste0('Number of used features: ',length(shared_feature_names)))
  return(shared_feature_names)
}
getSampleID <- function(file_name){
  df <- read.table(file = file_name, header = F, stringsAsFactors = F)
  job_list <- as.vector(t(df$V1))
  return(job_list)
}
sample_ids <- c('DP8480004854TR_L01_16','DP8480004853TR_L01_16','DP8480004854TR_L01_1','DP8480004853TR_L01_1')
feature_names <- list()
for (sidx in 1:length(sample_ids)){
  sid <- sample_ids[[sidx]]
  subdata <- subset(data.atac, subset=sample==sid)
  gene.activities <- GeneActivity(subdata, features = VariableFeatures(data.rna))  # only consider HVGs according to scRNA data
  print(paste0(sid, ': ',length(rownames(gene.activities))))
  feature_names[[length(feature_names)+1]] <- rownames(gene.activities)  # append list into a list
}
print(length(feature_names))
shared_feature_names <- Reduce(intersect,feature_names)
print(paste0('Number of used features: ',length(shared_feature_names)))
gene.activities <- GeneActivity(data.atac, features = shared_feature_names)
# gene.activities <- GeneActivity(data.atac, features = VariableFeatures(data.rna))
print('Gene Activity Done')

# -------- 2. Add gene activities as a new assay
data.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# -------- 3. Filter out cells with low ATAC counts
# data.atac <- subset(data.atac, subset = nCount_ATAC > 500)  # set number

# -------- 3. Normalize gene activities
DefaultAssay(data.atac) <- "ACTIVITY"
print('Normalizing and Scaling data...')
data.atac <- NormalizeData(data.atac)
# 2024-06-03
# data.atac <- readRDS('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/scripts/scmulti/data.atac.small.rds')
data.atac <- ScaleData(data.atac, features = rownames(data.atac))
# data.atac <- ScaleData(data.atac, features = VariableFeatures(data.rna))
# saveRDS(data.atac, '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/scripts/scmulti/data.atac.small_scaled.rds')
gc()

# -------- 4. Identify anchors
print('Finding Transfer Anchors...')
transfer.anchors <- FindTransferAnchors(reference = data.rna,
                                        query = data.atac,
                                        features = rownames(data.atac),  #VariableFeatures(object = data.rna),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "rpca")  # 2024-06-04, cca
print('-----Finished finding transfer anchors-----')
saveRDS(transfer.anchors, '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/scripts/scmulti/transfer.anchors.rds')
gc()

# -------- 5. Annotate scATAC-seq cells via label transfer
print('Transfering data...')
celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = data.rna@meta.data[[args$celltype]],
                                     weight.reduction = data.atac[["integrated_lsi"]],
                                     dims = 2:30)
print('-----Finished transferring data-----')
data.atac <- AddMetaData(data.atac, metadata = celltype.predictions)

saveRDS(data.atac, file.path(args$output, paste0(args$name, '.merged.atac.rds')))
print('-----Saved results to disk-----')


# -------- 6. Visualization
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
