library(Seurat)
library(Signac)
library(Matrix)
library(argparse)

# Create an argument parser object
parser <- argparse::ArgumentParser()

# Add the input directory argument
parser$add_argument("--rna", help = "scRNA rds file")
parser$add_argument("--peaks", required = TRUE, help = "Path to the directory containing peak.bed, barcodes.tsv and matrix.mtx files")
parser$add_argument("-m", "--metadata_fn", help = "Path to ATAC metadata.tsv file")
parser$add_argument("-f", "--fragments", help = "Path to ATAC fragments.tsv.gz file")
parser$add_argument("-g", "--gff", help = "Path to GFF gene annotation file")
parser$add_argument("--gene_name_col", help = "gene name columns in to GFF gene annotation file")
parser$add_argument("--name", help = "Project name")
# Parse the command-line arguments
args <- parser$parse_args()

# load scATAC data
load_peak_matrix <- function(directory_path){
    barcode_fn <- file.path(directory_path, "barcodes.tsv")
    peaks_fn <- file.path(directory_path, "peak.bed")
    matrix_fn<- file.path(directory_path, "matrix.mtx")

    # Load the data files
    barcodes <- read.table(barcode_fn, sep = "\t", header = FALSE)
    peaks <- read.table(peaks_fn, sep = "\t", check.names=FALSE)
    frag_matrix <- readMM(matrix_fn)

    # Remove the extra row from the matrix
    #frag_matrix <- frag_matrix[-nrow(frag_matrix), ]
    # Set the row names for the matrix from the peaks object
    rownames(frag_matrix) <- peaks$V1

    # Subset the barcodes to match the cells in the matrix
    #barcodes <- barcodes[barcodes$V1 %in% colnames(frag_matrix), , drop = FALSE]
    # Set the column names for the matrix
    colnames(frag_matrix) <- barcodes$V1
    
    return (frag_matrix)
}
# Create chromain assay
peak_dir <- args$peaks
peak_matrix <- load_peak_matrix(peak_dir)
print(dim(peak_matrix))
metadata <- read.table(args$metadata, sep='\t',header =TRUE)
frag_fn <- args$fragments

chrom_assay <- CreateChromatinAssay(
  counts = peak_matrix,
  sep = c(":", "-","_"),
  fragments = frag_fn,
  min.cells = 1,
  min.features = 2
)
# Create ATAC Seurat Object
data.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
# [optional] add other info
data.atac$tech <- 'atac'

# preprocessing
data.atac <- FindVariableFeatures(data.atac)#, assay = "ACTIVITY", selection.method = "vst", nfeatures = 2000)
data.atac <- NormalizeData(data.atac)
data.atac <- ScaleData(data.atac)

# ATAC analysis add gene annotation information
gff_fn <- args$gff
read_gff <- function(gff_fn, custom_column = NULL){
    # read file
    gff = rtracklayer::import.gff(gff_fn)

    # ensure 'gene_name', 'gene_biotype' are in the metadata columns
    column_names = colnames(GenomicRanges::mcols(gff))

    # only select gene and protein coding genes
    if ('type' %in% column_names){
        gff = gff[gff$type=='gene']
    }
    if ('gene_biotype' %in% column_names){
        gff = gff[gff$gene_biotype=='protein_coding']
    }else{
        gff$gene_biotype <- "protein_coding"
    }

    # Check if 'gene_name' is missing and update it based on other columns
    if (!'gene_name'  %in% column_names){
      if ("gene_id" %in% column_names) {
        gff$gene_name <- gff$gene_id
      } else if ("gene" %in% column_names) {
        gff$gene_name <- gff$gene
      } else if (!is.null(custom_column)) {
        if (custom_column %in% column_names) {
          gff$gene_name <- gff$custom_column
        } else {
          stop("Custom column '", custom_column, "' is not found in the metadata.")
        }
      } else {
        stop("Unable to assign gene_name. No suitable metadata column found.")
      }
    }
    return (gff)
}

gene_gff <- read_gff(gff_fn, custom_column=args$gene_name_col)
Annotation(data.atac) <- gene_gff

# Dimension Reduction of the peak matrix via Latent Semantic Indexing
# To study the internal structure of the scATAC-seq data, making anchor detection possible
VariableFeatures(data.atac) <- names(which(Matrix::rowSums(data.atac) > 100))
# We exclude the first dimension as this is typically correlated with sequencing depth
data.atac <- RunTFIDF(data.atac)
data.atac <- FindTopFeatures(data.atac, min.cutoff = "q0")
data.atac <- RunSVD(data.atac)
data.atac <- RunUMAP(data.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


# Load scRNA data
scrna_fn <- args$rna
data.rna <- readRDS(scrna_fn)

# Identifying anchors between scRNA-seq and scATAC-seq datasets
# transfer label
# quantify gene activity
gene.activities <- GeneActivity(data.atac, features = VariableFeatures(data.rna))

# add gene activities as a new assay
data.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(data.atac) <- "ACTIVITY"
data.atac <- NormalizeData(data.atac)
data.atac <- ScaleData(data.atac, features = rownames(data.atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = data.rna, query = data.atac, features = VariableFeatures(object = data.rna), reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
# data.atac = readRDS('mid.atac.rds')
print('-----Finished finding transfer anchors-----')

# Annotate scATAC-seq cells via label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = data.rna@meta.data$celltype, weight.reduction = data.atac[["lsi"]], dims = 2:30)
print('-----Finished transferring data-----')
data.atac <- AddMetaData(data.atac, metadata = celltype.predictions)

saveRDS(data.atac, paste0(args$name, '.merged.atac.rds'))
print('-----Saved results to disk-----')

# library(SeuratDisk)
# SaveH5Seurat(data.atac, filename = paste0(args$name, ".h5Seurat"))
# Convert(paste0(args$name, ".h5Seurat"), dest = "h5ad")
