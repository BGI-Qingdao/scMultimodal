# Load the necessary R packages
library(Seurat)
library(Matrix)
library(argparse)

# Create an argument parser object
parser <- argparse::ArgumentParser()

# Add the input directory argument
parser$add_argument("-d", "--directory", help = "Path to the directory containing peak.bed, barcodes.tsv and matrix.mtx files")
parser$add_argument("-m", "--metadata_fn", help = "Path to ATAC metadata file")

# Parse the command-line arguments
args <- parser$parse_args()

# Access the directory path
directory_path <- args$directory
#directory_path <- '/dellfsqd2/C_OCEAN/USERS/c-zhangjin/07.project/5.liver_ATAC/0.qc/QC/lsqsm-youti-gan-2_2/output/Peak_matrix'

# Assuming the file names are fixed
barcode_fn <- file.path(directory_path, "barcodes.tsv")
peaks_fn <- file.path(directory_path, "peak.bed")
matrix_fn<- file.path(directory_path, "matrix.mtx")


# Load the data files
barcodes <- read.table(barcode_fn, sep = "\t", header = FALSE)
peaks <- read.table(peaks_fn, sep = "\t")
frag_matrix <- readMM(matrix_fn)

# Remove the extra row from the matrix
#frag_matrix <- frag_matrix[-nrow(frag_matrix), ]

# Set the row names for the matrix from the peaks object
rownames(frag_matrix) <- peaks$V1

# Subset the barcodes to match the cells in the matrix
#barcodes <- barcodes[barcodes$V1 %in% colnames(frag_matrix), , drop = FALSE]

# Set the column names for the matrix
colnames(frag_matrix) <- barcodes$V1

# Create a Seurat object
#metadata <- read.table('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/atac/QC/lsqsm-youti-gan-2_1/output/lsqsm-youti-gan-2_1.Metadata.tsv',sep='\t',header =TRUE)
metadata_fn <- args$metadata_fn
metadata <- read.table(metadata_fn,sep='\t',header =TRUE)

seuratObj <- CreateSeuratObject(counts = frag_matrix, meta.data = metadata, features = peaks, assay = "peaks")
print(seuratObj)
#seuratObj <- FindVariableFeatures(seuratObj, assay = "RNA", selection.method = "vst", nfeatures = 2000)

# Save the Seurat object to an RDS file
saveRDS(seuratObj, "atac.rds")

