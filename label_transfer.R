library(Seurat)
library(Signac)
library(Matrix)


# load scRNA data
fn <- '/dellfsqd1/ST_OCEAN/C_OCEAN/USERS/c-zhangjin/07.project/2.atlas_lampery/5.rds_v3/a_liver.rds'
data.rna <- readRDS(fn)

# load scATAC data
load_atac <- function(directory_path){
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

    # Create a Seurat object
    #seuratObj <- CreateSeuratObject(counts = frag_matrix, meta.data = barcodes, features = peaks)
    #return (seuratObj)
    return (frag_matrix)
}

dir2 <- '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/atac/QC/lsqsm-youti-gan-2_1/output/Peak_matrix/'
peaks <- load_atac(dir2)

# 2. chromain assay
metadata <- read.table('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/atac/QC/lsqsm-youti-gan-2_1/output/lsqsm-youti-gan-2_1.Metadata.tsv',sep='\t',header =TRUE)

frag_fn <- '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/atac/QC/lsqsm-youti-gan-2_1/output/lsqsm-youti-gan-2_1.fragments.tsv.gz'

chrom_assay <- CreateChromatinAssay(
  counts = peaks,
  sep = c(":", "-","_"),
  fragments = frag_fn,
  min.cells = 10,
  min.features = 200
)

data.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

data.atac$tech <- 'atac'

# preprocessing
data.atac <- FindVariableFeatures(data.atac)#, assay = "ACTIVITY", selection.method = "vst", nfeatures = 2000)
data.atac <- NormalizeData(data.atac)
data.atac <- ScaleData(data.atac)

# ATAC analysis add gene annotation information
gff_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/atac/LRfinal.gff3'
gff_fn = '/dellfsqd1/ST_OCEAN/C_OCEAN/USERS/c-zhangjin/07.project/5.liver_ATAC/0.qc/1.input_gem/LRfinal.gff3'
gff = rtracklayer::import.gff(gff_fn)
gene_gff = gff[gff$type=='gene']
gene_gff$gene_id <- gene_gff$ID
gene_gff$gene_biotype <- "protein_coding"
gene_gff$gene_name <- gene_gff$gene_id
Annotation(data.atac) <- gene_gff


# Dimension Reduction of the peak matrix via Latent Semantic Indexing
# To study the internal structure of the scATAC-seq data, making anchor detection possible
VariableFeatures(data.atac) <- names(which(Matrix::rowSums(data.atac) > 100))
# We exclude the first dimension as this is typically correlated with sequencing depth
data.atac <- RunTFIDF(data.atac)
data.atac <- FindTopFeatures(data.atac, min.cutoff = "q0")
data.atac <- RunSVD(data.atac)
data.atac <- RunUMAP(data.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# plotting
#p1 <- DimPlot(data.rna, group.by = "celltype", label = TRUE) + NoLegend() + ggtitle("RNA")
#p2 <- DimPlot(data.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
#p1 + p2
#plot <- (p1 + p2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
#ggsave(filename = "../output/images/atacseq_integration_vignette.png", height = 7, width = 12, plot = plot, quality = 50)
#


#new_cell_names = colnames(data.rna)[1:8045]
#new.data.atac <- RenameCells(data.atac, new.names = new_cell_names)

data.atac <- readRDS('data.atac.rds')

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
transfer.anchors <- FindTransferAnchors(reference = data.rna, query = data.atac, features = VariableFeatures(object = data.rna),reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

#saveRDS(data.atac, 'merged.atac.rds')
#data.atac = readRDS('merged.atac.rds')

# Annotate scATAC-seq cells via label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = data.rna$seurat_annotations, weight.reduction = data.atac[["lsi"]], dims = 2:30)
data.atac <- AddMetaData(data.atac, metadata = celltype.predictions)

# Do not check is annotation is correct
#data.atac$annotation_correct <- data.atac$predicted.id == data.atac$seurat_annotations

saveRDS(data.atac, 'merged.atac.rds')

