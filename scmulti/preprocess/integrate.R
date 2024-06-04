suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomicRanges))
suppressMessages(library(future))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(argparse))
suppressMessages(library(data.table))


# Create an argument parser object
parser <- argparse::ArgumentParser()
# Add the input directory argument
parser$add_argument("--list", help = "library id list text file")
parser$add_argument("--data_path", required = TRUE, help = "Path to the directory containing peak.bed, barcodes.tsv and matrix.mtx files")
parser$add_argument("-n", "--name", help = "orgam/sample name")
# Parse the command-line arguments
args <- parser$parse_args()

print('------------------------------------------------------------------------------------------------------------')
print(args$name)
df <- read.table(file = args$list, header = F, stringsAsFactors = F)
job_list <- as.vector(t(df$V1))
print(job_list)

#### 1.create muti datasets combined peaks####
# read all peak sets 
datapath <- args$data_path
# Create an empty list to store the data
data_list <- list()
gr_list <- c()
for (job in job_list){
  # read peak
  peak <- read.table(file.path(datapath,job,"/Peak_matrix/peak.bed"),
                     col.names = c("chr", "start", "end"),
                     sep = c("-","-"))

  #convert to genomic ranges
  gr <- makeGRangesFromDataFrame(peak)

  # Append the line data to the list
  data_list <- c(data_list, peak)
  gr_list <- c(gr_list, gr)
}

# Create a unified set of peaks to quantify in each dataset
code <- paste("c(", paste0("gr_list[[", 1:length(job_list), "]]", collapse = ","), ")")
combined.peaks <- reduce(eval(parse(text = code)))
print(combined.peaks)

#annotation obj
read_gff <- function(gff_fn, choose_type=NULL, custom_column = NULL){
    # read file
    gff <- rtracklayer::import.gff(gff_fn)
    # ensure 'gene_name', 'gene_id', 'gene_biotype' are in the metadata columns
    column_names <- colnames(GenomicRanges::mcols(gff))
    # only select gene and protein coding genes
    if (!is.null(choose_type)){
        if ('type' %in% column_names){
            gff <- gff[gff$type == choose_type]
        }
    }
    if ('gene_biotype' %in% column_names){
        gff <- gff[gff$gene_biotype=='protein_coding']
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
    # Check if 'gene_id' is missing and update it based on other columns
    if (!'gene_id' %in% column_names){
      if ("gene_name" %in% column_names) {
        gff$gene_id <- gff$gene_name
      } else {
        print("Unable to assign gene_id. No suitable metadata column found.")
      }
    }
    return (gff)
}
print('Reading GTF/GFF gene annotation...')
gff_fn <- '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/03.Oryzias_melastigma/GCF_002922805.2_ASM292280v2_genomic.gtf'
gff <- read_gff(gff_fn, choose_type='gene')

seurat_obj_list <- list()
for (job in job_list){
  # read meta
  print('Reading fragments metadata...')
  md <- read.table(file.path(datapath,job,paste0(job,'.Metadata.tsv')),
                    stringsAsFactors = FALSE,
                    header = TRUE,
                    row.names = 1)

  # create fragment objects
  print('create fragment objects...')
  frags <- CreateFragmentObject(
    path = file.path(datapath,job, paste0(job, ".fragments.tsv.gz")),
    cells = rownames(md)
  )

  # 3.Quantify peaks in each dataset
  print('Create Feature Matrix...')
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    cells = rownames(md)
  )

  # Create ChromatinAssay
  print('Creating ChromatinAssay...')
  chromatin_assay <- CreateChromatinAssay(counts, fragments = frags)
  obj <- CreateSeuratObject(chromatin_assay,
                            assay = "ATAC",
                            meta.data=md,
                            min.cells = 1,
                            min.features = 2)
  Annotation(obj) <- gff
  obj$sample <- job

  #compute LSI
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = 10)
  obj <- RunSVD(obj)
  # save_path <- file.path(datapath, 'individual_rds')
  saveRDS(obj, file.path(datapath, paste0(job, "_atac.rds")))

  # Append the line data to the list
  seurat_obj_list <- c(seurat_obj_list, obj)
}

###merge obj####
# seurat_obj_list <- list()
# for (job in job_list){
#   d <- readRDS(file.path(datapath, paste0(job, '_atac.rds')))
#   print(dim(d))
#   seurat_obj_list <- c(seurat_obj_list, d)
# }
# print(seurat_obj_list[13])
# seurat_obj_list <- seurat_obj_list[-13]

for (d in seurat_obj_list){
  print(dim(d))
}
print('Merging Seurat Objects...')
combined <- merge(
  x = seurat_obj_list[[1]],
  y = seurat_obj_list[-1]
)

# combined <- readRDS('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/03.Oryzias_melastigma/01.merge_libraries/fin/combinded2.rds')
print('Processing merged Seurat Object...')
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
saveRDS(combined,file.path(datapath, "combinded.rds"))
# combined <- readRDS('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/03.Oryzias_melastigma/01.merge_libraries/fin/combinded3.rds')


####integrate atac####
# Create an empty list to store the VariableFeatures results
variable_features_list <- list()
# Apply VariableFeatures to each element in seurat_obj_list
variable_features_list <- lapply(seurat_obj_list, VariableFeatures)
# Pass the list of VariableFeatures results to Reduce
print('Finding intersecting peaks...')
peaks.use <- Reduce(intersect, variable_features_list)
saveRDS(peaks.use,file.path(datapath, "peaks.use.rds"))

print('Find Integration Anchors...')
integration.anchors <- FindIntegrationAnchors(
  object.list = seurat_obj_list,
  anchor.features = peaks.use,
  reduction = "rlsi",
  dims = 2:50
)
saveRDS(integration.anchors,file.path(datapath, "integration.anchors.rds"))

print('Integrate Embeddings...')
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)
saveRDS(integrated,file.path(datapath, "integrated.rds"))

print('Making UMAP plots...')
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p1 <- DimPlot(combined, group.by = "sample") + ggtitle("Merged")
p2 <- DimPlot(integrated, group.by = "sample") + ggtitle("Integrated")
(p1|p2)
ggsave(file.path(datapath, "integrate.pdf"),width = 16,height = 7)

print('DONE.')
print('------------------------------------------------------------------------------------------------------------')

# integrated <- readRDS(file.path(datapath, "integrated.rds"))
# print('Processing integrated data...')
# integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
# integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3)
#
# print('Creating gene activity matrix...')
# gene.activities <- GeneActivity(integrated)
# print('Creating assay object...')
# integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
# integrated <- NormalizeData(
#   object = integrated,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(integrated$nCount_RNA)
# )
#
# saveRDS(integrated, file.path(datapath, "integrated_atac.rds"))
# print('------------------------------------------------------------------------------------------------------------')
