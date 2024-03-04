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
parser$add_argument("-o", "--organ", help = "orgam/sample name")
# Parse the command-line arguments
args <- parser$parse_args()

print('------------------------------------------------------------------------------------------------------------')
print(args$organ)
# job_list <- c('DP8480004852TR_L01_3',
#               'DP8480004851TR_L01_3',
#               'DP8480004852TR_L01_2',
#               'DP8480004851TR_L01_2')
# /dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/03.Oryzias_melastigma/sample_ids
df <- read.table(file = args$list, header = F, stringsAsFactors = F)
job_list <- as.vector(t(df$V1))
print(job_list)

#### 1.create muti datasets combined peaks####
# read all peak sets 
datapath <- file.path('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/03.Oryzias_melastigma/01.merge_libraries/', args$organ)
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
#会有warning说有的scaffold在其他样本没有很正常，reduce不仅会合并相邻的还会保留特异的，这些特异的在别的时期没有
#https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/IRanges/html/inter-range-methods.html
combined.peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]], gr_list[[3]], gr_list[[4]]))
# combined.peaks <- reduce(x = gr_list)
# combined.peaks <- do.call(reduce, gr_list)
#peakwidths <- width(combined.peaks)
print(combined.peaks)

#annotation obj
read_gff <- function(gff_fn, choose_type=NULL, custom_column = NULL){
    # read file
    gff = rtracklayer::import.gff(gff_fn)

    # ensure 'gene_name', 'gene_id', 'gene_biotype' are in the metadata columns
    column_names = colnames(GenomicRanges::mcols(gff))

    # only select gene and protein coding genes
    if (!is.null(choose_type)){
        if ('type' %in% column_names){
            gff <- gff[gff$type == choose_type]
        }
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
print(gff)

seurat_obj_list <- list()
for (job in job_list){
  # read meta
  print('Reading fragments metadata...')
  md <- read.table(file.path(datapath,job,paste0(job,'.Metadata.tsv')),
                    stringsAsFactors = FALSE,
                    header = TRUE,
                    row.names = 1)
  print(head(rownames(md)))
  # create fragment objects
  print('create fragment objects...')
  frags <- CreateFragmentObject(
    path = file.path(datapath,job, paste0(job, ".fragments.tsv.gz")),
    cells = rownames(md)
  )
  print(frags)

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

####merge obj####
print('Merging Seurat Objects...')
combined <- merge(
  x = seurat_obj_list[[1]],
  y = seurat_obj_list[-1]
)

print('Processing merged Seurat Object...')
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
saveRDS(combined,file.path(datapath, "combined.rds"))

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
  dims = 2:30
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

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p1 <- DimPlot(combined, group.by = "sample") + ggtitle("Merged")
p2 <- DimPlot(integrated, group.by = "sample") + ggtitle("Integrated")
(p1|p2)

ggsave(file.path(datapath, "integrate.pdf"),width = 16,height = 7)

integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3)


gene.activities <- GeneActivity(integrated)

integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)

saveRDS(integrated,file.path(datapath, "integrated_atac.rds"))
print('------------------------------------------------------------------------------------------------------------')
