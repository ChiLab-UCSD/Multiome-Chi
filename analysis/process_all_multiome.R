# following https://satijalab.org/signac/articles/pbmc_multiomic.html

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

base_path<-'/data-2/'

file_list<-read.table('~/multiome_dir.txt')$V1

results_dir<-'/data/results/'

for (f in file_list)
{

# load the RNA and ATAC data

base_path<-paste0('/data-2/',f,'epigenomics.sdsc.edu/hjiao/bot_output/Neil_Chi_152161/counts/',f,'outs/')

counts <- Read10X_h5(paste0(base_path,f,"filtered_feature_bc_matrix.h5"))

fragpath <- paste0(base_path,"atac_fragments.tsv.gz")

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

run_name<-substr(f,1,nchar(f))

pdf(paste0(results_dir,run_name,'.pdf')

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

dev.off()

macs_path<-'/home/ubuntu/miniconda3/envs/Signac/bin/macs2'

peaks <- CallPeaks(pbmc, macs2.path = macs_path)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

save.image(paste0(results_dir,run_name,'.Rdata')

}