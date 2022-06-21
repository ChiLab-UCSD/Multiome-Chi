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

file_list<-read.table('~/multiome_dir_4.txt')$V1

results_dir<-'/data/results/'

#JB_745_1_2_JB_742_1_2 causing memory problems

for (f in file_list)
{

# load the RNA and ATAC data

  run_name<-substr(f,1,nchar(f)-1)
  print(run_name)
  load(paste0(results_dir,run_name,'.Rdata'))
  rm(list=c("peaks","annotation","macs2_counts","counts"))
  pbmc$dataset<-run_name

DefaultAssay(pbmc) <- "ATAC"

#pbmc <- NucleosomeSignal(pbmc)
#pbmc <- TSSEnrichment(pbmc)

run_name<-substr(f,1,nchar(f)-1)

#QC

pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

#clean up and save

save.image(paste0(results_dir,'qc/',run_name,'.Rdata'))
rm(pbmc)
gc(TRUE)
}