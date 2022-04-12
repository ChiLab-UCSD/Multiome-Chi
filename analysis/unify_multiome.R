library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 32000 * 1024^2)

file_list<-read.table('~/multiome_dir.txt')$V1

results_dir<-'/data/results/'

for (f in file_list)
{
	run_name<-substr(f,1,nchar(f)-1)
	print(run_name)
	load(paste0(results_dir,run_name,'.Rdata'))
	rm(list=c("peaks","annotation","macs2_counts","counts"))
	pbmc$dataset<-run_name
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
	if (!exists("merged_obj"))
	{
		merged_obj<-pbmc
		rm(pbmc)
	}else{
		merged_obj<-merge(x=merged_obj,y=pbmc)
		rm(pbmc)
	}

}

merged_obj