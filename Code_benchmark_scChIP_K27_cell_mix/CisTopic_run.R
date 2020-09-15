#CisTopic testing
library(cisTopic)
library(GenomicRanges)
library(tidyverse)

metadata = read.csv(file.path("Additional_Files/metadata.csv"),row.names = 1)

## From single-cell BAMs, uncomment and run lines below - or skip to loading cisTopicObject
# pathToBams <- file.path("CisTopic","single_cell_bams")
# bamFiles <- file.path(pathToBams, list.files(pathToBams,pattern = "*.bam$"))
# # if(length(grep("merged.bam", bamFiles))>0 ) bamFiles = bamFiles[-grep("merged.bam", bamFiles)]
# if(!dir.exists(file.path("CisTopic"))) dir.create(file.path("CisTopic"))
# regions <- file.path("Additional_Files/merged_peaks.5000.bed")
# 
# system.time({
#   cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions,
#                                                 project.name='Jurkat_Ramos_HBCx22_MM468'
#   )
# })

# cisTopicObject@cell.names = gsub("X.media.pacome.LaCie.InstitutCurie.Documents.Data.Tests.Benchmark_ChromSCape.Jurkat_Ramos_HBCx22_MM4681.CisTopic.single_cell_bams.","",cisTopicObject@cell.names)
# cisTopicObject@cell.names = gsub(".bam","",cisTopicObject@cell.names)
# colnames(cisTopicObject@count.matrix) = gsub("X.media.pacome.LaCie.InstitutCurie.Documents.Data.Tests.Benchmark_ChromSCape.Jurkat_Ramos_HBCx22_MM4681.CisTopic.single_cell_bams.","",
#                                              colnames(cisTopicObject@count.matrix))
# colnames(cisTopicObject@count.matrix) = gsub(".bam","",colnames(cisTopicObject@count.matrix)) 
# names = data.frame(barcode = gsub(".*_","",cisTopicObject@cell.names), ID =cisTopicObject@cell.names,
#                    label =  gsub("_BC.*","",cisTopicObject@cell.names))
# rownames(names)= names$ID
# rownames(cisTopicObject@cell.data) = cisTopicObject@cell.names
# cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = names)
# 
# cells_passing_QC = as.character(read.table("cells_passing_QC.txt")[,1])
# 
# # keep selected cells_passing_QC
# cisTopicObject@count.matrix = cisTopicObject@count.matrix[,which(cisTopicObject@cell.names %in% cells_passing_QC)]
# cisTopicObject@cell.data = cisTopicObject@cell.data[which(cisTopicObject@cell.names %in% cells_passing_QC),]
# cisTopicObject@cell.data$ID = factor(cisTopicObject@cell.data$ID)
# cisTopicObject@cell.data$barcode = factor(cisTopicObject@cell.data$barcode)
# cisTopicObject@binary.count.matrix =  cisTopicObject@binary.count.matrix[,which(cisTopicObject@cell.names %in% cells_passing_QC)]
# cisTopicObject@cell.names =  cisTopicObject@cell.names[which(cisTopicObject@cell.names %in% cells_passing_QC)]
# 
# CNV_regions = read.table("Additional_Files/MM468_identified_CNA.bed");
# 
# black_list.gr = GRanges(
#   CNV_regions[,1],
#   IRanges(CNV_regions[,2], CNV_regions[,3])
# );
# 
# idy = queryHits(findOverlaps(cisTopicObject@region.ranges, black_list.gr));
# if(length(idy) > 0){
#   cisTopicObject@count.matrix = cisTopicObject@count.matrix[-idy, ]
#   cisTopicObject@binary.count.matrix = cisTopicObject@binary.count.matrix[-idy, ]
#   cisTopicObject@region.names = cisTopicObject@region.names[-idy]
#   cisTopicObject@region.data = cisTopicObject@region.data[-idy,]
#   cisTopicObject@region.ranges = cisTopicObject@region.ranges[-idy,]
#   };
# 
# 
# system.time({
#   cisTopicObject <- cisTopic::runCGSModels( cisTopicObject, topic=c(10, 20, 25, 30, 35, 40), seed=987,
#                                             nCores=6, burnin = 120, iterations = 150, addModels=FALSE)
# })
# 
# cisTopicObject <- selectModel(cisTopicObject)
# logLikelihoodByIter(cisTopicObject, select=c(10, 20, 25, 30, 35, 40))

load("CisTopic/cisTopic_object.RData")
cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')

save(cellassign, file =file.path("CisTopic","reduced_dim_cisTopic.RData"))
