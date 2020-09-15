#SnapATAC 
library(SnapATAC)
library(GenomicRanges)
library(tidyverse)
library(RColorBrewer)
library(colorRamps)
library(colorspace)

## From Snap files, uncomment and run lines below - or skip to loading snapObject
# files = list.files(file.path("SnapATAC"),full.names = T, pattern = "*.snap$")
# samples = gsub(".snap","",basename(files))
# 
# x.sp = createSnap(file=files, sample = samples, do.par = T, num.cores = 6)
# x.sp = addBmatToSnap(x.sp, bin.size=50000, do.par = T, num.cores=6);

load("SnapATAC/SnapATAC_object.RData")

x.sp@barcode = paste0(x.sp@sample,"_",x.sp@barcode)
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(Bin Cov)", 
    col="lightblue", 
    breaks = 250
);

cells_passing_QC = as.character(read.table("cells_passing_QC.txt")[,1])

x.sp = x.sp[which(x.sp@barcode %in% cells_passing_QC),]

x.sp = makeBinary(x.sp, mat="bmat");

black_list = read.table("Additional_Files/hg38-blacklist.v2.bed",sep = "\t")
CNV_regions = read.table("Additional_Files/MM468_identified_CNA.bed");

exclude_regions = rbind(black_list[,1:3],CNV_regions)

black_list = exclude_regions
black_list.gr = GRanges(
    black_list[,1],
    IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};

idy2 = NULL
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
if(length(chr.exclude)>0 ) idy2 = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy2) > 0){x.sp = x.sp[,-idy2, mat="bmat"]};


options(repr.plot.width=4, repr.plot.height=4)
plotBinCoverage(
    x.sp,
    pdf.file.name=NULL,
    col="grey",
    border="grey",
    breaks=10,
    xlim=c(-6,6)
);

x.sp = filterBins(
    x.sp,
    low.threshold=-2,
    high.threshold=2,
    mat="bmat"
);

# filter out cells with less than 1000 reads
# low_covered_barcodes = x.sp@barcode[which(Matrix::rowSums(x.sp@bmat) <=1000)]
# length(low_covered_barcodes)
# x.sp = x.sp[which(!x.sp@barcode %in% low_covered_barcodes),]


x.sp = runDiffusionMaps(
    obj=x.sp,
    input.mat="bmat"
    ,num.eigs = 20
)

options(repr.plot.width=6, repr.plot.height=6)
plotDimReductElbow(
    obj=x.sp, 
    point.size=1.5,
    point.shape=19,
    point.color="red",
    point.alpha=1,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7,
    labs.title="PCA Elbow plot",
    labs.subtitle=NULL
);

plotDimReductPW(
    obj=x.sp, 
    eigs.dims=1:20,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
);

num.eigs = 6

x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:num.eigs,
    k=4
)

x.sp=runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    resolution = 1000,
    seed.use=10
)

x.sp@metaData$cluster = x.sp@cluster
x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:num.eigs, 
    method="Rtsne",
    seed.use=10
)

x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:num.eigs, 
    method="umap",
    seed.use=10
)

snap_reduced_dims = x.sp@smat@dmat
rownames(snap_reduced_dims) = x.sp@barcode
colnames(snap_reduced_dims) = paste0("PC",1:20)
save(snap_reduced_dims,file=file.path("SnapATAC","reduced_dim_SnapATAC.RData"))
