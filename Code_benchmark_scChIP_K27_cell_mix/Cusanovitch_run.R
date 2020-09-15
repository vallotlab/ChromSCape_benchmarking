#Cusanovitch 2018 testing
library(data.table)
library(Matrix)
library(proxy)
library(Rtsne)
library(densityClust)
library(data.table)
library(irlba)
library(umap)
library(ggplot2)

elbow_plot <- function(mat,num_pcs=50,scale=FALSE,center=FALSE,title='',width=3,height=3){
    set.seed(2019) 
    mat = data.matrix(mat)
    SVD = irlba(mat, num_pcs, num_pcs,scale=scale,center=center)
    options(repr.plot.width=width, repr.plot.height=height)
    df_plot = data.frame(PC=1:num_pcs, SD=SVD$d);
    #     print(SVD$d[1:num_pcs])
    p <- ggplot(df_plot, aes(x = PC, y = SD)) +
        geom_point(col="#cd5c5c",size = 1) + 
        ggtitle(title)
    return(p)
}

#Load peaks x cells matrix filtered by black regions & CNV, obtained by cisTopic
load(file.path("CisTopic","cisTopic_object.RData"))
names = cisTopicObject@cell.names
mat = cisTopicObject@count.matrix
colnames(mat) = names
dim(mat)
binary_mat = as.matrix((mat > 0) + 0)
binary_mat = Matrix(binary_mat, sparse = TRUE)

num_cells_ncounted = Matrix::rowSums(binary_mat)
ncounts = binary_mat[num_cells_ncounted >= dim(binary_mat)[2]*0.01,]
new_counts = Matrix::colSums(ncounts)
# ncounts = ncounts[,new_counts >= quantile(new_counts,probs=0.1)]
ncounts = ncounts[Matrix::rowSums(ncounts) > 0,]

png(file.path("Cusanovich2018","histogram_filtering.png"), width = 1000,height = 800)
par(mfrow=c(1,2))
hist(log10(num_cells_ncounted),main="No. of Cells Each Site is Observed In",breaks=50)
abline(v=log10(min(num_cells_ncounted[num_cells_ncounted >= dim(binary_mat)[2]*0.01])),lwd=2,col="indianred")
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
# abline(v=log10(quantile(new_counts,probs=0.1)),lwd=2,col="indianred")
dev.off()

sexsites = c(grep("chrY",rownames(ncounts)),grep("chrX",rownames(ncounts)))
ncounts.nosex = ncounts[-sexsites,]

nfreqs = t(t(ncounts.nosex) / Matrix::colSums(ncounts.nosex))
idf = as(log(1 + ncol(ncounts.nosex) / Matrix::rowSums(ncounts.nosex)), "sparseVector")
tf_idf_counts = as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs
dim(tf_idf_counts)

p_elbow_LSI <- elbow_plot(tf_idf_counts,num_pcs = 200, title = 'PCA on LSI')

png(file.path("Cusanovich2018","elbow_plot.png"), width = 800,height = 800)
p_elbow_LSI
dev.off()


set.seed(2019)
num_pcs = 6
SVDtsne = irlba(tf_idf_counts, num_pcs, num_pcs)
d_diagtsne = matrix(0, nrow=num_pcs, ncol=num_pcs)
diag(d_diagtsne) = SVDtsne$d
SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))

reduced_dim_Cusanovich2018 = t(SVDtsne_vd)
colnames(reduced_dim_Cusanovich2018) = names
rownames(reduced_dim_Cusanovich2018) = paste('PC',1:dim(SVDtsne_vd)[2])
dim(reduced_dim_Cusanovich2018)
save(reduced_dim_Cusanovich2018,file=file.path("Cusanovich2018","reduced_dim_Cusanovich2018.RData"))


