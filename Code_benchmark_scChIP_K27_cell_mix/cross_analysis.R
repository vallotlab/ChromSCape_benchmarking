# Cross analysis of Grosselin
library(tidyverse)
library(dplyr)
library(viridis)

rd = new.env()
# CisTopic
load(file.path("CisTopic","reduced_dim_cisTopic.RData"), rd)
rd_cisTopic = as.data.frame(t(rd$cellassign))
rd_cisTopic = rd_cisTopic[order(rownames(rd_cisTopic)),]

#Snap
load(file.path("SnapATAC","reduced_dim_SnapATAC.RData"), rd)
rd_snapATAC = as.data.frame(rd$snap_reduced_dims)
rd_snapATAC = rd_snapATAC[order(rownames(rd_snapATAC)),]

#ChromSCape
load(file.path("ChromSCape","reduced_dim_ChromSCape.RData"), rd)
rd_ChromSCape = as.data.frame(rd$reduced_dim_ChromSCape)
rd_ChromSCape = rd_ChromSCape[order(rownames(rd_ChromSCape)),]
#Cusanovich2018
load(file.path("Cusanovich2018","reduced_dim_Cusanovich2018.RData"), rd)
rd_Cusanovich2018 = as.data.frame(t(rd$reduced_dim_Cusanovich2018))
rd_Cusanovich2018 = rd_Cusanovich2018[order(rownames(rd_Cusanovich2018)),]

#epiScanpy
rd_epiScanpy = read.csv(file.path("EpiScanpy","reduced_dim_EpiScanpy.csv"),row.names = 1)
rd_epiScanpy = as.data.frame(t(rd_epiScanpy))
rd_epiScanpy = rd_epiScanpy[order(rownames(rd_epiScanpy)),]

affectation_cisTopic = affectation_ChromSCape = affectation_snapATAC = 
    affectation_Cusanovich2018 = affectation_epiScanpy = read.csv(
        "Additional_Files/metadata.csv",row.names = 1)

find_clusters_hc = function(reduced_dim, k = 4){
    dist_mat = as.dist(1-cor(t(reduced_dim)))
    hc_cor = stats::hclust(dist_mat, method = "ward.D")
    cell_clusters = stats::cutree(hc_cor, k = k)
    cell_clusters = paste0("C",cell_clusters)
    cell_clusters = as.factor(cell_clusters)
    return(cell_clusters)
}

affectation_cisTopic$cell_cluster = find_clusters_hc(rd_cisTopic)
affectation_ChromSCape$cell_cluster = find_clusters_hc(rd_ChromSCape)
affectation_snapATAC$cell_cluster = find_clusters_hc(rd_snapATAC)
affectation_Cusanovich2018$cell_cluster = find_clusters_hc(rd_Cusanovich2018)
affectation_epiScanpy$cell_cluster = find_clusters_hc(rd_epiScanpy)

table(affectation_cisTopic[,c("label","cell_cluster")])
table(affectation_snapATAC[,c("label","cell_cluster")])
table(affectation_ChromSCape[,c("label","cell_cluster")])
table(affectation_Cusanovich2018[,c("label","cell_cluster")])
table(affectation_epiScanpy[,c("label","cell_cluster")])

# Clusters are given arbitrary names (C1, C2, ...), this function tries to
# make names correspond to each other between methods
change_cluster_names = function(affectation,C1, C2, C3, C4){
    affectation$cell_cluster = as.character(affectation$cell_cluster)
    affectation$cell_cluster[which(affectation$cell_cluster=="C1")] = paste0(C1,"x")
    affectation$cell_cluster[which(affectation$cell_cluster=="C2")] = paste0(C2,"x")
    affectation$cell_cluster[which(affectation$cell_cluster=="C3")] = paste0(C3,"x")
    affectation$cell_cluster[which(affectation$cell_cluster=="C4")] = paste0(C4,"x")
    affectation$cell_cluster = as.factor(gsub("x","",affectation$cell_cluster))
    return(affectation)
}

# Manually assign each cluster to the closest cluster found in ChromSCape
table(affectation_snapATAC[,c("label","cell_cluster")])
affectation_snapATAC = change_cluster_names(affectation_snapATAC, "C1", "C4", "C3", "C2")
table(affectation_snapATAC[,c("label","cell_cluster")])

#### Plot UMAPs ####
set.seed(10)
library(umap)
# calculate UMAP
umap_cisTopic = umap(rd_cisTopic)
umap_cisTopic = as.data.frame(umap_cisTopic$layout)

umap_ChromSCape = umap(rd_ChromSCape)
umap_ChromSCape = as.data.frame(umap_ChromSCape$layout)

umap_snapATAC = umap(rd_snapATAC)
umap_snapATAC = as.data.frame(umap_snapATAC$layout)

umap_Cusanovich2018 = umap(rd_Cusanovich2018)
umap_Cusanovich2018 = as.data.frame(umap_Cusanovich2018$layout)

umap_epiScanpy = umap(rd_epiScanpy)
umap_epiScanpy = as.data.frame(umap_epiScanpy$layout)

plot_umap = function(umap, affectation, by = "cell_cluster", colors = c("#4285F4", "#F4B400", "#DB4437", "#0F9D58")
                     , leg.pos="none"){
    colnames(umap)[1:2] = paste0("Component_",1:2)
    umap %>% ggplot(aes(x=Component_1, y = Component_2, color = affectation[,by] )) +
        geom_point(size=7.5, alpha = 0.65) + theme_minimal() + scale_color_manual(values=colors) +
        theme(panel.grid = element_blank(), legend.position = leg.pos,text = element_blank()) 
}

# Print png
png(file.path("Cross_analysis","umap_cluster_ChromSCape_hc.png"), width = 800,height = 800)
plot_umap(umap_ChromSCape,affectation_ChromSCape)
dev.off()

png(file.path("Cross_analysis","umap_cluster_cisTopic_hc.png"), width = 800,height = 800)
plot_umap(umap_cisTopic,affectation_cisTopic)
dev.off()

png(file.path("Cross_analysis","umap_cluster_snapATAC_hc.png"), width = 800,height = 800)
plot_umap(umap_snapATAC,affectation_snapATAC)
dev.off()

png(file.path("Cross_analysis","umap_cluster_Cusanovich2018_hc.png"), width = 800,height = 800)
plot_umap(umap_Cusanovich2018,affectation_Cusanovich2018)
dev.off()

png(file.path("Cross_analysis","umap_cluster_epiScanpy_hc.png"), width = 800,height = 800)
plot_umap(umap_epiScanpy,affectation_epiScanpy)
dev.off()

# Print by label
png(file.path("Cross_analysis","umap_cluster_ChromSCape_label.png"), width = 800,height = 800)
plot_umap(umap_ChromSCape,affectation_ChromSCape,colors = c("#AE847E","#2C0E37","#690375","#CB429F"), by ="label")
dev.off()

png(file.path("Cross_analysis","umap_cluster_cisTopic_label.png"), width = 800,height = 800)
plot_umap(umap_cisTopic,affectation_cisTopic,colors = c("#AE847E","#2C0E37","#690375","#CB429F"), by ="label")
dev.off()

png(file.path("Cross_analysis","umap_cluster_snapATAC_label.png"), width = 800,height = 800)
plot_umap(umap_snapATAC,affectation_snapATAC,colors = c("#AE847E","#2C0E37","#690375","#CB429F"), by ="label")
dev.off()

png(file.path("Cross_analysis","umap_cluster_Cusanovich2018_label.png"), width = 800,height = 800)
plot_umap(umap_Cusanovich2018,affectation_Cusanovich2018,colors = c("#AE847E","#2C0E37","#690375","#CB429F"), by ="label")
dev.off()

png(file.path("Cross_analysis","umap_cluster_epiScanpy_label.png"), width = 800,height = 800)
plot_umap(umap_epiScanpy,affectation_epiScanpy,colors = c("#AE847E","#2C0E37","#690375","#CB429F"), by ="label")
dev.off()

# Adjusted Rand Index :
library(mclust)
true_clusters = as.numeric(as.factor(affectation_ChromSCape$label))
ari_ChromSCape = mclust::adjustedRandIndex(as.numeric(gsub("C","",affectation_ChromSCape$cell_cluster)),true_clusters)

ari_cisTopic = mclust::adjustedRandIndex(as.numeric(gsub("C","",affectation_cisTopic$cell_cluster)),true_clusters)

ari_Cusanovich2018 = mclust::adjustedRandIndex(as.numeric(gsub("C","",affectation_Cusanovich2018$cell_cluster)),true_clusters)

ari_snapATAC = mclust::adjustedRandIndex(as.numeric(gsub("C","",affectation_snapATAC$cell_cluster)),true_clusters)

ari_epiScanpy = mclust::adjustedRandIndex(as.numeric(gsub("C","",affectation_epiScanpy$cell_cluster)),true_clusters)

ARIs = c(ari_ChromSCape,ari_cisTopic,ari_Cusanovich2018,ari_snapATAC,ari_epiScanpy)
ggplot(data.frame(method =c("ChromSCape","cisTopic","Cusanovich2018","snapATAC","epiScanpy"), 
                  ARI = ARIs)) + geom_bar(aes(x=method,y=ARI), stat="identity",
                                          col="dark red", fill = "#D97604",width = 0.5) + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) + xlab("")
names(ARIs) =  c("ari_ChromSCape","ari_cisTopic","ari_Cusanovich2018","ari_snapATAC","ari_epiScanpy")
ARIs       
