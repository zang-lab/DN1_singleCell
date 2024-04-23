require("Matrix")
library(Seurat)
library(caret)
library(harmony)
library(gplots)
library(SeuratWrappers)
library(patchwork)
library(gridExtra)
library(monocle3)
library(slingshot)
library(SCENIC)
### p1 as example
## readADT

setwd("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/CITE/WTonly")

RNAprocess_DN1 <- readRDS(file="../allconds/RNAprocess_DN1merge.rds")
#RNAprocess_DN1_allgeneSCT <- readRDS(file="RNAprocess_DN1_allgeneSCT.rds")
#CTassign_DN1 <- read.table(file="../CTassign_DN1.txt",row.names=1,header=F)

RNAprocess <- function(inRNAobj){
	cbmc <- NormalizeData(inRNAobj)
	cbmc <- FindVariableFeatures(cbmc)
	cbmc <- ScaleData(cbmc)
	cbmc <- RunPCA(cbmc, verbose = FALSE)
	cbmc <- FindNeighbors(cbmc, dims = 1:30)
	cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
	cbmc <- RunUMAP(cbmc, dims = 1:30)
	return(cbmc)
}

WT_cell <- names(RNAprocess_DN1$ADT)[which(RNAprocess_DN1$ADT %in% c("WT","WTr2"))]
RNAprocess_DN1WT <- RNAprocess(RNAprocess_DN1[,WT_cell])


### UMAP + features
p1 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")
p2 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="ADT") + ggtitle("cluster")
#p3 <- FeaturePlot(RNAprocess_DN1WT,features="MTpercent")+ ggtitle("chrM%")

pdf(file="WTonly_UMAP_features.pdf",width=10,height=4)
grid.arrange(p1,p2, ncol=2,nrow=1)
dev.off()



#RNAprocess_DN1WT_res0.5 <- FindClusters(RNAprocess_DN1WT, resolution = 0.5, verbose = FALSE)
#RNAprocess_DN1WT_res0.6 <- FindClusters(RNAprocess_DN1WT, resolution = 0.6, verbose = FALSE)
#RNAprocess_DN1WT_res0.7 <- FindClusters(RNAprocess_DN1WT, resolution = 0.7, verbose = FALSE)
#RNAprocess_DN1WT_res0.8 <- FindClusters(RNAprocess_DN1WT, resolution = 0.8, verbose = FALSE)
#RNAprocess_DN1WT_res0.9 <- FindClusters(RNAprocess_DN1WT, resolution = 0.9, verbose = FALSE)
#RNAprocess_DN1WT_res1.0 <- FindClusters(RNAprocess_DN1WT, resolution = 1.0, verbose = FALSE)
#RNAprocess_DN1WT_res1.1 <- FindClusters(RNAprocess_DN1WT, resolution = 1.1, verbose = FALSE)
#RNAprocess_DN1WT_res1.2 <- FindClusters(RNAprocess_DN1WT, resolution = 1.2, verbose = FALSE)
#RNAprocess_DN1WT_res1.3 <- FindClusters(RNAprocess_DN1WT, resolution = 1.3, verbose = FALSE)
#RNAprocess_DN1WT_res1.4 <- FindClusters(RNAprocess_DN1WT, resolution = 1.4, verbose = FALSE)
#RNAprocess_DN1WT_res1.5 <- FindClusters(RNAprocess_DN1WT, resolution = 1.5, verbose = FALSE)
#RNAprocess_DN1WT_res1.6 <- FindClusters(RNAprocess_DN1WT, resolution = 1.6, verbose = FALSE)
#RNAprocess_DN1WT_res1.7 <- FindClusters(RNAprocess_DN1WT, resolution = 1.7, verbose = FALSE)
#RNAprocess_DN1WT_res1.8 <- FindClusters(RNAprocess_DN1WT, resolution = 1.8, verbose = FALSE)
#RNAprocess_DN1WT_res1.9 <- FindClusters(RNAprocess_DN1WT, resolution = 1.9, verbose = FALSE)
#RNAprocess_DN1WT_res2.0 <- FindClusters(RNAprocess_DN1WT, resolution = 2.0, verbose = FALSE)
#
#p0.5 <- DimPlot(RNAprocess_DN1WT_res0.5, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.5")
#p0.6 <- DimPlot(RNAprocess_DN1WT_res0.6, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.6")
#p0.7 <- DimPlot(RNAprocess_DN1WT_res0.7, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.7")
#p0.8 <- DimPlot(RNAprocess_DN1WT_res0.8, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.8")
#p0.9 <- DimPlot(RNAprocess_DN1WT_res0.9, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.9")
#p1.0 <- DimPlot(RNAprocess_DN1WT_res1.0, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.0")
#p1.1 <- DimPlot(RNAprocess_DN1WT_res1.1, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.1")
#p1.2 <- DimPlot(RNAprocess_DN1WT_res1.2, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.2")
#p1.3 <- DimPlot(RNAprocess_DN1WT_res1.3, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.3")
#p1.4 <- DimPlot(RNAprocess_DN1WT_res1.4, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.4")
#p1.5 <- DimPlot(RNAprocess_DN1WT_res1.5, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.5")
#p1.6 <- DimPlot(RNAprocess_DN1WT_res1.6, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.6")
#p1.7 <- DimPlot(RNAprocess_DN1WT_res1.7, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.7")
#p1.8 <- DimPlot(RNAprocess_DN1WT_res1.8, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.8")
#p1.9 <- DimPlot(RNAprocess_DN1WT_res1.9, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.9")
#p2.0 <- DimPlot(RNAprocess_DN1WT_res2.0, label = TRUE,group.by="seurat_clusters") + ggtitle("res2.0")

#pdf(file="UMAP_diffClusterResolution.pdf",width=30,height=30)
#p0.5|p0.6|p0.7|p0.8|p0.9|p1.0|p1.1|p1.2|p1.3|p1.4|p1.5|p1.6|p1.7|p1.8|p1.9|p2.0
#dev.off()

#pdf(file="WTonly_UMAP_diffClusterResolution.pdf",width=30,height=30)
#grid.arrange(p0.5,p0.6,p0.7,p0.8,p0.9,p1.0,p1.1,p1.2,p1.3,p1.4,p1.5,p1.6,p1.7,p1.8,p1.9,p2.0, ncol=4,nrow=4)
#dev.off()
#
#
#p1<-VlnPlot(RNAprocess_DN1WT, features="MTpercent",group.by="seurat_clusters")
#p2<-VlnPlot(RNAprocess_DN1WT, features="MTpercent",group.by="ADT")
#pdf(file="WTonly_boxplot_mtDNA.pdf",width=10,height=4)
#grid.arrange(p1,p2, ncol=2,nrow=1)
#dev.off()

#
#confmat <- function(inmat){
#	feature1 <- unique(inmat[,1])
#	feature2 <- unique(inmat[,2])
#	outmat <- matrix(rep(0, length(feature1)*length(feature2)),nrow=length(feature1),ncol=length(feature2))
#	rownames(outmat) <- feature1
#	colnames(outmat) <- feature2
#	for(i in 1:nrow(inmat)){
#		outmat[inmat[i,1], inmat[i,2]] <- outmat[inmat[i,1], inmat[i,2]] + 1
#	}
#	return(outmat)
#}
#
#confmat_ADT_cluster <- confmat(cbind(as.character(RNAprocess_DN1WT$seurat_clusters),RNAprocess_DN1WT$ADT))
#
#

### MKG
MKGcluster_DN1WT <- FindAllMarkers(RNAprocess_DN1WT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MKGclusterFilter_DN1WT <- MKGcluster_DN1WT[which(MKGcluster_DN1WT[,2] >= log2(2) & MKGcluster_DN1WT[,5] < 0.01),]
write.table(MKGclusterFilter_DN1WT, file="MKGclusterFilter_DN1WT.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(MKGcluster_DN1WT, file="MKGcluster_DN1WT.txt",row.names=T,col.names=T,sep="\t",quote=F)


outdata <- c()
for(C in sort(unique(MKGclusterFilter_DN1WT[,"cluster"])) ){
	tmpdata <- MKGclusterFilter_DN1WT[which(MKGclusterFilter_DN1WT[,"cluster"]==C),]
	outdata <- rbind(outdata, c("DN1",paste0("C",C),tmpdata[order(tmpdata[,"avg_log2FC"],decreasing=T),][1:5,"gene"]))
}

alltop5MKG <- c()
for(this_cluster in 0:17 ){
	#print(this_cluster)
	this_MKGdata <- MKGclusterFilter_DN1WT[which(MKGclusterFilter_DN1WT[,"cluster"] == this_cluster),]
	if(nrow(this_MKGdata)>0){
        this_MKGdata_top5 <- this_MKGdata[order(this_MKGdata[,"avg_log2FC"],decreasing=T),][1:min(5,nrow(this_MKGdata)),]
        tmpobj <- RNAprocess_DN1WT
        tmpobj$target_cluster <- rep(0, ncol(tmpobj))
        tmpobj$target_cluster[which(tmpobj$seurat_clusters == this_cluster)] <- 1
        
        #pdf(file=paste0("UMAP_keygeneExp_C",this_cluster,".pdf"),width=24,height=4)
        FeaturePlot(tmpobj,features=c("target_cluster",as.vector(this_MKGdata_top5[,"gene"])),ncol=6)
        #dev.off()
        ggsave(paste0("WTonly_UMAP_keygeneExp_C",this_cluster,".pdf"),width=24,height=4)
        alltop5MKG <- c(alltop5MKG, as.vector(this_MKGdata_top5[,"gene"]))		
	}
}




#RNAprocess_DN1WT_allgeneSCT <-  SCTransform(RNAprocess_DN1WT, assay = "RNA", verbose = FALSE, ,variable.features.n = 24538)
#
#usegene <- intersect(alltop5MKG, rownames(RNAprocess_DN1WT_allgeneSCT@assays$RNA@scale.data))
#
#expmat_scale_SCT <- RNAprocess_DN1WT_allgeneSCT@assays$SCT@scale.data
#expmat_scale_RNA <- RNAprocess_DN1WT@assays$RNA@scale.data
#
#aveExpMat_scale_RNA <-  matrix(rep(0, length(usegene)*length(unique(RNAprocess_DN1WT$seurat_clusters)) ),nrow=length(usegene))
#rownames(aveExpMat_scale_RNA) <- usegene
#colnames(aveExpMat_scale_RNA) <- paste0("C",0:(length(unique(RNAprocess_DN1WT$seurat_clusters)) -1))
#medianExpMat_scale_RNA <- aveExpMat_scale_RNA
#for(i in 0:(length(unique(RNAprocess_DN1WT$seurat_clusters)) -1)){
#	usecell <- names(RNAprocess_DN1WT$seurat_clusters)[which(RNAprocess_DN1WT$seurat_clusters==i)]
#	tempMat <- expmat_scale_RNA[usegene,usecell]
#	aveExpMat_scale_RNA[, paste0("C",i)] <- apply(tempMat, 1,mean)
#	medianExpMat_scale_RNA[, paste0("C",i)] <- apply(tempMat, 1,median)
#}
#
#aveExpMat_scale_SCT <-  matrix(rep(0, length(usegene)*length(unique(RNAprocess_DN1WT$seurat_clusters)) ),nrow=length(usegene))
#rownames(aveExpMat_scale_SCT) <- usegene
#colnames(aveExpMat_scale_SCT) <- paste0("C",0:(length(unique(RNAprocess_DN1WT$seurat_clusters)) -1))
#medianExpMat_scale_SCT <- aveExpMat_scale_SCT
#for(i in 0:(length(unique(RNAprocess_DN1WT$seurat_clusters)) -1)){
#	usecell <- names(RNAprocess_DN1WT$seurat_clusters)[which(RNAprocess_DN1WT$seurat_clusters==i)]
#	tempMat <- expmat_scale_SCT[usegene,usecell]
#	aveExpMat_scale_SCT[, paste0("C",i)] <- apply(tempMat, 1,mean)
#	medianExpMat_scale_SCT[, paste0("C",i)] <- apply(tempMat, 1,median)
#}
#
#
#usecolor <- c("blue","white","red")
#ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
#pdf(file="WTonly_keygeneExpRNA_cluster_heatmap.pdf")
#heatmap.2(aveExpMat_scale_RNA, trace="none", col=ColorRamp,main="aveExp scaleRNA",cexRow=0.6)
#dev.off()
#
#pdf(file="WTonly_keygeneExpSCT_cluster_heatmap.pdf")
#heatmap.2(aveExpMat_scale_SCT, trace="none", col=ColorRamp,main="aveExp scaleSCT",cexRow=0.6)
#dev.off()



######## old anno in r1WTonly

old_proj <- readRDS("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/CITE/old_WTonly/RNAprocess_DN1WT.rds")
r1CT4 <- (old_proj$CT_4subset)
names(r1CT4) <- paste0("r1_",names(r1CT4))

RNAprocess_DN1WT$r1CT4 <- rep("NA",length(RNAprocess_DN1WT$ADT))
names(RNAprocess_DN1WT$r1CT4) <- names(RNAprocess_DN1WT$ADT)
RNAprocess_DN1WT$r1CT4[names(r1CT4)] <- r1CT4


p1 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="seurat_clusters") + ggtitle("seurat_clusters")
p2 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="r1CT4") + ggtitle("r1WTonly CT4")
#p3 <- FeaturePlot(RNAprocess_DN1WT,features="MTpercent")+ ggtitle("chrM%")

pdf(file="WTonly_UMAP_r1CT4.pdf",width=10,height=4)
grid.arrange(p1, p2, ncol=2,nrow=1)
dev.off()




renameCluster <- rep("NA",length(RNAprocess_DN1WT$seurat_clusters))
names(renameCluster) <- names(RNAprocess_DN1WT$seurat_clusters) 
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 0  )] <- 1
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 1  )] <- 2
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 2  )] <- 4
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 3  )] <- 5
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 4  )] <- 11
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 5  )] <- 12
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 6  )] <- 6
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 7  )] <- 7
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 8  )] <- 13
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 9  )] <- 14
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 10 )] <- 3
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 11 )] <- 15
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 12 )] <- 8
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 13 )] <- 16
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 14 )] <- 17
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 15 )] <- 18
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 16 )] <- 9
renameCluster[which(RNAprocess_DN1WT$seurat_clusters == 17 )] <- 10

RNAprocess_DN1WT$renameCluster <- renameCluster
RNAprocess_DN1WT_allgeneSCT$renameCluster <- renameCluster



### old cluster cmp
#ETPs_quiescent: 0, 1, 10. 
#ETPs_proliferative: 2, 3, 6. 
#Non-T potentials: 7, 12, 16, 17
#T+ILCs: 4, 5, 8, 9, 11, 13, 14, 15. 

#CT_4subset <- rep("NA",length(RNAprocess_DN1WT$seurat_clusters))
#names(CT_4subset) <- names(RNAprocess_DN1WT$seurat_clusters) 
#CT_4subset[which(RNAprocess_DN1WT$seurat_clusters %in% c(0,1,10))] <- "ETPs_quiescent"
#CT_4subset[which(RNAprocess_DN1WT$seurat_clusters %in% c(2,3,6))] <- "ETPs_proliferative"
#CT_4subset[which(RNAprocess_DN1WT$seurat_clusters %in% c(7,12,16,17))] <- "nonT_potential"
#CT_4subset[which(RNAprocess_DN1WT$seurat_clusters %in% c(4, 5, 8, 9, 11, 13, 14, 15))] <- "T_ILCs"

### new cluster cmp
#ETPs_quiescent: 1, 2, 3. 
#ETPs_proliferative: 4, 5, 6. 
#Non-T potentials: 7-10
#T+ILCs: 11-18. Clusters 11, 15 and 16 have higher expression of Gata3, Rorc and Tox, called T_ILC_Gata3/Rorchi
# The rest has higher expression of Id2 and Il7r, called T_ILC_Id2/Il7rhi


CT_4subset <- rep("NA",length(RNAprocess_DN1WT$renameCluster))
names(CT_4subset) <- names(RNAprocess_DN1WT$renameCluster) 
CT_4subset[which(RNAprocess_DN1WT$renameCluster %in% c(1,2,3))] <- "ETPs_quiescent"
CT_4subset[which(RNAprocess_DN1WT$renameCluster %in% c(4,5,6))] <- "ETPs_proliferative"
CT_4subset[which(RNAprocess_DN1WT$renameCluster %in% c(7,8,9,10))] <- "nonT_potential"
CT_4subset[which(RNAprocess_DN1WT$renameCluster %in% c(11:18))] <- "T_ILCs"
RNAprocess_DN1WT$CT_4subset <- CT_4subset


CT_5subset <- rep("NA",length(RNAprocess_DN1WT$renameCluster))
names(CT_5subset) <- names(RNAprocess_DN1WT$renameCluster) 
CT_5subset[which(RNAprocess_DN1WT$renameCluster %in% c(1,2,3))] <- "ETPs_quiescent"
CT_5subset[which(RNAprocess_DN1WT$renameCluster %in% c(4,5,6))] <- "ETPs_proliferative"
CT_5subset[which(RNAprocess_DN1WT$renameCluster %in% c(7,8,9,10))] <- "nonT_potential"
CT_5subset[which(RNAprocess_DN1WT$renameCluster %in% c(11,15,16))] <- "T_ILCs_GRhi"
CT_5subset[which(RNAprocess_DN1WT$renameCluster %in% c(12,13,14,17,18))] <- "T_ILCs_IIhi"
RNAprocess_DN1WT$CT_5subset <- CT_5subset
RNAprocess_DN1WT_allgeneSCT$CT_5subset <- CT_5subset


### SCT data
RNAprocess_DN1WT_allgeneSCT <-  SCTransform(RNAprocess_DN1WT, assay = "RNA", verbose = FALSE, ,variable.features.n = 24538)


### UMAP + features
p1 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="ADT") + ggtitle("ADT")
p2 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="renameCluster") + ggtitle("rename cluster")
p3 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="CT_5subset") + ggtitle("cell type")
pdf(file="WTonly_UMAP_newFeatures.pdf",width=15,height=4)
grid.arrange(p1,p2,p3, ncol=3,nrow=1)
dev.off()

p1 <- DimPlot(RNAprocess_DN1WT, label = FALSE,group.by="ADT") + ggtitle("ADT")+ NoLegend()
p2 <- DimPlot(RNAprocess_DN1WT, label = FALSE,group.by="renameCluster") + ggtitle("rename cluster")+ NoLegend()
p3 <- DimPlot(RNAprocess_DN1WT, label = FALSE,group.by="CT_5subset") + ggtitle("cell type")+ NoLegend()
pdf(file="../WTonly/WTonly_UMAP_newFeatures_nolabel.pdf",width=12,height=4)
grid.arrange(p1,p2,p3, ncol=3,nrow=1)
dev.off()


#keygene_CT4subset <- c("Aurkb","Birc5","Ccna2","Ccnb1","Cdk1","Lig1","Mki67","Mcm3","Mcm6","Ezh2","Dnmt1","Bcl11a","Spi1","Lyl1","Cd24a","Kit","Flt3","Notch1","Hhex","Hes1","Tcf7","Cd19","Cd79a","Blnk","Irf8","Gata2","Elane","Mpo","Thy1","Cd3e","Zap70","Prkcq","Il7r","Ets1","Bcl11b","Id2","Tox","Gata3","Zbtb16","Rorc")
#neworder_CT4subset <- c("ETPs_proliferative","ETPs_quiescent","nonT_potential","T_ILCs")
#RNAprocess_DN1WT_allgeneSCT$CT4_subset_reorder <- factor(RNAprocess_DN1WT_allgeneSCT$CT_4subset, levels = neworder_CT4subset)
#pdf(file="SCTheatmap_WTonly_celltype.pdf",width=10,height=10)
#DoHeatmap(RNAprocess_DN1WT_allgeneSCT, features=keygene_CT4subset,group.by="CT4_subset_reorder")
#dev.off()


keygene_CT5subset <- c("Aurkb","Birc5","Ccna2","Ccnb1","Cdk1","Lig1","Mki67","Ezh2","Dnmt1","Mcm3","Mcm6","Bcl11a","Lyl1","Cd24a","Kit","Flt3","Notch1","Hhex","Hes1","Spi1","Mpo","Gata2","Elane","Irf8","Cd19","Cd79a","Blnk","Thy1","Cd3e","Zap70","Prkcq","Tcf7","Ets1","Bcl11b","Tox","Gata3","Rorc","Id2","Il7r","Zbtb16")
neworder_CT5subset <- c("ETPs_proliferative","ETPs_quiescent","nonT_potential","T_ILCs_GRhi","T_ILCs_IIhi")
RNAprocess_DN1WT_allgeneSCT$CT5_subset_reorder <- factor(RNAprocess_DN1WT_allgeneSCT$CT_5subset, levels = neworder_CT5subset)
pdf(file="SCTheatmap_WTonly_celltype5.pdf",width=10,height=10)
DoHeatmap(RNAprocess_DN1WT_allgeneSCT, features=keygene_CT5subset,group.by="CT5_subset_reorder")
dev.off()



######### CT5 MKG
Idents(RNAprocess_DN1WT) <- "CT_5subset"
MKGct5_DN1WT <- FindAllMarkers(RNAprocess_DN1WT, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)

MKGct5Filter_DN1WT <- MKGct5_DN1WT[which(MKGct5_DN1WT[,2] >= log2(2) & MKGct5_DN1WT[,5] < 0.01),]
write.table(MKGct5Filter_DN1WT, file="MKGct5Filter_DN1WT.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(MKGct5_DN1WT, file="MKGct5_DN1WT.txt",row.names=T,col.names=T,sep="\t",quote=F)


MKGct5_DN1WT <- read.table("MKGct5_DN1WT.txt",row.names=1,header=T,sep="\t")












########## monocle3
### monocle 3

cds <- as.cell_data_set(RNAprocess_DN1WT)
cds <- cluster_cells(cds, cluster_method = "louvain")
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 1]))
cds$renameCluster <- RNAprocess_DN1WT$renameCluster
cds$CT_5subset <- RNAprocess_DN1WT$CT_5subset

p1 <- plot_cells(cds, color_cells_by = "renameCluster", show_trajectory_graph = FALSE) + ggtitle("renameCluster")
p2 <- plot_cells(cds, color_cells_by = "CT_5subset", show_trajectory_graph = FALSE) + ggtitle("CT_5subset")
p3 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)+ ggtitle("partition")

p4 <- plot_cells(cds,
           color_cells_by = "renameCluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE) + ggtitle("cluster+trajectory")
p5 <- plot_cells(cds,
           color_cells_by = "CT_5subset",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,show_trajectory_graph = FALSE, 
           trajectory_graph_test_res = pr_graph_test_res(cds, neighbor_graph="knn"),
           label_branch_points=FALSE) + ggtitle("CT5+trajectory")

p6 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")+ggtitle("trajectory+pseudotime")

pdf(file="UMAP_monocle3Trajectory.pdf",width=15,height=10)
wrap_plots(p1, p2, p3,p4,p5,p6)
dev.off()







########### slingshot


data <- RNAprocess_DN1WT
dimred <- data@reductions$umap@cell.embeddings
clustering <- data$CT_5subset
counts <- as.matrix(data@assays$RNA@counts[data@assays$RNA@var.features, ])


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)
#pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
#pal <- display.brewer.pal(n = 22, name = 'Set1')
pal<- ggplotColours(n=5)

#par(mfrow = c(1, 2))
#plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
#for (i in levels(clustering)) {
#  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
#}
pdf(file="UMAP_slingshotTrajectory.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(dimred[, 1:2], col = pal[as.numeric(as.factor(clustering))], cex = 0.5, pch = 16,main="Linage")
lines(SlingshotDataSet(lineages), lwd = 1, col = "black",cex=0.5)

curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
plot(dimred[, 1:2], col = pal[as.numeric(as.factor(clustering))], cex = 0.5, pch = 16,main="Principal Curves")
lines(SlingshotDataSet(curves), lwd = 1, col = "black")
dev.off()
















######### SCENIC


#library(SCENIC)
#
#.openDev <- function(fileName, devType, ...)
#{
#  if(devType=="pdf")
#    pdf(paste0(fileName, ".pdf"), ...)
#  
#  if(devType=="png")
#    png(paste0(fileName, ".png", type="cairo"), ...)
#  
#  if(devType=="cairo_pfd") # similar to Cairo::CairoPDF?
#    grDevices::cairo_pdf(paste0(fileName, ".pdf"), ...)
#}
#
#.openDevHeatmap <- function(fileName, devType)
#{
#  if(devType!="pdf") 
#  {
#    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
#    if(devType!="png") .openDev(fileName=fileName, devType=devType)
#    fileName <- NA
#  }else{
#    fileName <- paste0(fileName,".pdf")
#  }
#  return(fileName)
#}
#
#.closeDevHeatmap <- function(devType)
#{
#  if(devType!="pdf") 
#  {
#    dev.off()
#  }
#}
#
#data <- RNAprocess_DN1WT
#exprMat <- as.matrix(data@assays$RNA@counts[data@assays$RNA@var.features, ])
#cellInfo <- FetchData(data,vars=c("CT_5subset","nFeature_RNA","nCount_RNA"))
##cellInfo[,1] <- paste0("C",cellInfo[,1])
#colnames(cellInfo) <- c("CellType","nGene","nUMI")
#cellInfo <- data.frame(cellInfo)
#dir.create("int")
#saveRDS(cellInfo, file="int/cellInfo.Rds")
#
#
#tmpdata <- ggplotColours(n=5)
#names(tmpdata) <- unique(data$CT_5subset)#paste0("C",seq(0,13))
#colVars <-list(CellType = tmpdata)
#colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
#saveRDS(colVars, file="int/colVars.Rds")
#plot.new(); legend("top",cex=0.4,bty="n", fill=colVars$CellType, legend=names(colVars$CellType))
#
#scenicOptions <- initializeScenic(org="mgi", dbDir="/nv/vol190/zanglab/sh8tv/Data/SCENIC/mm10/resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/", nCores=10)
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
#scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status
#
#genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
#                           minCountsPerGene=3*.01*ncol(exprMat),
#                           minSamples=ncol(exprMat)*.01)
#
#
#keygene_tmp <- c("Tcf7","Lef1","Tle1","Tle3","Tle4","Gata3","Hes1","Bcl11b","Runx1","Nfat","Il7r","Ctcf","Notch1","Nfatc1","Tcf3","Pbx1","Id3","Kit","Lck","Flt3")
#keygene_HSC <- c("Spi1","Lyl1","Hhex","Bcl11a","Hoxa9","Meis1")
#keygene_Tlineage <- c("Notch3","Erg","Tox","Eomes","Sox4","Sox13","Ltk","Prkcq","Tcf7","Bcl11b","Gata3")
#keygene_Blineage <- c("Ebf1","Pax5","Btk","Syk","Cd19","Cd79a","Cd79b","Blnk")
#keygene_Survival <- "Bax"
#keygene_Dendritic <- c("Irf8","Tcf4","Itgax","Batf3","Irf4","Tcf3","Nfil3","Zeb2","Klf4")
#keygene_Myeloid <- c("Cebpa","Cebpb","Cebpe","Rara","Mzf1","Hoxa10","Hoxb7","Hoxb8","Itgam","Gata2","Gfi1")
#keygene_InnateLymphoid <- c("Id2","Tox","Rorc","Ets1","Bcl11b","Gata3","Zbtb16")
#
#
#interestingGenes <- keygene_Tlineage#c("Sox9", "Sox10", "Dlx5")
## any missing?
#interestingGenes[which(!interestingGenes %in% genesKept)]
#
#exprMat_filtered <- exprMat[genesKept, ]
#dim(exprMat_filtered)
#
#runCorrelation(exprMat_filtered, scenicOptions)
#
#exprMat_filtered <- log2(exprMat_filtered+1) 
#runGenie3(exprMat_filtered, scenicOptions)
#
#exprMat_log <- log2(exprMat+1)
#
#
#
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status
#scenicOptions@settings$verbose <- TRUE
#scenicOptions@settings$nCores <- 10
#scenicOptions@settings$seed <- 123
#
## For a very quick run: 
## coexMethod=c("top5perTarget")
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
## save...
#
#scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
##
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status
#
#### step3
#regulons <- loadInt(scenicOptions, "regulons")
#regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
#regulons <- regulons[lengths(regulons)>=10]
#if(length(regulons) <2)  stop("Not enough regulons with at least 10 genes.")
#
#
## Add the TF to the regulon (keeping it only once) & rename regulon
#regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
#names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")
#saveRDS(regulons, file=getIntName(scenicOptions, "aucell_regulons"))
#
#msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 3. Analyzing the network activity in each individual cell")
#if(getSettings(scenicOptions, "verbose")) message(msg)
#
#msg <- paste0("\nNumber of regulons to evaluate on cells: ", length(regulons),
#              "\nBiggest (non-extended) regulons: \n",
#              paste("\t", grep("_extended",names(regulons),invert = T, value = T)[1:10], collapse="\n")) # TODO maxlen?
#if(getSettings(scenicOptions, "verbose")) message(msg)
#
#
## 3.1, create ranking
#nCores <- getSettings(scenicOptions, "nCores")
#
#if(is.data.frame(exprMat)) 
#{
#  supportedClasses <- paste(gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")), collapse=", ")
#  supportedClasses <- gsub("-method", "", supportedClasses)
#  
#  stop("'exprMat' should be one of the following classes: ", supportedClasses, 
#       "\n(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
#}
#
#set.seed(getSettings(scenicOptions,"seed"))
#tryCatch({
#  .openDev(fileName=getIntName(scenicOptions, "aucell_genesStatsPlot"),
#           devType=getSettings(scenicOptions, "devType"))
#  aucellRankings <- AUCell_buildRankings(exprMat, nCores=nCores, 
#                                         plotStats=TRUE, verbose=getSettings(scenicOptions, "verbose"))
#  abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3)
#  dev.off()
#},error = function(e) {
#  message("Catched error in AUCell_buildRankings() or in the histogram plot: ", e$message)
#})
#saveRDS(aucellRankings, file=getIntName(scenicOptions, "aucell_rankings"))
#
#
##scenicOptions <- readRDS(file="int/scenicOptions.Rds")
##regulons <- readRDS(file=getIntName(scenicOptions, "aucell_regulons"))
##aucellRankings <- readRDS(file=getIntName(scenicOptions, "aucell_rankings"))
#
#
## 3.2, calculate AUC
#regulonAUC <- AUCell_calcAUC(regulons, aucellRankings, 
#                             aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=nCores)
#
## Order the modules by similarity, for easier exploration in the upcoming steps & save
#regulonOrder <- orderAUC(regulonAUC) # added to AUCell 1.5.1
#regulonAUC <- regulonAUC[regulonOrder,]
#saveRDS(regulonAUC, file=getIntName(scenicOptions, "aucell_regulonAUC"))
#
##aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log) # default t-SNE
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status
#
#setwd("/Users/sh8tv/Dropbox/TCF1_Tcell/Result/DN1_analysis/CITE/combine4condition/trajectory_network")
#RNAprocess_DN1 <- readRDS(file="../RNAprocess_DN1.rds")
#data <- RNAprocess_DN1
#exprMat <- as.matrix(data@assays$RNA@counts[data@assays$RNA@var.features, ])
#exprMat_log <- log2(exprMat+1)
#scenicOptions <- readRDS(file="int/scenicOptions.Rds")
#regulons <- readRDS(file=getIntName(scenicOptions, "aucell_regulons"))
#regulonAUC <- readRDS(file=getIntName(scenicOptions, "aucell_regulonAUC"))
#aucellRankings <- readRDS(file=getIntName(scenicOptions, "aucell_rankings"))


#aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
#savedSelections <- shiny::runApp(aucellApp)
#
## Save the modified thresholds:
#newThresholds <- savedSelections$thresholds
#scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
#saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#
##step4
#scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
#
#fetchEXT <- function(inname){
#  tmp <-strsplit(inname, "_")[[1]][2]
#  return(strsplit(tmp," ")[[1]][1])
#}
#termEXT <- unlist(lapply(rownames(regulonAUC),fetchEXT))
#regulonAUC_filter <- regulonAUC[is.na(termEXT), ]
#
#
#### AUC heatmap
#cellCluster <- data$seurat_clusters[colnames(regulonAUC)]
#  
#NMF::aheatmap(getAUC(regulonAUC),#[,cells2plot],
#              annCol=cellInfo,
#              annColor=colVars,
#              main="regulonAUC primary+extend",
#              sub=paste("all cells"),
#              filename="output/regulonAUC_heatmap_primaryExtend.pdf")
##NMF::aheatmap(getAUC(regulonAUC[,order(as.numeric(cellCluster))]),#[,cells2plot],
##              annCol=cellInfo,
##              annColor=colVars,
##              Colv=F,
##              main="regulonAUC primary+extend noCellCluster",
##              sub=paste("all cells"),
##              filename="output/regulonAUC_heatmap_primaryExtend_noCellCluster.pdf")
##
#NMF::aheatmap(getAUC(regulonAUC_filter),#[,cells2plot],
#              annCol=cellInfo,
#              annColor=colVars,
#              main="regulonAUC primaryOnly",
#              sub=paste("all cells"),
#              filename="output/regulonAUC_heatmap_primaryOnly.pdf")
##NMF::aheatmap(getAUC(regulonAUC_filter)[,order(cellInfo[,"CellType"])],#[,cells2plot],
##              annCol=cellInfo,
##              annColor=colVars,
##              Colv=F,Rowv=F,
##              main="regulonAUC primaryOnly noCellCluster",
##              sub=paste("all cells"),
##              filename="output/regulonAUC_heatmap_primaryOnly_noCellCluster.pdf")
#
##my_list <- as.list(pal)
##names(my_list) <- paste0("C",seq(0,21))
#
## RSS score
#rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
#rssPlot <- plotRSS(rss)
#pdf(file="output/RSSheatmap.pdf",width=6,height=8)
#rssPlot$plot
#dev.off()
#
#rss_filter <- calcRSS(AUC=getAUC(regulonAUC_filter), cellAnnotation=cellInfo[colnames(regulonAUC_filter), "CellType"], )
#rssPlot_filter <- plotRSS(rss_filter)
#pdf(file="output/RSSheatmap_primaryOnly.pdf",width=6,height=8)
#rssPlot_filter$plot
#dev.off()
##pdf(file="output/RSS_")
##plotRSS_oneSet(rss, setName = "C0")
#
#
#dr_coords <- Embeddings(RNAprocess_DN1, reduction="umap")
#
#trimRN <- function(inname){
#  tmp <-strsplit(inname, "_")[[1]][1]
#  return(strsplit(tmp," ")[[1]][1])
#}
#selectedGenes <- unique(unlist(lapply(rownames(regulonAUC),trimRN)))
#
#for(usegene in selectedGenes){
#  print(usegene)
#  pdf(file=paste0("output/UMAP_SCENICauc",usegene,".pdf"),height=4,width=16)
#  par(mar=c(4,4,2,2),mfrow=c(1,4))
#  #layout(matrix(c(1,2,3,4,8,5,6,7),nrow=2,byrow=T))
#  plot(dimred[, 1:2], col = pal[as.numeric(clustering)], cex = 0.5, pch = 16,main=usegene)
#  AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, usegene,onlyNonDuplicatedExtended=T), exprMat=exprMat_log,plots = c("binaryAUC","AUC","expression"))
#  dev.off()
#}










### subset assign
#ETPs-quiescent: clusters 0, 1, 2
#ETPs-proliferative: cluster 3, 5, 7. 
#non-T potential" cluster 4
#T+ILCs: clusters 6, 8, 9, 10, 11, 12, 13. 
#CT_4subset <- rep("NA",length(RNAprocess_DN1WT$seurat_clusters))
#names(CT_4subset) <- names(RNAprocess_DN1WT$seurat_clusters) 
#CT_4subset[which(RNAprocess_DN1WT$seurat_clusters %in% c(0,1,2))] <- "ETPs_quiescent"
#CT_4subset[which(RNAprocess_DN1WT$seurat_clusters %in% c(3,5,7))] <- "ETPs_proliferative"
#CT_4subset[which(RNAprocess_DN1WT$seurat_clusters %in% c(4))] <- "nonT_potential"
#CT_4subset[which(RNAprocess_DN1WT$seurat_clusters %in% c(6, 8, 9, 10, 11, 12, 13))] <- "T_ILCs"
#
#RNAprocess_DN1WT$CT_4subset <- CT_4subset
#RNAprocess_DN1WT_allgeneSCT$CT_4subset <- CT_4subset
#
#
#p1 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")
#p2 <- DimPlot(RNAprocess_DN1WT, label = TRUE,group.by="CT_4subset") + ggtitle("CT_4subset")
#
#pdf(file="WTonly_UMAP_features_CT4subset.pdf",width=10,height=4)
#grid.arrange(p1,p2, ncol=2,nrow=1)
#dev.off()
#
#
#pdf(file="WTonly_keygeneExpSCT_CT4subset_heatmap.pdf")
#heatmap.2(aveExpMat_scale_SCT, trace="none", col=ColorRamp,main="aveExp scaleSCT",cexRow=0.6)
#dev.off()




##### find markers for CT_4subset

### MKG
#RNAprocess_DN1WT$CT <- as.factor(RNAprocess_DN1WT$CT_4subset)
#RNAprocess_DN1WT_forMarkerGene <- SetIdent(RNAprocess_DN1WT, value = "CT_4subset")
#MKG_CT4subset_DN1WT <- FindAllMarkers(RNAprocess_DN1WT_forMarkerGene,group.by="CT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#MKG_CT4subset_Filter_DN1WT <- MKG_CT4subset_DN1WT[which(MKG_CT4subset_DN1WT[,2] >= log2(2) & MKG_CT4subset_DN1WT[,5] < 0.01),]
#write.table(MKG_CT4subset_Filter_DN1WT, file="MKG_CT4subsetFilter_DN1WT.txt",row.names=T,col.names=T,sep="\t",quote=F)
#write.table(MKG_CT4subset_DN1WT, file="MKGCT4subset_DN1WT.txt",row.names=T,col.names=T,sep="\t",quote=F)
#
#
#outdata <- c()
#for(C in sort(unique(MKG_CT4subset_Filter_DN1WT[,"cluster"])) ){
#	tmpdata <- MKG_CT4subset_Filter_DN1WT[which(MKG_CT4subset_Filter_DN1WT[,"cluster"]==C),]
#	outdata <- rbind(outdata, c("DN1",C,tmpdata[order(tmpdata[,"avg_log2FC"],decreasing=T),][1:5,"gene"]))
#}
#write.table(outdata, file='MKG_CT4subsetFilter_DN1WT_top5names.txt',row.names=F,col.names=F,sep="\t",quote=F)
#alltop5MKG <- c()
#for(this_cluster in sort(unique(MKG_CT4subset_Filter_DN1WT[,"cluster"])) ){
#	#print(this_cluster)
#	this_MKGdata <- MKG_CT4subset_Filter_DN1WT[which(MKG_CT4subset_Filter_DN1WT[,"cluster"] == this_cluster),]
#	if(nrow(this_MKGdata)>0){
#        this_MKGdata_top5 <- this_MKGdata[order(this_MKGdata[,"avg_log2FC"],decreasing=T),][1:min(5,nrow(this_MKGdata)),]
#        tmpobj <- RNAprocess_DN1WT
#        tmpobj$target_cluster <- rep(0, ncol(tmpobj))
#        tmpobj$target_cluster[which(tmpobj$CT_4subset == this_cluster)] <- 1
#        
#        #pdf(file=paste0("UMAP_keygeneExp_C",this_cluster,".pdf"),width=24,height=4)
#        FeaturePlot(tmpobj,features=c("target_cluster",as.vector(this_MKGdata_top5[,"gene"])),ncol=6)
#        #dev.off()
#        ggsave(paste0("CT4subet_MKG/WTonly_UMAP_keygeneExp_CT4subet_",this_cluster,".pdf"),width=24,height=4)
#        alltop5MKG <- c(alltop5MKG, as.vector(this_MKGdata_top5[,"gene"]))		
#	}
#}
#





##### knowledge base marker
keygene <- c("Tcf7","Lef1","Tle1","Tle3","Tle4","Gata3","Hes1","Bcl11b","Runx1","Nfat","Il7r","Ctcf","Notch1","Nfatc1","Tcf3","Pbx1","Id3","Kit","Lck","Flt3")
usekeygene <- intersect(keygene, rownames(RNAprocess_DN1WT_allgeneSCT@assays$RNA@counts))


#allgenes <- rownames(RNAprocess_DN1WT_allgeneSCT@assays$RNA@counts)
#write.table(allgenes,file="allgenename.txt",row.names=F,col.names=F,sep="\t",quote=F)
keygene_HSC <- c("Spi1","Lyl1","Hhex","Bcl11a","Hoxa9","Meis1")
keygene_Tlineage <- c("Notch3","Erg","Tox","Eomes","Sox4","Sox13","Ltk","Prkcq","Tcf7","Bcl11b","Gata3")
keygene_Blineage <- c("Ebf1","Pax5","Btk","Syk","Cd19","Cd79a","Cd79b","Blnk")
keygene_Survival <- "Bax"
keygene_Dendritic <- c("Irf8","Tcf4","Itgax","Batf3","Irf4","Tcf3","Nfil3","Zeb2","Klf4")
keygene_Myeloid <- c("Cebpa","Cebpb","Cebpe","Rara","Mzf1","Hoxa10","Hoxb7","Hoxb8","Itgam","Gata2","Gfi1")
keygene_InnateLymphoid <- c("Id2","Tox","Rorc","Ets1","Bcl11b","Gata3","Zbtb16")
keygene_Cebp <- c("Cebpa","Cebpb","Cebpd","Cebpe","Cebpg","Cebpz","Cebpzos")
keygene_other <- c("Cd24a","Rag1","Rag2","Rbpj","Tbx21","Cd3e","Cd3d","Cd3g","Thy1","Satb1","Cxcr6","Mpo")
keygene_sp1 <- c("Notch1","Hes1","Hhex","Kit")
keygene_ETPproliferative <- c("Aurkb","Birc5","Bub1b","Ccna2","Ccnb1","Ccnb2","Cdk1","Dnmt1","Ezh2","Hmgb1","Lig1","Mki67")
keygene_ETPquiescent <- c("Cd24a","Notch1","Cdk6","Hes1","Hhex","Egr1","Etv6","Jun","Myc","Sox4","Spi1","Bcl11a")
keygene_nonTpotential <- c("Cebpa","Cebpb","Cebpd","Gapdh","Hk3","Fos","Lmo2","Ldha","Pgam1","Irf8","Spi1","Zeb2")
keygene_TILClike <- c("Cd3e","Cd3d","Ets1","Il2rg","Il7r","Thy1","Lck","Zap70","Gata3","Id2","Id3","Tox")
keygene_CT4subset <- c("Notch1","Hes1","Hhex","Cd24a","Cdk6","Etv6","Bcl11a","Egr1","Jun","Sox4","Myc","Spi1","Aurkb","Birc5","Bub1b","Ccna2","Ccnb1","Ccnb2","Hmgb1","Cdk1","Mki67","Dnmt1","Ezh2","Lig1","Gapdh","Ldha","Pgam1","Hk3 ","Cebpa","Cebpb","Cebpd","Fos","Lmo2","Irf8","Spi1","Zeb2","Cd3e","Cd3d","Lck","Id3","Tox","Ets1","Zap70","Il2rg","Thy1","Il7r","Id2","Gata3")
keygene_CT5subset <- c("Aurkb","Birc5","Ccna2","Ccnb1","Cdk1","Lig1","Mki67","Ezh2","Dnmt1","Mcm3","Mcm6","Bcl11a","Lyl1","Cd24a","Kit","Flt3","Notch1","Hhex","Hes1","Spi1","Mpo","Gata2","Elane","Irf8","Cd19","Cd79a","Blnk","Thy1","Cd3e","Zap70","Prkcq","Tcf7","Ets1","Bcl11b","Tox","Gata3","Rorc","Id2","Il7r","Zbtb16")
keygene_CT5subset <- c("Aurkb","Birc5","Ccna2","Ccnb1","Cdk1","Lig1","Mki67","Ezh2","Dnmt1","Mcm3","Mcm6","Bcl11a","Lyl1","Cd24a","Kit","Flt3","Notch1","Hhex","Hes1","Spi1","Mpo","Gata2","Elane","Irf8","Cd19","Cd79a","Blnk","Thy1","Cd3e","Zap70","Prkcq","Tcf7","Ets1","Bcl11b","Tox","Gata3","Rorc","Id2","Il7r","Zbtb16")
keygene_discussionRelated <- c("Klf12","Nfil3","Klrd1","Clnk","Cd160","Klrb1c","Itga2")
### key gene exp


pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_discussionRelated.pdf",width=16,height=8)
FeaturePlot(RNAprocess_DN1WT,features=keygene_discussionRelated,ncol=4)
dev.off()


pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_CT5subset.pdf",width=20,height=32)
FeaturePlot(RNAprocess_DN1WT,features=keygene_CT5subset,ncol=5)
dev.off()


pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_temp.pdf",width=20,height=16)
FeaturePlot(RNAprocess_DN1WT,features=usekeygene,ncol=5)
dev.off()


pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_HSC.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1WT,features=keygene_HSC,ncol=5)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_Tlineage.pdf",width=20,height=12)
FeaturePlot(RNAprocess_DN1WT,features=keygene_Tlineage,ncol=5)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_Blineage.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1WT,features=keygene_Blineage,ncol=5)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_Survival.pdf",width=20,height=4)
FeaturePlot(RNAprocess_DN1WT,features=keygene_Survival,ncol=5)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_Dendritic.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1WT,features=keygene_Dendritic,ncol=5)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_Myeloid.pdf",width=20,height=12)
FeaturePlot(RNAprocess_DN1WT,features=keygene_Myeloid,ncol=5)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_InnateLymphoid.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1WT,features=keygene_InnateLymphoid,ncol=5)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_other.pdf",width=20,height=12)
FeaturePlot(RNAprocess_DN1WT,features=keygene_other,ncol=5)
dev.off()


pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_sp1.pdf",width=20,height=4)
FeaturePlot(RNAprocess_DN1WT,features=keygene_sp1,ncol=5)
dev.off()


pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_ETPproliferative.pdf",width=16,height=12)
FeaturePlot(RNAprocess_DN1WT,features=keygene_ETPproliferative,ncol=4)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_ETPquiescent.pdf",width=16,height=12)
FeaturePlot(RNAprocess_DN1WT,features=keygene_ETPquiescent,ncol=4)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_nonTpotential.pdf",width=16,height=12)
FeaturePlot(RNAprocess_DN1WT,features=keygene_nonTpotential,ncol=4)
dev.off()

pdf(file="WTonly_UMAP_knowledgeKeyGeneExp_TILClike.pdf",width=16,height=12)
FeaturePlot(RNAprocess_DN1WT,features=keygene_TILClike,ncol=4)
dev.off()


#FeaturePlot(tmpobj,features=usekeygene,ncol=4)
#dev.off()
#ggsave(paste0("UMAP_keygeneExp_C",this_cluster,".pdf"),width=16,height=16)

#pdf(file="boxplot_knowledgeKeyGeneExp_ADTGroup.pdf",width=20,height=16)
#VlnPlot(RNAprocess_DN1WT, features=usekeygene,ncol=5,group.by="ADT")
#dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_ClusterGroup.pdf",width=25,height=16)
VlnPlot(RNAprocess_DN1WT, features=usekeygene,ncol=5,group.by="seurat_clusters")
dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_HSC.pdf",width=20,height=8)
VlnPlot(RNAprocess_DN1WT,features=keygene_HSC,ncol=5,group.by="seurat_clusters")
dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_Tlineage.pdf",width=20,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_Tlineage,ncol=5,group.by="seurat_clusters")
dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_Blineage.pdf",width=20,height=8)
VlnPlot(RNAprocess_DN1WT,features=keygene_Blineage,ncol=5,group.by="seurat_clusters")
dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_Survival.pdf",width=20,height=4)
VlnPlot(RNAprocess_DN1WT,features=keygene_Survival,ncol=5,group.by="seurat_clusters")
dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_Dendritic.pdf",width=20,height=8)
VlnPlot(RNAprocess_DN1WT,features=keygene_Dendritic,ncol=5,group.by="seurat_clusters")
dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_Myeloid.pdf",width=20,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_Myeloid,ncol=5,group.by="seurat_clusters")
dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_InnateLymphoid.pdf",width=20,height=8)
VlnPlot(RNAprocess_DN1WT,features=keygene_InnateLymphoid,ncol=5,group.by="seurat_clusters")
dev.off()

pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_other.pdf",width=20,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_other,ncol=5,group.by="seurat_clusters")
dev.off()


pdf(file="WTonly_boxplot_knowledgeKeyGeneExp_sp1.pdf",width=20,height=4)
VlnPlot(RNAprocess_DN1WT,features=keygene_sp1,ncol=5,group.by="seurat_clusters")
dev.off()







#### CT subset4
pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_ClusterGroup.pdf",width=25,height=16)
VlnPlot(RNAprocess_DN1WT, features=usekeygene,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_HSC.pdf",width=20,height=8)
VlnPlot(RNAprocess_DN1WT,features=keygene_HSC,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_Tlineage.pdf",width=20,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_Tlineage,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_Blineage.pdf",width=20,height=8)
VlnPlot(RNAprocess_DN1WT,features=keygene_Blineage,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_Survival.pdf",width=20,height=4)
VlnPlot(RNAprocess_DN1WT,features=keygene_Survival,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_Dendritic.pdf",width=20,height=8)
VlnPlot(RNAprocess_DN1WT,features=keygene_Dendritic,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_Myeloid.pdf",width=20,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_Myeloid,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_InnateLymphoid.pdf",width=20,height=8)
VlnPlot(RNAprocess_DN1WT,features=keygene_InnateLymphoid,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_other.pdf",width=20,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_other,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_sp1.pdf",width=20,height=4)
VlnPlot(RNAprocess_DN1WT,features=keygene_sp1,ncol=5,group.by="CT_4subset")
dev.off()

pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_ETPproliferative.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_ETPproliferative,ncol=4,group.by="CT_4subset")
dev.off()
pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_ETPquiescent.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_ETPquiescent,ncol=4,group.by="CT_4subset")
dev.off()
pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_nonTpotential.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_nonTpotential,ncol=4,group.by="CT_4subset")
dev.off()
pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_TILClike.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_TILClike,ncol=4,group.by="CT_4subset")
dev.off()




pdf(file="WTonly_boxplot_4subset_knowledgeKeyGeneExp_ETPproliferative.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_ETPproliferative,ncol=4,group.by="seurat_clusters")
dev.off()
pdf(file="WTonly_boxplot_4subset_knowledgeKeyGeneExp_ETPquiescent.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_ETPquiescent,ncol=4,group.by="seurat_clusters")
dev.off()
pdf(file="WTonly_boxplot_4subset_knowledgeKeyGeneExp_nonTpotential.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_nonTpotential,ncol=4,group.by="seurat_clusters")
dev.off()
pdf(file="WTonly_boxplot_4subset_knowledgeKeyGeneExp_TILClike.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_TILClike,ncol=4,group.by="seurat_clusters")
dev.off()
pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_ETPproliferative.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_ETPproliferative,ncol=4,group.by="CT_4subset")
dev.off()
pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_ETPquiescent.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_ETPquiescent,ncol=4,group.by="CT_4subset")
dev.off()
pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_nonTpotential.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_nonTpotential,ncol=4,group.by="CT_4subset")
dev.off()
pdf(file="WTonly_boxplotCT_4subset_knowledgeKeyGeneExp_TILClike.pdf",width=16,height=12)
VlnPlot(RNAprocess_DN1WT,features=keygene_TILClike,ncol=4,group.by="CT_4subset")
dev.off()



### keygene heatmap
usegene <- c(keygene_ETPproliferative,keygene_ETPquiescent,keygene_nonTpotential,keygene_TILClike)
usecell1 <- names(RNAprocess_DN1WT$CT_4subset)[which(RNAprocess_DN1WT$CT_4subset=="ETPs_proliferative")]
usecell2 <- names(RNAprocess_DN1WT$CT_4subset)[which(RNAprocess_DN1WT$CT_4subset=="ETPs_quiescent")]
usecell3 <- names(RNAprocess_DN1WT$CT_4subset)[which(RNAprocess_DN1WT$CT_4subset=="nonT_potential")]
usecell4 <- names(RNAprocess_DN1WT$CT_4subset)[which(RNAprocess_DN1WT$CT_4subset=="T_ILCs")]

#expmat_scale_SCT <- RNAprocess_DN1WT_allgeneSCT@assays$SCT@scale.data
#expmat_scale_RNA <- RNAprocess_DN1WT@assays$RNA@scale.data
#
#
#SCTmat <- expmat_scale_SCT[usegene, c(usecell1,usecell2,usecell3,usecell4)]
pdf(file="SCTheatmap_CT4subsetKeygenes_clusterGrouping.pdf",width=10,height=10)
DoHeatmap(RNAprocess_DN1WT_allgeneSCT, features=usegene, cells=c(usecell1,usecell2,usecell3,usecell4),group.by="seurat_clusters")
dev.off()

pdf(file="SCTheatmap_CT4subsetKeygenes_CT4subsetGrouping.pdf",width=10,height=10)
DoHeatmap(RNAprocess_DN1WT_allgeneSCT, features=usegene, cells=c(usecell1,usecell2,usecell3,usecell4),group.by="CT_4subset")
dev.off()

uesgene1 <- c("Notch1","Hes1","Hhex","Cd24a","Cdk6","Etv6","Bcl11a","Egr1","Jun","Sox4","Myc","Spi1")
pdf(file="SCTheatmap_newsubset1Keygenes_clusterGrouping.pdf",width=10,height=10)
DoHeatmap(RNAprocess_DN1WT_allgeneSCT, features=uesgene1, cells=c(usecell1,usecell2,usecell3,usecell4),group.by="seurat_clusters")
dev.off()

pdf(file="SCTheatmap_newsubset1Keygenes_CT4subsetGrouping.pdf",width=10,height=10)
DoHeatmap(RNAprocess_DN1WT_allgeneSCT, features=uesgene1, cells=c(usecell1,usecell2,usecell3,usecell4),group.by="CT_4subset")
dev.off()

keygene_CT4subset <- c("Notch1","Hes1","Hhex","Cd24a","Cdk6","Etv6","Bcl11a","Egr1","Jun","Sox4","Myc","Spi1","Aurkb","Birc5","Bub1b","Ccna2","Ccnb1","Ccnb2","Hmgb1","Cdk1","Mki67","Dnmt1","Ezh2","Lig1","Gapdh","Ldha","Pgam1","Hk3","Cebpa","Cebpb","Cebpd","Fos","Lmo2","Irf8","Spi1","Zeb2","Cd3e","Cd3d","Lck","Id3","Tox","Ets1","Zap70","Il2rg","Thy1","Il7r","Id2","Gata3")

neworder <- c(0,1,2,3,5,7,4,9,6,8,10,11,12,13)
RNAprocess_DN1WT_allgeneSCT$seurat_clusters_reorder <- factor(RNAprocess_DN1WT_allgeneSCT$seurat_clusters, levels = neworder)
pdf(file="SCTheatmap_newsubset1Keygenes_doubleLabel_material.pdf",width=10,height=10)
DoHeatmap(RNAprocess_DN1WT_allgeneSCT, features=keygene_CT4subset,group.by="seurat_clusters_reorder")
dev.off()

neworder_CT4subset <- c("ETPs_quiescent","ETPs_proliferative","nonT_potential","T_ILCs")
RNAprocess_DN1WT_allgeneSCT$CT4_subset_reorder <- factor(RNAprocess_DN1WT_allgeneSCT$CT_4subset, levels = neworder_CT4subset)
pdf(file="SCTheatmap_newsubset1Keygenes_doubleLabel_material2.pdf",width=10,height=10)
DoHeatmap(RNAprocess_DN1WT_allgeneSCT, features=keygene_CT4subset,group.by="CT4_subset_reorder")
dev.off()

### saveRDS


#saveRDS(RNAprocess_DN1WT, file="RNAprocess_DN1WT.rds")
#saveRDS(RNAprocess_DN1WT_allgeneSCT, file="RNAprocess_DN1WT_allgeneSCT.rds")
RNAprocess_DN1WT <- readRDS(file="RNAprocess_DN1WT.rds")
RNAprocess_DN1WT_allgeneSCT <- readRDS(file="RNAprocess_DN1WT_allgeneSCT.rds")


ETPp_exp <- apply(RNAprocess_DN1WT@assays$RNA@data[,which(RNAprocess_DN1WT$CT_5subset == "ETPs_proliferative")],1,mean)
ETPq_exp <- apply(RNAprocess_DN1WT@assays$RNA@data[,which(RNAprocess_DN1WT$CT_5subset == "ETPs_quiescent")],1,mean)
nonT_exp <- apply(RNAprocess_DN1WT@assays$RNA@data[,which(RNAprocess_DN1WT$CT_5subset == "nonT_potential")],1,mean)
TILCGR_exp <- apply(RNAprocess_DN1WT@assays$RNA@data[,which(RNAprocess_DN1WT$CT_5subset == "T_ILCs_GRhi")],1,mean)
TILCII_exp <- apply(RNAprocess_DN1WT@assays$RNA@data[,which(RNAprocess_DN1WT$CT_5subset == "T_ILCs_IIhi")],1,mean)

CTexp <- cbind(ETPp_exp,ETPq_exp,nonT_exp,TILCGR_exp,TILCII_exp)
meanCT <- apply(CTexp,1,mean) 
sdCT <-  apply(CTexp,1,sd) 
cvCT <- sdCT/meanCT

MKGct5 <- read.table("MKGct5_DN1WT.txt",row.names=1,header=T)
MKGene <- unique(MKGct5[,"gene"])

pdf("CITE_WTonly_mean_cv_CTlevelEXP_scatter.pdf")
plot(meanCT,cvCT,pch=".",main="WT CTlevel exp")
points(meanCT[MKGene],cvCT[MKGene],pch=".",col="red")
abline(h=0.1, col="blue")
abline(v=0.5,col="blue")
dev.off()

stable_gene <- setdiff(names(meanCT)[which(meanCT >= 0.5 & cvCT < 0.1)],MKGene)
write.table(stable_gene,file="CITE_WTonly_stable_gene.txt",row.names=F,col.names=F,sep="\t",quote=F)


ETPq_cells <- names(RNAprocess_DN1WT$CT_4subset)[which(RNAprocess_DN1WT$CT_4subset ==  "ETPs_quiescent")]

ETPq_exp <- apply(RNAprocess_DN1WT@assays$RNA[,ETPq_cells],1,mean)
write.table(ETPq_exp, file="ETPq_RNAexp.txt",row.names=T,col.names=F,sep="\t",quote=F)









### Cxcr5 on/off

tmpdata <- RNAprocess_D2D3@assays$RNA@scale.data
RNAprocess_D2D3_CXCR5 <- rep("off",ncol(tmpdata))
names(RNAprocess_D2D3_CXCR5) <- colnames(tmpdata)
RNAprocess_D2D3_CXCR5[which(tmpdata["Cxcr5",]>0)] <- "on"
RNAprocess_D2D3$CXCR5on <- RNAprocess_D2D3_CXCR5

confmat_D2D3_CXCR5on <- confmat(cbind(as.character(RNAprocess_D2D3$seurat_clusters),RNAprocess_D2D3$CXCR5on))

tmpdata <- confmat_D2D3_CXCR5on[,"on"] / (confmat_D2D3_CXCR5on[,"on"] + confmat_D2D3_CXCR5on[,"off"])
idxTmp <- names(RNAprocess_D2D3$seurat_clusters)
CXCR5onoff <- rep("NA",length(idxTmp))
names(CXCR5onoff) <- idxTmp
for(i in 0:13){
	CXCR5onoff[which(RNAprocess_D2D3$seurat_clusters == i)] <- tmpdata[as.character(i)]
}
RNAprocess_D2D3$CXCR5onPercent <- as.numeric(CXCR5onoff)

pdf(file="D2D3_UMAP_features_CXCR5onPercent.pdf",width=4,height=4)
p1 <- FeaturePlot(RNAprocess_D2D3,features="CXCR5onPercent",ncol=1)+ ggtitle("D2D3 all")
p1
dev.off()

########## QC plots


p1 <- DimPlot(RNAprocess_d2, label = TRUE,group.by="ADT") + ggtitle("D2")
p2 <- DimPlot(RNAprocess_d2clean, label = TRUE,group.by="ADT") + ggtitle("D2 no naive")
p3 <- DimPlot(RNAprocess_d3, label = TRUE,group.by="ADT") + ggtitle("D3")
p4 <- DimPlot(RNAprocess_d5, label = TRUE,group.by="ADT") + ggtitle("D5")
p5 <- DimPlot(RNAprocess_d30, label = TRUE,group.by="ADT") + ggtitle("D30")

pdf(file="UMAP_ADTcolor.pdf",width=20,height=4)
p1|p2|p3|p4|p5
#DimPlot(RNAprocess_d2, label = TRUE,group.by="ADT") + ggtitle("D2")
#DimPlot(RNAprocess_d2clean, label = TRUE,group.by="ADT") + ggtitle("D2 no naive")
#DimPlot(RNAprocess_d3, label = TRUE,group.by="ADT") + ggtitle("D3")
#DimPlot(RNAprocess_d5, label = TRUE,group.by="ADT") + ggtitle("D5")
#DimPlot(RNAprocess_d30, label = TRUE,group.by="ADT") + ggtitle("D30")
dev.off()


p1 <- DimPlot(RNAprocess_d2, label = TRUE,group.by="seurat_clusters") + ggtitle("D2")
p2 <- DimPlot(RNAprocess_d2clean, label = TRUE,group.by="seurat_clusters") + ggtitle("D2 no naive")
p3 <- DimPlot(RNAprocess_d3, label = TRUE,group.by="seurat_clusters") + ggtitle("D3")
p4 <- DimPlot(RNAprocess_d5, label = TRUE,group.by="seurat_clusters") + ggtitle("D5")
p5 <- DimPlot(RNAprocess_d30, label = TRUE,group.by="seurat_clusters") + ggtitle("D30")

pdf(file="UMAP_Clustercolor.pdf",width=20,height=4)
p1|p2|p3|p4|p5
dev.off()

pdf(file="UMAP_keygeneExp.pdf",width=4,height=12)
FeaturePlot(RNAprocess_d2,features=c("Tcf7","Bcl6","Cxcr5"),ncol=1)
FeaturePlot(RNAprocess_d2clean,features=c("Tcf7","Bcl6","Cxcr5"),ncol=1)
FeaturePlot(RNAprocess_d3,features=c("Tcf7","Bcl6","Cxcr5"),ncol=1)
FeaturePlot(RNAprocess_d5,features=c("Tcf7","Bcl6","Cxcr5"),ncol=1)
FeaturePlot(RNAprocess_d30,features=c("Tcf7","Bcl6","Cxcr5"),ncol=1)
dev.off()

pdf(file="barplot_keygeneExp_ADTgroup.pdf",width=4,height=16)
VlnPlot(RNAprocess_d2, features=c("Tcf7","Bcl6","Cxcr5","Pdcd1"),ncol=1,group.by="ADT")
VlnPlot(RNAprocess_d2clean, features=c("Tcf7","Bcl6","Cxcr5","Pdcd1"),ncol=1,group.by="ADT")
VlnPlot(RNAprocess_d3, features=c("Tcf7","Bcl6","Cxcr5","Pdcd1"),ncol=1,group.by="ADT")
VlnPlot(RNAprocess_d5, features=c("Tcf7","Bcl6","Cxcr5","Pdcd1"),ncol=1,group.by="ADT")
VlnPlot(RNAprocess_d30, features=c("Tcf7","Bcl6","Cxcr5","Pdcd1"),ncol=1,group.by="ADT")
dev.off()

keygene <- "Tcf7"
p1 <- VlnPlot(RNAprocess_d2, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D2 ",keygene))
p2 <- VlnPlot(RNAprocess_d2clean, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D2 no naive ",keygene))
p3 <- VlnPlot(RNAprocess_d3, features=c(keygene),group.by="seurat_clusters") + ggtitle(paste0("D3 ",keygene))
p4 <- VlnPlot(RNAprocess_d5, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D5 ",keygene))
p5 <- VlnPlot(RNAprocess_d30, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D30 ",keygene))
pdf(file=paste0("barplot_",keygene,"Exp_ClusterGroup.pdf"),width=20,height=4)
p1|p2|p3|p4|p5
dev.off()	

keygene <- "Bcl6"
p1 <- VlnPlot(RNAprocess_d2, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D2 ",keygene))
p2 <- VlnPlot(RNAprocess_d2clean, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D2 no naive ",keygene))
p3 <- VlnPlot(RNAprocess_d3, features=c(keygene),group.by="seurat_clusters") + ggtitle(paste0("D3 ",keygene))
p4 <- VlnPlot(RNAprocess_d5, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D5 ",keygene))
p5 <- VlnPlot(RNAprocess_d30, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D30 ",keygene))
pdf(file=paste0("barplot_",keygene,"Exp_ClusterGroup.pdf"),width=20,height=4)
p1|p2|p3|p4|p5
dev.off()	

keygene <- "Cxcr5"
p1 <- VlnPlot(RNAprocess_d2, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D2 ",keygene))
p2 <- VlnPlot(RNAprocess_d2clean, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D2 no naive ",keygene))
p3 <- VlnPlot(RNAprocess_d3, features=c(keygene),group.by="seurat_clusters") + ggtitle(paste0("D3 ",keygene))
p4 <- VlnPlot(RNAprocess_d5, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D5 ",keygene))
p5 <- VlnPlot(RNAprocess_d30, features=c(keygene),group.by="seurat_clusters")+ ggtitle(paste0("D30 ",keygene))
pdf(file=paste0("barplot_",keygene,"Exp_ClusterGroup.pdf"),width=20,height=4)
p1|p2|p3|p4|p5
dev.off()	



# Get unique list of genes and cells across all matrices
#all_genes <- sort(unique(c(rownames(RNA_DN1r1), rownames(RNA_DN1r2), rownames(RNA_DN1r3), rownames(RNA_DN1r4))))
##all_ADTs <- 
#all_cells <- sort(unique(c(colnames(matrix1), colnames(matrix2), colnames(matrix3), colnames(matrix4))))
#
## Function to align a matrix to the master list of genes and cells
#align_matrix <- function(mat, all_genes, all_cells) {
#  missing_genes <- setdiff(all_genes, rownames(mat))
#  missing_cells <- setdiff(all_cells, colnames(mat))
#  
#  if(length(missing_genes) > 0) {
#    mat <- rbind(mat, matrix(0, nrow=length(missing_genes), ncol=ncol(mat), dimnames=list(missing_genes, colnames(mat))))
#  }
#  
#  if(length(missing_cells) > 0) {
#    mat <- cbind(mat, matrix(0, nrow=nrow(mat), ncol=length(missing_cells), dimnames=list(rownames(mat), missing_cells)))
#  }
#  
#  return(mat[all_genes, all_cells])
#}
#
## Align all matrices
#matrix1 <- align_matrix(matrix1, all_genes, all_cells)
#matrix2 <- align_matrix(matrix2, all_genes, all_cells)
#matrix3 <- align_matrix(matrix3, all_genes, all_cells)
#matrix4 <- align_matrix(matrix4, all_genes, all_cells)
#
## Now you can sum them
#combinedMatrix <- matrix1 + matrix2 + matrix3 + matrix4


