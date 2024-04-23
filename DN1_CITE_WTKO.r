require("Matrix")
library(Seurat)
library(caret)
library(harmony)
library(gplots)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(gridExtra)
library(ggplot2)

#library(monocle3)
#library(slingshot)
#library(SCENIC)
### p1 as example
setwd("/sfs/ceph/standard/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/CITE/allconds")
## readADT
readADT <- function(inname){
	ADTfolder <- paste0("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/CITE/preprocess/DN1/results/",inname,"_ADT/")
	ADT <- Matrix::readMM(paste0(ADTfolder,"matrix.mtx.gz"))
	cellinfo <- read.table(paste0(ADTfolder,"barcodes.tsv.gz"))
	Antibody <- read.table(paste0(ADTfolder,"features.tsv.gz"))
	row.names(ADT) <- Antibody[,1]#row.names(peakinfo)
	colnames(ADT) <- paste0(cellinfo[,1],"-1")
	return(ADT)
}
ADT_DN1 <- readADT("DN1")

#ADT_DN1r1 <- readADT("DN1r1")
#ADT_DN1r2 <- readADT("DN1r2")
#ADT_DN1r3 <- readADT("DN1r3")
#ADT_DN1r4 <- readADT("DN1r4")




readRNA <- function(inname){
	RNAfolder <- paste0("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/CITE/preprocess/DN1/results/",inname,"_RNA/")
	RNA.data <- Read10X(data.dir = paste0(RNAfolder))
	RNA <- CreateSeuratObject(counts = RNA.data, project = inname)
	return(RNA)
}

RNA_DN1 <- readRNA("DN1")
#RNA_DN1r1 <- readRNA("DN1r1")
#RNA_DN1r2 <- readRNA("DN1r2")
#RNA_DN1r3 <- readRNA("DN1r3")
#RNA_DN1r4 <- readRNA("DN1r4")




#pdf(file="ADTcount_scatter.pdf",width=9,height=9)
#par(mfrow=c(2,2),mar=c(4,4,2,2))
#plot(ADT_d2["WT-ACCCACCAGTAAGAC",coCell_d2], ADT_d2["Y181F-GGTCGAGAGCATTCA",coCell_d2], pch=".",main="D2 ADTcount", xlab="WT",ylab="Y181F" )
#plot(ADT_d3["WT-AAAGCATTCTTCACG",coCell_d3], ADT_d3["Y181F-CTTTGTCTTTGTGAG",coCell_d3], pch=".", main="D3 ADTcount", xlab="WT",ylab="Y181F" )
#plot(ADT_d5["Inf-ACCCACCAGTAAGAC",coCell_d5], ADT_d5["Imm-CTTGCCGCATGTCAT",coCell_d5], pch=".", main="D5 ADTcount", xlab="Inf",ylab="Imm" )
#plot(ADT_d30["Inf-ACCCACCAGTAAGAC",coCell_d30], ADT_d30["Imm-GGTCGAGAGCATTCA",coCell_d30], pch=".", main="D30 ADTcount", xlab="Inf",ylab="Imm" )
#dev.off()
## readRNA

line_gt0 <- function(inline){
	return(length(which(inline > 0)))
}
line_maxIdx <- function(inline){
	return(which(inline == max(inline)))
}
trimRN <- function(inname){
	return(strsplit(inname, "-")[[1]][1])
}
QCpass_Cell <- function(ADT1,RNA1, totalADT_cutoff, unmapRate_cutoff, maxP_cutoff){
	totalADT1 <- apply(ADT1,2,sum)
	unmapRate1 <- ADT1["unmapped",]/apply(ADT1,2,sum)
	maxPercent1 <- apply(ADT1[setdiff(rownames(ADT1),"unmapped"),],2,max)/apply(ADT1[setdiff(rownames(ADT1),"unmapped"),],2,sum)
	coverGene <- apply(RNA1@assays$RNA@counts, 2, line_gt0)

    useCell <-  intersect(intersect(intersect(names(totalADT1)[which(totalADT1 >= 1000)],
					                		  names(unmapRate1)[which(unmapRate1 <= unmapRate_cutoff)]),
				          			names(maxPercent1)[which(maxPercent1 >= maxP_cutoff)]),
						  names(coverGene)[which(coverGene >= 1000)])
    useADT <- ADT1[,useCell]
    RNshort <- unlist(lapply(rownames(useADT),trimRN))
    CTassign <- RNshort[apply(useADT,2,line_maxIdx)]
    names(CTassign) <- useCell
	return(CTassign)
}

#QC_Cell <- function(ADT1,RNA1){
#	totalADT1 <- apply(ADT1,2,sum)
#	unmapRate1 <- ADT1["unmapped",]/apply(ADT1,2,sum)
#	maxPercent1 <- apply(ADT1[setdiff(rownames(ADT1),"unmapped"),],2,max)/apply(ADT1[setdiff(rownames(ADT1),"unmapped"),],2,sum)
#	cbdata <- cbind(totalADT1, unmapRate1, maxPercent1, rep(0, length(totalADT1)))
#	colnames(cbdata) <- c("totalADT","unmapRate","maxPercent","coverGene")
#
#	coverGene <- apply(RNA1@assays$RNA@counts, 2, line_gt0)
#
#	commonCell <- intersect(rownames(cbdata), names(coverGene))
#	cbdata[commonCell, "coverGene"] <- coverGene[commonCell]
#
#    useCell <-  intersect(intersect(intersect(names(totalADT1)[which(totalADT1 >= 1000)],
#					                		  names(unmapRate1)[which(unmapRate1 <= 0.2)]),
#				          			names(maxPercent1)[which(maxPercent1 >= 0.8)]),
#						  names(coverGene)[which(coverGene >= 1000)])
#    useADT <- ADT1[,useCell]
#    RNshort <- unlist(lapply(rownames(useADT),trimRN))
#    CTassign <- RNshort[apply(useADT,2,line_maxIdx)]
#    names(CTassign) <- useCell
#	return(list(CTassign,cbdata))
#}
#
#


CTassign_DN1 <- QCpass_Cell(ADT_DN1, RNA_DN1, 1000, 0.1, 0.9)
#CTassign_DN1_v2 <- QCpass_Cell(ADT_DN1, RNA_DN1, 1000, 0.2, 0.8)
#CTassign_DN1_v3 <- QCpass_Cell(ADT_DN1, RNA_DN1, 1000, 0.3, 0.7)

write.table(CTassign_DN1,file="CTassign_DN1.txt",row.names=T,col.names=F,sep="\t",quote=F)
CTassign_DN1 <- read.table("CTassign_DN1.txt",row.names=1,header=F)
#a <- read.table("CTassign_DN1.txt",row.names=1,header=F)

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

RNAprocess_DN1 <- RNAprocess(RNA_DN1[,rownames(CTassign_DN1)])
RNAprocess_DN1$ADT <- CTassign_DN1


WT_cell <- rownames(CTassign_DN1)[which(CTassign_DN1[,1] == "WT")]
Tcf1KO_cell <- rownames(CTassign_DN1)[which(CTassign_DN1[,1] == "Tcf1KO")]
Tcf1Lef1KO_cell <- rownames(CTassign_DN1)[which(CTassign_DN1[,1] == "Tcf1Lef1KO")]
Tle134KO_cell <- rownames(CTassign_DN1)[which(CTassign_DN1[,1] == "Tle134KO")]

### UMAP + features
p1 <- DimPlot(RNAprocess_DN1, label = TRUE,group.by="ADT") + ggtitle("ADT")
p2 <- DimPlot(RNAprocess_DN1, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")

pdf(file="UMAP_features.pdf",width=10,height=4)
grid.arrange(p1,p2, ncol=2,nrow=1)
dev.off()


RNAprocess_DN1_res0.5 <- FindClusters(RNAprocess_DN1, resolution = 0.5, verbose = FALSE)
RNAprocess_DN1_res0.6 <- FindClusters(RNAprocess_DN1, resolution = 0.6, verbose = FALSE)
RNAprocess_DN1_res0.7 <- FindClusters(RNAprocess_DN1, resolution = 0.7, verbose = FALSE)
RNAprocess_DN1_res0.8 <- FindClusters(RNAprocess_DN1, resolution = 0.8, verbose = FALSE)
RNAprocess_DN1_res0.9 <- FindClusters(RNAprocess_DN1, resolution = 0.9, verbose = FALSE)
RNAprocess_DN1_res1.0 <- FindClusters(RNAprocess_DN1, resolution = 1.0, verbose = FALSE)
RNAprocess_DN1_res1.1 <- FindClusters(RNAprocess_DN1, resolution = 1.1, verbose = FALSE)
RNAprocess_DN1_res1.2 <- FindClusters(RNAprocess_DN1, resolution = 1.2, verbose = FALSE)
RNAprocess_DN1_res1.3 <- FindClusters(RNAprocess_DN1, resolution = 1.3, verbose = FALSE)
RNAprocess_DN1_res1.4 <- FindClusters(RNAprocess_DN1, resolution = 1.4, verbose = FALSE)
RNAprocess_DN1_res1.5 <- FindClusters(RNAprocess_DN1, resolution = 1.5, verbose = FALSE)
RNAprocess_DN1_res1.6 <- FindClusters(RNAprocess_DN1, resolution = 1.6, verbose = FALSE)
RNAprocess_DN1_res1.7 <- FindClusters(RNAprocess_DN1, resolution = 1.7, verbose = FALSE)
RNAprocess_DN1_res1.8 <- FindClusters(RNAprocess_DN1, resolution = 1.8, verbose = FALSE)
RNAprocess_DN1_res1.9 <- FindClusters(RNAprocess_DN1, resolution = 1.9, verbose = FALSE)
RNAprocess_DN1_res2.0 <- FindClusters(RNAprocess_DN1, resolution = 2.0, verbose = FALSE)

p0.5 <- DimPlot(RNAprocess_DN1_res0.5, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.5")
p0.6 <- DimPlot(RNAprocess_DN1_res0.6, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.6")
p0.7 <- DimPlot(RNAprocess_DN1_res0.7, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.7")
p0.8 <- DimPlot(RNAprocess_DN1_res0.8, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.8")
p0.9 <- DimPlot(RNAprocess_DN1_res0.9, label = TRUE,group.by="seurat_clusters") + ggtitle("res0.9")
p1.0 <- DimPlot(RNAprocess_DN1_res1.0, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.0")
p1.1 <- DimPlot(RNAprocess_DN1_res1.1, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.1")
p1.2 <- DimPlot(RNAprocess_DN1_res1.2, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.2")
p1.3 <- DimPlot(RNAprocess_DN1_res1.3, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.3")
p1.4 <- DimPlot(RNAprocess_DN1_res1.4, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.4")
p1.5 <- DimPlot(RNAprocess_DN1_res1.5, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.5")
p1.6 <- DimPlot(RNAprocess_DN1_res1.6, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.6")
p1.7 <- DimPlot(RNAprocess_DN1_res1.7, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.7")
p1.8 <- DimPlot(RNAprocess_DN1_res1.8, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.8")
p1.9 <- DimPlot(RNAprocess_DN1_res1.9, label = TRUE,group.by="seurat_clusters") + ggtitle("res1.9")
p2.0 <- DimPlot(RNAprocess_DN1_res2.0, label = TRUE,group.by="seurat_clusters") + ggtitle("res2.0")

#pdf(file="UMAP_diffClusterResolution.pdf",width=30,height=30)
#p0.5|p0.6|p0.7|p0.8|p0.9|p1.0|p1.1|p1.2|p1.3|p1.4|p1.5|p1.6|p1.7|p1.8|p1.9|p2.0
#dev.off()

pdf(file="UMAP_diffClusterResolution.pdf",width=30,height=30)
grid.arrange(p0.5,p0.6,p0.7,p0.8,p0.9,p1.0,p1.1,p1.2,p1.3,p1.4,p1.5,p1.6,p1.7,p1.8,p1.9,p2.0, ncol=4,nrow=4)
dev.off()

### MTDNA, to be done

mtDNA <- read.table("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/CITE/preprocess/DN1/RNA/mtDNA/DN1_mtDNAreadnum.txt",row.names=1)

mtDNA_use <- as.numeric(mtDNA[names(RNAprocess_DN1$nCount_RNA),1])
names(mtDNA_use) <- names(RNAprocess_DN1$nCount_RNA)
RNAprocess_DN1$mtDNAreads <- mtDNA_use

idx <- names(mtDNA_use)
mtDNApercent <- mtDNA_use[idx] / (RNAprocess_DN1$nCount_RNA[idx]+mtDNA_use[idx])
names(mtDNApercent) <- idx
RNAprocess_DN1$MTpercent <- mtDNApercent
RNAprocess_DN1_allgeneSCT$MTpercent <- RNAprocess_DN1$MTpercent 


p1 <- DimPlot(RNAprocess_DN1, label = TRUE,group.by="ADT") + ggtitle("ADT")
p2 <- DimPlot(RNAprocess_DN1, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")
p3 <- FeaturePlot(RNAprocess_DN1,features="MTpercent")+ ggtitle("chrM%")

pdf(file="UMAP_features_mtDNA.pdf",width=15,height=4)
grid.arrange(p1,p2,p3, ncol=3,nrow=1)
dev.off()


p1<-VlnPlot(RNAprocess_DN1, features="MTpercent",group.by="ADT")
p2<-VlnPlot(RNAprocess_DN1, features="MTpercent",group.by="seurat_clusters")

pdf(file="boxplot_mtDNA.pdf",width=10,height=4)
grid.arrange(p1,p2, ncol=2,nrow=1)
dev.off()


confmat <- function(inmat){
	feature1 <- unique(inmat[,1])
	feature2 <- unique(inmat[,2])
	outmat <- matrix(rep(0, length(feature1)*length(feature2)),nrow=length(feature1),ncol=length(feature2))
	rownames(outmat) <- feature1
	colnames(outmat) <- feature2
	for(i in 1:nrow(inmat)){
		outmat[inmat[i,1], inmat[i,2]] <- outmat[inmat[i,1], inmat[i,2]] + 1
	}
	return(outmat)
}

confmat_ADT_cluster <- confmat(cbind(as.character(RNAprocess_DN1$seurat_clusters),RNAprocess_DN1$ADT))




### MKG
MKGcluster_DN1 <- FindAllMarkers(RNAprocess_DN1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MKGclusterFilter_DN1 <- MKGcluster_DN1[which(MKGcluster_DN1[,2] >= log2(2) & MKGcluster_DN1[,5] < 0.01),]
write.table(MKGclusterFilter_DN1, file="MKGclusterFilter_DN1.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(MKGcluster_DN1, file="MKGcluster_DN1.txt",row.names=T,col.names=T,sep="\t",quote=F)


outdata <- c()
for(C in sort(unique(MKGclusterFilter_DN1[,"cluster"])) ){
	tmpdata <- MKGclusterFilter_DN1[which(MKGclusterFilter_DN1[,"cluster"]==C),]
	outdata <- rbind(outdata, c("DN1",paste0("C",C),tmpdata[order(tmpdata[,"avg_log2FC"],decreasing=T),][1:5,"gene"]))
}

alltop5MKG <- c()
for(this_cluster in 0:21 ){
	#print(this_cluster)
	this_MKGdata <- MKGclusterFilter_DN1[which(MKGclusterFilter_DN1[,"cluster"] == this_cluster),]
	if(nrow(this_MKGdata)>0){
        this_MKGdata_top5 <- this_MKGdata[order(this_MKGdata[,"avg_log2FC"],decreasing=T),][1:min(5,nrow(this_MKGdata)),]
        tmpobj <- RNAprocess_DN1
        tmpobj$target_cluster <- rep(0, ncol(tmpobj))
        tmpobj$target_cluster[which(tmpobj$seurat_clusters == this_cluster)] <- 1
        
        #pdf(file=paste0("UMAP_keygeneExp_C",this_cluster,".pdf"),width=24,height=4)
        FeaturePlot(tmpobj,features=c("target_cluster",as.vector(this_MKGdata_top5[,"gene"])),ncol=6)
        #dev.off()
        ggsave(paste0("UMAP_keygeneExp_C",this_cluster,".pdf"),width=24,height=4)
        alltop5MKG <- c(alltop5MKG, as.vector(this_MKGdata_top5[,"gene"]))		
	}
}




RNAprocess_DN1_allgeneSCT <-  SCTransform(RNAprocess_DN1, assay = "RNA", verbose = FALSE, ,variable.features.n = 24538)

usegene <- intersect(alltop5MKG, rownames(RNAprocess_DN1_allgeneSCT@assays$RNA@scale.data))

expmat_scale_SCT <- RNAprocess_DN1_allgeneSCT@assays$SCT@scale.data
expmat_scale_RNA <- RNAprocess_DN1@assays$RNA@scale.data

aveExpMat_scale_RNA <-  matrix(rep(0, length(usegene)*length(unique(RNAprocess_DN1$seurat_clusters)) ),nrow=length(usegene))
rownames(aveExpMat_scale_RNA) <- usegene
colnames(aveExpMat_scale_RNA) <- paste0("C",0:(length(unique(RNAprocess_DN1$seurat_clusters)) -1))
medianExpMat_scale_RNA <- aveExpMat_scale_RNA
for(i in 0:(length(unique(RNAprocess_DN1$seurat_clusters)) -1)){
	usecell <- names(RNAprocess_DN1$seurat_clusters)[which(RNAprocess_DN1$seurat_clusters==i)]
	tempMat <- expmat_scale_RNA[usegene,usecell]
	aveExpMat_scale_RNA[, paste0("C",i)] <- apply(tempMat, 1,mean)
	medianExpMat_scale_RNA[, paste0("C",i)] <- apply(tempMat, 1,median)
}

aveExpMat_scale_SCT <-  matrix(rep(0, length(usegene)*length(unique(RNAprocess_DN1$seurat_clusters)) ),nrow=length(usegene))
rownames(aveExpMat_scale_SCT) <- usegene
colnames(aveExpMat_scale_SCT) <- paste0("C",0:(length(unique(RNAprocess_DN1$seurat_clusters)) -1))
medianExpMat_scale_SCT <- aveExpMat_scale_SCT
for(i in 0:(length(unique(RNAprocess_DN1$seurat_clusters)) -1)){
	usecell <- names(RNAprocess_DN1$seurat_clusters)[which(RNAprocess_DN1$seurat_clusters==i)]
	tempMat <- expmat_scale_SCT[usegene,usecell]
	aveExpMat_scale_SCT[, paste0("C",i)] <- apply(tempMat, 1,mean)
	medianExpMat_scale_SCT[, paste0("C",i)] <- apply(tempMat, 1,median)
}


usecolor <- c("blue","white","red")
ColorRamp <- colorRampPalette(usecolor, bias=1)(101)   #color list
pdf(file="keygeneExpRNA_cluster_heatmap.pdf")
heatmap.2(aveExpMat_scale_RNA, trace="none", col=ColorRamp,main="aveExp scaleRNA",cexRow=0.6)
dev.off()

pdf(file="keygeneExpSCT_cluster_heatmap.pdf")
heatmap.2(aveExpMat_scale_SCT, trace="none", col=ColorRamp,main="aveExp scaleSCT",cexRow=0.6)
dev.off()



######### KO vs WT DEG in every cluster
outdata <- c()
outdata_raw <- c()
for(this_cluster in unique(RNAprocess_DN1$seurat_clusters)){
	thisdata <- RNAprocess_DN1[,which(RNAprocess_DN1$seurat_clusters == this_cluster)]
	WT_cell <- as.numeric(table(thisdata$ADT)["WT"])
	if(WT_cell < 20 || is.na(WT_cell)){
		next
	}
	for(KOtype in c("Tcf1KO","Tcf1Lef1KO", "Tle134KO")){
		KO_cell <- as.numeric(table(thisdata$ADT)[KOtype])
		if(KO_cell < 20 || is.na(KO_cell)){
			next
		}
		print(paste0("c",this_cluster,"_",KOtype,"vsWT"))
		this_diff <- FindMarkers(thisdata, ident.1=KOtype,ident.2="WT",group.by="ADT")
		this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
		outdata <- rbind(outdata, cbind(this_diff_filter, rownames(this_diff_filter),rep(paste0("c",this_cluster,"_",KOtype,"vsWT"),nrow(this_diff_filter))))
		outdata_raw <- rbind(outdata_raw, cbind(this_diff, rownames(this_diff),rep(paste0("c",this_cluster,"_",KOtype,"vsWT"),nrow(this_diff))))
#		write.table(this_diff_filter[order(this_diff_filter[,2]),],file=paste0("clusterLevel_ADTdiff/c",this_clsuter,"_",KOtype,"vsWT_DEGfc2q01.txt"),row.names=T,col.names=T,sep="\t",quote=F)
	}

}
colnames(outdata) <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","gene","cmp")
colnames(outdata_raw) <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","gene","cmp")
write.table(outdata, file="clusterLevel_ADTdiff_DEGfc2q01.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(outdata, file="clusterLevel_ADTdiff_DEGraw.txt",row.names=F,col.names=T,sep="\t",quote=F)





### saveRDS


#saveRDS(RNAprocess_DN1, file="RNAprocess_DN1.rds")
#saveRDS(RNAprocess_DN1_allgeneSCT, file="RNAprocess_DN1_allgeneSCT.rds")
RNAprocess_DN1 <- readRDS(file="RNAprocess_DN1.rds")
RNAprocess_DN1_allgeneSCT <- readRDS(file="RNAprocess_DN1_allgeneSCT.rds")



########## DN1+DN1r2 analysis (add Runx3ko)
set.seed(1)
ADT_DN1r2 <- readADT("DN1r2")
RNA_DN1r2 <- readRNA("DN1r2")
CTassign_DN1r2 <- QCpass_Cell(ADT_DN1r2, RNA_DN1r2, 1000, 0.1, 0.9)

write.table(CTassign_DN1r2,file="CTassign_DN1r2.txt",row.names=T,col.names=F,sep="\t",quote=F)
CTassign_DN1r2 <- read.table("CTassign_DN1r2.txt",row.names=1,header=F)

RNAprocess_DN1r2 <- RNAprocess(RNA_DN1r2[,rownames(CTassign_DN1r2)])
RNAprocess_DN1r2$ADT <- CTassign_DN1r2


p1 <- DimPlot(RNAprocess_DN1r2, label = TRUE,group.by="ADT") + ggtitle("ADT")
p2 <- DimPlot(RNAprocess_DN1r2, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")

pdf(file="UMAP_features_DN1r2.pdf",width=10,height=4)
grid.arrange(p1,p2, ncol=2,nrow=1)
dev.off()

saveRDS(RNAprocess_DN1r2, file="RNAprocess_DN1r2.rds")
RNAprocess_DN1r2 <- readRDS(file="RNAprocess_DN1r2.rds")


#### r2 mtdNA

mtDNA <- read.table("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/CITE/preprocess/DN1/RNA/mtDNA/DN1r2_mtDNAreadnum.txt",row.names=1)

mtDNA_use <- as.numeric(mtDNA[names(RNAprocess_DN1r2$nCount_RNA),1])
names(mtDNA_use) <- names(RNAprocess_DN1r2$nCount_RNA)
RNAprocess_DN1r2$mtDNAreads <- mtDNA_use

idx <- names(mtDNA_use)
mtDNApercent <- mtDNA_use[idx] / (RNAprocess_DN1r2$nCount_RNA[idx]+mtDNA_use[idx])
names(mtDNApercent) <- idx
RNAprocess_DN1r2$MTpercent <- mtDNApercent
#RNAprocess_DN1_allgeneSCT$MTpercent <- RNAprocess_DN1r2$MTpercent 


#p1 <- DimPlot(RNAprocess_DN1, label = TRUE,group.by="ADT") + ggtitle("ADT")
#p2 <- DimPlot(RNAprocess_DN1, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")
#p3 <- FeaturePlot(RNAprocess_DN1,features="MTpercent")+ ggtitle("chrM%")
#
#pdf(file="UMAP_features_mtDNA.pdf",width=15,height=4)
#grid.arrange(p1,p2,p3, ncol=3,nrow=1)
#dev.off()
#
#
#p1<-VlnPlot(RNAprocess_DN1, features="MTpercent",group.by="ADT")
#p2<-VlnPlot(RNAprocess_DN1, features="MTpercent",group.by="seurat_clusters")
#
#pdf(file="boxplot_mtDNA.pdf",width=10,height=4)
#grid.arrange(p1,p2, ncol=2,nrow=1)
#dev.off()


### merge r1, r2
DN1merge <- merge(RNAprocess_DN1, y = c(RNAprocess_DN1r2), add.cell.ids = c("r1", "r2"), project = "DN1merge")


#DN1_r1_norm <- NormalizeData(RNAprocess_DN1)
#DN1_r2_norm <- NormalizeData(RNAprocess_DN1r2)
#DN1normmerge <- merge(DN1_r1_norm, y = DN1_r2_norm, add.cell.ids = c("r1", "r2"), project = "DN1normmerge",
#    merge.data = TRUE)
#GetAssayData(pbmc.combined)[1:10, 1:15]


RNAprocess_DN1merge <- RNAprocess(DN1merge)

## combine mtDNA
idx_r2 <- colnames(RNAprocess_DN1r2)
idx_r2_formerge <- paste0("r2_",idx_r2)

RNAprocess_DN1merge$MTpercent[idx_r2_formerge] <- RNAprocess_DN1r2$MTpercent[idx_r2]
RNAprocess_DN1merge$mtDNAreads[idx_r2_formerge] <- RNAprocess_DN1r2$mtDNAreads[idx_r2]


confmat_ADT_cluster_DN1merge <- confmat(cbind(as.character(RNAprocess_DN1merge$seurat_clusters),RNAprocess_DN1merge$ADT))
rownames(confmat_ADT_cluster_DN1merge) <- paste0("C", rownames(confmat_ADT_cluster_DN1merge))
confmat_ADT_cluster_DN1merge_percent <- t(t(confmat_ADT_cluster_DN1merge)/ apply(confmat_ADT_cluster_DN1merge,2,sum))

write.table(confmat_ADT_cluster_DN1merge, file="scRNA_DN1merge_ADTvsCluster_confmat.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(confmat_ADT_cluster_DN1merge_percent, file="scRNA_DN1merge_ADTvsCluster_confmat_percent.txt",row.names=T,col.names=T,sep="\t",quote=F)

heatmap.2(confmat_ADT_cluster_DN1merge_percent,col=blues9,trace="none",margins=c(10,5))
dev.off()

p1 <- DimPlot(RNAprocess_DN1merge, label = TRUE,group.by="ADT") + ggtitle("ADT")
p2 <- DimPlot(RNAprocess_DN1merge, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")

pdf(file="UMAP_features_DN1merge.pdf",width=10,height=4)
grid.arrange(p1,p2, ncol=2,nrow=1)
dev.off()


########## QC
p1 <- FeaturePlot(RNAprocess_DN1merge,features="MTpercent")+ ggtitle("chrM%")
p2 <- FeaturePlot(RNAprocess_DN1merge,features="nCount_RNA")+ ggtitle("#reads")
p3 <- FeaturePlot(RNAprocess_DN1merge,features="nFeature_RNA")+ ggtitle("#coverGene")

pdf(file="UMAP_QC_DN1merge.pdf",width=15,height=4)
grid.arrange(p1,p2,p3, ncol=3,nrow=1)
dev.off()
#
#
p1<-VlnPlot(RNAprocess_DN1merge, features="MTpercent",group.by="ADT")
p2<-VlnPlot(RNAprocess_DN1merge, features="nCount_RNA",group.by="ADT")
p3<-VlnPlot(RNAprocess_DN1merge, features="nFeature_RNA",group.by="ADT")
p4<-VlnPlot(RNAprocess_DN1merge, features="MTpercent",group.by="orig.ident")
p5<-VlnPlot(RNAprocess_DN1merge, features="nCount_RNA",group.by="orig.ident")
p6<-VlnPlot(RNAprocess_DN1merge, features="nFeature_RNA",group.by="orig.ident")
#p2<-VlnPlot(RNAprocess_DN1, features="MTpercent",group.by="seurat_clusters")
#
pdf(file="boxplot_QC_DN1merge.pdf",width=15,height=8)
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3,nrow=2)
dev.off()


#### rename cluster


renameCluster <- rep("NA",length(RNAprocess_DN1merge$seurat_clusters))
names(renameCluster) <- names(RNAprocess_DN1merge$seurat_clusters) 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 2  )] <- 1 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 3  )] <- 2 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 1  )] <- 3 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 11 )] <- 4 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 7  )] <- 5 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 13 )] <- 6 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 10 )] <- 7 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 15 )] <- 8 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 4  )] <- 9 
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 20 )] <- 10
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 19 )] <- 11
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 23 )] <- 12
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 18 )] <- 13
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 6  )] <- 14
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 16 )] <- 15
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 17 )] <- 16
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 8  )] <- 17
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 12 )] <- 18
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 0  )] <- 19
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 9  )] <- 20
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 5  )] <- 21
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 14 )] <- 22
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 21 )] <- 23
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 22 )] <- 24
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 24 )] <- 25
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 25 )] <- 26
renameCluster[which(RNAprocess_DN1merge$seurat_clusters == 26 )] <- 27

RNAprocess_DN1merge$renameCluster <- renameCluster
#RNAprocess_DN1merge_allgeneSCT$renameCluster <- renameCluster


CT <- rep("NA",length(RNAprocess_DN1merge$renameCluster))
names(CT) <- names(RNAprocess_DN1merge$renameCluster) 
CT[which(RNAprocess_DN1merge$renameCluster %in% c(1,2,3))] <- "ETPs_quiescent"
CT[which(RNAprocess_DN1merge$renameCluster %in% c(4,5,6))] <- "ETPs_proliferative"
CT[which(RNAprocess_DN1merge$renameCluster %in% c(7,8,9,10,11,12))] <- "nonT_potential"
CT[which(RNAprocess_DN1merge$renameCluster %in% c(13,14,15,16,17))] <- "T_ILC"
CT[which(RNAprocess_DN1merge$renameCluster %in% c(18))] <- "TcfKO_up"
CT[which(RNAprocess_DN1merge$renameCluster %in% c(19,20))] <- "RunxKO_up"
CT[which(RNAprocess_DN1merge$renameCluster %in% c(21,22,23))] <- "TleKO_up"
CT[which(RNAprocess_DN1merge$renameCluster %in% c(24,25,26,27))] <- "nonT_minor"
RNAprocess_DN1merge$CT <- CT

RNAprocess_DN1merge_allgeneSCT <-  SCTransform(RNAprocess_DN1merge, assay = "RNA", verbose = FALSE, ,variable.features.n = 24538)
RNAprocess_DN1merge_allgeneSCT$CT <- RNAprocess_DN1merge$CT
RNAprocess_DN1merge_allgeneSCT$renameCluster <- RNAprocess_DN1merge$renameCluster

### UMAP + features
p1 <- DimPlot(RNAprocess_DN1merge, label = TRUE,group.by="ADTcb") + ggtitle("condition")
p2 <- DimPlot(RNAprocess_DN1merge, label = TRUE,group.by="renameCluster") + ggtitle("rename cluster")
p3 <- DimPlot(RNAprocess_DN1merge, label = TRUE,group.by="CT") + ggtitle("cell type")
pdf(file="merge_UMAP_newFeatures.pdf",width=15,height=4)
grid.arrange(p1,p2,p3, ncol=3,nrow=1)
dev.off()

p1 <- DimPlot(RNAprocess_DN1merge, label = FALSE,group.by="ADTcb") + ggtitle("condition")+ NoLegend()
p2 <- DimPlot(RNAprocess_DN1merge, label = FALSE,group.by="renameCluster") + ggtitle("rename cluster")+ NoLegend()
p3 <- DimPlot(RNAprocess_DN1merge, label = FALSE,group.by="CT") + ggtitle("cell type")+ NoLegend()
pdf(file="merge_UMAP_newFeatures_nolabel.pdf",width=12,height=4)
grid.arrange(p1,p2,p3, ncol=3,nrow=1)
dev.off()


confmat_ADT_CT <- confmat(cbind(as.character(RNAprocess_DN1merge$ADTcb),RNAprocess_DN1merge$CT))
confmat_ADT_renameClsuter <- confmat(cbind(as.character(RNAprocess_DN1merge$ADTcb),RNAprocess_DN1merge$renameCluster))

write.table(t(confmat_ADT_CT),file="confmat_ADT_CT.txt",row.names=T,col.names=T,sep="\t")
write.table(t(confmat_ADT_renameClsuter),file="confmat_ADT_renameClsuter.txt",row.names=T,col.names=T,sep="\t")





########### try colors

library(RColorBrewer)
library(viridisLite)
library(viridis)


#barplot(rep(1, 5), col= brewer.pal(5, "Set1"))
#barplot(rep(1, 5), col= brewer.pal(5, "Set2"))
#barplot(rep(1, 5), col= brewer.pal(5, "Set3"))
#barplot(rep(1, 5), col= viridis::viridis_pal()(5))
#barplot(rep(1, 5), col=  terrain.colors(5))
#barplot(rep(1, 5), col= viridisLite::inferno(5))
#barplot(rep(1, 5), col= viridisLite::cividis(5))
#barplot(rep(1, 5), col= viridisLite::magma(5))
#barplot(rep(1, 5), col= viridisLite::mako(5))
#barplot(rep(1, 5), col= viridisLite::plasma(5))
#barplot(rep(1, 5), col= viridisLite::rocket(5))
#barplot(rep(1, 5), col= viridisLite::turbo(5))


pdf(file="merge_UMAP_tryColors.pdf",width=48,height=8)
p1 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = brewer.pal(5, "Set1"))  + ggtitle("brewer.pal.set1")
p2 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = brewer.pal(5, "Set2"))  + ggtitle("brewer.pal.set2")
p3 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = brewer.pal(5, "Set3"))  + ggtitle("brewer.pal.set3")
p4 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = viridis::viridis_pal()(5))  + ggtitle("viridis::viridis_pal")
p5 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols =  terrain.colors(5))  + ggtitle("terrain")
p6 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = viridisLite::inferno(5))  + ggtitle("viridisLite::inferno")
p7 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = viridisLite::cividis(5))  + ggtitle("viridisLite::cividis")
p8 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = viridisLite::magma(5))  + ggtitle("viridisLite::magma")
p9 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = viridisLite::mako(5))  + ggtitle("viridisLite::mako")
p10 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = viridisLite::plasma(5))  + ggtitle("viridisLite::plasma")
p11 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = viridisLite::rocket(5))  + ggtitle("viridisLite::rocket")
p12 <- DimPlot(RNAprocess_DN1merge, group.by = "ADTcb", cols = viridisLite::turbo(5))  + ggtitle("viridisLite::turbo")

p13 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = brewer.pal(8, "Set1"))  + ggtitle("brewer.pal.set1")
p14 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = brewer.pal(8, "Set2"))  + ggtitle("brewer.pal.set2")
p15 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = brewer.pal(8, "Set3"))  + ggtitle("brewer.pal.set3")
p16 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = viridis::viridis_pal()(8))  + ggtitle("viridis::viridis_pal")
p17 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols =  terrain.colors(8))  + ggtitle("terrain")
p18 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = viridisLite::inferno(8))  + ggtitle("viridisLite::inferno")
p19 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = viridisLite::cividis(8))  + ggtitle("viridisLite::cividis")
p20 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = viridisLite::magma(8))  + ggtitle("viridisLite::magma")
p21 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = viridisLite::mako(8))  + ggtitle("viridisLite::mako")
p22 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = viridisLite::plasma(8))  + ggtitle("viridisLite::plasma")
p23 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = viridisLite::rocket(8))  + ggtitle("viridisLite::rocket")
p24 <- DimPlot(RNAprocess_DN1merge, group.by = "CT", cols = viridisLite::turbo(8))  + ggtitle("viridisLite::turbo")
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24, ncol=12,nrow=2)
dev.off()

tmpobj <- RNAprocess_DN1merge
Runx3KO <- rep("NA",length(RNAprocess_DN1merge$ADTcb))
Runx3KO[which(RNAprocess_DN1merge$ADTcb == "Runx3KO")] <- "Runx3KO"
Tcf1KO <- rep("NA",length(RNAprocess_DN1merge$ADTcb))
Tcf1KO[which(RNAprocess_DN1merge$ADTcb == "Tcf1KO")] <- "Tcf1KO"
Tcf1Lef1KO <- rep("NA",length(RNAprocess_DN1merge$ADTcb))
Tcf1Lef1KO[which(RNAprocess_DN1merge$ADTcb == "Tcf1Lef1KO")] <- "Tcf1Lef1KO"
Tle134KO <- rep("NA",length(RNAprocess_DN1merge$ADTcb))
Tle134KO[which(RNAprocess_DN1merge$ADTcb == "Tle134KO")] <- "Tle134KO"
WT <- rep("NA",length(RNAprocess_DN1merge$ADTcb))
WT[which(RNAprocess_DN1merge$ADTcb == "WT")] <- "WT"
tmpobj$Runx3KO <- Runx3KO
tmpobj$Tcf1KO <- Tcf1KO
tmpobj$Tcf1Lef1KO <- Tcf1Lef1KO
tmpobj$Tle134KO <- Tle134KO
tmpobj$WT <- WT

pdf(file="merge_UMAP_ADTsingle_brewerSet1.pdf",width=20,height=4)
p1 <- DimPlot(tmpobj, group.by = "Runx3KO", cols = c(brewer.pal(5, "Set1")[1],"grey")[2:1])  + ggtitle("Runx3KO")
p2 <- DimPlot(tmpobj, group.by = "Tcf1KO", cols = c(brewer.pal(5, "Set1")[2],"grey")[2:1])  + ggtitle("Tcf1KO")
p3 <- DimPlot(tmpobj, group.by = "Tcf1Lef1KO", cols = c(brewer.pal(5, "Set1")[3],"grey")[2:1])  + ggtitle("Tcf1Lef1KO")
p4 <- DimPlot(tmpobj, group.by = "Tle134KO", cols = c(brewer.pal(5, "Set1")[4],"grey")[2:1])  + ggtitle("Tle134KO")
p5 <- DimPlot(tmpobj, group.by = "WT", cols = c(brewer.pal(5, "Set1")[5],"grey")[2:1])  + ggtitle("WT")
grid.arrange(p1,p2,p3,p4,p5,nrow=1,ncol=5)
dev.off()


pdf(file="merge_UMAP_ADTsingle_brewerSet3.pdf",width=20,height=4)
p1 <- DimPlot(tmpobj, group.by = "Runx3KO", cols = c(brewer.pal(5, "Set3")[1],"grey")[2:1])  + ggtitle("Runx3KO")
p2 <- DimPlot(tmpobj, group.by = "Tcf1KO", cols = c(brewer.pal(5, "Set3")[2],"grey")[2:1])  + ggtitle("Tcf1KO")
p3 <- DimPlot(tmpobj, group.by = "Tcf1Lef1KO", cols = c(brewer.pal(5, "Set3")[3],"grey")[2:1])  + ggtitle("Tcf1Lef1KO")
p4 <- DimPlot(tmpobj, group.by = "Tle134KO", cols = c(brewer.pal(5, "Set3")[4],"grey")[2:1])  + ggtitle("Tle134KO")
p5 <- DimPlot(tmpobj, group.by = "WT", cols = c(brewer.pal(5, "Set3")[5],"grey")[2:1])  + ggtitle("WT")
grid.arrange(p1,p2,p3,p4,p5,nrow=1,ncol=5)
dev.off()








##### knowledge base marker
keygene <- c("Tcf7","Lef1","Tle1","Tle3","Tle4","Gata3","Hes1","Bcl11b","Runx1","Nfat","Il7r","Ctcf","Notch1","Nfatc1","Tcf3","Pbx1","Id3","Kit","Lck","Flt3")
usekeygene <- intersect(keygene, rownames(RNAprocess_DN1merge@assays$RNA@counts))


allgenes <- rownames(RNAprocess_DN1merge@assays$RNA@counts)
write.table(allgenes,file="allgenename.txt",row.names=F,col.names=F,sep="\t",quote=F)

allgenes <- read.table("allgenename.txt")[,1]
keygene_HSC <- c("Spi1","Lyl1","Hhex","Bcl11a","Hoxa9","Meis1")
keygene_Tlineage <- c("Notch3","Erg","Tox","Eomes","Sox4","Sox13","Ltk","Prkcq","Tcf7","Bcl11b","Gata3")
keygene_Blineage <- c("Ebf1","Pax5","Btk","Syk","Cd19","Cd79a","Cd79b","Blnk")
keygene_Survival <- "Bax"
keygene_Dendritic <- c("Irf8","Tcf4","Itgax","Batf3","Irf4","Tcf3","Nfil3","Zeb2","Klf4")
keygene_Myeloid <- c("Cebpa","Cebpb","Cebpe","Rara","Mzf1","Hoxa10","Hoxb7","Hoxb8","Itgam","Gata2","Gfi1")
keygene_InnateLymphoid <- c("Id2","Tox","Rorc","Ets1","Bcl11b","Gata3","Zbtb16")
keygene_Cebp <- c("Cebpa","Cebpb","Cebpd","Cebpe","Cebpg","Cebpz","Cebpzos")
keygene_other <- c("Cd24a","Rag1","Rag2","Rbpj","Tbx21","Cd3e","Cd3d","Cd3g","Thy1","Satb1","Cxcr6","Mpo","Cd34","Maml1","Maml2","Maml3","Cd74")
keygene_CT5subset <- c("Aurkb","Birc5","Ccna2","Ccnb1","Cdk1","Lig1","Mki67","Ezh2","Dnmt1","Mcm3","Mcm6","Bcl11a","Lyl1","Cd24a","Kit","Flt3","Notch1","Hhex","Hes1","Spi1","Mpo","Gata2","Elane","Irf8","Cd19","Cd79a","Blnk","Thy1","Cd3e","Zap70","Prkcq","Tcf7","Ets1","Bcl11b","Tox","Gata3","Rorc","Id2","Il7r","Zbtb16")

keygeneNEW_ETPquiescent <- c("Cd34","Plac8","Jun","Fos","Mif")
keygeneNEW_nonT <- c("Lyz2","Gfi1b","Hmgb1","Hmgb2","Hmgb3","Klf1","Lmo2","Top2a","Ccr7","Apoe","Apobec3","Ly6a","Ly6d","Pou2af1","Spib","Cd7","Irf5","Irf7","Zeb2","Cx3cr1","Aff3","Aff4","Atp1a1","Atp2a2","Etv6","Sox13","Tcf12","Fosb","Ikzf1")
keygeneNEW_Tcf1KOup <- c("Ifngr1","Il17re","Nfkb1","Ahr","Rora","Maf")
keygeneNEW_RunxKOup <- c("Cd4","Dgka","Lat","Lck","Nfatc3","Socs1","Socs3","Cish","Themis","Ctla4","Il2ra","Zap70","Sh2d1a")
keygeneNEW_TleKO <- c("Foxo1","Gzma","Id3","Tox2","Klf2","Il23r","Pdcd1","Bhlhe40","Cxcr5","Klf4","Klrd1","Lag3","Sell","Tcf4")
keygeneNEW_TILC <- c("Foxo1","Ikzf2","Ccr7","Cd28","Nr4a1","Nr4a2","Klrb1b","Klrd1","Klrk1","Nkg7","Xcl1","Gzmb","Klrb1","Klrc2","Prf1","Tbx21","Etv5","Sox13","Tcf12")


### key gene exp

pdf(file="UMAP_knowledgeKeyGeneExp_WTCT5subset.pdf",width=20,height=32)
FeaturePlot(RNAprocess_DN1merge,features=keygene_CT5subset,ncol=5)
dev.off()




pdf(file="UMAP_knowledgeKeyGeneExp_WTCT5subset.pdf",width=20,height=32)
FeaturePlot(RNAprocess_DN1merge,features=keygene_CT5subset,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp.pdf",width=20,height=16)
FeaturePlot(RNAprocess_DN1merge,features=usekeygene,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_HSC.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1merge,features=keygene_HSC,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_Tlineage.pdf",width=20,height=12)
FeaturePlot(RNAprocess_DN1merge,features=keygene_Tlineage,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_Blineage.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1merge,features=keygene_Blineage,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_Survival.pdf",width=20,height=4)
FeaturePlot(RNAprocess_DN1merge,features=keygene_Survival,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_Dendritic.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1merge,features=keygene_Dendritic,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_Myeloid.pdf",width=20,height=12)
FeaturePlot(RNAprocess_DN1merge,features=keygene_Myeloid,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_InnateLymphoid.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1merge,features=keygene_InnateLymphoid,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_Cebp.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1merge,features=keygene_Cebp,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_other.pdf",width=20,height=16)
FeaturePlot(RNAprocess_DN1merge,features=keygene_other,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExp_SNAI.pdf",width=20,height=4)
FeaturePlot(RNAprocess_DN1merge,features=c("Snai1","Snai2"),ncol=5)
dev.off()
pdf(file="UMAP_knowledgeKeyGeneExp_TCF.pdf",width=20,height=4)
FeaturePlot(RNAprocess_DN1merge,features=c("Tcf12","Tcf3","Tcf4","Tcf7"),ncol=5)
dev.off()
keygeneNEW_ETPquiescent <- c("Cd34","Plac8","Jun","Fos","Mif")
keygeneNEW_nonT <- c("Lyz2","Gfi1b","Hmgb1","Hmgb2","Hmgb3","Klf1","Lmo2","Top2a","Ccr7","Apoe","Apobec3","Ly6a","Ly6d","Pou2af1","Spib","Cd7","Irf5","Irf7","Zeb2","Cx3cr1","Aff3","Aff4","Atp1a1","Atp2a2","Etv6","Sox13","Tcf12","Fosb","Ikzf1")
keygeneNEW_Tcf1KOup <- c("Ifngr1","Il17re","Nfkb1","Ahr","Rora","Maf")
keygeneNEW_RunxKOup <- c("Cd4","Dgka","Lat","Lck","Nfatc3","Socs1","Socs3","Cish","Themis","Ctla4","Il2ra","Zap70","Sh2d1a")
keygeneNEW_TleKOup <- c("Foxo1","Gzma","Id3","Tox2","Klf2","Il23r","Pdcd1","Bhlhe40","Cxcr5","Klf4","Klrd1","Lag3","Sell","Tcf4")
keygeneNEW_TILC <- c("Foxo1","Ikzf2","Ccr7","Cd28","Nr4a1","Nr4a2","Klrb1b","Klrd1","Klrk1","Nkg7","Xcl1","Gzmb","Klrb1","Klrc2","Prf1","Tbx21","Etv5","Sox13","Tcf12")
keygeneNEW_other <- c("Plac8","Jun","Fos","Mif","Psen1","Psenen","Ncstn","Aph1a")


pdf(file="UMAP_knowledgeKeyGeneExpNEW_other.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1merge,features=keygeneNEW_other,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExpNEW_ETPquiescent.pdf",width=20,height=4)
FeaturePlot(RNAprocess_DN1merge,features=keygeneNEW_ETPquiescent,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExpNEW_nonT.pdf",width=20,height=24)
FeaturePlot(RNAprocess_DN1merge,features=keygeneNEW_nonT,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExpNEW_Tcf1KOup.pdf",width=20,height=8)
FeaturePlot(RNAprocess_DN1merge,features=keygeneNEW_Tcf1KOup,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExpNEW_RunxKOup.pdf",width=20,height=12)
FeaturePlot(RNAprocess_DN1merge,features=keygeneNEW_RunxKOup,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExpNEW_TleKOup.pdf",width=20,height=12)
FeaturePlot(RNAprocess_DN1merge,features=keygeneNEW_TleKOup,ncol=5)
dev.off()

pdf(file="UMAP_knowledgeKeyGeneExpNEW_TILC.pdf",width=20,height=16)
FeaturePlot(RNAprocess_DN1merge,features=keygeneNEW_TILC,ncol=5)
dev.off()


#### compare to WT cell type

RNAprocess_DN1WT <- readRDS(file="../WTonly/RNAprocess_DN1WT.rds")
#RNAprocess_DN1WT_allgeneSCT <- readRDS(file="../WTonly/RNAprocess_DN1WT_allgeneSCT.rds")

length(names(RNAprocess_DN1merge$orig.ident))
length(names(RNAprocess_DN1WT$orig.ident))
length(intersect(names(RNAprocess_DN1WT$orig.ident), names(RNAprocess_DN1merge$orig.ident)))

WTonlyCTname <- rep("NA",length(names(RNAprocess_DN1merge$orig.ident)))
names(WTonlyCTname) <- names(RNAprocess_DN1merge$orig.ident)
WTonlyCTname[names(RNAprocess_DN1WT$orig.ident)] <- RNAprocess_DN1WT$CT_5subset[names(RNAprocess_DN1WT$orig.ident)]
RNAprocess_DN1merge$WTonlyCTname <- WTonlyCTname

RNAprocess_DN1merge$ADTcb <- RNAprocess_DN1merge$ADT
RNAprocess_DN1merge$ADTcb[which(RNAprocess_DN1merge$ADT %in% c("Runx3KOr1","Runx3KOr2"))] <- "Runx3KO"
RNAprocess_DN1merge$ADTcb[which(RNAprocess_DN1merge$ADT %in% c("WTr1","WTr2"))] <- "WT"


p1 <- DimPlot(RNAprocess_DN1merge, label = TRUE,group.by="ADTcb") + ggtitle("ADTcb")
p2 <- DimPlot(RNAprocess_DN1merge, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")
p3 <- DimPlot(RNAprocess_DN1merge, label = TRUE,group.by="WTonlyCTname") + ggtitle("WTonlyCTname")
pdf(file="UMAP_features_DN1merge_WTcelltype.pdf",width=15,height=4)
grid.arrange(p1,p2,p3, ncol=3,nrow=1)
dev.off()


p1 <- DimPlot(RNAprocess_DN1merge, label = FALSE,group.by="ADTcb") + ggtitle("ADTcb")+ NoLegend()
p2 <- DimPlot(RNAprocess_DN1merge, label = FALSE,group.by="seurat_clusters") + ggtitle("cluster")+ NoLegend()
p3 <- DimPlot(RNAprocess_DN1merge, label = FALSE,group.by="WTonlyCTname") + ggtitle("WTonlyCTname")+ NoLegend()
pdf(file="UMAP_features_DN1merge_WTcelltype_nolabel.pdf",width=12,height=4)
grid.arrange(p1,p2,p3, ncol=3,nrow=1)
dev.off()


confmat_ADT_cluster <- confmat(cbind(as.character(RNAprocess_DN1merge$seurat_clusters),RNAprocess_DN1merge$ADTcb))
confmat_WTcelltype_cluster <- confmat(cbind(as.character(RNAprocess_DN1merge$seurat_clusters),RNAprocess_DN1merge$WTonlyCTname))

confmat_WTct_ADT_cluster <- cbind(confmat_ADT_cluster, confmat_WTcelltype_cluster[rownames(confmat_ADT_cluster),])
write.table(confmat_WTct_ADT_cluster, file="confmat_WTct_ADT_cluster.txt",row.names=T,col.names=T,sep="\t",quote=F)

pdf(file="confmat_WTct_ADT_cluster_heatmap.pdf")
heatmap.2(log10(confmat_WTct_ADT_cluster[,c("Tcf1KO","Tle134KO","Tcf1Lef1KO","Runx3KO","ETPs_quiescent","ETPs_proliferative","nonT_potential","T_ILCs_IIhi","T_ILCs_GRhi")]+1),col=blues9,trace="none",margins=c(10,5))
dev.off()


#### ETPq exp
ETPq_cell <- names(RNAprocess_DN1merge$CT)[which(RNAprocess_DN1merge$CT=="ETPs_quiescent")]
ETPq_aveExp <- apply(RNAprocess_DN1merge@assays$RNA@data[, ETPq_cell],1,mean)
write.table(ETPq_aveExp, file="ETPq_aveExp.txt",row.names=T,col.names=F,sep="\t",quote=F)

RunxKO_cell <- names(RNAprocess_DN1merge$CT)[which(RNAprocess_DN1merge$CT=="RunxKO_up")]
RunxKO_aveExp <- apply(RNAprocess_DN1merge@assays$RNA@data[, RunxKO_cell],1,mean)
write.table(RunxKO_aveExp, file="RunxKO_up_aveExp.txt",row.names=T,col.names=F,sep="\t",quote=F)

TcfKO_cell <- names(RNAprocess_DN1merge$CT)[which(RNAprocess_DN1merge$CT=="TcfKO_up")]
TcfKO_aveExp <- apply(RNAprocess_DN1merge@assays$RNA@data[, TcfKO_cell],1,mean)
write.table(TcfKO_aveExp, file="TcfKO_up_aveExp.txt",row.names=T,col.names=F,sep="\t",quote=F)

TleKO_cell <- names(RNAprocess_DN1merge$CT)[which(RNAprocess_DN1merge$CT=="TleKO_up")]
TleKO_aveExp <- apply(RNAprocess_DN1merge@assays$RNA@data[, TleKO_cell],1,mean)
write.table(TleKO_aveExp, file="TleKO_up_aveExp.txt",row.names=T,col.names=F,sep="\t",quote=F)

idx <- names(ETPq_aveExp)
aveExpAll <- cbind(ETPq_aveExp[idx],
				   RunxKO_aveExp[idx],
				   TcfKO_aveExp[idx],
				   TleKO_aveExp[idx])
colnames(aveExpAll) <- c("ETPq_aveExp","RunxKO_aveExp","TcfKO_aveExp","TleKO_aveExp")
##### ETP vs KO DEG
#ETPq_vs_RunxKO_DEG
this_diff <- FindMarkers(RNAprocess_DN1merge, ident.1="ETPs_quiescent", ident.2="RunxKO_up",group.by="CT")
this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
this_diff_order <- this_diff[order(this_diff[,2]),]
ETPq_vs_RunxKO_DEG<-this_diff_order
write.table(ETPq_vs_RunxKO_DEG,file=paste0("ETPq_vs_RunxKO_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)


#ETPq_vs_RunxKO_DEG
this_diff <- FindMarkers(RNAprocess_DN1merge, ident.1="ETPs_quiescent", ident.2="TcfKO_up",group.by="CT")
this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
this_diff_order <- this_diff[order(this_diff[,2]),]
ETPq_vs_RunxKO_DEG<-this_diff_order
write.table(ETPq_vs_RunxKO_DEG,file=paste0("ETPq_vs_TcfKO_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)


this_diff <- FindMarkers(RNAprocess_DN1merge, ident.1="ETPs_quiescent", ident.2="TleKO_up",group.by="CT")
this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
this_diff_order <- this_diff[order(this_diff[,2]),]
ETPq_vs_RunxKO_DEG<-this_diff_order
write.table(ETPq_vs_RunxKO_DEG,file=paste0("ETPq_vs_TleKO_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)




#############33 non-diff gene

##### ETP vs KO DEG
#ETPq_vs_RunxKO_DEG
diff_ETP_RunxKO <- FindMarkers(RNAprocess_DN1merge, ident.1="ETPs_quiescent", ident.2="RunxKO_up",group.by="CT")
diff_ETP_TcfKO <- FindMarkers(RNAprocess_DN1merge, ident.1="ETPs_quiescent", ident.2="TcfKO_up",group.by="CT")
diff_ETP_TleKO <- FindMarkers(RNAprocess_DN1merge, ident.1="ETPs_quiescent", ident.2="TleKO_up",group.by="CT")

summary_diff <- matrix(rep(0,nrow(aveExpAll)*3),ncol=3)
rownames(summary_diff) <- rownames(aveExpAll)
colnames(summary_diff) <- c("ETP_RunxKO","ETP_TcfKO","ETP_TleKO")
summary_diff[rownames(diff_ETP_RunxKO), "ETP_RunxKO"] <- 1
summary_diff[rownames(diff_ETP_TcfKO), "ETP_TcfKO"] <- 1
summary_diff[rownames(diff_ETP_TleKO), "ETP_TleKO"] <- 1

pdf(file="aveExp_hist.pdf",width=12,height=12)
par(mfrow=c(2,2),mar=c(4,4,2,2))
hist(outdata[,4],n=200,xlab="ave Exp", main=colnames(outdata)[4])
hist(outdata[,5],n=200,xlab="ave Exp", main=colnames(outdata)[5])
hist(outdata[,6],n=200,xlab="ave Exp", main=colnames(outdata)[6])
hist(outdata[,7],n=200,xlab="ave Exp", main=colnames(outdata)[7])
dev.off()
outdata <- cbind(summary_diff, aveExpAll)
write.table(outdata,file="allcond_DEGstatus_aveExp.txt",row.names=T,col.names=T,sep="\t",quote=F)




this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
this_diff_order <- this_diff[order(this_diff[,2]),]
ETPq_vs_RunxKO_DEG<-this_diff_order
write.table(ETPq_vs_RunxKO_DEG,file=paste0("ETPq_vs_RunxKO_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)


#ETPq_vs_RunxKO_DEG
this_diff <- FindMarkers(RNAprocess_DN1merge, ident.1="ETPs_quiescent", ident.2="TcfKO_up",group.by="CT")
this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
this_diff_order <- this_diff[order(this_diff[,2]),]
ETPq_vs_RunxKO_DEG<-this_diff_order
write.table(ETPq_vs_RunxKO_DEG,file=paste0("ETPq_vs_TcfKO_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)


this_diff <- FindMarkers(RNAprocess_DN1merge, ident.1="ETPs_quiescent", ident.2="TleKO_up",group.by="CT")
this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
this_diff_order <- this_diff[order(this_diff[,2]),]
ETPq_vs_RunxKO_DEG<-this_diff_order
write.table(ETPq_vs_RunxKO_DEG,file=paste0("ETPq_vs_TleKO_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)











############# target pair DEG
g1name <- "ETPs_quiescent"
g2name <- "TcfKO_up"

thisdata <- RNAprocess_DN1merge
C123cond <- rep("NA", length(thisdata$ADTcb))
C123cond[which(thisdata$CT == "ETPs_quiescent" & thisdata$WTonlyCTname == "ETPs_quiescent")] <- "WT"
C123cond[which(thisdata$CT == "ETPs_quiescent" & thisdata$ADTcb == "Tcf1KO")] <- "Tcf1KO"
thisdata$C123cond <- C123cond
this_diff <- FindMarkers(thisdata, ident.1="WT",ident.2="Tcf1KO",group.by="C123cond")
this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
this_diff_order <- this_diff[order(this_diff[,2]),]
C123condDEG<-this_diff_order
write.table(this_diff_order,file=paste0("C123cond_WT_Tcf1KO_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)



markers <- this_diff
markers$negLog10_p_val_adj <- -log10(markers$p_val_adj)
significant_genes <- ifelse(abs(markers[,"avg_log2FC"])>=log2(1.5) & markers[,"p_val_adj"] < 0.001, "Significant", "Not Significant")

# Subset significant genes for labeling
label_genes <- markers[which(significant_genes=="Significant"),]#subset(markers, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
volcano <- ggplot(markers, aes(x = avg_log2FC, y = negLog10_p_val_adj, color = significant_genes)) + 
  geom_point(alpha = 0.6, size = 1.5) + 
  scale_color_manual(values = c("gray", "red")) + 
  ggtitle("Volcano plot") + 
  theme_minimal() +
  xlab("Average log2 fold change") + 
  ylab("-log10 adjusted p-value")
volcano + 
  geom_text(data = label_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(label_genes)), vjust = 1.5, hjust = 0.5, size = 3, angle = 45, check_overlap = TRUE, inherit.aes = FALSE)
ggsave("C123cond_WT_Tcf1KO_DEG_volcano.pdf")


thisdata <- RNAprocess_DN1merge
C456cond <- rep("NA", length(thisdata$ADTcb))
C456cond[which(thisdata$CT == "ETPs_proliferative" & thisdata$WTonlyCTname == "ETPs_proliferative")] <- "WT"
C456cond[which(thisdata$CT == "ETPs_proliferative" & thisdata$ADTcb == "Tcf1KO")] <- "Tcf1KO"
thisdata$C456cond <- C456cond
this_diff <- FindMarkers(thisdata, ident.1="WT",ident.2="Tcf1KO",group.by="C456cond")
this_diff_filter <- this_diff[which(abs(this_diff[,"avg_log2FC"])>=1 & this_diff[,"p_val_adj"] < 0.01),]
this_diff_order <- this_diff[order(this_diff[,2]),]
C456condDEG<-this_diff_order
write.table(this_diff_order,file=paste0("C456cond_WT_Tcf1KO_DEG.txt"),row.names=T,col.names=T,sep="\t",quote=F)



markers <- this_diff
markers$negLog10_p_val_adj <- -log10(markers$p_val_adj)
significant_genes <- ifelse(abs(markers[,"avg_log2FC"])>=log2(1.5) & markers[,"p_val_adj"] < 0.001, "Significant", "Not Significant")

# Subset significant genes for labeling
label_genes <- markers[which(significant_genes=="Significant"),]#subset(markers, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
volcano <- ggplot(markers, aes(x = avg_log2FC, y = negLog10_p_val_adj, color = significant_genes)) + 
  geom_point(alpha = 0.6, size = 1.5) + 
  scale_color_manual(values = c("gray", "red")) + 
  ggtitle("Volcano plot") + 
  theme_minimal() +
  xlab("Average log2 fold change") + 
  ylab("-log10 adjusted p-value")
volcano + 
  geom_text(data = label_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), label = rownames(label_genes)), vjust = 1.5, hjust = 0.5, size = 3, angle = 45, check_overlap = TRUE, inherit.aes = FALSE)
ggsave("C456cond_WT_Tcf1KO_DEG_volcano.pdf")




idx <- intersect(rownames(C123condDEG),rownames(C456condDEG))
cmp_condDEG <- cbind(C123condDEG[idx,'avg_log2FC'],C456condDEG[idx,'avg_log2FC'])
rownames(cmp_condDEG) <- idx
colnames(cmp_condDEG) <- c("C123log2fc","C456log2fc")
pdf(file="C123log2fc_vs_C456log2fc_scatter.pdf")
plot(cmp_condDEG,pch=16,xlim=c(-2,2),ylim=c(-2,2))
abline(a=0,b=1,col="red")
abline(h=0,v=0,col="red")
dev.off()
idx[which(cmp_condDEG[,1]<0 & cmp_condDEG[,2]>0)]











######## cluster mkg
MKGcluster_DN1 <- FindAllMarkers(RNAprocess_DN1merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MKGclusterFilter_DN1 <- MKGcluster_DN1[which(MKGcluster_DN1[,2] >= log2(2) & MKGcluster_DN1[,5] < 0.01),]
write.table(MKGclusterFilter_DN1, file="MKGclusterFilter_DN1merge.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(MKGcluster_DN1, file="MKGcluster_DN1merge.txt",row.names=T,col.names=T,sep="\t",quote=F)


outdata <- c()
for(C in sort(unique(MKGclusterFilter_DN1[,"cluster"])) ){
	tmpdata <- MKGclusterFilter_DN1[which(MKGclusterFilter_DN1[,"cluster"]==C),]
	outdata <- rbind(outdata, c("DN1",paste0("C",C),tmpdata[order(tmpdata[,"avg_log2FC"],decreasing=T),][1:5,"gene"]))
}
write.table(outdata,file="MKGclusterTop5_DN1merge.txt",row.names=T,col.names=T,sep="\t",quote=F)

alltop5MKG <- c()
for(this_cluster in 0:26 ){
	#print(this_cluster)
	this_MKGdata <- MKGclusterFilter_DN1[which(MKGclusterFilter_DN1[,"cluster"] == this_cluster),]
	if(nrow(this_MKGdata)>0){
        this_MKGdata_top5 <- this_MKGdata[order(this_MKGdata[,"avg_log2FC"],decreasing=T),][1:min(5,nrow(this_MKGdata)),]
        tmpobj <- RNAprocess_DN1merge
        tmpobj$target_cluster <- rep(0, ncol(tmpobj))
        tmpobj$target_cluster[which(tmpobj$seurat_clusters == this_cluster)] <- 1
        
        #pdf(file=paste0("UMAP_keygeneExp_C",this_cluster,".pdf"),width=24,height=4)
        FeaturePlot(tmpobj,features=c("target_cluster",as.vector(this_MKGdata_top5[,"gene"])),ncol=6)
        #dev.off()
        ggsave(paste0("UMAP_keygeneExp_C",this_cluster,".pdf"),width=24,height=4)
        alltop5MKG <- c(alltop5MKG, as.vector(this_MKGdata_top5[,"gene"]))		
	}
}








#saveRDS(RNAprocess_DN1merge, file="RNAprocess_DN1merge.rds")
RNAprocess_DN1merge <- readRDS(file="RNAprocess_DN1merge.rds")

#saveRDS(RNAprocess_DN1merge_allgeneSCT, file="RNAprocess_DN1merge_allgeneSCT.rds")
RNAprocess_DN1merge_allgeneSCT <- readRDS(file="RNAprocess_DN1merge_allgeneSCT.rds")






## print cell type reads
cell_types <- unique(RNAprocess_DN1merge$CT)
for (ct in cell_types) {
  # Subset Seurat object for cells of the current type
  cells_of_type <- WhichCells(RNAprocess_DN1merge, expression = CT == ct)
  subset_seurat <- subset(RNAprocess_DN1merge, cells = cells_of_type)

  # Extract read counts
  # Adjust the slot (e.g., "counts" or "data") based on where your read counts are stored
  read_counts <- subset_seurat@assays$RNA@counts 

  # Process read_counts as needed, e.g., convert to matrix, sum per gene, etc.
  # You can print, save to a file, or perform further analysis as required

  print(paste("Read counts for cell type:", ct))
  print(read_counts)
}











### WT + dko
RNAprocess_DN1WTKO <- RNAprocess(RNAprocess_DN1[,c(WT_cell,Tcf1Lef1KO_cell)])

confmat_ADT_cluster_WTKO <- confmat(cbind(as.character(RNAprocess_DN1WTKO$seurat_clusters),RNAprocess_DN1WTKO$ADT))
write.table(confmat_ADT_cluster_WTKO, file="scRNA_WTdKO_ADTvsCluster_confmat.txt",row.names=T,col.names=T,sep="\t",quote=F)
p1 <- DimPlot(RNAprocess_DN1WTKO, label = TRUE,group.by="ADT") + ggtitle("ADT")
p2 <- DimPlot(RNAprocess_DN1WTKO, label = TRUE,group.by="seurat_clusters") + ggtitle("cluster")

pdf(file="UMAP_features_WTKO.pdf",width=10,height=4)
grid.arrange(p1,p2, ncol=2,nrow=1)
dev.off()






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


