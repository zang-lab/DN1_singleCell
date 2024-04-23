library(ArchR) # version 
set.seed(1)
addArchRGenome("mm10")
addArchRThreads(1)

inputFiles <- c(
    "/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/cellrangerATAC/WT/DN1_WT/outs/fragments.tsv.gz",
    "/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/cellrangerATAC/Runx1KO/DN1_Runx1KO/outs/fragments.tsv.gz",
    "/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/cellrangerATAC/Tcf1Lef1dko/DN1_Tcf1Lef1dko/outs/fragments.tsv.gz",    
    "/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/cellrangerATAC/Tle34KO/DN1_Tle34KO/outs/fragments.tsv.gz")
#names(inputFiles)<-c("WT","Tcf1Lef1dko")
names(inputFiles)<-c("WT","Runx1KO","Tcf1Lef1dKO","Tle34KO")





#doubScores2 <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
#  LSIMethod = 1,
#  force=TRUE
#)
#inputFiles <- c("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/ArchR/combine4/WT.arrow",
#                "/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/ArchR/combine4/Runx1KO.arrow",
#                "/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/ArchR/combine4/Tcf1Lef1dKO.arrow",
#                "/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/ArchR/combine4/Tle34KO.arrow")

# Create an ArchRProject
#proj <- createArchRProject(
#  inputFiles = inputFiles,
#  outputDirectory = "NameOfYourProject",
#  projectName = "NameOfYourProject",
#  genome = "hg38" # Adjust this to the genome you're working with
#)
# create ArchR object
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  force=TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force=TRUE
)

projDN1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "DN1",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)


# basic QC 
proj_DN1_1 <- projDN1
df <- getCellColData(proj_DN1_1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj_DN1_1, addDOC = FALSE)

p1 <- plotGroups(
    ArchRProj = proj_DN1_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p2 <- plotGroups(
    ArchRProj = proj_DN1_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p3 <- plotGroups(
    ArchRProj = proj_DN1_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj_DN1_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj_DN1_1, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = proj_DN1_1)
p2 <- plotTSSEnrichment(ArchRProj = proj_DN1_1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj_DN1_1, addDOC = FALSE, width = 5, height = 5)





# filter cells
proj_DN1_2 <- filterDoublets(proj_DN1_1,filterRatio=1.5)
idxPass <- which(proj_DN1_2$TSSEnrichment >= 7 & proj_DN1_2$nFrags >= 10000)
df2 <- getCellColData(proj_DN1_2,select = c("log10(nFrags)", "TSSEnrichment"))
cellsPass <- proj_DN1_2$cellNames[idxPass]
proj_DN1_2 <- proj_DN1_2[cellsPass, ]

allbarcode_measure <- getCellColData(proj_DN1_1)
write.table(allbarcode_measure, file="allbarcode_measure.txt",row.names=T,col.names=T,sep="\t",quote=F)
highQcell_measure <- getCellColData(proj_DN1_2)
write.table(highQcell_measure, file="highQcell_measure.txt",row.names=T,col.names=T,sep="\t",quote=F)


get_tag1 <- function(inname){
    return(strsplit(inname,"#")[[1]][1])
}


sampleTags <- unlist(lapply(proj_DN1_2$cellNames,get_tag1))
WTcells <- proj_DN1_2$cellNames[which(sampleTags == "WT")]
Runx1KOcells <- proj_DN1_2$cellNames[which(sampleTags == "Runx1KO")]
Tcf1Lef1dKOcells <- proj_DN1_2$cellNames[which(sampleTags == "Tcf1Lef1dKO")]
Tle34KOcells <- proj_DN1_2$cellNames[which(sampleTags == "Tle34KO")]
proj_DN1_2_WT <- proj_DN1_2[WTcells,]
proj_DN1_2_Runx1KO <- proj_DN1_2[Runx1KOcells,]
proj_DN1_2_Tcf1Lef1dKO <- proj_DN1_2[Tcf1Lef1dKOcells,]
proj_DN1_2_Tle34KO <- proj_DN1_2[Tle34KOcells,]



p1 <- plotFragmentSizes(ArchRProj = proj_DN1_2_WT)
p2 <- plotTSSEnrichment(ArchRProj = proj_DN1_2_WT)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_WT.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)

p1 <- plotFragmentSizes(ArchRProj = proj_DN1_2_Runx1KO)
p2 <- plotTSSEnrichment(ArchRProj = proj_DN1_2_Runx1KO)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_Runx1KO.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)

p1 <- plotFragmentSizes(ArchRProj = proj_DN1_2_Tcf1Lef1dKO)
p2 <- plotTSSEnrichment(ArchRProj = proj_DN1_2_Tcf1Lef1dKO)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_Tcf1Lef1dKO.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)

p1 <- plotFragmentSizes(ArchRProj = proj_DN1_2_Tle34KO)
p2 <- plotTSSEnrichment(ArchRProj = proj_DN1_2_Tle34KO)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_Tle34KO.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)


p <- ggPoint(
    x = df2[,"log10(nFrags)"], 
    y = df2[,"TSSEnrichment"], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(3, quantile(df2[,"log10(nFrags)"], probs = 0.99)),
    ylim = c(4, quantile(df2[,"TSSEnrichment"], probs = 0.99))
) + geom_hline(yintercept = 7, lty = "dashed") + geom_vline(xintercept = 4, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags_cutoff.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE)

# dimensional reduction
proj_DN1_2 <- addIterativeLSI(
    ArchRProj = proj_DN1_2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    seed=1,force=T
)

# basic clustering 
proj_DN1_2 <- addClusters(
    input = proj_DN1_2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force=T,seed=1
)


# UMAP embedding
proj_DN1_2 <- addUMAP(
    ArchRProj = proj_DN1_2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",force=T
)

p1 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_LSI.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)

# QC score projected on UMAP
p1 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "PromoterRatio", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "NucleosomeRatio", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-QC_LSI.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)








# gene score with ArchR default method
proj_DN1_2 <- addGeneScoreMatrix(proj_DN1_2,force=TRUE)
proj_DN1_2 <- addImputeWeights(proj_DN1_2,seed=1)

#saveArchRProject(ArchRProj = proj_DN1_2, outputDirectory = "SaveDN1_2", load = FALSE)
#saveArchRProject(ArchRProj = proj_DN1_2_HAR, outputDirectory = "SaveDN1_2_HAR", load = FALSE)

proj_DN1_2 <- loadArchRProject(path = "SaveDN1_2")



### marker gene
markersGS <- getMarkerFeatures(
    ArchRProj = proj_DN1_2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# cluster specific genes
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
#
saveRDS(markersGS, file="DN1_2_markerGS.rds")
#saveRDS(markerList, file="DN1_2_markerList.rds")
markerList <- readRDS("DN1_2_markerList.rds")


printMKG <- function(MKGdata){
  outdata <- matrix(rep(0,nrow(MKGdata)*ncol(MKGdata)),nrow=nrow(MKGdata))
  colnames(outdata) <- colnames(MKGdata)
  for(i in 1:ncol(outdata)){
    outdata[,i] <- as.vector(MKGdata[,i])
  }
  return(outdata)
}
alltop5 <- c()
for(Group in names(markerList)){
    thisdata <- printMKG(markerList[Group][[1]])
    alltop5 <- rbind(alltop5, cbind(rep(Group, 5),thisdata[order(thisdata[,"Log2FC"],decreasing=T),][1:5,]))
    write.table(thisdata[order(thisdata[,"Log2FC"]),], file=paste0("MKG_geneScore_clusters/MKG_geneScore_cluster",Group,".txt"),row.names=F,col.names=T,sep="\t",quote=F)
}
write.table(alltop5[,c(1,6,8,9)], file="MKG_geneScore_clusters/MKG_geneScore_clusterTop5summary.txt",row.names=F,col.names=T,sep="\t",quote=F)



### mark key genes
allgenes <- as.vector(proj_DN1_2@geneAnnotation$genes$symbol)

keygene <- c("Tcf7","Lef1","Tle1","Tle3","Tle4","Gata3","Hes1","Bcl11b","Runx1","Nfat","Il7r","Ctcf","Notch1","Nfatc1","Tcf3","Pbx1","Id3","Kit","Lck","Flt3")
keygene_HSC <- c("Spi1","Lyl1","Hhex","Bcl11a","Hoxa9","Meis1")
keygene_Tlineage <- c("Notch3","Erg","Tox","Eomes","Sox4","Sox13","Ltk","Prkcq","Tcf7","Bcl11b","Gata3")
keygene_Blineage <- c("Ebf1","Pax5","Btk","Syk","Cd19","Cd79a","Cd79b","Blnk")
keygene_Survival <- "Bax"
keygene_Dendritic <- c("Irf8","Tcf4","Itgax","Batf3","Irf4","Tcf3","Nfil3","Zeb2","Klf4")
keygene_Myeloid <- c("Cebpa","Cebpb","Cebpe","Rara","Mzf1","Hoxa10","Hoxb7","Hoxb8","Itgam","Gata2","Gfi1")
keygene_InnateLymphoid <- c("Id2","Tox","Rorc","Ets1","Bcl11b","Gata3","Zbtb16")
keygene_Cebp <- c("Cebpa","Cebpb","Cebpd","Cebpe","Cebpg","Cebpz","Cebpzos")
keygene_other <- c("Cd24a","Rag1","Rag2","Rbpj","Tbx21","Cd3e","Cd3d","Cd3g","Thy1","Satb1","Cxcr6","Mpo")
cbkeygene <- intersect(allgenes, c(keygene, keygene_HSC,keygene_Tlineage,keygene_Blineage,
                                   keygene_Survival,keygene_Dendritic,keygene_Myeloid,keygene_InnateLymphoid,
                                   keygene_Cebp,keygene_other))

pdf(file="UMAP_knowledgeKeyGeneExp.pdf",width=20,height=16)
FeaturePlot(RNAprocess_DN1,features=usekeygene,ncol=5)
dev.off()

keygeneExpUmap <- function(useMKG0){
    useMKG <- intersect(useMKG0,allgenes )
    pImp <- plotEmbedding(
    ArchRProj = proj_DN1_2, 
    colorBy = "GeneScoreMatrix", 
    name = useMKG, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_DN1_2)
    )
#    p <- plotEmbedding(
#        ArchRProj = proj_DN1_2, 
#        colorBy = "GeneScoreMatrix", 
#        name = useMKG, 
#        embedding = "UMAP",
#        quantCut = c(0.01, 0.95),
#        imputeWeights = NULL
#    )
    p2Imp <- lapply(pImp, function(x){
        x + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x=element_blank(), 
            axis.ticks.x=element_blank(), 
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank()
        )
    })
    return(p2Imp)
    #pdf(file=paste0("UMAP_geneScoreImpute_",outname,".pdf"),width=12,height=6)
    #cowplot::plot_grid(plotlist=p2Imp[1:8],nrow=2,ncol=4)
    #if(length(p2Imp)>8){
    #    cowplot::plot_grid(plotlist=p2Imp[9:16],nrow=2,ncol=4)
    #}
    #if(length(p2Imp)>8){
    #    cowplot::plot_grid(plotlist=p2Imp[17:24],nrow=2,ncol=4)
    #}
    #cowplot::plot_grid(plotlist=p2Imp[25:32],nrow=2,ncol=4)
    #cowplot::plot_grid(plotlist=p2Imp[33:40],nrow=2,ncol=4)
    #cowplot::plot_grid(plotlist=p2Imp[41:48],nrow=2,ncol=4)
    #cowplot::plot_grid(plotlist=p2Imp[49:54],nrow=2,ncol=4)
    #dev.off()
}

keygene_Img <- keygeneExpUmap(keygene)
keygene_HSC_Img <- keygeneExpUmap(keygene_HSC)
keygene_Tlineage_Img <- keygeneExpUmap(keygene_Tlineage)
keygene_Blineage_Img <- keygeneExpUmap(keygene_Blineage)
keygene_Survival_Img <- keygeneExpUmap(keygene_Survival)
keygene_Dendritic_Img <- keygeneExpUmap(keygene_Dendritic)
keygene_Myeloid_Img <- keygeneExpUmap(keygene_Myeloid)
keygene_InnateLymphoid_Img <- keygeneExpUmap(keygene_InnateLymphoid)
keygene_Cebp_Img <- keygeneExpUmap(keygene_Cebp)
keygene_other_Img <- keygeneExpUmap(keygene_other)

pdf(file=paste0("UMAP_geneScore_keygene.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_Img[1:8],nrow=2,ncol=4)
cowplot::plot_grid(plotlist=keygene_Img[9:16],nrow=2,ncol=4)
cowplot::plot_grid(plotlist=keygene_Img[17:24],nrow=2,ncol=4)
dev.off()

pdf(file=paste0("UMAP_geneScore_HSC.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_HSC_Img[1:8],nrow=2,ncol=4)
dev.off()
pdf(file=paste0("UMAP_geneScore_Tlineage.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_Tlineage_Img[1:8],nrow=2,ncol=4)
cowplot::plot_grid(plotlist=keygene_Tlineage_Img[9:16],nrow=2,ncol=4)
dev.off()
pdf(file=paste0("UMAP_geneScore_Blineage.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_Blineage_Img[1:8],nrow=2,ncol=4)
dev.off()
pdf(file=paste0("UMAP_geneScore_Survival.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_Survival_Img[1:8],nrow=2,ncol=4)
dev.off()
pdf(file=paste0("UMAP_geneScore_Dendritic.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_Dendritic_Img[1:8],nrow=2,ncol=4)
cowplot::plot_grid(plotlist=keygene_Dendritic_Img[9:16],nrow=2,ncol=4)
dev.off()
pdf(file=paste0("UMAP_geneScore_Myeloid.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_Myeloid_Img[1:8],nrow=2,ncol=4)
cowplot::plot_grid(plotlist=keygene_Myeloid_Img[9:16],nrow=2,ncol=4)
dev.off()
pdf(file=paste0("UMAP_geneScore_InnateLymphoid.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_InnateLymphoid_Img[1:8],nrow=2,ncol=4)
dev.off()
pdf(file=paste0("UMAP_geneScore_Cebp.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_Cebp_Img[1:8],nrow=2,ncol=4)
dev.off()
pdf(file=paste0("UMAP_geneScore_other.pdf"),width=12,height=6)
cowplot::plot_grid(plotlist=keygene_other_Img[1:8],nrow=2,ncol=4)
cowplot::plot_grid(plotlist=keygene_other_Img[9:16],nrow=2,ncol=4)
dev.off()


p1heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p2heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_HSC),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p3heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Tlineage),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p4heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Blineage),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p5heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Survival),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p6heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Dendritic),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p7heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Myeloid),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p8heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_InnateLymphoid),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p9heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Cebp),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p10heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_other),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p1 <- ComplexHeatmap::draw(p1heatmapGS, column_title = "keygene",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p2 <- ComplexHeatmap::draw(p2heatmapGS, column_title = "HSC",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p3 <- ComplexHeatmap::draw(p3heatmapGS, column_title = "T lineage",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p4 <- ComplexHeatmap::draw(p4heatmapGS, column_title = "B lineage",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p5 <- ComplexHeatmap::draw(p5heatmapGS, column_title = "Survival",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p6 <- ComplexHeatmap::draw(p6heatmapGS, column_title = "Dendritic",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p7 <- ComplexHeatmap::draw(p7heatmapGS, column_title = "Myeloid",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p8 <- ComplexHeatmap::draw(p8heatmapGS, column_title = "InnateLymphoid",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p9 <- ComplexHeatmap::draw(p9heatmapGS, column_title = "Cebp",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p10 <- ComplexHeatmap::draw(p10heatmapGS, column_title = "other",heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(p1, name = "MKGexp_heatmap_keygene.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2, name = "MKGexp_heatmap_HSC.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p3, name = "MKGexp_heatmap_Tlineage.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p4, name = "MKGexp_heatmap_Blineage.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p5, name = "MKGexp_heatmap_Survival.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p6, name = "MKGexp_heatmap_Dendritic.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p7, name = "MKGexp_heatmap_Myeloid.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p8, name = "MKGexp_heatmap_InnateLymphoid.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p9, name = "MKGexp_heatmap_Cebp.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p10, name = "MKGexp_heatmap_other.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)





############# RNA projection
seRNA <- readRDS("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/CITE/allconds/RNAprocess_DN1merge.rds")

proj_DN1_3 <- addGeneIntegrationMatrix(
    ArchRProj = proj_DN1_2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    #groupList = groupList,
    groupRNA = "CT",
    nameCell = "predCell_CT",
    nameGroup = "predGroup_CT",
    nameScore = "predScore_CT"
)



ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

pal_CT <- c(ggplotColours(n=8))
#pal_CT <- brewer.pal(8, "Set1"))
names(pal_CT) <- c("ETPs_proliferative","ETPs_quiescent","nonT_minor","nonT_potential","RunxKO_up","T_ILC","TcfKO_up","TleKO_up")
proj_DN1_3$CT <- proj_DN1_3$predGroup_CT

p1 <- plotEmbedding(proj_DN1_3, colorBy = "cellColData", name = "CT", pal = pal_CT, baseSize=2)
plotPDF(p1, name = "Plot-UMAP-RNAproject_CT.pdf", ArchRProj = proj_DN1_3, addDOC = FALSE, width = 5, height = 5)

plot(proj_DN1_3@embeddings$UMAP$df)


saveArchRProject(ArchRProj = proj_DN1_3, outputDirectory = "SaveDN1_3", load = FALSE)
proj_DN1_3 <- loadArchRProject(path = "SaveDN1_3")





pathToMacs2 <- findMacs2()
#pathToMacs2 <- "macs3"
#pathToMacs2 <- "~/anaconda3/envs/condaPY3/bin/macs3"
proj_DN1_4 <- addGroupCoverages(ArchRProj = proj_DN1_3, groupBy = "CT",force=T)
#proj_DN1_4tmp <- proj_DN1_4
#### now here

proj_DN1_4 <- addReproduciblePeakSet(
    ArchRProj = proj_DN1_4, 
    groupBy = "CT", 
    pathToMacs2 = "macs3",
    extsize=100,
    cutOff=0.01,
    shift=0,
    extendSummits=200,
    promoterRegion=c(2000,2000),
    genomeSize="mm",
    reproducibility = "(n+1)/2",
    threads = getArchRThreads()
)
proj_DN1_4 <- addPeakMatrix(proj_DN1_4)
allpeaks <- getPeakSet(proj_DN1_4)
proj_DN1_4 <- addImputeWeights(proj_DN1_4,seed=1)

saveRDS(allpeaks, file="allpeaks_v2.rds")
### change name:
#proj_DN1_4$CT5 <- proj_DN1_3$CT5#[which(proj_DN1_4$CT5 == "T_ILCs_sub")] <- "T_ILCs_GRhi"


######33 marker peak detection
#dir.create("CT_MKpeak/")
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_DN1_4, 
    useMatrix = "PeakMatrix", 
    groupBy = "CT",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, file="CT_MKpeak/CT_MKpeak.rds")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

output_diffpeak <- function(this_data,cellname){
    out_data <- cbind(as.vector(this_data@seqnames),
                      this_data@ranges@start,
                      this_data@ranges@start + this_data@ranges@width,
                      paste0(cellname,"_",seq(length(this_data@seqnames))),
                      this_data$Log2FC,
                      this_data$FDR)
    write.table(out_data,file=paste0("CT_MKpeak/",cellname,"_markerPeak.bed"),quote=F,row.names=F,col.names=F,sep="\t")
}
output_diffpeak(markerList$ETPs_quiescent,"ETPs_quiescent")
output_diffpeak(markerList$nonT_potential,"nonT_potential")
output_diffpeak(markerList$T_ILC,"T_ILC")
output_diffpeak(markerList$nonT_minor,"nonT_minor")
output_diffpeak(markerList$TleKO_up,"TleKO_up")
output_diffpeak(markerList$TcfKO_up,"TcfKO_up")
output_diffpeak(markerList$RunxKO_up,"RunxKO_up")


saveArchRProject(ArchRProj = proj_DN1_4, outputDirectory = "SaveDN1_4", load = FALSE)
proj_DN1_4 <- loadArchRProject(path = "SaveDN1_4")



############### nonT potential sub-analysis
proj_nonT <- proj_DN1_4[which(proj_DN1_4$CT=="nonT_potential")]

set.seed(1)
# dimensional reduction
proj_nonT <- addIterativeLSI(
    ArchRProj = proj_nonT,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    seed=1,force=T
)

# basic clustering 
proj_nonT <- addClusters(
    input = proj_nonT,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force=T,seed=1
)


# UMAP embedding
proj_nonT <- addUMAP(
    ArchRProj = proj_nonT, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",force=T
)



proj_nonT$nonTsub <- rep("NA",length(proj_nonT$Sample))
proj_nonT$nonTsub[which(proj_nonT$Clusters %in% paste0("C",c(1,3,4,5,7,8,9)))] <- "nonTsub1"
proj_nonT$nonTsub[which(proj_nonT$Clusters %in% paste0("C",c(2,6)))] <- "nonTsub2"
proj_nonT$nonTsub[which(proj_nonT$Clusters %in% paste0("C",c(10:15)))] <- "nonTsub3"

proj_nonT$BvsSample <- rep("NA",length(proj_nonT$Sample))
proj_nonT$BvsSample[which(proj_nonT$Clusters %in% paste0("C",c(1,2,3,4,5,6,7,8,9)) & proj_nonT$Sample == "WT")] <- "B_WT"
proj_nonT$BvsSample[which(proj_nonT$Clusters %in% paste0("C",c(1,2,3,4,5,6,7,8,9)) & proj_nonT$Sample == "Tcf1Lef1dKO")] <- "B_Tcf1KO"
proj_nonT$BvsSample[which(proj_nonT$Clusters %in% paste0("C",c(1,2,3,4,5,6,7,8,9)) & proj_nonT$Sample == "Tle34KO")] <- "B_TleKO"
proj_nonT$BvsSample[which(proj_nonT$Clusters %in% paste0("C",c(10:15)) & proj_nonT$Sample == "WT")] <- "nonB_WT"
proj_nonT$BvsSample[which(proj_nonT$Clusters %in% paste0("C",c(10:15)) & proj_nonT$Sample == "Tcf1Lef1dKO")] <- "nonB_Tcf1KO"
proj_nonT$BvsSample[which(proj_nonT$Clusters %in% paste0("C",c(10:15)) & proj_nonT$Sample == "Tle34KO")] <- "nonB_TleKO"



p1 <- plotEmbedding(ArchRProj = proj_nonT, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_nonT, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj_nonT, colorBy = "cellColData", name = "nonTsub", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj_nonT, colorBy = "cellColData", name = "BvsSample", embedding = "UMAP")
#plotPDF(p1,p2,p3, name = "Plot-UMAP-Sample-Clusters_LSI_nonTcellOnly.pdf", ArchRProj = proj_DN1_4, addDOC = FALSE, width = 5, height = 5)
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP-Sample-Clusters_LSI_nonTcellOnly_new.pdf", ArchRProj = proj_nonT, addDOC = FALSE, width = 5, height = 5)




peakdata <- read.table("/standard/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/scATAC_ArchR/combine4/CT_reads_new/combine4_mergePeak_CT.bed",row.names=4)
peakGR <- GRanges(seqnames=peakdata[,1],ranges=IRanges(peakdata[,2],peakdata[,3]))
proj_nonT <- addPeakSet(ArchRProj=proj_nonT, peakSet=peakGR,force=TRUE)
proj_nonT <- addPeakMatrix(proj_nonT)


saveArchRProject(ArchRProj = proj_nonT, outputDirectory = "proj_nonT", load = FALSE)
proj_nonT <- loadArchRProject(path = "proj_nonT")


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

a<-confmat(cbind(proj_DN1_4$Sample, proj_DN1_4$CT))
write.table(a,file="combine4scATAC_CTvsSample.txt",row.names=T,col.names=T,sep="\t",quote=F)



p1 <- plotEmbedding(ArchRProj = proj_DN1_4, colorBy = "cellColData", name = "nonTsub", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-nonTsubLabel_LSI.pdf", ArchRProj = proj_DN1_4, addDOC = FALSE, width = 5, height = 5)


nonTsub1_cell <- proj_nonT$cellNames[which(proj_nonT$nonTsub == "nonTsub1")]
nonTsub2_cell <- proj_nonT$cellNames[which(proj_nonT$nonTsub == "nonTsub2")]
nonTsub3_cell <- proj_nonT$cellNames[which(proj_nonT$nonTsub == "nonTsub3")]

proj_DN1_4$nonTsub <- rep("NA",length(proj_DN1_4$Sample))
proj_DN1_4$nonTsub[which(proj_DN1_4$cellNames %in% nonTsub1_cell)] <- "nonTsub1"
proj_DN1_4$nonTsub[which(proj_DN1_4$cellNames %in% nonTsub2_cell)] <- "nonTsub2"
proj_DN1_4$nonTsub[which(proj_DN1_4$cellNames %in% nonTsub3_cell)] <- "nonTsub3"



B_WT_cell <- proj_nonT$cellNames[which(proj_nonT$BvsSample == "B_WT")]
B_TcfKO_cell <- proj_nonT$cellNames[which(proj_nonT$BvsSample == "B_Tcf1KO")]
B_TleKO_cell <- proj_nonT$cellNames[which(proj_nonT$BvsSample == "B_TleKO")]
nonB_WT_cell <- proj_nonT$cellNames[which(proj_nonT$BvsSample == "nonB_WT")]
nonB_TcfKO_cell <- proj_nonT$cellNames[which(proj_nonT$BvsSample == "nonB_Tcf1KO")]
nonB_TleKO_cell <- proj_nonT$cellNames[which(proj_nonT$BvsSample == "nonB_TleKO")]

proj_DN1_4$BvsSample <- rep("NA",length(proj_DN1_4$Sample))
proj_DN1_4$BvsSample[which(proj_DN1_4$cellNames %in% B_WT_cell)] <- "B_WT"
proj_DN1_4$BvsSample[which(proj_DN1_4$cellNames %in% B_TcfKO_cell)] <- "B_TcfKO_cell"
proj_DN1_4$BvsSample[which(proj_DN1_4$cellNames %in% B_TleKO_cell)] <- "B_TleKO_cell"
proj_DN1_4$BvsSample[which(proj_DN1_4$cellNames %in% nonB_WT_cell)] <- "nonB_WT"
proj_DN1_4$BvsSample[which(proj_DN1_4$cellNames %in% nonB_TcfKO_cell)] <- "nonB_TcfKO_cell"
proj_DN1_4$BvsSample[which(proj_DN1_4$cellNames %in% nonB_TleKO_cell)] <- "nonB_TleKO_cell"



p <- plotEmbedding(
    ArchRProj = proj_DN1_4, 
    colorBy = "GeneScoreMatrix", 
    name = c("Pax5", "Elane", "Mpo","Irf8"), 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-ALLcells.pdf", 
    ArchRProj = proj_DN1_4, 
    addDOC = FALSE, width = 10, height = 10)

p <- plotEmbedding(
    ArchRProj = proj_nonT, 
    colorBy = "GeneScoreMatrix", 
    name = c("Pax5", "Elane", "Mpo","Irf8","Blnk","Cd19","Cd79a","Cd79b","Lyz2","Plac8","Cebpa","Cebpb"), 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-nonTcells.pdf", 
    ArchRProj = proj_nonT, 
    addDOC = FALSE, width = 10, height = 10)

markersGSnonT <- getMarkerFeatures(
    ArchRProj = proj_nonT, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "nonTsub",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)


markerList <- getMarkers(markersGSnonT, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

write.table(markerList$nonTsub1, file = "nonTsub1_MKG.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(markerList$nonTsub2, file = "nonTsub2_MKG.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(markerList$nonTsub3, file = "nonTsub3_MKG.txt", sep = "\t", row.names = FALSE, quote = FALSE)




combine2mklist <- function(list1, list2,name1,name2){
    a <- list1[[1]][,7:9]
    rownames(a) <- list1[[1]][,5]
    b <- list2[[1]][,7:9]
    rownames(b) <- list2[[1]][,5]
    idx <- intersect(rownames(a),rownames(b))
    cbdata <- cbind(a[idx,],b[idx,])
    colnames(cbdata) <- paste0(rep(c(name1,name2),each=3),"_",colnames(cbdata))
    return(cbdata)
}




MKG_RunxKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "GeneScoreMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "ETPs_quiescent")
MKG_RunxKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "GeneScoreMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "T_ILC")
MKG_RunxKOup_ETP_list <- getMarkers(MKG_RunxKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
MKG_RunxKOup_TILC_list <- getMarkers(MKG_RunxKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
MKG_RunxKOup <-combine2mklist(MKG_RunxKOup_ETP_list,MKG_RunxKOup_TILC_list,"vsETP","vsTILC")
write.table(MKG_RunxKOup,file="MKG_RunxKOup.txt",row.names=T,col.names=T,sep="\t",quote=F)

MKG_TcfKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "GeneScoreMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "ETPs_quiescent")
MKG_TcfKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "GeneScoreMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "T_ILC")
MKG_TcfKOup_ETP_list <- getMarkers(MKG_TcfKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
MKG_TcfKOup_TILC_list <- getMarkers(MKG_TcfKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
MKG_TcfKOup <-combine2mklist(MKG_TcfKOup_ETP_list,MKG_TcfKOup_TILC_list,"vsETP","vsTILC")
write.table(MKG_TcfKOup,file="MKG_TcfKOup.txt",row.names=T,col.names=T,sep="\t",quote=F)

MKG_TleKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "GeneScoreMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "ETPs_quiescent")
MKG_TleKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "GeneScoreMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "T_ILC")
MKG_TleKOup_ETP_list <- getMarkers(MKG_TleKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
MKG_TleKOup_TILC_list <- getMarkers(MKG_TleKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
MKG_TleKOup <-combine2mklist(MKG_TleKOup_ETP_list,MKG_TleKOup_TILC_list,"vsETP","vsTILC")
write.table(MKG_TleKOup,file="MKG_TleKOup.txt",row.names=T,col.names=T,sep="\t",quote=F)

MKG_nonT_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "GeneScoreMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "nonT_potential",   bgdGroups = "ETPs_quiescent")
MKG_nonT_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "GeneScoreMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "nonT_potential",   bgdGroups = "T_ILC")
MKG_nonT_ETP_list <- getMarkers(MKG_nonT_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
MKG_nonT_TILC_list <- getMarkers(MKG_nonT_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
MKG_nonT <-combine2mklist(MKG_nonT_ETP_list,MKG_nonT_TILC_list,"vsETP","vsTILC")
write.table(MKG_nonT,file="MKG_nonT.txt",row.names=T,col.names=T,sep="\t",quote=F)



######33 marker peak detection
#dir.create("markerPeak_CT/")
#markersPeaks <- getMarkerFeatures(
#    ArchRProj = proj_DN1_4, 
#    useMatrix = "PeakMatrix", 
#    groupBy = "CT",
#  bias = c("TSSEnrichment", "log10(nFrags)"),
#  testMethod = "wilcoxon"
#)
#saveRDS(markersPeaks, file="markerPeak_CT/markerPeak_CT.rds")
#markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
#markerList$


combine2mklist_peak <- function(list1, list2,name1,name2){
    a <- mcols(list1[[1]])
    rownames(a) <- paste0(as.character(seqnames(list1[[1]])),":",
          start(list1[[1]]),"-",
          end(list1[[1]]))
    
    b <- mcols(list2[[1]])
    rownames(b) <- paste0(as.character(seqnames(list2[[1]])),":",
          start(list2[[1]]),"-",
          end(list2[[1]]))
    idx <- intersect(rownames(a),rownames(b))
    cbdata <- cbind(a[idx,],b[idx,])
    colnames(cbdata) <- paste0(rep(c(name1,name2),each=3),"_",colnames(cbdata))
    return(cbdata)
}


MKpeak_RunxKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_RunxKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "T_ILC")
MKpeak_RunxKOup_ETP_list <- getMarkers(MKpeak_RunxKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
MKpeak_RunxKOup_TILC_list <- getMarkers(MKpeak_RunxKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")

MKpeak_TcfKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_TcfKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "T_ILC")
MKpeak_TcfKOup_ETP_list <- getMarkers(MKpeak_TcfKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
MKpeak_TcfKOup_TILC_list <- getMarkers(MKpeak_TcfKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")

MKpeak_TleKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_TleKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "T_ILC")
MKpeak_TleKOup_ETP_list <- getMarkers(MKpeak_TleKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
MKpeak_TleKOup_TILC_list <- getMarkers(MKpeak_TleKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")

# down
MKpeak_ETP_RunxKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "RunxKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_RunxKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "RunxKO_up",   useGroups = "T_ILC")
MKpeak_ETP_RunxKOup_list <- getMarkers(MKpeak_ETP_RunxKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
MKpeak_TILC_RunxKOup_list <- getMarkers(MKpeak_TILC_RunxKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")

MKpeak_ETP_TcfKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TcfKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_TcfKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TcfKO_up",   useGroups = "T_ILC")
MKpeak_ETP_TcfKOup_list <- getMarkers(MKpeak_ETP_TcfKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
MKpeak_TILC_TcfKOup_list <- getMarkers(MKpeak_TILC_TcfKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")

MKpeak_ETP_TleKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TleKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_TleKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TleKO_up",   useGroups = "T_ILC")
MKpeak_ETP_TleKOup_list <- getMarkers(MKpeak_ETP_TleKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
MKpeak_TILC_TleKOup_list <- getMarkers(MKpeak_TILC_TleKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")


transformBED <- function(indata){
    #regions <- rownames(indata)
    #split_regions <- strsplit(regions, ":|-")
    # Extract chromosomes, start, and end positions
    chromosomes <- as.vector(indata[,"seqnames"])#sapply(split_regions, `[`, 1)
    starts <- as.vector(indata[,"start"])#as.integer(sapply(split_regions, `[`, 2))
    ends <- as.vector(indata[,"end"])#as.integer(sapply(split_regions, `[`, 3))
    newdata <-  cbind(chromosomes,starts,ends)
    return(newdata)
}

write.table(as.matrix(transformBED(MKpeak_RunxKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_RunxKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_RunxKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_RunxKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TcfKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TcfKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TcfKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TcfKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TleKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TleKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TleKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TleKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)

write.table(as.matrix(transformBED(MKpeak_ETP_RunxKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_ETP_RunxKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_RunxKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TILC_RunxKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_ETP_TcfKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_ETP_TcfKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_TcfKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TILC_TcfKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_ETP_TleKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_ETP_TleKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_TleKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TILC_TleKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)


#MKpeak_nonT_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "nonT_potential",   bgdGroups = "ETPs_quiescent")
#MKpeak_nonT_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_4,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "nonT_potential",   bgdGroups = "T_ILC")
#MKpeak_nonT_ETP_list <- getMarkers(MKpeak_nonT_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
#MKpeak_nonT_TILC_list <- getMarkers(MKpeak_nonT_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
#MKpeak_nonT <-combine2mklist_peak(MKpeak_nonT_ETP_list,MKpeak_nonT_TILC_list,"vsETP","vsTILC")
#write.table(MKpeak_nonT,file="MKpeak_nonT.txt",row.names=T,col.names=T,sep="\t",quote=F)







#proj_nonT <- addGroupCoverages(ArchRProj = proj_nonT, groupBy = "BvsSample",force=T)
#
#proj_nonT <- addReproduciblePeakSet(
#    ArchRProj = proj_nonT, 
#    groupBy = "BvsSample", 
#    pathToMacs2 = "macs3",
#    extsize=100,
#    cutOff=0.01,
#    shift=0,
#    extendSummits=200,
#    promoterRegion=c(2000,2000),
#    genomeSize="mm",
#    reproducibility = "(n+1)/2",
#    threads = getArchRThreads()
#)
#proj_nonT <- addPeakMatrix(proj_nonT)
#allpeaks <- getPeakSet(proj_nonT)
#proj_nonT <- addImputeWeights(proj_nonT,seed=1)

#saveRDS(allpeaks, file="allpeaks_v2.rds")

MKpeak_B_WTvsTcf1KO     <- getMarkerFeatures(ArchRProj = proj_nonT,useGroups = "B_WT",bgdGroups = "B_Tcf1KO", useMatrix = "PeakMatrix",groupBy = "BvsSample",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"))
MKpeak_B_WTvsTleKO      <- getMarkerFeatures(ArchRProj = proj_nonT,useGroups = "B_WT",bgdGroups = "B_TleKO", useMatrix = "PeakMatrix",groupBy = "BvsSample",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"))
MKpeak_nonB_WTvsTcf1KO  <- getMarkerFeatures(ArchRProj = proj_nonT,useGroups = "nonB_WT",bgdGroups = "nonB_Tcf1KO", useMatrix = "PeakMatrix",groupBy = "BvsSample",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"))
MKpeak_nonB_WTvsTleKO   <- getMarkerFeatures(ArchRProj = proj_nonT,useGroups = "nonB_WT",bgdGroups = "nonB_TleKO", useMatrix = "PeakMatrix",groupBy = "BvsSample",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"))
MKpeak_B_Tcf1KOvsWT     <- getMarkerFeatures(ArchRProj = proj_nonT,bgdGroups = "B_WT",useGroups = "B_Tcf1KO", useMatrix = "PeakMatrix",groupBy = "BvsSample",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"))
MKpeak_B_TleKOvsWT      <- getMarkerFeatures(ArchRProj = proj_nonT,bgdGroups = "B_WT",useGroups = "B_TleKO", useMatrix = "PeakMatrix",groupBy = "BvsSample",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"))
MKpeak_nonB_Tcf1KOvsWT  <- getMarkerFeatures(ArchRProj = proj_nonT,bgdGroups = "nonB_WT",useGroups = "nonB_Tcf1KO", useMatrix = "PeakMatrix",groupBy = "BvsSample",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"))
MKpeak_nonB_TleKOvsWT   <- getMarkerFeatures(ArchRProj = proj_nonT,bgdGroups = "nonB_WT",useGroups = "nonB_TleKO", useMatrix = "PeakMatrix",groupBy = "BvsSample",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"))

MKpeak_B_WTvsTcf1KO_list <- getMarkers(MKpeak_B_WTvsTcf1KO    , cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_B_WTvsTleKO_list <- getMarkers(MKpeak_B_WTvsTleKO     , cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_nonB_WTvsTcf1KO_list <- getMarkers(MKpeak_nonB_WTvsTcf1KO , cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_nonB_WTvsTleKO_list <- getMarkers(MKpeak_nonB_WTvsTleKO  , cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_B_Tcf1KOvsWT_list <- getMarkers(MKpeak_B_Tcf1KOvsWT    , cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_B_TleKOvsWT_list <- getMarkers(MKpeak_B_TleKOvsWT     , cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_nonB_Tcf1KOvsWT_list <- getMarkers(MKpeak_nonB_Tcf1KOvsWT , cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_nonB_TleKOvsWT_list <- getMarkers(MKpeak_nonB_TleKOvsWT  , cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")



write.table(as.matrix(transformBED(MKpeak_B_WTvsTcf1KO_list[[1]])),file="BvsSample_MKpeak/MKpeak_B_WTvsTcf1KO.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_B_WTvsTleKO_list[[1]])),file="BvsSample_MKpeak/MKpeak_B_WTvsTleKO.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_nonB_WTvsTcf1KO_list[[1]])),file="BvsSample_MKpeak/MKpeak_nonB_WTvsTcf1KO.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_nonB_WTvsTleKO_list[[1]])),file="BvsSample_MKpeak/MKpeak_nonB_WTvsTleKO.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_B_Tcf1KOvsWT_list[[1]])),file="BvsSample_MKpeak/MKpeak_B_Tcf1KOvsWT.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_B_TleKOvsWT_list[[1]])),file="BvsSample_MKpeak/MKpeak_B_TleKOvsWT.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_nonB_Tcf1KOvsWT_list[[1]])),file="BvsSample_MKpeak/MKpeak_nonB_Tcf1KOvsWT.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_nonB_TleKOvsWT_list[[1]])),file="BvsSample_MKpeak/MKpeak_nonB_TleKOvsWT.bed",row.names=F,col.names=F,sep="\t",quote=F)




proj_nonT <- addMotifAnnotations(ArchRProj = proj_nonT, motifSet = "homer", name = "Motif",force = TRUE)

MotifMKpeak_B_WTvsTcf1KO <- peakAnnoEnrichment(seMarker = MKpeak_B_WTvsTcf1KO ,ArchRProj = proj_nonT,peakAnnotation = "Motif",cutOff = "FDR <= 0.01 & Log2FC >= 1")
MotifMKpeak_B_WTvsTleKO <- peakAnnoEnrichment(seMarker = MKpeak_B_WTvsTleKO ,ArchRProj = proj_nonT,peakAnnotation = "Motif",cutOff = "FDR <= 0.01 & Log2FC >= 1")
MotifMKpeak_nonB_WTvsTcf1KO <- peakAnnoEnrichment(seMarker = MKpeak_nonB_WTvsTcf1KO ,ArchRProj = proj_nonT,peakAnnotation = "Motif",cutOff = "FDR <= 0.01 & Log2FC >= 1")
MotifMKpeak_nonB_WTvsTleKO <- peakAnnoEnrichment(seMarker = MKpeak_nonB_WTvsTleKO ,ArchRProj = proj_nonT,peakAnnotation = "Motif",cutOff = "FDR <= 0.01 & Log2FC >= 1")
MotifMKpeak_B_Tcf1KOvsWT <- peakAnnoEnrichment(seMarker = MKpeak_B_Tcf1KOvsWT ,ArchRProj = proj_nonT,peakAnnotation = "Motif",cutOff = "FDR <= 0.01 & Log2FC >= 1")
MotifMKpeak_B_TleKOvsWT <- peakAnnoEnrichment(seMarker = MKpeak_B_TleKOvsWT ,ArchRProj = proj_nonT,peakAnnotation = "Motif",cutOff = "FDR <= 0.01 & Log2FC >= 1")
MotifMKpeak_nonB_Tcf1KOvsWT <- peakAnnoEnrichment(seMarker = MKpeak_nonB_Tcf1KOvsWT ,ArchRProj = proj_nonT,peakAnnotation = "Motif",cutOff = "FDR <= 0.01 & Log2FC >= 1")
MotifMKpeak_nonB_TleKOvsWT <- peakAnnoEnrichment(seMarker = MKpeak_nonB_TleKOvsWT ,ArchRProj = proj_nonT,peakAnnotation = "Motif",cutOff = "FDR <= 0.01 & Log2FC >= 1")

idx <- rownames(MotifMKpeak_B_WTvsTcf1KO@assays@data[[1]])
outdata <- cbind(MotifMKpeak_B_WTvsTcf1KO@assays@data[[1]][idx,1],
MotifMKpeak_B_WTvsTleKO@assays@data[[1]][idx,1],
MotifMKpeak_nonB_WTvsTcf1KO@assays@data[[1]][idx,1],
MotifMKpeak_nonB_WTvsTleKO@assays@data[[1]][idx,1],
MotifMKpeak_B_Tcf1KOvsWT@assays@data[[1]][idx,1],
MotifMKpeak_B_TleKOvsWT@assays@data[[1]][idx,1],
MotifMKpeak_nonB_Tcf1KOvsWT@assays@data[[1]][idx,1],
MotifMKpeak_nonB_TleKOvsWT@assays@data[[1]][idx,1])

rownames(outdata) <- idx
colnames(outdata) <- c("B_WTvsTcf1KO","B_WTvsTleKO","nonB_WTvsTcf1KO","nonB_WTvsTleKO","B_Tcf1KOvsWT","B_TleKOvsWT","nonB_Tcf1KOvsWT","nonB_TleKOvsWT")

write.table(outdata,file="BvsSample_motifEnrichment.txt",row.names=T,col.names=T,sep="\t",quote=F)
filterdata= outdata[which(apply(outdata,1,max) >= 10),]
write.table(filterdata,file="BvsSample_motifEnrichment_filter10.txt",row.names=T,col.names=T,sep="\t",quote=F)
filterdata[,c(1,2,4,5,6,8)]

library(ggplot2)
library(dplyr)
library(pheatmap)
your_matrix <- as.matrix(filterdata[,c(1,2,4,5,6,8)])
row_clusters <- hclust(dist(your_matrix))
ordered_matrix <- your_matrix[row_clusters$order, ]
data_for_plot <- as.data.frame(as.table(as.matrix(ordered_matrix)))
colnames(data_for_plot) <- c("RowName", "ColName", "Value")
heatmap_plot <- ggplot(data_for_plot, aes(x = ColName, y = RowName, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x labels for clarity
        axis.title = element_blank()) # Remove axis titles
ggsave(heatmap_plot,file="BvsSample_motifEnrichment_filter10_heatmap.pdf")
# P
########## PRINT READS
options(scipen = 999)
Cdata <- getCellColData(ArchRProj = proj_DN1_4, select = NULL, drop = FALSE)

# all CT
# "nonT_potential"
# "TleKO_up"      
# "ETPs_quiescent"
# "nonT_minor"    
# "T_ILC"         
# "RunxKO_up"     
# "TcfKO_up"      

### ETPq reads & onlyWT
cellnames_WT<-rownames(Cdata[which(Cdata$CT== "ETPs_quiescent" & Cdata$Sample=="WT"),])
frag_WT<-getFragmentsFromArrow(
    ArrowFile = "WT.arrow",
    cellNames = cellnames_WT,
    verbose = TRUE,
)
fragOut <- c(frag_WT)
df <- data.frame(seqnames=seqnames(fragOut),
starts=start(fragOut)-1,
ends=end(fragOut))
options(scipen = 999)
write.table(df, file="CT_reads_new/ETPs_quiescent_read_v2.bed", quote=F, sep="\t", row.names=F, col.names=F)

### TleKO_up reads & only TleKO
cellnames_Tle34KO<-rownames(Cdata[which(Cdata$CT== "TleKO_up" & Cdata$Sample=="Tle34KO"),])
frag_Tle34KO<-getFragmentsFromArrow(
    ArrowFile = "Tle34KO.arrow",
    cellNames = cellnames_Tle34KO,
    verbose = TRUE,
)
fragOut <- c(frag_Tle34KO)
df <- data.frame(seqnames=seqnames(fragOut),
starts=start(fragOut)-1,
ends=end(fragOut))
options(scipen = 999)
write.table(df, file="CT_reads_new/TleKO_up_read_v2.bed", quote=F, sep="\t", row.names=F, col.names=F)


for(ct in c("nonT_potential","T_ILC","RunxKO_up","TcfKO_up")){
    print(ct)
    cellnames_WT<-rownames(Cdata[which(Cdata$CT==ct & Cdata$Sample=="WT"),])
    frag_WT<-getFragmentsFromArrow(
        ArrowFile = "WT.arrow",
        cellNames = cellnames_WT,
        verbose = TRUE,
    )
    cellnames_Tcf1Lef1dKO<-rownames(Cdata[which(Cdata$CT==ct & Cdata$Sample=="Tcf1Lef1dKO"),])
    frag_Tcf1Lef1dKO<-getFragmentsFromArrow(
        ArrowFile = "Tcf1Lef1dKO.arrow",
        cellNames = cellnames_Tcf1Lef1dKO,
        verbose = TRUE,
    )
    cellnames_Tle34KO<-rownames(Cdata[which(Cdata$CT==ct & Cdata$Sample=="Tle34KO"),])
    frag_Tle34KO<-getFragmentsFromArrow(
        ArrowFile = "Tle34KO.arrow",
        cellNames = cellnames_Tle34KO,
        verbose = TRUE,
    )
    cellnames_Runx1KO<-rownames(Cdata[which(Cdata$CT==ct & Cdata$Sample=="Runx1KO"),])
    frag_Runx1KO<-getFragmentsFromArrow(
        ArrowFile = "Runx1KO.arrow",
        cellNames = cellnames_Runx1KO,
        verbose = TRUE,
    )
    fragOut <- c(frag_WT, frag_Tcf1Lef1dKO,frag_Tle34KO,frag_Runx1KO)
    
    df <- data.frame(seqnames=seqnames(fragOut),
    starts=start(fragOut)-1,
    ends=end(fragOut))
    write.table(df, file=paste0("CT_reads_new/",ct,"_read.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}




########## PRINT READS for nonT groups
options(scipen = 999)

Cdata <- getCellColData(ArchRProj = proj_nonT, select = NULL, drop = FALSE)
for(c in unique(Cdata$BvsSample)){
    print(c)
    cellnames_WT<-rownames(Cdata[which(Cdata$BvsSample==c & Cdata$Sample=="WT"),])
    frag_WT<-getFragmentsFromArrow(
        ArrowFile = "WT.arrow",
        cellNames = cellnames_WT,
        verbose = TRUE,
    )
    cellnames_Tcf1Lef1dKO<-rownames(Cdata[which(Cdata$BvsSample==c & Cdata$Sample=="Tcf1Lef1dKO"),])
    frag_Tcf1Lef1dKO<-getFragmentsFromArrow(
        ArrowFile = "Tcf1Lef1dKO.arrow",
        cellNames = cellnames_Tcf1Lef1dKO,
        verbose = TRUE,
    )
    cellnames_Tle34KO<-rownames(Cdata[which(Cdata$BvsSample==c & Cdata$Sample=="Tle34KO"),])
    frag_Tle34KO<-getFragmentsFromArrow(
        ArrowFile = "Tle34KO.arrow",
        cellNames = cellnames_Tle34KO,
        verbose = TRUE,
    )
    cellnames_Runx1KO<-rownames(Cdata[which(Cdata$BvsSample==c & Cdata$Sample=="Runx1KO"),])
    frag_Runx1KO<-getFragmentsFromArrow(
        ArrowFile = "Runx1KO.arrow",
        cellNames = cellnames_Runx1KO,
        verbose = TRUE,
    )
    fragOut <- c(frag_WT, frag_Tcf1Lef1dKO,frag_Tle34KO,frag_Runx1KO)
    
    df <- data.frame(seqnames=seqnames(fragOut),
    starts=start(fragOut)-1,
    ends=end(fragOut))
    write.table(df, file=paste0("nonT_BvsSample_reads/",c,"_read.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}




####################### add cb reads macs2 peaks as peakset
peakdata <- read.table("/standard/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/scATAC_ArchR/combine4/CT_reads_new/combine4_mergePeak_CT.bed",row.names=4)
peakGR <- GRanges(seqnames=peakdata[,1],ranges=IRanges(peakdata[,2],peakdata[,3]))
proj_DN1_5 <- addPeakSet(ArchRProj=proj_DN1_4, peakSet=peakGR,force=TRUE)
proj_DN1_5 <- addPeakMatrix(proj_DN1_5)



########### wait to decide whether use TleKOup celltype & TleKO sample for cells in marker peak detection
clean4CT <- rep("NA", length(proj_DN1_5$CT))
clean4CT[which(proj_DN1_5$CT == "nonT_potential")] <- "nonT_potential"
clean4CT[which(proj_DN1_5$CT == "T_ILC")] <- "T_ILC"
clean4CT[which(proj_DN1_5$CT == "RunxKO_up")] <- "RunxKO_up"
clean4CT[which(proj_DN1_5$CT == "TcfKO_up")] <- "TcfKO_up"
clean4CT[which(proj_DN1_5$CT == "TleKO_up" & proj_DN1_5$Sample == "Tle34KO")] <- "TleKO_up"
clean4CT[which(proj_DN1_5$CT == "ETPs_quiescent" & proj_DN1_5$Sample == "WT")] <- "ETPs_quiescent"

proj_DN1_5$clean4CT <- clean4CT


saveArchRProject(ArchRProj = proj_DN1_5, outputDirectory = "SaveDN1_5", load = FALSE)
proj_DN1_5 <- loadArchRProject(path = "SaveDN1_5")


MKpeak_RunxKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_RunxKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "T_ILC")
MKpeak_RunxKOup_ETP_list <- getMarkers(MKpeak_RunxKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_RunxKOup_TILC_list <- getMarkers(MKpeak_RunxKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")

MKpeak_TcfKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_TcfKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "T_ILC")
MKpeak_TcfKOup_ETP_list <- getMarkers(MKpeak_TcfKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_TcfKOup_TILC_list <- getMarkers(MKpeak_TcfKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")

MKpeak_TleKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_TleKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "T_ILC")
MKpeak_TleKOup_ETP_list <- getMarkers(MKpeak_TleKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_TleKOup_TILC_list <- getMarkers(MKpeak_TleKOup_TILC, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")

# down
MKpeak_ETP_RunxKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "RunxKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_RunxKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "RunxKO_up",   useGroups = "T_ILC")
MKpeak_ETP_RunxKOup_list <- getMarkers(MKpeak_ETP_RunxKOup, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_TILC_RunxKOup_list <- getMarkers(MKpeak_TILC_RunxKOup, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")

MKpeak_ETP_TcfKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TcfKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_TcfKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TcfKO_up",   useGroups = "T_ILC")
MKpeak_ETP_TcfKOup_list <- getMarkers(MKpeak_ETP_TcfKOup, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_TILC_TcfKOup_list <- getMarkers(MKpeak_TILC_TcfKOup, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")

MKpeak_ETP_TleKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TleKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_TleKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TleKO_up",   useGroups = "T_ILC")
MKpeak_ETP_TleKOup_list <- getMarkers(MKpeak_ETP_TleKOup, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")
MKpeak_TILC_TleKOup_list <- getMarkers(MKpeak_TILC_TleKOup, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625")



#MKpeak_RunxKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "ETPs_quiescent")
#MKpeak_TcfKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "ETPs_quiescent")
#MKpeak_TleKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "ETPs_quiescent")
#MKpeak_ETP_RunxKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "RunxKO_up",   useGroups = "ETPs_quiescent")
#MKpeak_ETP_TcfKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TcfKO_up",   useGroups = "ETPs_quiescent")
#MKpeak_ETP_TleKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TleKO_up",   useGroups = "ETPs_quiescent")
#
#MKpeak_RunxKOup_ETP_list <- getMarkers(MKpeak_RunxKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
#MKpeak_TcfKOup_ETP_list <- getMarkers(MKpeak_TcfKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
#MKpeak_TleKOup_ETP_list <- getMarkers(MKpeak_TleKOup_ETP, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
#MKpeak_ETP_RunxKOup_list <- getMarkers(MKpeak_ETP_RunxKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
#MKpeak_ETP_TcfKOup_list <- getMarkers(MKpeak_ETP_TcfKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
#MKpeak_ETP_TleKOup_list <- getMarkers(MKpeak_ETP_TleKOup, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")

peakAnno <- read.table("/standard/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/scATAC_ArchR/combine4/CT_reads_new/combine4_mergePeak_CT.bed",row.names=4)
peak1kb <- cbind(peakAnno[,1],round((peakAnno[,2]+peakAnno[,3])/2)-1000,round((peakAnno[,2]+peakAnno[,3])/2)+1000,rownames(peakAnno))
peakCenter <- cbind(peakAnno[,1],round((peakAnno[,2]+peakAnno[,3])/2),round((peakAnno[,2]+peakAnno[,3])/2)+1,rownames(peakAnno))

ETP_vs_RunxKO_status <- rep(0, nrow(peak1kb))
ETP_vs_TleKO_status <- rep(0, nrow(peak1kb))
ETP_vs_TcfKO_status <- rep(0, nrow(peak1kb))
#ETP_vs_RunxKO_lfc <- rep(0, nrow(peak1kb))
#ETP_vs_TleKO_lfc <- rep(0, nrow(peak1kb))
#ETP_vs_TcfKO_lfc <- rep(0, nrow(peak1kb))
ETP_vs_RunxKO_status[as.numeric(rownames(MKpeak_RunxKOup_ETP_list$RunxKO_up))] <- -1
ETP_vs_RunxKO_status[as.numeric(rownames(MKpeak_ETP_RunxKOup_list$ETPs_quiescent))] <- 1
ETP_vs_TleKO_status[as.numeric(rownames(MKpeak_TleKOup_ETP_list$TleKO_up))] <- -1
ETP_vs_TleKO_status[as.numeric(rownames(MKpeak_ETP_TleKOup_list$ETPs_quiescent))] <- 1
ETP_vs_TcfKO_status[as.numeric(rownames(MKpeak_TcfKOup_ETP_list$TcfKO_up))] <- -1
ETP_vs_TcfKO_status[as.numeric(rownames(MKpeak_ETP_TcfKOup_list$ETPs_quiescent))] <- 1


outdata <- cbind(
      ETP_vs_RunxKO_status, 
      ETP_vs_TleKO_status, 
      ETP_vs_TcfKO_status, 
      MKpeak_ETP_RunxKOup@assays@data$Log2FC,
      MKpeak_ETP_TleKOup@assays@data$Log2FC,
      MKpeak_ETP_TcfKOup@assays@data$Log2FC)
colnames(outdata) <- c("ETP_vs_RunxKO_status","ETP_vs_TleKO_status","ETP_vs_TcfKO_status",
                       "ETP_vs_RunxKO_lfc","ETP_vs_TleKO_lfc","ETP_vs_TcfKO_lfc")
rownames(outdata) <- rownames(peakAnno)

options(scipen = 999)
write.table(outdata,file="mergePeak_CT_ETPvsKOdiff.txt",row.names=T,col.names=T,sep="\t",quote=F)
write.table(peak1kb, file="mergePeak_CT_center1kb.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(peakCenter, file="mergePeak_CT_center.bed",row.names=F,col.names=F,sep="\t",quote=F)


transformBED <- function(indata){
    #regions <- rownames(indata)
    #split_regions <- strsplit(regions, ":|-")
    # Extract chromosomes, start, and end positions
    chromosomes <- as.vector(indata[,"seqnames"])#sapply(split_regions, `[`, 1)
    starts <- as.vector(indata[,"start"])#as.integer(sapply(split_regions, `[`, 2))
    ends <- as.vector(indata[,"end"])#as.integer(sapply(split_regions, `[`, 3))
    newdata <-  cbind(chromosomes,starts,ends)
    return(newdata)
}

write.table(as.matrix(transformBED(MKpeak_RunxKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_RunxKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_RunxKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_RunxKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TcfKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TcfKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TcfKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TcfKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TleKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TleKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TleKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TleKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)

write.table(as.matrix(transformBED(MKpeak_ETP_RunxKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_ETP_RunxKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_RunxKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TILC_RunxKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_ETP_TcfKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_ETP_TcfKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_TcfKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TILC_TcfKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_ETP_TleKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_ETP_TleKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_TleKOup_list[[1]])),file="CT_MKpeak_pairwise/MKpeak_TILC_TleKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)





############ FC2 ver

MKpeak_RunxKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_RunxKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "RunxKO_up",   bgdGroups = "T_ILC")
MKpeak_RunxKOup_ETP_list <- getMarkers(MKpeak_RunxKOup_ETP, cutOff = "FDR <= 0.001 & Log2FC >= 1")
MKpeak_RunxKOup_TILC_list <- getMarkers(MKpeak_RunxKOup_TILC, cutOff = "FDR <= 0.001 & Log2FC >= 1")

MKpeak_TcfKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_TcfKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TcfKO_up",   bgdGroups = "T_ILC")
MKpeak_TcfKOup_ETP_list <- getMarkers(MKpeak_TcfKOup_ETP, cutOff = "FDR <= 0.001 & Log2FC >= 1")
MKpeak_TcfKOup_TILC_list <- getMarkers(MKpeak_TcfKOup_TILC, cutOff = "FDR <= 0.001 & Log2FC >= 1")

MKpeak_TleKOup_ETP <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "ETPs_quiescent")
MKpeak_TleKOup_TILC <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   useGroups = "TleKO_up",   bgdGroups = "T_ILC")
MKpeak_TleKOup_ETP_list <- getMarkers(MKpeak_TleKOup_ETP, cutOff = "FDR <= 0.001 & Log2FC >= 1")
MKpeak_TleKOup_TILC_list <- getMarkers(MKpeak_TleKOup_TILC, cutOff = "FDR <= 0.001 & Log2FC >= 1")

# down
MKpeak_ETP_RunxKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "RunxKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_RunxKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "RunxKO_up",   useGroups = "T_ILC")
MKpeak_ETP_RunxKOup_list <- getMarkers(MKpeak_ETP_RunxKOup, cutOff = "FDR <= 0.001 & Log2FC >= 1")
MKpeak_TILC_RunxKOup_list <- getMarkers(MKpeak_TILC_RunxKOup, cutOff = "FDR <= 0.001 & Log2FC >= 1")

MKpeak_ETP_TcfKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TcfKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_TcfKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TcfKO_up",   useGroups = "T_ILC")
MKpeak_ETP_TcfKOup_list <- getMarkers(MKpeak_ETP_TcfKOup, cutOff = "FDR <= 0.001 & Log2FC >= 1")
MKpeak_TILC_TcfKOup_list <- getMarkers(MKpeak_TILC_TcfKOup, cutOff = "FDR <= 0.001 & Log2FC >= 1")

MKpeak_ETP_TleKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TleKO_up",   useGroups = "ETPs_quiescent")
MKpeak_TILC_TleKOup <- getMarkerFeatures(  ArchRProj = proj_DN1_5,   useMatrix = "PeakMatrix",   groupBy = "clean4CT",   testMethod = "wilcoxon",   bias = c("TSSEnrichment","log10(nFrags)"),   bgdGroups = "TleKO_up",   useGroups = "T_ILC")
MKpeak_ETP_TleKOup_list <- getMarkers(MKpeak_ETP_TleKOup, cutOff = "FDR <= 0.001 & Log2FC >= 1")
MKpeak_TILC_TleKOup_list <- getMarkers(MKpeak_TILC_TleKOup, cutOff = "FDR <= 0.001 & Log2FC >= 1")



peakAnno <- read.table("/standard/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/scATAC_ArchR/combine4/CT_reads_new/combine4_mergePeak_CT.bed",row.names=4)
peak1kb <- cbind(peakAnno[,1],round((peakAnno[,2]+peakAnno[,3])/2)-1000,round((peakAnno[,2]+peakAnno[,3])/2)+1000,rownames(peakAnno))
peakCenter <- cbind(peakAnno[,1],round((peakAnno[,2]+peakAnno[,3])/2),round((peakAnno[,2]+peakAnno[,3])/2)+1,rownames(peakAnno))

ETP_vs_RunxKO_status <- rep(0, nrow(peak1kb))
ETP_vs_TleKO_status <- rep(0, nrow(peak1kb))
ETP_vs_TcfKO_status <- rep(0, nrow(peak1kb))
#ETP_vs_RunxKO_lfc <- rep(0, nrow(peak1kb))
#ETP_vs_TleKO_lfc <- rep(0, nrow(peak1kb))
#ETP_vs_TcfKO_lfc <- rep(0, nrow(peak1kb))
ETP_vs_RunxKO_status[as.numeric(rownames(MKpeak_RunxKOup_ETP_list$RunxKO_up))] <- -1
ETP_vs_RunxKO_status[as.numeric(rownames(MKpeak_ETP_RunxKOup_list$ETPs_quiescent))] <- 1
ETP_vs_TleKO_status[as.numeric(rownames(MKpeak_TleKOup_ETP_list$TleKO_up))] <- -1
ETP_vs_TleKO_status[as.numeric(rownames(MKpeak_ETP_TleKOup_list$ETPs_quiescent))] <- 1
ETP_vs_TcfKO_status[as.numeric(rownames(MKpeak_TcfKOup_ETP_list$TcfKO_up))] <- -1
ETP_vs_TcfKO_status[as.numeric(rownames(MKpeak_ETP_TcfKOup_list$ETPs_quiescent))] <- 1


outdata <- cbind(
      ETP_vs_RunxKO_status, 
      ETP_vs_TleKO_status, 
      ETP_vs_TcfKO_status, 
      MKpeak_ETP_RunxKOup@assays@data$Log2FC,
      MKpeak_ETP_TleKOup@assays@data$Log2FC,
      MKpeak_ETP_TcfKOup@assays@data$Log2FC)
colnames(outdata) <- c("ETP_vs_RunxKO_status","ETP_vs_TleKO_status","ETP_vs_TcfKO_status",
                       "ETP_vs_RunxKO_lfc","ETP_vs_TleKO_lfc","ETP_vs_TcfKO_lfc")
rownames(outdata) <- rownames(peakAnno)

options(scipen = 999)
write.table(outdata,file="mergePeak_CT_ETPvsKOdiff_FC2.txt",row.names=T,col.names=T,sep="\t",quote=F)



transformBED <- function(indata){
    #regions <- rownames(indata)
    #split_regions <- strsplit(regions, ":|-")
    # Extract chromosomes, start, and end positions
    chromosomes <- as.vector(indata[,"seqnames"])#sapply(split_regions, `[`, 1)
    starts <- as.vector(indata[,"start"])#as.integer(sapply(split_regions, `[`, 2))
    ends <- as.vector(indata[,"end"])#as.integer(sapply(split_regions, `[`, 3))
    newdata <-  cbind(chromosomes,starts,ends)
    return(newdata)
}

write.table(as.matrix(transformBED(MKpeak_RunxKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_RunxKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_RunxKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_RunxKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TcfKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_TcfKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TcfKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_TcfKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TleKOup_ETP_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_TleKOup_ETP.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TleKOup_TILC_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_TleKOup_TILC.bed",row.names=F,col.names=F,sep="\t",quote=F)

write.table(as.matrix(transformBED(MKpeak_ETP_RunxKOup_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_ETP_RunxKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_RunxKOup_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_TILC_RunxKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_ETP_TcfKOup_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_ETP_TcfKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_TcfKOup_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_TILC_TcfKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_ETP_TleKOup_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_ETP_TleKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeak_TILC_TleKOup_list[[1]])),file="CT_MKpeak_pairwise_FC2/MKpeak_TILC_TleKOup.bed",row.names=F,col.names=F,sep="\t",quote=F)










a <- cbind(MKpeak_RunxKOup_ETP@assays@data$Log2FC,MKpeak_ETP_RunxKOup@assays@data$Log2FC )


dim(MKpeak_RunxKOup_ETP@assays@data$Log2FC)
dim(MKpeak_ETP_RunxKOup@assays@data$Log2FC)

proj_DN1_5@peakSet
plot(MKpeak_RunxKOup_ETP@assays@data$Log2FC[,1],
     MKpeak_ETP_RunxKOup@assays@data$Log2FC[,1],pch="."
     )
























