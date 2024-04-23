library(ArchR) # version 
set.seed(1)
addArchRGenome("mm10")
addArchRThreads(1)

inputFiles <- c(
    "/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/cellrangerATAC/WT/DN1_WT/outs/fragments.tsv.gz")#,
    #"/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Data/process/scATAC/DN1/cellrangerATAC/Tcf1Lef1dko/DN1_Tcf1Lef1dko/outs/fragments.tsv.gz")    
#names(inputFiles)<-c("WT","Tcf1Lef1dko")
names(inputFiles)<-c("WT")


#doubScores2 <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
#  LSIMethod = 1,
#  force=TRUE
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
#KOcells <- proj_DN1_2$cellNames[which(sampleTags == "Tcf1Lef1dko")]
proj_DN1_2_WT <- proj_DN1_2[WTcells,]
#proj_DN1_2_KO <- proj_DN1_2[KOcells,]


#p1 <- plotFragmentSizes(ArchRProj = proj_DN1_2_WT)
#p2 <- plotTSSEnrichment(ArchRProj = proj_DN1_2_WT)
#plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_WTonly.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)

#p1 <- plotFragmentSizes(ArchRProj = proj_DN1_2_KO)
#p2 <- plotTSSEnrichment(ArchRProj = proj_DN1_2_KO)
#plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_KOonly.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)


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

#p1 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p1 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-Clusters_LSI.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)

# QC score projected on UMAP
p1 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "PromoterRatio", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "NucleosomeRatio", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = proj_DN1_2, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-QC_LSI.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)


proj_DN1_2 <- addGeneScoreMatrix(proj_DN1_2,force=TRUE)
proj_DN1_2 <- addImputeWeights(proj_DN1_2,seed=1)

saveArchRProject(ArchRProj = proj_DN1_2, outputDirectory = "SaveDN1_2", load = FALSE)
#saveArchRProject(ArchRProj = proj_DN1_2_HAR, outputDirectory = "SaveDN1_2_HAR", load = FALSE)

proj_DN1_2 <- loadArchRProject(path = "SaveDN1_2")




############# RNA projection
seRNA <- readRDS("/nv/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/CITE/WTonly/RNAprocess_DN1WT.rds")

proj_DN1_3 <- addGeneIntegrationMatrix(
    ArchRProj = proj_DN1_2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    #groupList = groupList,
    groupRNA = "CT_5subset",
    nameCell = "predCell_CT5",
    nameGroup = "predGroup_CT5",
    nameScore = "predScore_CT5"
)


#proj_DN1_3 <- addGeneIntegrationMatrix(
#    ArchRProj = proj_DN1_3, 
#    useMatrix = "GeneScoreMatrix",
#    matrixName = "GeneIntegrationMatrix",
#    reducedDims = "IterativeLSI",
#    seRNA = seRNA,
#    addToArrow = TRUE,
#    force= TRUE,
#    #groupList = groupList,
#    groupRNA = "renameCluster",
#    nameCell = "predCell_renameCluster",
#    nameGroup = "predGroup_renameCluster",
#    nameScore = "predScore_renameCluster"
#)


#palCT <- paletteDiscrete(values = unique(seRNA$CT_4subset))
#pal_RNAcelltype[c("Fibroblast","Endothelial","Macrophage","Fibro","T_NK","SMC","Pericyte1","unknown1",
#                  "Pericyte2","B","Plasma","unknown2","Neuron","unknown3","Mast")] <- rainbow(15)
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

#color_list <- ggplotColours(n=12)

#pal_CT4 <- ggplotColours(n=4)
#names(pal_CT4) <- c("ETPs_proliferative","ETPs_quiescent","nonT_potential","T_ILCs")
#pal_renameCluster <- ggplotColours(n=18)
#names(pal_renameCluster) <- c(1,seq(10,18),seq(2,9))

pal_CT5 <- c(ggplotColours(n=5))
names(pal_CT5) <- c("ETPs_proliferative","ETPs_quiescent","nonT_potential","T_ILCs_GRhi","T_ILCs_IIhi")

proj_DN1_3$CT5 <- proj_DN1_3$predGroup_CT5
#proj_DN1_3$CT5[which(proj_DN1_3$predGroup_renameCluster %in% c(11,15,16))] <- "T_ILCs_GRhi"
#proj_DN1_3$CT5[which(proj_DN1_3$predGroup_renameCluster %in% c(12,13,14,17,18))] <- "T_ILCs_IIhi"
#
#proj_DN1_3$predGroup_renameCluster[which(proj_DN1_3$CT5=="T_ILCs")]
#

#p1 <- plotEmbedding(proj_DN1_3, colorBy = "cellColData", name = "predGroup_CT4", pal = pal_CT4, baseSize=2)
#p2 <- plotEmbedding(proj_DN1_3, colorBy = "cellColData", name = "predScore_CT4")
#p3 <- plotEmbedding(proj_DN1_3, colorBy = "cellColData", name = "predGroup_renameCluster", pal = pal_renameCluster,baseSize=2)
#p4 <- plotEmbedding(proj_DN1_3, colorBy = "cellColData", name = "predScore_renameCluster")
#p5 <- plotEmbedding(proj_DN1_3, colorBy = "cellColData", name = "CT5", pal = pal_CT5, baseSize=2)
#plotPDF(p1,p2,p3,p4,p5, name = "Plot-UMAP-RNAproject.pdf", ArchRProj = proj_DN1_3, addDOC = FALSE, width = 10, height = 10)

#p1 <- plotEmbedding(proj_DN1_3, colorBy = "cellColData", name = "predGroup_renameCluster", pal = pal_renameCluster,baseSize=2)
p1 <- plotEmbedding(proj_DN1_3, colorBy = "cellColData", name = "CT5", pal = pal_CT5, baseSize=2)
plotPDF(p1, name = "Plot-UMAP-RNAproject_CT5.pdf", ArchRProj = proj_DN1_3, addDOC = FALSE, width = 5, height = 5)

plot(proj_DN1_3@embeddings$UMAP$df)



saveArchRProject(ArchRProj = proj_DN1_3, outputDirectory = "SaveDN1_3", load = FALSE)
proj_DN1_3 <- loadArchRProject(path = "SaveDN1_3")





###########3 peak calling

# peak calling with macs2 v2.1.2
pathToMacs2 <- findMacs2()
pathToMacs2 <- "~/anaconda3/envs/condaPY3/bin/macs3"
proj_DN1_4 <- addGroupCoverages(ArchRProj = proj_DN1_3, groupBy = "CT5")
#proj_DN1_4tmp <- proj_DN1_4

proj_DN1_4 <- addReproduciblePeakSet(
    ArchRProj = proj_DN1_4, 
    groupBy = "CT5", 
    pathToMacs2 = "pathToMacs2",
    extsize=100,
    cutOff=0.01,
    shift=0,
    extendSummits=200,
    promoterRegion=c(2000,2000),
    genomeSize="mm10",
    reproducibility = "(n+1)/2",
    threads = getArchRThreads()
)
proj_DN1_4 <- addPeakMatrix(proj_DN1_4)
allpeaks <- getPeakSet(proj_DN1_4)
proj_DN1_4 <- addImputeWeights(proj_DN1_4,seed=1)

saveRDS(allpeaks, file="allpeaks.rds")
### change name:
proj_DN1_4$CT5 <- proj_DN1_3$CT5#[which(proj_DN1_4$CT5 == "T_ILCs_sub")] <- "T_ILCs_GRhi"
saveArchRProject(ArchRProj = proj_DN1_4, outputDirectory = "SaveDN1_4", load = FALSE)
proj_DN1_4 <- loadArchRProject(path = "SaveDN1_4")


#######33 marker gene detection

######33 marker gene

markersGS <- getMarkerFeatures(
    ArchRProj = proj_DN1_4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "CT5",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# cluster specific genes
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
#
saveRDS(markersGS, file="CT5_markerGS.rds")
#saveRDS(markerList, file="DN1_2_markerList.rds")
#markerList <- readRDS("CT5_markerGS.rds")

dir.create("markerGene_CT5")
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
    alltop5 <- rbind(alltop5, cbind(rep(Group, 10),thisdata[order(thisdata[,"Log2FC"],decreasing=T),][1:10,]))
    write.table(thisdata[order(thisdata[,"Log2FC"]),], file=paste0("markerGene_CT5/",Group,"_markerGene.txt"),row.names=F,col.names=T,sep="\t",quote=F)
}
write.table(alltop5[,c(1,6,8,9)], file="markerGene_CT5/markerGene_CT5_Top10summary.txt",row.names=F,col.names=T,sep="\t",quote=F)

markerGenes <- c(
"Aurkb","Birc5","Ccna2","Ccnb1","Cdk1","Lig1","Mki67","Ezh2","Dnmt1","Mcm3","Mcm6",
"Bcl11a","Lyl1","Cd24a","Kit","Flt3","Notch1","Hhex","Hes1","Spi1",
"Mpo","Gata2","Elane","Irf8","Cd19","Cd79a","Blnk",
"Thy1","Cd3e","Zap70","Prkcq","Tcf7","Ets1","Bcl11b","Tox","Gata3","Rorc","Id2","Il7r","Zbtb16")


heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot", use_raster = TRUE)
plotPDF(heatmapGS, name = "GeneScoresMarker_CT5_SIGHeatmap.pdf", width = 8, height = 6, ArchRProj = proj_DN1_4, addDOC = FALSE)



#heatmapPeaks <- plotMarkerHeatmap(
#  seMarker = markersPeaks, 
#  cutOff = "FDR <= 0.01 & Log2FC >= 1",
#  labelMarkers="",
#  transpose = TRUE
#)
#library(ComplexHeatmap)
##ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot", use_raster = TRUE)
#plotPDF(heatmapPeaks, name = "markerPeak_CT5_SIGheatmap.pdf", width = 8, height = 6, ArchRProj = proj_DN1_4, addDOC = FALSE)
#




######33 marker peak detection
dir.create("markerPeak_CT5/")
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_DN1_4, 
    useMatrix = "PeakMatrix", 
    groupBy = "CT5",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, file="markerPeak_CT5/markerPeak_CT5.rds")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

output_diffpeak <- function(this_data,cellname){
    out_data <- cbind(as.vector(this_data@seqnames),
                      this_data@ranges@start,
                      this_data@ranges@start + this_data@ranges@width,
                      paste0(cellname,"_",seq(length(this_data@seqnames))),
                      this_data$Log2FC,
                      this_data$FDR)
    write.table(out_data,file=paste0("markerPeak_CT5/",cellname,"_markerPeak.bed"),quote=F,row.names=F,col.names=F,sep="\t")
}
output_diffpeak(markerList$ETPs_proliferative,"ETPs_proliferative")
output_diffpeak(markerList$ETPs_quiescent,"ETPs_quiescent")
output_diffpeak(markerList$nonT_potential,"nonT_potential")
output_diffpeak(markerList$T_ILCs_GRhi,"T_ILCs_GRhi")
output_diffpeak(markerList$T_ILCs_IIhi,"T_ILCs_IIhi")


library(ComplexHeatmap)
library(grDevices)

# Define a custom color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))
# Set the color palette in the current R session
options("ArchR.Colors" = my_palette(100))



heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE,
)
#ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot", use_raster = TRUE)
plotPDF(heatmapPeaks, name = "markerPeak_CT5_SIGheatmap.pdf", width = 8, height = 6, ArchRProj = proj_DN1_4, addDOC = FALSE)

#for(celltype in c("ETPs_proliferative","ETPs_quiescent","nonT_potential","T_ILCs_GRhi","T_ILCs_IIhi")){
#    pma <- plotMarkers(seMarker = markersPeaks, name =celltype, cutOff = "FDR <= 0.01 & Log2FC >= 1", plotAs = "MA")
#    pv <- plotMarkers(seMarker = markersPeaks, name =celltype, cutOff = "FDR <= 0.01 & Log2FC >= 1", plotAs = "Volcano")
#    plotPDF(pma, pv, name = paste0("markerPeak_CT5_scatterMaVolcano_",celltype,".pdf"), width = 5, height = 5, ArchRProj = proj_DN1_4, addDOC = FALSE)
#}
proj_DN1_4$poolCT <- proj_DN1_4$CT5 
proj_DN1_4$poolCT[which(proj_DN1_4$CT5  %in% c("ETPs_proliferative","ETPs_quiescent"))] <- "ETPs"
proj_DN1_4$poolCT[which(proj_DN1_4$CT5  %in% c("T_ILCs_GRhi","T_ILCs_IIhi"))] <- "T_ILCs"

markerPeak_ETP_nonT <- getMarkerFeatures(
  ArchRProj = proj_DN1_4, 
  useMatrix = "PeakMatrix",
  groupBy = "poolCT",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ETPs",
  bgdGroups = "nonT_potential"
)
markerPeak_TILC_nonT <- getMarkerFeatures(
  ArchRProj = proj_DN1_4, 
  useMatrix = "PeakMatrix",
  groupBy = "poolCT",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_ILCs",
  bgdGroups = "nonT_potential"
)
markerPeak_ETP_TILC <- getMarkerFeatures(
  ArchRProj = proj_DN1_4, 
  useMatrix = "PeakMatrix",
  groupBy = "poolCT",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ETPs",
  bgdGroups = "T_ILCs"
)
markerPeak_TILCsubset <- getMarkerFeatures(
  ArchRProj = proj_DN1_4, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_ILCs_GRhi",
  bgdGroups = "T_ILCs_IIhi"
)

p1 <- markerPlot(seMarker = markerPeak_ETP_nonT, name = "ETPs", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
p2 <- markerPlot(seMarker = markerPeak_TILC_nonT, name = "T_ILCs", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
p3 <- markerPlot(seMarker = markerPeak_ETP_TILC, name = "ETPs", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
p4 <- markerPlot(seMarker = markerPeak_TILCsubset, name = "T_ILCs_GRhi", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(p1,p2,p3,p4, name = "pairwise_markerPeakVolcano", width = 5, height = 5, ArchRProj = proj_DN1_4, addDOC = FALSE)

markerList1 <- getMarkers(markerPeak_ETP_nonT, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)



ETP_nonT_enrichMotifs1 <- peakAnnoEnrichment(
    seMarker = markerPeak_ETP_nonT,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )
ETP_nonT_enrichMotifs2 <- peakAnnoEnrichment(
    seMarker = markerPeak_ETP_nonT,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -1"
  )

ETP_nonT_padj_mat1 <- ETP_nonT_enrichMotifs1@assays@data$mlog10Padj#[usemotif,]
ETP_nonT_padj_mat2 <- ETP_nonT_enrichMotifs2@assays@data$mlog10Padj#[usemotif,]

TILC_nonT_enrichMotifs1 <- peakAnnoEnrichment(
    seMarker = markerPeak_TILC_nonT,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )
TILC_nonT_enrichMotifs2 <- peakAnnoEnrichment(
    seMarker = markerPeak_TILC_nonT,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -1"
  )

TILC_nonT_padj_mat1 <- TILC_nonT_enrichMotifs1@assays@data$mlog10Padj#[usemotif,]
TILC_nonT_padj_mat2 <- TILC_nonT_enrichMotifs2@assays@data$mlog10Padj#[usemotif,]

ETP_TILC_enrichMotifs1 <- peakAnnoEnrichment(
    seMarker = markerPeak_ETP_TILC,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )
ETP_TILC_enrichMotifs2 <- peakAnnoEnrichment(
    seMarker = markerPeak_ETP_TILC,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -1"
  )

ETP_TILC_padj_mat1 <- ETP_TILC_enrichMotifs1@assays@data$mlog10Padj#[usemotif,]
ETP_TILC_padj_mat2 <- ETP_TILC_enrichMotifs2@assays@data$mlog10Padj#[usemotif,]

TILCsubset_enrichMotifs1 <- peakAnnoEnrichment(
    seMarker = markerPeak_TILCsubset,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )
TILCsubset_enrichMotifs2 <- peakAnnoEnrichment(
    seMarker = markerPeak_TILCsubset,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -1"
  )

TILCsubset_padj_mat1 <- TILCsubset_enrichMotifs1@assays@data$mlog10Padj#[usemotif,]
TILCsubset_padj_mat2 <- TILCsubset_enrichMotifs2@assays@data$mlog10Padj#[usemotif,]

pairmotif<-cbind(ETP_nonT_padj_mat1,
ETP_nonT_padj_mat2,
TILC_nonT_padj_mat1,
TILC_nonT_padj_mat2,
ETP_TILC_padj_mat1 ,
ETP_TILC_padj_mat2 ,
TILCsubset_padj_mat1, 
TILCsubset_padj_mat2 
)
filtermotif <- pairmotif[which(apply(pairmotif,1,max)>=10),]
write.table(filtermotif,file="pairmotif_filter.txt",row.names=T,col.names=T,sep="\t",quote=F)



###### motif dot plot
filtermotif <- read.table("pairmotif_filter.txt",row.names=1,header=T)











#######33 motif anno

proj_DN1_5 <- addMotifAnnotations(ArchRProj = proj_DN1_4, motifSet = "homer", name = "Motif",force = TRUE)
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_DN1_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "enrichMotifs_CT5_Heatmap.pdf", width = 8, height = 6, ArchRProj = proj_DN1_5, addDOC = FALSE)
usemotif <- unlist(lapply(colnames(heatmapEM@matrix), function(x) return(unlist(strsplit(x," ")[[1]][1]))))#[,1]
raw_padj_mat <- enrichMotifs@assays@data$mlog10Padj[usemotif,rownames(heatmapEM@matrix)]
write.table(enrichMotifs@assays@data$mlog10Padj,file="markerPeak_CT5/enrichMotifs_CT5.txt", row.names=T,col.names=T,sep="\t",quote=F)

motifMat <- enrichMotifs@assays@data$mlog10Padj
motifMat_filter <- motifMat[which(apply(motifMat,1,max)> 100),]

motifMat_filter_tmp <- motifMat_filter
motifMat_filter_tmp[motifMat_filter_tmp>100] <- 100
custom_palette <- colorRampPalette(c("white", "blue", "black"))(100)
heatmap_obj <- Heatmap(
  t(motifMat_filter_tmp),
  name = "Motif -log10(adjust.pvalue)",
  col = custom_palette,#colorRamp2(c(min(motifMat_filter), max(motifMat_filter)), colors = c("white", "red")),
  show_row_names = TRUE,   # Show row names
  show_column_names = TRUE,  # Show column names
  row_names_gp = gpar(fontsize = 10),  # Customize row name appearance
  column_names_gp = gpar(fontsize = 10)  # Customize column name appearance
)
draw(heatmap_obj, heatmap_legend_side = "bot", annotation_legend_side = "bot", use_raster = TRUE)
plotPDF(heatmap_obj, name = "markerPeakMotifadjp_CT5_SIGheatmap.pdf", width = 10, height = 8, ArchRProj = proj_DN1_4, addDOC = FALSE)


saveRDS(enrichMotifs, file="DN1_scATAC_CT5_enrichMotifs.rds")
#saveRDS(markerList, file="DN1_2_markerList.rds")
#enrichMotifs <- readRDS("DN1_scATAC_CT5_enrichMotifs.rds")
saveArchRProject(ArchRProj = proj_DN1_5, outputDirectory = "SaveDN1_5", load = FALSE)
proj_DN1_5 <- loadArchRProject(path = "SaveDN1_5")


#######33add coA

proj_DN1_6 <- addCoAccessibility(
    ArchRProj = proj_DN1_5,
    reducedDims = "IterativeLSI"
)
cA <- getCoAccessibility(
    ArchRProj = proj_DN1_6,
    corCutOff = 0.5,
    resolution = 10000,
    returnLoops = TRUE
)
proj_DN1_6 <- addPeak2GeneLinks(
    ArchRProj = proj_DN1_6,
    reducedDims = "IterativeLSI"
)
p2g <- getPeak2GeneLinks(
    ArchRProj = proj_DN1_6,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)
saveArchRProject(ArchRProj = proj_DN1_6, outputDirectory = "SaveDN1_6", load = FALSE)
proj_DN1_6 <- loadArchRProject(path = "SaveDN1_6")

write.table(p2g$Peak2GeneLinks,file="peak2geneLinks_final.txt",quote=F,sep="\t",row.names=F)
write.table(cA$CoAccessibility,file="CoAccessibility_final.txt",quote=F,sep="\t",row.names=F)



markerGenes  <- c(
"Cd19","Cd79a","Blnk","Ebf1",
"Irf8","Tcf4","Itgax",
"Gata2","Mpo","Cebpa","Elane",
"Cd3e","Cd3d","Cd3g","Thy1",
"Prkcq","Lck","Zap70","Erg","Bcl11b","Ets1",
"Id2","Tox","Gata3","Tcf7","Il7r","Rorc","Zbtb16",
"Spi1","Lyl1","Hhex","Bcl11a",
"Flt3","Hhex","Kit","Kit","Notch1","Hes1",
"Ezh2","Dnmt1","Mki67","Lig1","Cdk1","Ccnb1","Ccnb2","Ccna2","Mcm3","Mcm6","Mcm5","Aurkb","Birc5","Bub1b"
 )


p <- plotBrowserTrack(
    ArchRProj = proj_DN1_6, 
    groupBy = "CT5", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = cA
)

plotPDF(plotList = p, 
    name = "multiTracks_CT5.pdf", 
    ArchRProj = proj_DN1_6, 
    addDOC = FALSE, width = 5, height = 5)






#### pairwise diffpeak, combine sub type
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
pal_CT5 <- c(ggplotColours(n=5))
names(pal_CT5) <- c("ETPs_proliferative","ETPs_quiescent","nonT_potential","T_ILCs_sub","T_ILCs")

p1 <- plotEmbedding(proj_DN1_6, colorBy = "cellColData", name = "CT5", pal = pal_CT5, baseSize=2)
plotPDF(p1, name = "Plot-UMAP-RNAproject_CT5new.pdf", ArchRProj = proj_DN1_6, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(proj_DN1_6, colorBy = "cellColData", name = "CT5", baseSize=2)
plotPDF(p1, name = "Plot-UMAP-RNAproject_CT5new.pdf", ArchRProj = proj_DN1_6, addDOC = FALSE, width = 5, height = 5)


proj_DN1_6$CT5cb3 <- rep("nonT_potential",length(proj_DN1_6$CT5))
proj_DN1_6$CT5cb3[which(proj_DN1_6$CT5 %in% c("ETPs_proliferative","ETPs_quiescent"))] <- "ETPs"
proj_DN1_6$CT5cb3[which(proj_DN1_6$CT5 %in% c("T_ILCs","T_ILCs_sub"))] <- "T_ILCs"

MKpeakCT5cb3_nonT_ETPs <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "nonT_potential",
  bgdGroups = "ETPs"
)

MKpeakCT5cb3_nonT_TILCs <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "nonT_potential",
  bgdGroups = "T_ILCs"
)

MKpeakCT5cb3_ETPs_nonT <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ETPs",
  bgdGroups = "nonT_potential"
)

MKpeakCT5cb3_ETPs_TILCs <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ETPs",
  bgdGroups = "T_ILCs"
)

MKpeakCT5cb3_TILCs_nonT <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_ILCs",
  bgdGroups = "nonT_potential"
)

MKpeakCT5cb3_TILCs_ETPs <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_ILCs",
  bgdGroups = "ETPs"
)

MKpeakCT5_TILCs_E <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_ILCs",
  bgdGroups = "ETPs"
)

markerPeak_TILCsubset_TILChigh <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_ILCs",
  bgdGroups = "T_ILCs_sub"
)

markerPeak_TILCsubset_TILCsubhigh <- getMarkerFeatures(
  ArchRProj = proj_DN1_6, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = "T_ILCs",
  useGroups = "T_ILCs_sub"
)

MKpeakCT5cb3_nonT_ETPs_markerList <- getMarkers(MKpeakCT5cb3_nonT_ETPs, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_nonT_TILCs_markerList <- getMarkers(MKpeakCT5cb3_nonT_TILCs, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_ETPs_nonT_markerList <- getMarkers(MKpeakCT5cb3_ETPs_nonT, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_ETPs_TILCs_markerList <- getMarkers(MKpeakCT5cb3_ETPs_TILCs, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_TILCs_nonT_markerList <- getMarkers(MKpeakCT5cb3_TILCs_nonT, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_TILCs_ETPs_markerList <- getMarkers(MKpeakCT5cb3_TILCs_ETPs, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerPeak_TILCsubset_TILChigh_markerList <- getMarkers(markerPeak_TILCsubset_TILChigh, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerPeak_TILCsubset_TILCsubhigh_markerList <- getMarkers(markerPeak_TILCsubset_TILCsubhigh, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

combine2mklist <- function(list1, list2,name1,name2){
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
singlemklist <- function(list1){
    a <- mcols(list1[[1]])
    rownames(a) <- paste0(as.character(seqnames(list1[[1]])),":",
          start(list1[[1]]),"-",
          end(list1[[1]]))
    return(a)
}

MKpeakCT5cb3_nonT<-combine2mklist(MKpeakCT5cb3_nonT_ETPs_markerList,MKpeakCT5cb3_nonT_TILCs_markerList,"NonTgETP","NonTgTILC")
MKpeakCT5cb3_ETPs<-combine2mklist(MKpeakCT5cb3_ETPs_nonT_markerList,MKpeakCT5cb3_ETPs_TILCs_markerList,"ETPgNonT","ETPgTILC")
MKpeakCT5cb3_TILCs<-combine2mklist(MKpeakCT5cb3_TILCs_nonT_markerList,MKpeakCT5cb3_TILCs_ETPs_markerList,"TILCgNonT","TILCgETP")
MKpeakCT5_TILChigh <- singlemklist(markerPeak_TILCsubset_TILChigh_markerList)
MKpeakCT5_TILCsubhigh <- singlemklist(markerPeak_TILCsubset_TILCsubhigh_markerList)



transformBED <- function(indata){
    regions <- rownames(indata)
    split_regions <- strsplit(regions, ":|-")
    # Extract chromosomes, start, and end positions
    chromosomes <- sapply(split_regions, `[`, 1)
    starts <- as.integer(sapply(split_regions, `[`, 2))
    ends <- as.integer(sapply(split_regions, `[`, 3))
    newdata <-  cbind(chromosomes,starts,ends,indata)
    return(newdata)
}

write.table(as.matrix(transformBED(MKpeakCT5cb3_nonT)),file="newMKpeak_CT5cb3/MKpeakCT5cb3_nonT.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_ETPs)),file="newMKpeak_CT5cb3/MKpeakCT5cb3_ETPs.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_TILCs)),file="newMKpeak_CT5cb3/MKpeakCT5cb3_TILCs.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5_TILChigh)),file="newMKpeak_CT5cb3/MKpeakCT5_TILChigh.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5_TILCsubhigh)),file="newMKpeak_CT5cb3/MKpeakCT5_TILCsubhigh.txt",row.names=F,col.names=T,sep="\t",quote=F)


write.table(as.matrix(transformBED(MKpeakCT5cb3_nonT))[,1:3],file="newMKpeak_CT5cb3/MKpeakCT5cb3_nonT.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_ETPs))[,1:3],file="newMKpeak_CT5cb3/MKpeakCT5cb3_ETPs.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_TILCs))[,1:3],file="newMKpeak_CT5cb3/MKpeakCT5cb3_TILCs.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5_TILChigh))[,1:3],file="newMKpeak_CT5cb3/MKpeakCT5_TILChigh.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5_TILCsubhigh))[,1:3],file="newMKpeak_CT5cb3/MKpeakCT5_TILCsubhigh.bed",row.names=F,col.names=F,sep="\t",quote=F)










######33 print reads
Cdata <- getCellColData(ArchRProj = proj_DN1_6, select = NULL, drop = FALSE)

for(c in unique(Cdata$CT5)){
    print(c)
    cellnames<-rownames(Cdata[which(Cdata$CT5==c),])
    frag1<-getFragmentsFromArrow(
        ArrowFile = "WT.arrow",
        cellNames = cellnames,
        verbose = TRUE,
    )
    #cellnames<-rownames(Cdata[which(Cdata$superClusters==c & Cdata$Sample=="E9.5_r2"),])
    #frag2<-getFragmentsFromArrow(
    #    ArrowFile = "E9.5_r2.arrow",
    #    cellNames = cellnames,
    #    verbose = TRUE,
    #)
    fragOut<-c(frag1)#,frag2)

    df <- data.frame(seqnames=seqnames(fragOut),
    starts=start(fragOut)-1,
    ends=end(fragOut))

    write.table(df, file=paste0("CT5_reads/",c,"_read.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}


for(c in c("T_ILCs","T_ILCs_sub")){
    print(c)
    cellnames<-rownames(Cdata[which(Cdata$CT5==c),])
    frag1<-getFragmentsFromArrow(
        ArrowFile = "WT.arrow",
        cellNames = cellnames,
        verbose = TRUE,
    )
    #cellnames<-rownames(Cdata[which(Cdata$superClusters==c & Cdata$Sample=="E9.5_r2"),])
    #frag2<-getFragmentsFromArrow(
    #    ArrowFile = "E9.5_r2.arrow",
    #    cellNames = cellnames,
    #    verbose = TRUE,
    #)
    fragOut<-c(frag1)#,frag2)

    df <- data.frame(seqnames=seqnames(fragOut),
    starts=start(fragOut)-1,
    ends=end(fragOut))

    write.table(df, file=paste0("CT5_reads/",c,"_read.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}


####################### add cb reads macs2 peaks as peakset
peakdata <- read.table("/sfs/ceph/standard/vol190/zanglab/sh8tv/Project/TCF1_Tcell/Result/DN1/scATAC_ArchR/WTonly/CT5_reads/WTonly_mergePeak_CT.bed",row.names=4)
peakGR <- GRanges(seqnames=peakdata[,1],ranges=IRanges(peakdata[,2],peakdata[,3]))
proj_DN1_7 <- addPeakSet(ArchRProj=proj_DN1_6, peakSet=peakGR,force=TRUE)
proj_DN1_7 <- addPeakMatrix(proj_DN1_7)



proj_DN1_7$CT5cb3 <- rep("nonT_potential",length(proj_DN1_7$CT5))
proj_DN1_7$CT5cb3[which(proj_DN1_7$CT5 %in% c("ETPs_proliferative","ETPs_quiescent"))] <- "ETPs"
proj_DN1_7$CT5cb3[which(proj_DN1_7$CT5 %in% c("T_ILCs","T_ILCs_sub"))] <- "T_ILCs"




proj_DN1_7 <- addPeak2GeneLinks(
    ArchRProj = proj_DN1_7,
    reducedDims = "IterativeLSI"
)
p2g <- getPeak2GeneLinks(
    ArchRProj = proj_DN1_7,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)

write.table(p2g$Peak2GeneLinks,file="peak2geneLinks_newMergePeak.txt",quote=F,sep="\t",row.names=F)


p2geneDF <- metadata(proj_DN1_7@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]

usep2g <- p2geneDF[which(p2geneDF$Correlation > 0.45 & p2geneDF$FDR < 1e-4),]
outdat <- cbind(usep2g$peakName, usep2g$geneName, usep2g$Correlation, usep2g$FDR) 
colnames(outdat) <- c("peakname","genename","cor","fdr")
write.table(outdat,file="peak2geneLinks_newMergePeak_precise_pre.txt",row.names=F,col.name=F,sep="\t",quote=F) 

proj_DN1_8 <- addCoAccessibility(
    ArchRProj = proj_DN1_7,
    useMatrix = "GeneIntegrationMatrix", 
    reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
    ArchRProj = proj_DN1_6,
    corCutOff = 0.5,
    resolution = 10000,
    returnLoops = TRUE
)
proj_DN1_8 <- addPeak2GeneLinks(
    ArchRProj = proj_DN1_7,
    useMatrix = "GeneIntegrationMatrix", 
    reducedDims = "IterativeLSI"

)
p2g <- getPeak2GeneLinks(
    ArchRProj = proj_DN1_7,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)

write.table(p2g$Peak2GeneLinks,file="peak2geneLinks_newMergePeak.txt",quote=F,sep="\t",row.names=F)















saveArchRProject(ArchRProj = proj_DN1_7, outputDirectory = "SaveDN1_7", load = FALSE)
proj_DN1_7 <- loadArchRProject(path = "SaveDN1_7")


cellnames_WT_ETP  <- proj_DN1_7$cellNames[ which(proj_DN1_7$CT5cb3 == "ETPs")]
cellnames_WT_nonT <- proj_DN1_7$cellNames[ which(proj_DN1_7$CT5cb3 == "nonT_potential")]
cellnames_WT_TILC <- proj_DN1_7$cellNames[ which(proj_DN1_7$CT5cb3 == "T_ILCs")]
 
proj_WT_ETP  <- proj_DN1_7[cellnames_WT_ETP,]
proj_WT_nonT <- proj_DN1_7[cellnames_WT_nonT,]
proj_WT_TILC <- proj_DN1_7[cellnames_WT_TILC,]


proj_WT_ETP <- addPeak2GeneLinks(
    ArchRProj = proj_WT_ETP,
    reducedDims = "IterativeLSI"
)
p2g_ETP <- getPeak2GeneLinks(
    ArchRProj = proj_WT_ETP,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)


proj_WT_nonT <- addPeak2GeneLinks(
    ArchRProj = proj_WT_nonT,
    reducedDims = "IterativeLSI"
)
p2g_nonT <- getPeak2GeneLinks(
    ArchRProj = proj_WT_nonT,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)

proj_WT_TILC <- addPeak2GeneLinks(
    ArchRProj = proj_WT_TILC,
    reducedDims = "IterativeLSI"
)
p2g_TILC <- getPeak2GeneLinks(
    ArchRProj = proj_WT_TILC,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)

write.table(p2g_ETP$Peak2GeneLinks,file="peak2geneLinks_newMergePeak_ETP.txt",quote=F,sep="\t",row.names=F)
write.table(p2g_nonT$Peak2GeneLinks,file="peak2geneLinks_newMergePeak_nonT.txt",quote=F,sep="\t",row.names=F)
write.table(p2g_TILC$Peak2GeneLinks,file="peak2geneLinks_newMergePeak_TILC.txt",quote=F,sep="\t",row.names=F)



keygeneNEW<-c("Runx3","Irf7","Il2ra","Ly6c1","Il5","Il13","Cd34","Klrd1","Klrb1b","Itga2","Klrb1c","Il1rl1","Ly6c1","Cd63")
c("Notch1","Kit","Mki67","Lig1","Cepba","Cd19","Ets1","Tcf7","Ly6c1","Cd63","Klrd1","Klrb1b","Itga2","Klrb1c","Il1rl1","Il2ra","Il5","Il13","Runx3","Irf7","Flt3","Cd34","Hes1","Hhex","Spi1","Lyl1","Hoxa9","Bcl11a","")

markerGenes  <- c(
"Cd19","Cd79a","Blnk","Ebf1")#,
"Irf8","Tcf4","Itgax",
"Gata2","Mpo","Cebpa","Elane",
"Cd3e","Cd3d","Cd3g","Thy1",
"Prkcq","Lck","Zap70","Erg","Bcl11b","Ets1",
"Id2","Tox","Gata3","Tcf7","Il7r","Rorc","Zbtb16",
"Spi1","Lyl1","Hhex","Bcl11a",
"Flt3","Hhex","Kit","Kit","Notch1","Hes1",
"Ezh2","Dnmt1","Mki67","Lig1","Cdk1","Ccnb1","Ccnb2","Ccna2","Mcm3","Mcm6","Mcm5","Aurkb","Birc5","Bub1b"
 )


p <- plotBrowserTrack(
    ArchRProj = proj_DN1_7, 
    groupBy = "CT5cb3", 
    geneSymbol = keygeneNEW, 
    features =  MKpeakCT5cb3_ETPs_nonT_markerList,#getMarkers(markersPeaks_CT5cb3, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"],
    upstream = 50000,
    downstream = 50000,
    loops = p2g
)

plotPDF(plotList = p, 
    name = "multiTracks_CT5cb3.pdf", 
    ArchRProj = proj_DN1_7, 
    addDOC = FALSE, width = 5, height = 5)






#markersPeaks_CT5cb3 <- getMarkerFeatures(
#    ArchRProj = proj_DN1_7, 
#    useMatrix = "PeakMatrix", 
#    groupBy = "CT5cb3",
#  bias = c("TSSEnrichment", "log10(nFrags)"),
#  testMethod = "wilcoxon"
#)
MKpeakCT5cb3_nonT_ETPs <- getMarkerFeatures(
  ArchRProj = proj_DN1_7, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "nonT_potential",
  bgdGroups = "ETPs"
)

MKpeakCT5cb3_nonT_TILCs <- getMarkerFeatures(
  ArchRProj = proj_DN1_7, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "nonT_potential",
  bgdGroups = "T_ILCs"
)

MKpeakCT5cb3_ETPs_nonT <- getMarkerFeatures(
  ArchRProj = proj_DN1_7, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ETPs",
  bgdGroups = "nonT_potential"
)

MKpeakCT5cb3_ETPs_TILCs <- getMarkerFeatures(
  ArchRProj = proj_DN1_7, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ETPs",
  bgdGroups = "T_ILCs"
)

MKpeakCT5cb3_TILCs_nonT <- getMarkerFeatures(
  ArchRProj = proj_DN1_7, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_ILCs",
  bgdGroups = "nonT_potential"
)

MKpeakCT5cb3_TILCs_ETPs <- getMarkerFeatures(
  ArchRProj = proj_DN1_7, 
  useMatrix = "PeakMatrix",
  groupBy = "CT5cb3",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_ILCs",
  bgdGroups = "ETPs"
)

#MKpeakCT5_TILCs_E <- getMarkerFeatures(
#  ArchRProj = proj_DN1_7, 
#  useMatrix = "PeakMatrix",
#  groupBy = "CT5cb3",
#  testMethod = "wilcoxon",
#  bias = c("TSSEnrichment", "log10(nFrags)"),
#  useGroups = "T_ILCs",
#  bgdGroups = "ETPs"
#)
#
#markerPeak_TILCsubset_TILChigh <- getMarkerFeatures(
#  ArchRProj = proj_DN1_7, 
#  useMatrix = "PeakMatrix",
#  groupBy = "CT5",
#  testMethod = "wilcoxon",
#  bias = c("TSSEnrichment", "log10(nFrags)"),
#  useGroups = "T_ILCs",
#  bgdGroups = "T_ILCs_sub"
#)
#
#markerPeak_TILCsubset_TILCsubhigh <- getMarkerFeatures(
#  ArchRProj = proj_DN1_7, 
#  useMatrix = "PeakMatrix",
#  groupBy = "CT5",
#  testMethod = "wilcoxon",
#  bias = c("TSSEnrichment", "log10(nFrags)"),
#  bgdGroups = "T_ILCs",
#  useGroups = "T_ILCs_sub"
#)
#
#MKpeakCT5cb3_all_markerList <- getMarkers(markersPeaks_CT5cb3, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)


MKpeakCT5cb3_nonT_ETPs_markerList <- getMarkers(MKpeakCT5cb3_nonT_ETPs, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625", returnGR = TRUE)
MKpeakCT5cb3_nonT_TILCs_markerList <- getMarkers(MKpeakCT5cb3_nonT_TILCs, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625", returnGR = TRUE)
MKpeakCT5cb3_ETPs_nonT_markerList <- getMarkers(MKpeakCT5cb3_ETPs_nonT, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625", returnGR = TRUE)
MKpeakCT5cb3_ETPs_TILCs_markerList <- getMarkers(MKpeakCT5cb3_ETPs_TILCs, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625", returnGR = TRUE)
MKpeakCT5cb3_TILCs_nonT_markerList <- getMarkers(MKpeakCT5cb3_TILCs_nonT, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625", returnGR = TRUE)
MKpeakCT5cb3_TILCs_ETPs_markerList <- getMarkers(MKpeakCT5cb3_TILCs_ETPs, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625", returnGR = TRUE)
#markerPeak_TILCsubset_TILChigh_markerList <- getMarkers(markerPeak_TILCsubset_TILChigh, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
#markerPeak_TILCsubset_TILCsubhigh_markerList <- getMarkers(markerPeak_TILCsubset_TILCsubhigh, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)




combine2mklist <- function(list1, list2,name1,name2){
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
singlemklist <- function(list1){
    a <- mcols(list1[[1]])
    rownames(a) <- paste0(as.character(seqnames(list1[[1]])),":",
          start(list1[[1]]),"-",
          end(list1[[1]]))
    return(a)
}

MKpeakCT5cb3_nonT<-combine2mklist(MKpeakCT5cb3_nonT_ETPs_markerList,MKpeakCT5cb3_nonT_TILCs_markerList,"NonTgETP","NonTgTILC")
MKpeakCT5cb3_ETPs<-combine2mklist(MKpeakCT5cb3_ETPs_nonT_markerList,MKpeakCT5cb3_ETPs_TILCs_markerList,"ETPgNonT","ETPgTILC")
MKpeakCT5cb3_TILCs<-combine2mklist(MKpeakCT5cb3_TILCs_nonT_markerList,MKpeakCT5cb3_TILCs_ETPs_markerList,"TILCgNonT","TILCgETP")
#MKpeakCT5_TILChigh <- singlemklist(markerPeak_TILCsubset_TILChigh_markerList)
#MKpeakCT5_TILCsubhigh <- singlemklist(markerPeak_TILCsubset_TILCsubhigh_markerList)


#intersect(MKpeakCT5cb3_nonT_ETPs_markerList[[1]], MKpeakCT5cb3_nonT_TILCs_markerList[[1]] )

transformBED <- function(indata){
    regions <- rownames(indata)
    split_regions <- strsplit(regions, ":|-")
    # Extract chromosomes, start, and end positions
    chromosomes <- sapply(split_regions, `[`, 1)
    starts <- as.integer(sapply(split_regions, `[`, 2))
    ends <- as.integer(sapply(split_regions, `[`, 3))
    newdata <-  cbind(chromosomes,starts,ends,indata)
    return(newdata)
}



write.table(as.matrix(transformBED(MKpeakCT5cb3_nonT)),file="MKpeak_CT5readsMacs2/MKpeakCT5cb3_nonT.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_ETPs)),file="MKpeak_CT5readsMacs2/MKpeakCT5cb3_ETPs.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_TILCs)),file="MKpeak_CT5readsMacs2/MKpeakCT5cb3_TILCs.txt",row.names=F,col.names=T,sep="\t",quote=F)
#write.table(as.matrix(transformBED(MKpeakCT5_TILChigh)),file="MKpeak_CT5readsMacs2/MKpeakCT5_TILChigh.txt",row.names=F,col.names=T,sep="\t",quote=F)
#write.table(as.matrix(transformBED(MKpeakCT5_TILCsubhigh)),file="MKpeak_CT5readsMacs2/MKpeakCT5_TILCsubhigh.txt",row.names=F,col.names=T,sep="\t",quote=F)

write.table(as.matrix(transformBED(MKpeakCT5cb3_nonT))[,1:3],file="MKpeak_CT5readsMacs2/MKpeakCT5cb3_nonT.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_ETPs))[,1:3],file="MKpeak_CT5readsMacs2/MKpeakCT5cb3_ETPs.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_TILCs))[,1:3],file="MKpeak_CT5readsMacs2/MKpeakCT5cb3_TILCs.bed",row.names=F,col.names=F,sep="\t",quote=F)
#write.table(as.matrix(transformBED(MKpeakCT5_TILChigh))[,1:3],file="MKpeak_CT5readsMacs2/MKpeakCT5_TILChigh.bed",row.names=F,col.names=F,sep="\t",quote=F)
#write.table(as.matrix(transformBED(MKpeakCT5_TILCsubhigh))[,1:3],file="MKpeak_CT5readsMacs2/MKpeakCT5_TILCsubhigh.bed",row.names=F,col.names=F,sep="\t",quote=F)




### FC2
MKpeakCT5cb3_nonT_ETPs_markerList <- getMarkers(MKpeakCT5cb3_nonT_ETPs, cutOff = "FDR <= 0.001 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_nonT_TILCs_markerList <- getMarkers(MKpeakCT5cb3_nonT_TILCs, cutOff = "FDR <= 0.001 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_ETPs_nonT_markerList <- getMarkers(MKpeakCT5cb3_ETPs_nonT, cutOff = "FDR <= 0.001 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_ETPs_TILCs_markerList <- getMarkers(MKpeakCT5cb3_ETPs_TILCs, cutOff = "FDR <= 0.001 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_TILCs_nonT_markerList <- getMarkers(MKpeakCT5cb3_TILCs_nonT, cutOff = "FDR <= 0.001 & Log2FC >= 1", returnGR = TRUE)
MKpeakCT5cb3_TILCs_ETPs_markerList <- getMarkers(MKpeakCT5cb3_TILCs_ETPs, cutOff = "FDR <= 0.001 & Log2FC >= 1", returnGR = TRUE)

MKpeakCT5cb3_nonT<-combine2mklist(MKpeakCT5cb3_nonT_ETPs_markerList,MKpeakCT5cb3_nonT_TILCs_markerList,"NonTgETP","NonTgTILC")
MKpeakCT5cb3_ETPs<-combine2mklist(MKpeakCT5cb3_ETPs_nonT_markerList,MKpeakCT5cb3_ETPs_TILCs_markerList,"ETPgNonT","ETPgTILC")
MKpeakCT5cb3_TILCs<-combine2mklist(MKpeakCT5cb3_TILCs_nonT_markerList,MKpeakCT5cb3_TILCs_ETPs_markerList,"TILCgNonT","TILCgETP")



write.table(as.matrix(transformBED(MKpeakCT5cb3_nonT)),file="MKpeak_CT5readsMacs2_FC2/MKpeakCT5cb3_nonT.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_ETPs)),file="MKpeak_CT5readsMacs2_FC2/MKpeakCT5cb3_ETPs.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_TILCs)),file="MKpeak_CT5readsMacs2_FC2/MKpeakCT5cb3_TILCs.txt",row.names=F,col.names=T,sep="\t",quote=F)

write.table(as.matrix(transformBED(MKpeakCT5cb3_nonT))[,1:3],file="MKpeak_CT5readsMacs2_FC2/MKpeakCT5cb3_nonT.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_ETPs))[,1:3],file="MKpeak_CT5readsMacs2_FC2/MKpeakCT5cb3_ETPs.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(as.matrix(transformBED(MKpeakCT5cb3_TILCs))[,1:3],file="MKpeak_CT5readsMacs2_FC2/MKpeakCT5cb3_TILCs.bed",row.names=F,col.names=F,sep="\t",quote=F)







##### 5 CT marker peak

output_diffpeak <- function(this_data,outname){
    cellname <- strsplit(strsplit(outname,"/")[[1]][3],"_markerPeak.bed")[[1]][1]
    out_data <- cbind(as.vector(this_data@seqnames),
                      this_data@ranges@start,
                      this_data@ranges@start + this_data@ranges@width,
                      paste0(cellname,"_",seq(length(this_data@seqnames))),
                      this_data$Log2FC,
                      this_data$FDR)
    write.table(out_data,file=outname,quote=F,row.names=F,col.names=F,sep="\t")
}

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_DN1_7, 
    useMatrix = "PeakMatrix", 
    groupBy = "CT5",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, file="markerPeak_CT5/markerPeak_CT5.rds")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.001 & Log2FC >= 1", returnGR = TRUE)
output_diffpeak(markerList$ETPs_proliferative,"markerPeak_CT5/FC2/ETPs_proliferative_markerPeak.bed")
output_diffpeak(markerList$ETPs_quiescent,"markerPeak_CT5/FC2/ETPs_quiescent_markerPeak.bed")
output_diffpeak(markerList$nonT_potential,"markerPeak_CT5/FC2/nonT_potential_markerPeak.bed")
output_diffpeak(markerList$T_ILCs,"markerPeak_CT5/FC2/T_ILCs_markerPeak.bed")
output_diffpeak(markerList$T_ILCs_sub,"markerPeak_CT5/FC2/T_ILCsSub_markerPeak.bed")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5849625", returnGR = TRUE)
output_diffpeak(markerList$ETPs_proliferative,"markerPeak_CT5/FC1.5/ETPs_proliferative_markerPeak.bed")
output_diffpeak(markerList$ETPs_quiescent,"markerPeak_CT5/FC1.5/ETPs_quiescent_markerPeak.bed")
output_diffpeak(markerList$nonT_potential,"markerPeak_CT5/FC1.5/nonT_potential_markerPeak.bed")
output_diffpeak(markerList$T_ILCs,"markerPeak_CT5/FC1.5/T_ILCs_markerPeak.bed")
output_diffpeak(markerList$T_ILCs_sub,"markerPeak_CT5/FC1.5/T_ILCsSub_markerPeak.bed")










######## marker peak heatmap
MKpeakCT5cb3_nonT_GR <-intersect(MKpeakCT5cb3_nonT_ETPs_markerList[[1]],MKpeakCT5cb3_nonT_TILCs_markerList[[1]])#,"NonTgETP","NonTgTILC")
MKpeakCT5cb3_ETPs_GR<-intersect(MKpeakCT5cb3_ETPs_nonT_markerList[[1]],MKpeakCT5cb3_ETPs_TILCs_markerList[[1]])#,"ETPgNonT","ETPgTILC")
MKpeakCT5cb3_TILCs_GR<-intersect(MKpeakCT5cb3_TILCs_nonT_markerList[[1]],MKpeakCT5cb3_TILCs_ETPs_markerList[[1]])#,"TILCgNonT","TILCgETP")

combineGR <- c(MKpeakCT5cb3_nonT_GR,MKpeakCT5cb3_ETPs_GR,MKpeakCT5cb3_TILCs_GR)
useidx <- which(MKpeakCT5cb3_all_markerList[[1]] %in% MKpeakCT5cb3_nonT_GR)# | MKpeakCT5cb3_all_markerList[[1]] %in%  MKpeakCT5cb3_ETPs_GR | MKpeakCT5cb3_all_markerList[[1]] %in%  MKpeakCT5cb3_TILCs_GR)



MKpeakCT5cb3_obj <- markersPeaks_CT5cb3[useidx,]

heatmapPeaks <- markerHeatmap(
  seMarker = MKpeakCT5cb3_obj,
  #plotLog2FC=TRUE, 
  labelMarkers="",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot", use_raster = TRUE)
plotPDF(heatmapPeaks, name = "test.pdf", width = 8, height = 6, ArchRProj = proj_DN1_7, addDOC = FALSE)


heatmapPeaks <- markerHeatmap(
  seMarker = MKpeakCT5cb3_obj,
  #plotLog2FC=TRUE, 
  labelMarkers="",
  cutOff = "FDR <= 1 & Log2FC >= 0",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot", use_raster = TRUE)
plotPDF(heatmapPeaks, name = "test2.pdf", width = 8, height = 6, ArchRProj = proj_DN1_7, addDOC = FALSE)





heatmapPeaks <- markerHeatmap(
  seMarker = MKpeakCT5cb3_nonT_ETPs, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

peakByCell <- getMatrixFromProject(
  ArchRProj = proj_DN1_7,
  useMatrix = "PeakMatrix" # or "InsertionMatrix" for read counts
)

peakMat <- peakByCell@assays@data$PeakMatrix
cellNames <- proj_DN1_7$cellNames
CT5cb3Names <- proj_DN1_7$CT5cb3
peakNames_tmp <- getPeakSet(proj_DN1_7)
peakNames <- paste0(as.character(seqnames(peakNames)),":",
      start(peakNames),"-",
      end(peakNames))
rownames(peakMat) <- peakNames
colnames(peakMat) <- cellNames

peakCB <- cbind(apply(peakMat[, cellNames[which(CT5cb3Names == "ETPs")]],1,sum),
                apply(peakMat[, cellNames[which(CT5cb3Names == "nonT_potential")]],1,sum),
                apply(peakMat[, cellNames[which(CT5cb3Names == "T_ILCs")]],1,sum)
                )
colnames(peakCB) <- c("ETPs","nonT_potential","T_ILCs")

cellSums <- Matrix::colSums(peakCB)
normFactor <- 1e6 / cellSums
peakCB_tpm <- peakCB %*% Diagonal(x = normFactor)
colnames(peakCB_tpm) <- colnames(peakCB)
CT5cb3MKpeak_tpm <- peakCB_tpm[c(rownames(MKpeakCT5cb3_ETPs),rownames(MKpeakCT5cb3_nonT),rownames(MKpeakCT5cb3_TILCs)),]

useMat <- scale(CT5cb3MKpeak_tpm)
km<-kmeans(useMat,7)

### function
library(gplots)
bi_heatmap<-function(data0,usecolor,M){
data<-data0
data <- c(as.matrix(data))
data <- sort(data)
temp<-data[round(c(0.010000,0.5,0.99)*length(data))]
p20<-temp[1]
p50<-temp[2]
p80<-temp[3]
zmin=p20
zmax=p80
ColorRamp <- colorRampPalette(usecolor, bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
data0[data0<zmin] <- zmin
data0[data0>zmax] <- zmax
ColorRamp_ex <- ColorRamp[round( (min(data0)-zmin)*10000/(zmax-zmin) ) : round( (max(data0)-zmin)*10000/(zmax-zmin) )]
image(1:ncol(data0), 1:nrow(data0), t(data0), axes=F, col=ColorRamp_ex, xlab="", ylab="",main=M,useRaster=T,cex.main=1.5)
#image(1:ncol(data0), 1:nrow(data0), t(data0), axes=F, col=ColorRamp_ex, xlab=paste(zmin,zmax,sep=":"), ylab="",main=M,useRaster=T,cex.main=2)
#axis(side=2)
#axis(side=1,at=c(1,100,200),labels=c("-5k","0","5k"))
#box()
}
bi_heatmap(useMat[order(km$cluster),],c("blue","white","red"),"MKpeak")
#abline(km$size)


library(pheatmap)
tpmMatrix <- as.matrix(CT5cb3MKpeak_tpm)
tpmMatrix[is.na(tpmMatrix) | is.infinite(tpmMatrix)] <- 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Generate the heatmap without clustering and keeping the input order
pheatmap(
  tpmMatrix,
  scale = "row",              # No scaling, assuming input is already z-scores
  cluster_rows = FALSE,        # Do not cluster rows
  cluster_cols = FALSE,        # Do not cluster columns
  color = color_palette,       # Use the defined color palette
  fontsize = 10,               # Font size for the text in the heatmap
  fontsize_row = 8,            # Font size for row names
  fontsize_col = 12,           # Font size for column names (cell types)
  show_rownames = FALSE,       # Set to TRUE if you want to show row names (peaks)
  show_colnames = TRUE         # Set to TRUE to show column names (cell types)
)

ggsave("testheatmap.pdf", plot = last_plot(), width = 10, height = 8)

intersect(rownames(MKpeakCT5cb3_ETPs),rownames(MKpeakCT5cb3_nonT))



length(intersect(rownames(MKpeakCT5cb3_ETPs),rownames(MKpeakCT5cb3_nonT)))




# check batcheffect
#proj_DN1_2 <- addHarmony(
#    ArchRProj = proj_DN1_2,
#    reducedDims = "IterativeLSI",
#    name = "Harmony",
#    groupBy = "Sample",force=T
#)
#proj_DN1_2_HAR <- addClusters(
#    input = proj_DN1_2,
#    reducedDims = "Harmony",
#    method = "Seurat",
#    name = "Clusters",
#    resolution = 0.8,
#    force=T,seed=1
#)
#
#proj_DN1_2_HAR <- addUMAP(
#    ArchRProj = proj_DN1_2_HAR, 
#    reducedDims = "Harmony", 
#    name = "UMAP", 
#    nNeighbors = 30, 
#    minDist = 0.5, 
#    metric = "cosine",force=T
#)
#
#cM <- confusionMatrix(paste0(proj_DN1_2$Clusters), paste0(proj_DN1_2$Sample))
#write.table(as.matrix(cM), file="cluster_sample_confusionMap.txt",row.names=T,col.names=T,sep="\t",quote=F)
#
#cM <- cM / Matrix::rowSums(cM)
#p <- pheatmap::pheatmap(
#    mat = as.matrix(cM), 
#    color = paletteContinuous("whiteBlue"), 
#    border_color = "black"
#)
#plotPDF(p, name = "confusionMap_heatmap_LSI.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE)
#
#cM_HAR <- confusionMatrix(paste0(proj_DN1_2_HAR$Clusters), paste0(proj_DN1_2_HAR$Sample))
#cM_HAR <- cM_HAR / Matrix::rowSums(cM_HAR)
#p2 <- pheatmap::pheatmap(
#    mat = as.matrix(cM_HAR), 
#    color = paletteContinuous("whiteBlue"), 
#    border_color = "black"
#)
#
#plotPDF(p2, name = "confusionMap_heatmap_HAR.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE)
#write.table(cM_HAR, file="cluster_sample_confusionMap_HAR.txt",row.names=T,col.names=T,sep="\t",quote=F)




#p1 <- plotEmbedding(ArchRProj = proj_DN1_2_HAR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
#p2 <- plotEmbedding(ArchRProj = proj_DN1_2_HAR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
#plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_HAR.pdf", ArchRProj = proj_DN1_2_HAR, addDOC = FALSE, width = 5, height = 5)
#


# gene score with ArchR default method
proj_DN1_2 <- addGeneScoreMatrix(proj_DN1_2,force=TRUE)
proj_DN1_2 <- addImputeWeights(proj_DN1_2,seed=1)

saveArchRProject(ArchRProj = proj_DN1_2, outputDirectory = "SaveDN1_2", load = FALSE)
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


#pdf(file="UMAP_knowledgeKeyGeneExp.pdf",width=20,height=16)
#FeaturePlot(proj_DN1_2,features=usekeygene,ncol=5)
#dev.off()
#
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
#pdf(file=paste0("UMAP_geneScore_Survival.pdf"),width=12,height=6)
#cowplot::plot_grid(plotlist=keygene_Survival_Img[1:8],nrow=2,ncol=4)
#dev.off()
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
#p5heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Survival),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p6heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Dendritic),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p7heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Myeloid),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p8heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_InnateLymphoid),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p9heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_Cebp),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p10heatmapGS <- plotMarkerHeatmap(labelMarkers = c(keygene_other),seMarker = markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
p1 <- ComplexHeatmap::draw(p1heatmapGS, column_title = "keygene",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p2 <- ComplexHeatmap::draw(p2heatmapGS, column_title = "HSC",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p3 <- ComplexHeatmap::draw(p3heatmapGS, column_title = "T lineage",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p4 <- ComplexHeatmap::draw(p4heatmapGS, column_title = "B lineage",heatmap_legend_side = "bot", annotation_legend_side = "bot")
#p5 <- ComplexHeatmap::draw(p5heatmapGS, column_title = "Survival",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p6 <- ComplexHeatmap::draw(p6heatmapGS, column_title = "Dendritic",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p7 <- ComplexHeatmap::draw(p7heatmapGS, column_title = "Myeloid",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p8 <- ComplexHeatmap::draw(p8heatmapGS, column_title = "InnateLymphoid",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p9 <- ComplexHeatmap::draw(p9heatmapGS, column_title = "Cebp",heatmap_legend_side = "bot", annotation_legend_side = "bot")
p10 <- ComplexHeatmap::draw(p10heatmapGS, column_title = "other",heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(p1, name = "MKGexp_heatmap_keygene.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2, name = "MKGexp_heatmap_HSC.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p3, name = "MKGexp_heatmap_Tlineage.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p4, name = "MKGexp_heatmap_Blineage.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
#plotPDF(p5, name = "MKGexp_heatmap_Survival.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p6, name = "MKGexp_heatmap_Dendritic.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p7, name = "MKGexp_heatmap_Myeloid.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p8, name = "MKGexp_heatmap_InnateLymphoid.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p9, name = "MKGexp_heatmap_Cebp.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p10, name = "MKGexp_heatmap_other.pdf", ArchRProj = proj_DN1_2, addDOC = FALSE, width = 5, height = 5)


