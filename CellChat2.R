# do both objects separately
# repeat these following lines for control cells
library("NMF")
library("circlize")
library("ComplexHeatmap")
library("CellChat")
library(patchwork)
library(Seurat)
library(future)
options(stringsAsFactors = FALSE)

PROJECT="/path/"
setwd(PROJECT)
seuratObj=readRDS("PD.rds")  # subset your complete objects to the PD cells only, and to the clusters you wish to include in the analysis

Idents(seuratObj)= "Subpop"
head(seuratObj@meta.data)
table(seuratObj@active.ident)

# remove any space, plus or minus signs
seuratObj@meta.data$Subpop <- gsub(' ', '_', seuratObj@meta.data$Subpop)
seuratObj@meta.data$Subpop <- gsub('\\+', '_Pos', seuratObj@meta.data$Subpop)
seuratObj@meta.data$Subpop <- gsub('-', '_', seuratObj@meta.data$Subpop)
Idents(seuratObj)= "Subpop"
table(seuratObj@active.ident)
# you need the sample column to be called "samples"
names(seuratObj@meta.data)[names(seuratObj@meta.data) == "SubPool"] <- "samples"
head(seuratObj@meta.data,2)

# Seurat 4 : @
#data.input <- seuratObj[["RNA"]]@data # normalized data matrix

# For Seurat version >= “5.0.0”,: $
# get the normalized data via `seurat_object[["RNA"]]$data`
data.input <- seuratObj[["RNA"]]$data # normalized data matrix

Idents(seuratObj)="Subpop"
labels <- Idents(seuratObj)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
head(meta)
head(seuratObj@meta.data,1)
cellChat <- createCellChat(object = seuratObj, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
pdf("CellChatDB.pdf", height=5, width=8)
showDatabaseCategory(CellChatDB)
dev.off()

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB # I want to use them all
cellChat@DB <- CellChatDB.use
cellchat=cellChat
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 12) # do parallel
options(future.globals.maxSize = 400000 * 1024^2) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#  For a faster implementation of the Wilcoxon Test, please install the presto package
#--------------------------------------------
# devtools::install_github('immunogenomics/presto')
#--------------------------------------------
# Otherwise, plase set `do.fast = FALSE` for running the standard Wilcoxon Test!
# if you still want to use the standard ones:
#cellchat <- identifyOverExpressedGenes(cellchat,do.fast = FALSE)
#cellchat <- identifyOverExpressedInteractions(cellchat)

ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,population.size = TRUE, nboot=20) # allows permutations
cellchat <- filterCommunication(cellchat, min.cells = 5) # standard=10
cellchat
df.net <- subsetCommunication(cellchat) 
head(df.net)
write.table(df.net, "dfnet.txt",col.names=NA, sep="\t")
library(readr)
write_csv(df.net, "dfnet.csv")
df.netPathways <- subsetCommunication(cellchat, slot.name = "netP") 
head(df.netPathways)
write.table(df.netPathways, "dfnetPathways.txt",col.names=NA, sep="\t")
write_csv(df.netPathways, "dfnetPathways.csv")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
groupSize <- as.numeric(table(cellchat@idents))
pdf("All_interactions.pdf", height=8, width=14)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
pdf("Interactions_per_cluster.pdf", height=20, width=25)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
saveRDS(cellchat, "cellchat_2.rds")
#cellchat=readRDS("cellchat_2.rds")
print(cellchat@netP$pathways)
Idents(seuratObj)="Subpop"
table(seuratObj@active.ident)

for (i in 1:length(cellchat@netP$pathways)){
pathways.show <- cellchat@netP$pathways[i]
print(pathways.show) 

pdf(paste0("Significant_pathway_",pathways.show,"_4pages.pdf"), height=15, width=15)
par(mfrow=c(1,1))
print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")) 
print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord"))
print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds"))
print(netAnalysis_contribution(cellchat, signaling = pathways.show))
dev.off()
}

df=read.table("dfnet.txt", sep="\t", header=T)
head(df) # 
a=df$source
b=df$target
c=c(a,b)
d=unique(c)
d

for (i in (1:length(d))){                                                                                                     
SOURCE=d[i]
print(SOURCE)
TARG=d[-i]
print(TARG)

for (j in(1:length(TARG))){
TARGET=TARG[j]
print(TARGET)

tryCatch({
pdf(paste0("chord_all_interactions_source_",SOURCE,"_and_target_",TARGET,".pdf"), height=15, width=15)
netVisual_chord_gene(cellchat, sources.use = SOURCE, targets.use = TARGET,  title.name = paste0("Signaling from ",SOURCE," to ",TARGET), legend.pos.x = 10)
dev.off()
},error=function(e){})
}}


# once these have been done (until the cellchat_2.rds and dfnet.txt objects)
# do the comparison of the 2 objects:

library("NMF")
library("circlize")
library("ComplexHeatmap")
library("CellChat")
packageVersion("CellChat") # 2.1.2
library(patchwork)
library(Seurat)
ptm = Sys.time()
options(stringsAsFactors = FALSE)

data.dir <- '/path/comparison'
setwd(data.dir)

data_control=readRDS("/path/to/control/cellchat_2.rds")
df_control=read.table("/path/to/control/dfnet.txt", sep="\t", header=T)
x=unique(df_control$pathway_name)
x

data_disease=readRDS("/path/to/PD/cellchat_2.rds")
df_disease=read.table("/path/to/PD/dfnet.txt", sep="\t", header=T)

levels(data_control@idents) 
levels(data_disease@idents) 

y=unique(df_disease$pathway_name)
y

join=c(x,y)
JOIN=unique(join)
JOIN

plan("multicore", workers = 12) 
options(future.globals.maxSize = 400000 * 1024^2) 

data_control=netAnalysis_computeCentrality(data_control)
data_disease=netAnalysis_computeCentrality(data_disease)

# same indexes?
levels(data_control@idents) 
levels(data_disease@idents)  

    # if not the same order:
    data_disease=updateClusterLabels(
      data_disease,
      old.cluster.name = NULL,
      new.cluster.name = NULL,
      new.order = c("All_neurons","HSP_Glia","APOE_Pos_Glia","OTOGL_Pos_Glia","IGSF21_Pos_MF","APOE_Pos_MF","NAMPT_Pos_MF"),
      new.cluster.metaname = "new.labels"
    )
    # is it the same as control now?
    levels(data_disease@idents) # YES <3
    
object.list <- list(disease = data_disease, CTRL = data_control)
object.list 
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
save(object.list, file = "cellchat_diseaseandCtrl_list.RData")
save(cellchat, file = "cellchat_diseaseandCtrl.RData")

colors <- c("disease" = "red", "CTRL" = "blue")

for (CLUSTER in c("All_neurons","HSP_Glia","APOE_Pos_Glia","OTOGL_Pos_Glia","IGSF21_Pos_MF","APOE_Pos_MF","NAMPT_Pos_MF")){
pdf(paste0("compare_communications_bubbleplot_",CLUSTER,".pdf"), height=30, width=10)
print(netVisual_bubble(cellchat, sources.use = CLUSTER, targets.use =NULL,  comparison = c(1, 2), angle.x = 45,color.text.use = TRUE, color.text = colors))
print(netVisual_bubble(cellchat, sources.use = NULL, targets.use =CLUSTER,  comparison = c(1, 2), angle.x = 45,,color.text.use = TRUE, color.text = colors))
dev.off()
}
for (CLUSTER in c("All_neurons","APOE_Pos_Glia")){
pdf(paste0("compare_communications_bubbleplot_",CLUSTER,".pdf"), height=50, width=10)
print(netVisual_bubble(cellchat, sources.use = CLUSTER, targets.use =NULL,  comparison = c(1, 2), angle.x = 45,color.text.use = TRUE, color.text = colors))
print(netVisual_bubble(cellchat, sources.use = NULL, targets.use =CLUSTER,  comparison = c(1, 2), angle.x = 45,,color.text.use = TRUE, color.text = colors))
dev.off()
}

ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

pdf("comparison_barplot.pdf", height=6, width=10)
gg1 + gg2
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
gg
class(gg)
saveRDS(gg, file="gg.RData")
pdf("incoming_and_outgoing_interaction_strength.pdf", height=6, width=8)
patchwork::wrap_plots(plots = gg)
dev.off()

head(df_control) # same clusters as the other, only take one
# take unique of source and target, just in case one cluster is not in 1 column
a=df_control$source
b=df_control$target
c=c(a,b)
d=unique(c)
d

for (i in (1:length(d))){                                                                                                     
CLUSTER=d[i]
print(CLUSTER)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = CLUSTER, signaling.exclude = "")
pdf(paste0("sign_changes_in_",CLUSTER,".pdf"), height=6, width=10)
print(gg1)
dev.off()
}

colors <- c("disease" = "red", "CTRL" = "blue")
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, color.use=colors)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE, color.use=colors)
pdf("compare_signaling_pathways_barplot.pdf", height=15, width=12)
gg1 + gg2
dev.off()

# dataset 1 is disease
# dataset 2 is control
pdf("compare_communications_bubbleplot.pdf", height=80, width=24)
netVisual_bubble(cellchat, sources.use = NULL, targets.use =NULL,  comparison = c(1, 2), angle.x = 45)
dev.off()


gg1 <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL,  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ",names(object.list[2])), angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL,  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in ",names(object.list[2])), angle.x = 45, remove.isolate = T)
pdf(paste0("up_down_in_",names(object.list[2]),"_communications_bubbleplot.pdf"), height=60, width=25)
gg1 + gg2
dev.off()

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("disease", "CTRL")) # set factor level
for (i in 1:length(JOIN)) {
    pathways.show=JOIN[[i]]
    print(pathways.show)
tryCatch({                                                                                                                                                                                                  
pdf(paste0("specific_pathway_",pathways.show,".pdf"), height=10, width=15)
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
 print(netVisual_aggregate(object.list[[i]],sources.use = NULL, targets.use = NULL,remove.isolate = TRUE, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i])))
}
print(netVisual_heatmap(object.list[[1]],sources.use = NULL, targets.use = NULL, remove.isolate = TRUE,signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[1])))
print(netVisual_heatmap(object.list[[2]],sources.use = NULL,targets.use = NULL, remove.isolate = TRUE,signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[2])))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  print(netVisual_aggregate(object.list[[i]],sources.use = NULL, targets.use = NULL,remove.isolate = TRUE, signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i])))
}
colors <- c("disease" = "red", "CTRL" = "blue")
print(plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T, type = "violin", color.use=colors))
dev.off()
dev.list()
for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
   dev.off()
}
dev.list()
},error=function(e){}) 

}

for (i in (1:length(d))){                                                                                                     
SOURCE=d[i]
print(SOURCE)
TARG=d[-i]
print(TARG)

for (j in(1:length(TARG))){
TARGET=TARG[j]
print(TARGET)

tryCatch({
pdf(paste0("chord_all_interactions_source_",SOURCE,"_and_",TARGET,".pdf"), height=10, width=15)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = SOURCE, targets.use = TARGET,  title.name = paste0("Signaling from ",SOURCE," to ",TARGET," ", names(object.list)[i]), legend.pos.x = 10)
}
dev.off()
},error=function(e){})
}}
