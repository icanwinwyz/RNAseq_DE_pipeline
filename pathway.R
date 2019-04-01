source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
biocLite("gage")
biocLite("gageData")
#install.packages("pathfindR")
library(pathfindR)
data("RA_input")
head(RA_input)
data <- read.csv("M_vs_MC_DEGs_padj0.05.csv")
gene_list <- gsub("\\..*","",data$Row.names) #remove all after "."
data$Row.names <- gene_list
colnames(data)[1]<-"Gene"
colnames(data)[3]<-"log2FC"
colnames(data)[5]<-"p_adj"
head(data)
data$symbol = mapIds(org.Hs.eg.db,
                     keys=data$Gene, 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
data$entrez = mapIds(org.Hs.eg.db,
                    keys=data$Gene, 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
data$name =   mapIds(org.Hs.eg.db,
                    keys=data$Gene, 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")

pathway_input <- data[, c("Gene","log2FC","p_adj", "symbol", "entrez", "name")]
pathway_input <- data[, c("symbol","log2FC","p_adj")]
RA_output <- run_pathfindR(RA_input,gene_sets="GO-BP")
KEGG_output <- run_pathfindR(pathway_input,human_genes = TRUE, gene_sets="KEGG", output_dir=paste0(name,"_KEGG")) # KEGG is the default gene_sets
KEGG_output$Count <- (lengths(gregexpr(",", KEGG_output$Up_regulated)) + 1) + (lengths(gregexpr(",", KEGG_output$Down_regulated)) + 1)
# add Count col to show how many genes (both upR and downR) in each pathway
pdf(paste(name,"DE_gene_KEGG.pdf",sep = "_"),12,20)
#ggplot(test,aes(x=Fold_Enrichment_Score,y=Term,size=Gene_Number,colour=adjusted_pvalue))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment Score")+ylab("Term")+ggtitle("Enriched KEGG pathway")
ggplot(KEGG_output,aes(x=Fold_Enrichment,y=Pathway,size=Count,colour=lowest_p))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment")+ylab("Pathway")+ggtitle("Enriched KEGG pathway")
dev.off()


pdf(paste(name,"_KEGG.pdf",sep=""),16,12)
pathway_output_1 <- run_pathfindR(pathway_input,human_genes = TRUE, gene_sets="GO-BP", output_dir=paste0(name,"_GO_BP")) 
dev.off()
dev.off()
write.csv(pathway_output_1, paste0(name,"_KEGG.csv"))
library(gage)
library(gageData)
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
data(kegg.sets.hs)
data(sigmet.idx.hs)
#The main gage() function requires a named vector of fold changes, where the names of the values are the Entrez gene IDs.
foldchanges = pathway_input$log2FC
names(foldchanges) = pathway_input$entrez
head(foldchanges)

data(go.sets.hs)
data(go.subs.hs)
data(kegg.gs)
head(kegg.gs[[1]]) 

gobpsets = go.sets.hs[go.subs.hs$BP]
goccsets = go.sets.hs[go.subs.hs$CC]
gomfsets = go.sets.hs[go.subs.hs$MF]

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)
head(gobpres$greater[,1:5])
gobpres_up <- as.data.frame(gobpres$greater)
gobpres_down <- as.data.frame(gobpres$less)
#gobpres_up <- as.data.frame(gobpres_up)
gobpres_up <- subset(gobpres_up, p.val < 0.1)
gobpres_down <- subset(gobpres_down, p.val < 0.1)
#GO_terms <- gsub(".* ","",rownames(gobpres_up)) #remove all and before up to " "
gobpres_up$Term <- GO_terms
pdf(paste(name,"DE_gene_GO_BP.pdf",sep = "_"),12,20)
ggplot(as.data.frame(gobpres_up),aes(x=stat.mean,y=rownames(gobpres_up),size=set.size,colour=p.val))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment")+ylab("GO terms")+ggtitle("Enriched GO terms (upR genes)")
ggplot(as.data.frame(gobpres_down),aes(x=stat.mean,y=rownames(gobpres_down),size=set.size,colour=p.val))+geom_point()+scale_colour_gradient(low="red",high="blue")+xlab("Fold Enrichment")+ylab("GO terms")+ggtitle("Enriched GO terms (downR genes)")
dev.off()
dev.off()


kegg = gage(foldchanges, gsets=kegg.gs, same.dir=FALSE, use.fold = F)
kegg <- as.data.frame(kegg$greater)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")
devtools::install_github( c("guangchuangyu/enrichplot", "guangchuangyu/DOSE", "guangchuangyu/clusterProfiler"))
library(enrichplot)
library(clusterProfiler)



source("https://bioconductor.org/biocLite.R")
biocLite("fgsea")
library(fgsea)
Sys.getenv()['R_HOME']
data(examplePathways)
data(exampleRanks)

# pathways are stored in a list
class(examplePathways)
# List

# the example list has 1,457 "pathways"
length(examplePathways)
#[1] 1457

# this is the first pathway, which contains
# all the Entrez gene IDs of the Meiotic_Synapsis pathway
examplePathways[1]


data(exampleRanks)

# ranks are stored as a numeric vector
class(exampleRanks)
#[1] "numeric"

# the names of the vector are the Entrez gene IDs
head(exampleRanks)
#170942    109711     18124     12775     72148     16010 
#-63.33703 -49.74779 -43.63878 -41.51889 -33.26039 -32.77626

# the values are the rank metric
head(unname(exampleRanks))
#[1] -63.33703 -49.74779 -43.63878 -41.51889 -33.26039 -32.77626

# ranks are from lowest to highest
tail(unname(exampleRanks))


fgseaRes <- fgsea(pathways = examplePathways, 
                  stats = exampleRanks,
                  minSize=15,
                  maxSize=500,
                  nperm=100000)

# see https://github.com/Rdatatable/data.table/wiki
# if you are unfamiliar with the data.table package
class(fgseaRes)
#[1] "data.table" "data.frame"

# top 6 enriched pathways
head(fgseaRes[order(pval), ])
#pathway         pval         padj        ES      NES nMoreExtreme
#1:                  5990980_Cell_Cycle 1.226723e-05 0.0002701041 0.5388497 2.680900            0
#2:         5990979_Cell_Cycle,_Mitotic 1.255556e-05 0.0002701041 0.5594755 2.743471            0
#3:    5991210_Signaling_by_Rho_GTPases 1.317141e-05 0.0002701041 0.4238512 2.010545            0
#4:                     5991454_M_Phase 1.375667e-05 0.0002701041 0.5576247 2.552055            0
#5: 5991023_Metabolism_of_carbohydrates 1.391188e-05 0.0002701041 0.4944766 2.240926            0
#6:        5991209_RHO_GTPase_Effectors 1.394739e-05 0.0002701041 0.5248796 2.373183            0
#size                             leadingEdge
#1:  369   66336,66977,12442,107995,66442,19361,
#2:  317   66336,66977,12442,107995,66442,12571,
#3:  231 66336,66977,20430,104215,233406,107995,
#4:  173   66336,66977,12442,107995,66442,52276,
#5:  160    11676,21991,15366,58250,12505,20527,
#6:  157 66336,66977,20430,104215,233406,107995,

# number of significant pathways at padj < 0.01
sum(fgseaRes[, padj < 0.01])
#[1] 77

# plot the most significantly enriched pathway
plotEnrichment(examplePathways[[head(fgseaRes[order(pval), ], 1)$pathway]],
               exampleRanks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)



data <- read.csv("M_vs_MC_DEGs_padj0.05.csv")
gene_list <- gsub("\\..*","",data$Row.names) #remove all after "."
data$Row.names <- gene_list
colnames(data)[1]<-"Gene"
colnames(data)[3]<-"log2FC"
colnames(data)[5]<-"p_adj"
head(data)
data$symbol = mapIds(org.Hs.eg.db,
                     keys=data$Gene, 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
data$entrez = mapIds(org.Hs.eg.db,
                     keys=data$Gene, 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
data$name =   mapIds(org.Hs.eg.db,
                     keys=data$Gene, 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

input <- data[, c("entrez","log2FC")]
input.ord <- input[ order(-input[,"log2FC"]), ]
