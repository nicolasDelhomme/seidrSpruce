#' ---
#' title: "InfoMap analysis to Picea abies metanetwork"
#' author: "Aaron Ayllon Benitez and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' * Libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(pander)
  library(pracma)
  library(jsonlite)
  library(parallel)
  library(tidyverse)
  library(KEGGREST)
})
#' Data
load("/mnt/picea/home/aabenitez/Git/seidrSpruce/R/gopher.RData")
clusterResolve <- read.delim("/mnt/picea/projects/spruce/facility/seidr/infomap-bb-p1/clusterResolve.txt", 
                             header=FALSE,stringsAsFactors=FALSE)
# Add the names to the data frame
numberLevels <- length(clusterResolve) -3 - 1
levels <- sapply(rep(1:numberLevels),function(x){return(paste("Level",x,sep = ""))})
colnames(clusterResolve) <- c('Path', levels, 'Flow', 'Index', 'Gene')
# create a new results data frame
res <- data.frame(Gene = clusterResolve[,ncol(clusterResolve)],
                  Path = clusterResolve[,1],
                  stringsAsFactors = F)

#' Infomap data always have level 1 and none of the data is NA, we can add it directly
res$Level1 <- clusterResolve[,2]

#' loop though the rest of the levels and attach the name from the previous ones
for (level in 2:numberLevels) {
  
  currentLevel <- paste0("Level",level)
  prevLevel <- paste0("Level",(level-1))
  
  # join names
  res[[currentLevel]] <- paste0(res[[prevLevel]], ":", clusterResolve[[currentLevel]])
  
  # if there is an NA inside the current name, that gene doesn't belong to a cluster in that level, it turns into NA
  res[[currentLevel]] <- base::ifelse(is.na(res[[currentLevel]]), NA, res[[currentLevel]])  
}

pander(data.frame("ncl_L1"=length(unique(res$Level1)), 
           "ncl_L2"=length(unique(res$Level2)),
           "ncl_L3"=length(unique(res$Level3))))



gene.freq <- table(res$Level1)
gene.cumfreq <- cumsum(gene.freq)

cluster.frequency.per = data.frame("cluster"=as.numeric(names(gene.freq)),
           "gene.freq.per"=as.numeric(gene.freq)/length(res$Gene)*100,
           "gene.cumfreq.per"=gene.cumfreq/length(res$Gene)*100)

cluster.freq.per.melt <- reshape2::melt(cluster.frequency.per,id=1)

chosenX = 6
#+ fig.width=10, fig.height=4
ggplot(subset(cluster.freq.per.melt, cluster<51), aes(x=cluster, y=value, group=variable))+
  geom_line(aes(linetype=variable))+ 
  scale_y_continuous(name="gene percentage\n (gene in Cluster/total genes)", 
                     limit=c(0,100), breaks= seq(0,100,by=10)) +
  scale_x_continuous(name="clusters", breaks =  seq(2, 50, by = 2)) +
  geom_vline(xintercept=chosenX, color='red')+
  annotate("text", x=chosenX+1.8, y=cluster.frequency.per$gene.cumfreq.per[chosenX]+4,
           label= paste(round(cluster.frequency.per$gene.cumfreq.per[chosenX],2),"%",sep=""),
           color='red')+ 
  theme_classic()+theme(legend.title=element_blank())


#' ## Level 2

res.filtered <- res[res$Level1<chosenX+1,]

pander(data.frame("ncl_L1"=length(unique(res.filtered$Level1)), 
                  "ncl_L2"=length(unique(res.filtered$Level2)),
                  "ncl_L3"=length(unique(res.filtered$Level3))))


geneSet.test <- res.filtered[res.filtered$Level2=="1:1",1]

cluster.enrichment.annotation<-mclapply(unique(res.filtered$Level2),
         function(x){
           return(length(res.filtered[res.filtered$Level2==x,1]))
           }, mc.cores=1)

names(cluster.enrichment.annotation) = unique(res.filtered$Level2)

result <- gopher(geneSet.test, task= list("go","pfam","kegg","mapman"),url="pabies")

alpha=.1
#write.table(enrichment$kegg,file="analysis/DESeq2/autumn-vs-winter_KEGG.txt",
#            sep="\t",quote = FALSE, row.names=FALSE)
#' Find the pathway of the significantly enriched enzymes
pathways <- keggLink("pathway",enrichment$kegg$id[enrichment$kegg$padj<=alpha])
#' Pathways are duplicated, clean up
stopifnot(length(grep("map",pathways)) == length(grep("ec",pathways)))
pathways <- pathways[grep("map",pathways)]
#' Proportion of the pathways
barplot(sort(table(sub("path:","",pathways)),decreasing = TRUE),las=2)
#' Get the pathway info
pinfo <- lapply(split(unique(pathways), ceiling(seq_along(unique(pathways))/10)),keggGet)
pathway.df <- do.call(rbind,lapply(pinfo,function(x){data.frame(ENTRY=sapply(x,"[[","ENTRY"),
                                                NAME=sapply(x,"[[","NAME"),
                                                CLASS=unlist(lapply(x,function(y){
                                                                          return(ifelse(!is.null(y$CLASS),
                                                                                        y$CLASS,NA)
                                                                                 )})))}))
tab <- table(sub("path:","",pathways))
pathway.df$OCCURENCE <- tab[match(pathway.df$ENTRY,names(tab))]
write.table(pathway.df,file="analysis/DESeq2/autumn-vs-winter_KEGG-pathway.txt",
            sep="\t",quote = FALSE, row.names=FALSE)
#' Create a wordcloud of the pathway names
wordcloud(pathway.df$NAME,pathway.df$OCCURENCE/sum(pathway.df$OCCURENCE),
          colors = hpal,scale = c(2,.5),rot.per = 0)
