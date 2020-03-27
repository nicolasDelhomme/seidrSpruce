#' ---
#' title: "Seidr ROC curves"
#' author: "Aaron Ayllon Benitez and Nicolas Delhome"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' * Libraries
suppressPackageStartupMessages({
  library(gplots)
  library(here)
  library(magrittr)
  library(matrixStats)
  library(pander)
  library(pracma)
  library(tidyverse)
})

#' * Colors
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Functions
plotRoc <- function(f){
  
  dat <- read_tsv(f,
                  col_names = c("TPR","FPR","PR","ALGORITHM"),
                  col_types = cols(.default=col_double(),
                                   ALGORITHM = col_character()),
                  comment="#") %>% 
    mutate(TNR=1-FPR,FNR=1-TPR,balanced_accuracy=(TPR+TNR)/2) %>% 
    group_by(ALGORITHM) 
  
  vals <- read_lines(f) %>% subset(grepl("#",.))
  lvals <- length(vals)/3
  tabs <- tibble(ALGO=sapply(strsplit(vals[lvals + seq(1,length.out=lvals,by=2)],"\t"),"[",2),
                 positiveEdges=as.integer(sapply(strsplit(vals[1:lvals],"\t"),"[",2)),
                 negativeEdges=as.integer(sapply(strsplit(vals[1:lvals],"\t"),"[",3)),
                 AUC=as.double(sub(".* ","",sapply(strsplit(vals[lvals + seq(1,length.out=lvals,by=2)],"\t"),"[",1))),
                 AUPR=as.double(sub(".* ","",sapply(strsplit(vals[lvals + seq(2,length.out=lvals,by=2)],"\t"),"[",1))),
  ) 
  
  dat %<>% group_modify(~mutate(.,
                             p_edge=tabs[tabs$ALGO==.y[[1]],"positiveEdges",drop=TRUE],
                             n_edge=tabs[tabs$ALGO==.y[[1]],"negativeEdges",drop=TRUE]
                             ))
  dat %<>% mutate(TP=TPR * p_edge,
                  TN=(1-FPR) * n_edge,
                  FP=FPR * n_edge,
                  accuracy=(TP+TN)/(p_edge+n_edge),
                  PPCR=(TP+FP)/(p_edge+n_edge))
  
  aucs <- unlist(dat %>% group_by(ALGORITHM) %>% group_map(~round(trapz(.x$FPR,.x$TPR),digits=3)))
  names(aucs) <- attr(dat,"groups")[,1,drop=TRUE]
  
  p <- ggplot(dat,aes(x=FPR,y=TPR,col=ALGORITHM,group=ALGORITHM)) +
    geom_line() + 
    scale_x_continuous(name="1-specificity (FPR)") + 
    scale_y_continuous(name="sensitivity (TPR)") +
    ggtitle(label=paste(sub("\\.roc","",basename(f)), " ROC curve"))+
    theme_classic()
  
  suppressMessages(suppressWarnings(plot(p)))
  
  p2 <- ggplot(dat %>% mutate(R=TPR,P=TPR/(TPR+FPR)),
               aes(x=R,y=P,col=ALGORITHM,group=ALGORITHM)) +
    geom_line() + 
    scale_x_continuous(name="Recall") + 
    scale_y_continuous(name="Precision") +
    ggtitle(label=paste(sub("\\.roc","",basename(f)), " PR curve"))+
    theme_classic()
  
  suppressMessages(suppressWarnings(plot(p2)))
  
  sequence=rep(seq(0,1,length.out=nrow(dat)/nrow(attr(dat,"groups"))),nrow(attr(dat,"groups")))
  
  p3 <- ggplot(dat,aes(x=sequence,y=accuracy,col=ALGORITHM,group=ALGORITHM)) +
    geom_line() + 
    scale_x_continuous(name="Normalized rank") + 
    scale_y_continuous(name="Accuracy") +
    ggtitle(label=paste(sub("\\.roc","",basename(f)), " Accuracy curve"))+
    theme_classic()
  
  suppressMessages(suppressWarnings(plot(p3)))
  
  p4 <- ggplot(dat,aes(x=sequence,y=PPCR,col=ALGORITHM,group=ALGORITHM)) +
    geom_line() + 
    scale_x_continuous(name="Normalized rank") + 
    scale_y_continuous(name="Predictive Positive Condition Rate") +
    ggtitle(label=paste(sub("\\.roc","",basename(f)), " PPCR curve"))+
    theme_classic()
  
  suppressMessages(suppressWarnings(plot(p4)))
  
  
  return(list(stats=tabs,auc=aucs))
}

#' # ROC
#' ## Aggregated
#' ```{R CHANGEME1, echo=FALSE, eval=FALSE}
#' Change the path to the aggregated ROC results, if required. That ROC file must have been created using seidr roc -a option
#' ```
res <- plotRoc(here("data/seidr/roc/aggregated.roc"))

#' ### Stats of the gold standard analysis
pander(res$stats)

#' ## Backbone
#' ```{R CHANGEME2, echo=FALSE, eval=FALSE}
#' Change the path to the backbone ROC results directory as well as the file matching patter,if required. 
#' These ROC files must have been created using seidr roc -a option
#' ```
files <- dir(here("data/seidr/roc"),pattern="backbone.*\\.roc",full.names=TRUE)
names(files) <- gsub(".roc","",basename(files))
files <- files[order(as.integer(sub("aggregated-backbone-p","",names(files))))]
resb <- lapply(files,plotRoc)

#' ### Stats of the gold standard (GS) analysis
sts <- lapply(names(files),function(n,resb){
  resb[[n]]$stats
},resb)

names(sts) <- names(files)

pander(sts)

#' # Summary
#' Report all AUCs
aucs <- cbind(sapply(resb,"[[","auc"),aggregated=res$auc)

pander(aucs)

#' ## Heatmaps
#' ### AUCs
heatmap.2(aucs,trace="none",col=hpal,margins=c(7.1,7.1))

heatmap.2(aucs,trace="none",col=hpal,margins=c(7.1,7.1),Colv=FALSE,dendrogram="row")

#' ### AUCs penalised by the number of GS edges used
resb$aggregated=res

gsNum <- sapply(lapply(lapply(resb,"[[","stats"),"[",2:3),rowSums)

paucs <- aucs * t(t(gsNum) / colMaxs(gsNum))

pander(paucs)

heatmap.2(paucs,trace="none",col=hpal,margins=c(7.1,7.1),Colv=FALSE,dendrogram="row")

#' # AUPR
auprs <- cbind(sapply(sts,"[[","AUPR"),aggregated=res$stats$AUPR)
rownames(auprs) <- rownames(aucs)
pander(auprs)

#' The difference is striking to that of the AUCs
heatmap.2(auprs,trace="none",col=hpal,margins=c(7.1,7.1),Colv=FALSE,dendrogram="row")

#' Penalising the AUCs by the AUPRs is actually interesting
heatmap.2(aucs * auprs,trace="none",col=hpal,margins=c(7.1,7.1),Colv=FALSE,dendrogram="row")

#' # Conclusion
#' From the AUC and AUPR results a backbone cutoff at 1% seems adequate. Taking a cutoff at 1% 
#' (which maximises the `irp` AUC and AUPR) might be too conservative in the number of
#' edges retained. 
#' 
#' 
#' ```{r empty, eval=FALSE, echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```