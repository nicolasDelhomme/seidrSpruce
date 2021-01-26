#' ---
#' title: "Annotation QA to Picea abies"
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
  library(curl)
  library(stringr)
  library(scales)
  library(gridExtra)
  library(grid)
  library(ontologyReader)
})
#' # Functions
level1Visu <- function(go,GOA,t){
  c <- go$children(t)
  df.Annotationlevel1 <- tibble()
  for(x in 1:length(c)){
    
    genes <- unique(GOA$gene_id[GOA$go_id%in%go$descendants(c[x])])
    genes <- unique(c(genes, GOA$gene_id[GOA$go_id%in%c[x]]))
    row <- c(c[x],go$name(c[x]), length(genes))
    df.Annotationlevel1 <- rbind(df.Annotationlevel1, row)
    
  }
  colnames(df.Annotationlevel1) = c("GOID","GOName","value")
  df.Annotationlevel1$value = as.numeric(df.Annotationlevel1$value)
  df.Annotationlevel1<-df.Annotationlevel1 %>% arrange(desc(value))
  df.Annotationlevel1$GOName = factor(df.Annotationlevel1$GOName,levels = df.Annotationlevel1$GOName)
  df.Annotationlevel1$GOID = factor(df.Annotationlevel1$GOID,levels = df.Annotationlevel1$GOID)
  
  df.Annotationlevel1 <- df.Annotationlevel1 %>% filter(value!=0)
  
  genes<-c()
  cumulative <- c()
  for(x in 1:length(df.Annotationlevel1$GOID)){
    genes <- unique(c(genes,unique(GOA$gene_id[GOA$go_id%in%go$descendants(as.character(df.Annotationlevel1$GOID[x]))])))
    genes <- unique(c(genes,unique(GOA$gene_id[GOA$go_id%in%as.character(df.Annotationlevel1$GOID[x])])))
    cumulative <- c(cumulative,length(genes))
  }
  
  df.Annotationlevel1 <- df.Annotationlevel1 %>% mutate("Cumulative" = cumulative)
  
  p <- ggplot(data=df.Annotationlevel1, aes(x=GOName, y=value)) +
    geom_bar(stat="identity")+
    geom_point(aes(y=Cumulative))+
    geom_line(aes(y=Cumulative,group=1))+
    ylab("Number of genes")+
    xlab("GO terms - level 1")+
    ggtitle(go$name(t))+
    theme_classic()+ theme(axis.text.x = element_text(angle = 45,hjust=1))+
    annotate("segment", x = 0, xend =which(df.Annotationlevel1$Cumulative==length(genes))[1], y = length(genes), yend = length(genes), colour = "red",linetype="longdash")+
    annotate("segment", x = which(df.Annotationlevel1$Cumulative==length(genes))[1], xend =which(df.Annotationlevel1$Cumulative==length(genes))[1], y = 0, yend = length(genes), colour = "red",linetype="longdash")+
    annotate("text", x = which(df.Annotationlevel1$Cumulative==length(genes))[1]/2, y = length(genes)+500, colour = "red",label=paste("max. annotated genes", length(genes)))
  return(p)
}

#' # Load data
data(go)
url = "ftp://beta.plantgenie.org/Data/ConGenIE/Picea_abies/v1.0/Annotation/"
h = new_handle(dirlistonly=TRUE)
con = curl(url, "r", h)
tbl = read.table(con, stringsAsFactors=TRUE, fill=TRUE)
names(tbl) = c("Avail.")
close(con)
pander(tbl)
#' ## Load the GOA annotation
urlFile <- paste(url,"Pabies1.0_go_ids_description_name_space.tsv.gz",sep="")
GOA <- read_tsv(urlFile,col_names = TRUE,col_types=cols(
  col_character(),
  col_character(),
  col_character(),
  col_character()
))
#' ## Load the gene description (id,description)
urlFile <- paste(url,"Pabies1.0-b2g_ID-Desc.tsv.gz",sep="")
gene_des <- read_tsv(urlFile,col_names = TRUE,col_types=cols(
  col_character(),
  col_character()
)) %>%
  mutate(SeqName, gene_id=str_replace(SeqName,"\\.[123456789]",""))%>%
  dplyr::select("gene_id","Description")
#' # How many genes are annotated and well annotated?
#' The next bar chart represents the genes with ("known") and without annotation (unknown).
#' Additionally, considering the paper of Tomczak et al. 2018 (ScienceReport journal), a well
#' characterized genes should have more than 10 GO annotations. 
wellknown <- length(names(table(GOA$gene_id)[(table(GOA$gene_id)>10)==TRUE]))
known <- length(unique(GOA$gene_id)) - wellknown
unknown <- length(unique(gene_des$gene_id))-(known+wellknown)

df <- data.frame("category"=factor(c("unknown","annotated","well-annotated"),levels = c("unknown","annotated","well-annotated")), "value"=c(unknown,known,wellknown))
bp<- ggplot(df, aes(x="", y=value, fill=category))+
  geom_bar(width = 1, stat = "identity")+
  geom_text(aes(y=value, label=percent((value/sum(value)))), vjust=1.6, 
            color="black", size=3.5)
bar <- bp +  theme_classic()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) +scale_fill_brewer("category")+
  theme(axis.text.x=element_blank()) 
bar
#' # Gene Ontology analysis
#' ## Information content evaluation
#' The information content (IC) represents a quantified value represented how specific/informative a 
#' Gene Ontology (GO) term is. There are many ways to compute that scores, be using the gene annotated or
#' considering the GO graph structure. For this analysis, the IC based on the GO graph is more suitable since
#' the annotation based IC are really organism dependent, that require a high percentage of well characterize 
#' genes. 
#' First, the IC considering the GO graph is computed for every GO term. Then a IC distribution is used to extract the
#' different percentiles: 1% ,5%, 10%, 25%, 50%, 75%. 

GOINFO <- go$all_ic()
perc <- quantile(as.numeric(unlist(GOINFO$ics)),probs=c(0.01,0.05,0.1,0.25,0.50,0.75))

#' Second, the IC distribution is used as filter for the gene annotation. The idea is to check what is the
#' the impact of the specificity in the annotation. For example, if all the annotation GO terms are involved in the 1% or 5%, 
#' the annotation provide very general terms, otherwise, if the annotation GO terms are localized in the 50% or 75% means
#' that the gene provided really specific terms.
P0 <- unique(GOA$go_id)
P1 <- unique(GOA$go_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[1]])
P5 <- unique(GOA$go_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[2]])
P10<- unique(GOA$go_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[3]])
P25<- unique(GOA$go_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[4]])
P50<- unique(GOA$go_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[5]])
P75<- unique(GOA$go_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[6]])

df.annotation = tibble("P0"=length(P0),"P1"=length(P1),
                       "P5"=length(P5),"P10"=length(P10),
                       "P25"=length(P25),"P50"=length(P50),
                       "P75"=length(P75))

df.annotation <- tibble("Probability"=colnames(df.annotation),"size" = t(df.annotation)[,1])
df.annotation$Probability <- factor(x = df.annotation$Probability, levels=c("P0","P1","P5","P10","P25","P50","P75"))

P0 <- unique(GOA$gene_id)
P1 <- unique(GOA$gene_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[1]])
P5 <- unique(GOA$gene_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[2]])
P10<- unique(GOA$gene_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[3]])
P25<- unique(GOA$gene_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[4]])
P50<- unique(GOA$gene_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[5]])
P75<- unique(GOA$gene_id[GOINFO[GOINFO$goid%in%unique(GOA$go_id)==TRUE,2]>perc[6]])

df.annotatedgenes = tibble("P0"=length(P0),"P1"=length(P1),
                           "P5"=length(P5),"P10"=length(P10),
                           "P25"=length(P25),"P50"=length(P50),
                           "P75"=length(P75))
df.annotatedgenes <- tibble("Probability"=colnames(df.annotatedgenes),"size" = t(df.annotatedgenes)[,1])
df.annotatedgenes$Probability <- factor(x = df.annotatedgenes$Probability, levels=c("P0","P1","P5","P10","P25","P50","P75"))
#' ## Visualization impact of IC for the number of gene and the number of annotation 
p<-ggplot(data=df.annotation, aes(x=Probability, y=size, group=1)) +
  geom_line()+
  geom_point()+
  ylim(c(0,4000))+
  ylab("Number of GO term annotating genes")+
  xlab("IC probability")+
  theme_classic()


t<-ggplot(data=df.annotatedgenes, aes(x=Probability, y=size, group=1)) +
  geom_line()+
  geom_point()+
  #ylab("number of GO term annotating genes")+
  # Custom the Y scales:
  scale_y_continuous(
    
    # Features of the first axis
    name = "Number of annotated genes",
    limits=c(0,15000),
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~./length(gene_des$gene_id), name="Percentage of annotated genes")
  )+
  xlab("IC probability")+
  theme_classic()


n <- 2000
df22 <- data.frame(x = n:1, val = seq(0, 0.5, length.out = n), type = 1)
m<-ggplot(df22, aes(x = x, y = val)) +
  geom_ribbon(aes(ymax = val, ymin = 0, group = type)) +
  annotate("text", x=250, y=-0.03, label= "Including all GO terms") + 
  annotate("text", x=1700, y=-0.03, label="Including specific GO terms") + 
  geom_col(aes(fg=val)) +
  theme_classic()+ guides(fill = guide_legend(override.aes = list(color = NA)), 
                          color = FALSE, 
                          shape = FALSE)  +theme(axis.title.x=element_blank(),
                                                 axis.text.x=element_blank(),
                                                 axis.ticks.x=element_blank(),
                                                 axis.title.y=element_blank(),
                                                 axis.text.y=element_blank(),
                                                 axis.ticks.y=element_blank(),
                                                 axis.line = element_blank(),
                                                 legend.text = element_blank(),
                                                 legend.title = element_blank(),
                                                 legend.background  = element_blank())+
  theme(plot.margin = unit(c(1,1,1,2), "cm"))


grid.arrange(p,m, heights = c(2, 1.5),
             ncol = 1, nrow = 2,
             top=textGrob("GO terms number based on IC probability",gp=gpar(fontsize=20,font=3)))

grid.arrange(t,m, heights = c(2, 1.5),
             ncol = 1, nrow = 2,
             top=textGrob("Annotated gene number based on IC probability",gp=gpar(fontsize=20,font=3)))
#' ## Annotation summary
#' The idea here is to categorize the GO term that annotated gene into the first GO term level, that means
#' the children of the root label ("biological process", "cellular component" and "molecular function"). 
tops <- c("GO:0008150","GO:0003674","GO:0005575")
root <- tibble()
toplist = list()

for(x in 1:length(tops)){
  l = list()
  name <- go$name(tops[x])
  godesc <- go$descendants(tops[x])
  godescAn <- unique(GOA$go_id[GOA$go_id%in%godesc])
  genesAn <-  length(unique(GOA$gene_id[GOA$go_id%in%godesc]))
  l$name = name
  l$Andesc = godescAn
  toplist[[tops[x]]] = l
  root<-rbind(root, c(name,length(godesc),length( godescAn),genesAn))
}
colnames(root) = c("GOName","GO_desc","GO_desc_in_GOA","Annotated_genes")
pander(root)
for(x in 1:length(tops)){
  print(level1Visu(go, GOA, tops[x]))
}
#' # Network evaluation
#' The idea here is to charge the seidr network result and apply semantic similarity to evaluate the
#' whole network and after apply clustering, how similar are the gene cluster. Thus, one can have a quick
#' overview that what happens with the annotation and potentially see the limitation by using GO for the 
#' given organism. 
#'  

gene_cluster <- read_tsv("/mnt/picea/home/aabenitez/Git/seidrSpruce/data/seidr/infomap-bb-p1/clusterResolve.txt",
                         col_names = c("cluster","level1","level2","level3","level4","level5","score","indx","gene_id"),
                         col_types = cols(
                           col_character(),
                           col_integer(),
                           col_integer(),
                           col_integer(),
                           col_integer(),
                           col_integer(),
                           col_double(),
                           col_integer(),
                           col_character()
                         ))

levels <- colnames(gene_cluster)[grepl("[lL][Ee][Vv][Ee][Ll]*",colnames(gene_cluster))]

gn <- as.numeric(table(gene_cluster[levels[1]]))
total <- sum(gn)
cgn = cumsum(gn)
tgn <- tibble("geneNumber"=as.integer(gn), "cumsum(geneNumber)"=as.integer(cgn))
tgn$idx = as.integer(rownames(tgn))
mtgn <- melt(tgn,c("idx"))
sc = 5
ggplot(data=mtgn[mtgn$idx%in%1:100,], aes(x =idx, y = value,group=variable,color=variable))+
  geom_line()+
  geom_point()+
  geom_vline(xintercept=sc)+
  annotation_custom(grob = textGrob(paste("the first",sc,"clusters involve the",round(cgn[sc]/total*100,1),"% of genes")),  xmin =sc, xmax = 100, ymin = 0, ymax = total) +
  ylab("number of genes")+
  xlab("Cluster")+
  ggtitle("Representation of the gene number for the first 100 clusters")+
  theme_classic()+
  theme(axis.text.x = element_blank())

#' # Level 2
annotation <- GOA[,1:2]
names <- annotation %>%
  group_by(gene_id) %>% group_keys()
names <- names$gene_id
annotation<- annotation %>%
  group_by(gene_id) %>%
  group_map(~ head(.x$go_id))
names(annotation) <- names

infoDf <- tibble()
for(x in 1:5){
  gc1 <- gene_cluster[gene_cluster$level1==x,]
  lev <- unique(gc1$level2)
  for(y in 1:length(lev)){
    
  gc11 <- gc1[gc1$level2==lev[y],]
  gs <- gc11$gene_id
  
  
  annotationGS <- GOA[GOA$gene_id%in%gs,]
  GOterms <- unique(annotationGS$go_id)
  bp.GOterms <- GOterms[GOterms%in%toplist[["GO:0008150"]]$Andesc]
  mf.GOterms <- GOterms[GOterms%in%toplist[["GO:0003674"]]$Andesc]
  cc.GOterms <- GOterms[GOterms%in%toplist[["GO:0005575"]]$Andesc]
  
  annotatedGenes <- unique(annotationGS$gene_id)
  tab <- as.numeric(table(annotationGS$gene_id))
  mean.tab <- mean(tab)
  sd.tab <- sd(tab)
  median.tab <- median(tab)
  gs2 <- 0
  if(length(gs)>1){
  gs2 <-round(ontologyReader::GS2(gs,annotation,go),3)
  }
  gs2.onlyannotated <- 0
  if(gs2>0.){
  gs2.onlyannotated <-round(ontologyReader::GS2(annotatedGenes,annotation,go),3)
  }
  row <- c(x,y,length(gs),length(bp.GOterms),
           length(mf.GOterms),length(cc.GOterms),
           round(length(annotatedGenes)/length(gs)*100,2),mean.tab,median.tab,sd.tab,gs2,gs2.onlyannotated)
  infoDf<-rbind(infoDf,row)
  }
}
colnames(infoDf) <- c("l1","l2","Gene_number","BP_terms","MF_terms","CC_terms",
                      "AnnotatedGenePer","meanAnnotationPerGene",
                      "median_annotation_per_gene","sd","gs2","gs2_only_AG")
infoDf <- as_tibble(infoDf)


p<-ggplot(infoDf, aes(x=l1, y=AnnotatedGenePer,group=l1)) + 
  geom_violin(trim=TRUE)+
  ylim(0,100)+
  geom_hline(yintercept=50, color='red',linetype = 2) +
 # geom_point(shape=16, position=position_jitter(0.005))+
  stat_summary(fun.data=mean_sdl, 
               geom="point", color="red")+
  ylab("Percentage of annotated genes")+
  xlab("Clusters")+
  theme_classic()
p



for(x in 1:5){
infoDf.x = infoDf[infoDf$l1==x,]
print(ggplot(infoDf.x,
       aes(x = AnnotatedGenePer, y = gs2)) +
  geom_point(aes(size=Gene_number)) + 
  ylim(0,1)+
  scale_fill_continuous(type = "viridis") +
  geom_hline(yintercept=0.5, color='red',linetype = 2) +
  geom_vline(xintercept=50, color='red',linetype = 2) +
  scale_size(name = "Gene number",breaks = c(25,seq(50,1000,by=200),1000), range = c(1, 15))+
  ylab("GS2")+
  xlab("Percentage of annotated genes")+
  ggtitle(paste("Level 2 clusters - Level 1 cluster:",x))+
  theme_classic())

print(ggplot(infoDf.x,
       aes(x = Gene_number, y = gs2_only_AG )) +
  geom_point(aes(size=AnnotatedGenePer)) + 
  ylim(0,1)+
  scale_fill_continuous(type = "viridis") +
  scale_x_log10()+
  geom_hline(yintercept=0.5, color='red',linetype = 2) +
  scale_size(name = "Annotated Gene %",breaks = c(seq(0,100,by=20)), range = c(1, 15))+
  ylab("GS2")+
  xlab("log10(Gene number)")+
  ggtitle(paste("Level 2 clusters - Level 1 cluster:",x))+
  theme_classic())
}

#' # Conclusion
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' 



