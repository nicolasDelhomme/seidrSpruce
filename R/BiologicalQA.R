#' ---
#' title: "Sortmerna_trimmonatico Biological QA - Norway spruce"
#' author: "Aaron Ayllon Benitez and Nicolas Delhome"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(ggdendro)
  library(here)
  library(hyperSpec)
  library(magrittr)
  library(parallel)
  library(pander)
  library(plotly)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))
source(here("UPSCb-common/src/R/percentile.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Metadata
#' Norway spruce samples information: atlas_sample, seasonal1_sample, seasonal2_sample, flakaliden_samples, diurnal1_sample, and diurnal2_sample
#' ```{r CHANGEME1,eval=FALSE,echo=FALSE}
#' # The csv file should contain the sample information, including the sequencing file name, 
#' # any relevant identifier, and the metadata of importance to the study design
#' # as columns, e.g. the SamplingTime for a time series experiment
#'  ```
atlas<- read_csv(here("doc/atlas_samples.csv"),
                      col_types=cols(
                        col_character(),
                        col_character(),
                        col_factor(),
                        col_integer()
                      )) %>% filter(sample.type == "needles")%>%
  mutate(Description=paste(sample.ID,sample.desc,sep="_")) %>% 
  select(sample.ID,Description)

seasonal <- bind_rows(read_csv(here("doc/seasonal_2011_2012_samples_set1.csv"),
                                       col_types=cols(
                                         .default=col_character(),
                                         Sampling_ID=col_number()
                                       )),
                              read_csv(here("doc/seasonal_2012-2013_samples_set2.csv"),
                                       col_types=cols(
                                         .default=col_character(),
                                         Sampling_ID=col_number()
                                       ))) %>% transmute(sample.ID=Sequencing_ID,Description = paste(Sampling_Date,abs(Sampling_ID),sep="_"))

flakaliden<- suppressMessages(read_csv2(here("doc/flakaliden_samples.csv"),
                               col_types=cols(
                                 .default=col_character()
                               )))%>% mutate(Description=paste(treatment,date,plot,sep = "/")) %>%
  select(SciLifeID, Description) %>%
  rename(sample.ID=SciLifeID) %>% 
  filter(sample.ID != "P1169_167")
diurnal <-  bind_rows(suppressWarnings(read_csv(here("doc/diurnal_samples_set1.csv"), 
                              col_types=cols(
                                .default=col_character()))) %>%
  transmute(sample.ID=SciLifeID,Description=paste(SciLifeID,SamplingTime,Experiment,sep = "_")),
  read_csv(here("doc/diurnal_samples_set2.csv"),
                              col_types=cols(
                                .default=col_character()
                              ))%>% transmute(sample.ID=SciLifeID,
                                              Description =SamplingTime))


 

#' # * Raw data
atlas.filelist <- list.files(here("data/seidr/salmon/atlas"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)
flakaliden.filelist <- list.files(here("data/seidr/salmon/flakaliden"), 
                             recursive = TRUE, 
                             pattern = "quant.sf",
                             full.names = TRUE)
seasonal.filelist <- list.files(here("data/seidr/salmon/seasonal"), 
                             recursive = TRUE, 
                             pattern = "quant.sf",
                             full.names = TRUE)
diurnal.filelist <- list.files(here("data/seidr/salmon/diurnal"), 
                             recursive = TRUE, 
                             pattern = "quant.sf",
                             full.names = TRUE)

#' # Checking step
#' Sanity check to ensure that the data is sorted according to the sample info
#' ```{r CHANGEME3,eval=FALSE,echo=FALSE}
#' # This step is to validate that the salmon files are inthe same order as 
#' # described in the samples object. If not, then they need to be sorted
#' ````
#' ## * atlas dataset
atlas.filelist <- atlas.filelist[match(atlas$sample.ID,sub("_.*","",basename(dirname(atlas.filelist))))]
stopifnot(all(sub("_.*","",basename(dirname(atlas.filelist))) == atlas$sample.ID))

#' ## * flakaliden dataset
flakaliden %<>% arrange(sample.ID) #just in case to arrange 
flakaliden.filelist <- flakaliden.filelist[match(flakaliden$sample.ID,gsub(".*X_|_s.*","",
                            basename(dirname(flakaliden.filelist))))]

#' ## * seasonal dataset
seasonal.filelist <- seasonal.filelist[match(seasonal$sample.ID,
                                             sub("_s.*","",
                                                 basename(dirname(seasonal.filelist))))]

#' ## * diurnal dataset
#' The diurnal dataset is store in a same folder but it contains two differents works: diurnal1 and diurnal2. Moreover, in diurnal1 there is
#' repeated samples files and it is hard to match with the sample metadata. For this reason, the stopifnot function is not used in that case
#' and the matching was manually. 
diurnal <- diurnal[match(gsub(".*CXX_|_d.*|_1_.*","",
                              basename(dirname(diurnal.filelist))),diurnal$sample.ID),]
#' ## * Combining datasets
dataSet <- bind_cols(reduce(list(atlas,flakaliden,seasonal,diurnal),bind_rows),
          tibble(DataSet=factor(rep(c("atlas","flakaliden","seasonal","diurnal"),
                              sapply(list(atlas,flakaliden,seasonal,diurnal),nrow))),
                 Path=c(atlas.filelist,flakaliden.filelist,seasonal.filelist,diurnal.filelist)))

filelist <- dataSet$Path 
names(filelist) <- dataSet$sample.ID

#' Read the expression at the gene level
#' ```{r CHANGEME4,eval=FALSE,echo=FALSE}
#' If the species has only one transcript per gene, replace with the following
#' counts <- suppressMessages(round(tximport(files = filelist, type = "salmon",txOut=TRUE)$counts))
#' ```
counts <- suppressMessages(round(tximport(files = filelist, 
                                  type = "salmon",
                                  txOut=TRUE)$counts))

#' ## Merge technical replicates
#' seasonal and diurnal
counts <- sapply(split.data.frame(t(counts),dataSet$Description),colSums)

dataSet %<>% filter(!duplicated(dataSet$Description))

counts <- counts[,match(dataSet$Description,colnames(counts))]

stopifnot(all(colnames(counts) == dataSet$Description))

colnames(counts) <- dataSet$sample.ID

#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, colouring by DataSet
#' ```{r CHANGEME5,eval=FALSE,echo=FALSE}
#' # In the following most often you need to replace CHANGEME by your
#' # variable of interest, i.e. the metadata represented as column in
#' # your samples object, e.g. SamplingTime
#' ```
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(dataSet)

ggplot(dat,aes(x,y,fill=DataSet)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")+
  theme_classic()

#' The same is done for the individual samples colored by CHANGEME. 
#' ```{r CHANGEME6,eval=FALSE,echo=FALSE}
#' # In the following, the second mutate also needs changing, I kept it 
#' # as an example to illustrate the first line. SampleID would be 
#' # a column in the samples object (the metadata) that uniquely indentify
#' # the samples.
#' # If you have only a single metadata, then remove the second mutate call
#' # If you have more, add them as needed.
#' ```
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(DataSet=dataSet$DataSet[match(ind,dataSet$sample.ID)]) #%>% 
#  mutate(SamplingTime=samples$SamplingTime[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=DataSet)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")+
  theme_classic()

#' One sample is obviously rubbish, remove it
range(colSums(counts))
percentile(colSums(counts))

#sel <- which(colSums(counts)==min(colSums(counts)))
#counts <- counts[,-sel]
#dataSet <- dataSet[-sel,]

#' ## Export
dir.create(here("data/seidr/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/seidr/analysis/salmon/raw-unormalised-gene-expression_data.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#' 
#' The used design was the different kind of dataset associated to samples
#'  
#'  ```{r CHANGEME7,eval=FALSE,echo=FALSE}
#'  # In the following, we provide the expected expression model, based on the study design.
#'  # It is technically irrelevant here, as we are only doing the quality assessment of the data, 
#'  # but it does not harm setting it correctly for the differential expression analyses that may follow.
#'  ```
dds <- suppressMessages(DESeqDataSetFromMatrix(
  countData = counts,
  colData = dataSet,
  design = ~ DataSet))

save(dds,file=here("data/analysis/salmon/dds.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
#pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)

#' The following operation is to get zero for the minimal value in vst.
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model. We used only one variable in the design.
nvar=1

#' An the number of possible combinations
#' ```{r CHANGEME8,eval=FALSE,echo=FALSE}
#' This needs to be adapted to your study design. Add or drop variables as needed.
#' ```
#nlevel=nlevels(dds$MDay) * nlevels(dds$MGenotype)
nlevel=nlevels(dds$DataSet)
#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)+
  theme_classic()

#' focus on 0-20
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component",limits = c(0,20)) + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)+
  theme_classic()


#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    dataSet)

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=DataSet,shape=DataSet,text=Description)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")+
  theme_classic()
 
ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))



#' ### Heatmap
#' 
#' Filter for noise
#' 
conds <- factor(dataSet$DataSet)
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 5

#' * Heatmap of "all" genes
#' 
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

hc <- as.hclust(hm$colDendrogram)
ggdendrogram(hc, rotate=FALSE)
#' ## Conclusion
#' CHANGEME
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
