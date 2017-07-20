---
title: "Differential Gene Expression Analysis"
teaching: 0
exercises: 0
questions:
- "Key question"
objectives:
- "To check for sample mixup in data."
- "To identify genes differentially expressed between treatment groups."
- "To describe the effect of sample size on differential gene expression analysis."
keypoints:
- "Non-coding RNA can be used to check for sample mixups in RNA sequencing data."
---

Once we have aligned sequence reads and have quantified expression abundances, we can continue the pipeline with differential expression analysis. One of the most common applications of RNA sequencing technology is to identify genes that are differentially expressed between sample groups, for example, between wild type and mutant, or between tumor and normal samples. Count data report the number of sequence reads (fragments) assigned to each gene, which describes the expression abundance of a gene. Similar data can be found in ChIP-Seq, HiC, shRNA screening, or mass spectrometry experiments.

We will use read counts at the gene level and the R package [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html). 

R Libraries and Data Import
------------------------------------
Load the R packages and libraries.

~~~
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
library(ggplot2)
library(dplyr)
~~~
{: .r}

Load the data file from the URL. See the documentation for *gzcon* for an explanation of uncompressing a zipped (.gz) file through a connection, in this case, a URL.

~~~
z <- gzcon(con = url("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72759&format=file&file=GSE72759_DO192_RNAseq_UpperQuartileNormalized_n21454genes_forGEOSubmission.txt.gz")) # uncompress the file
raw <- textConnection(readLines(z)) # prepare for the read.table function
data <- read.table(raw, header = TRUE)
dim(data)
names(data)
data[1:5, 1:5] # a quick look at the data
~~~
{: .r}


~~~
                   X    F326   F327    F328    F329
1 ENSMUSG00000090025  0.0145 0.0000  0.0667  0.0149
2 ENSMUSG00000051951  0.0000 0.0000  0.0267  0.0000
3 ENSMUSG00000025902  0.4928 0.2930  0.5067  0.3134
4 ENSMUSG00000098104  0.0145 0.0000  0.0133  0.0000
5 ENSMUSG00000033845 11.3478 8.2895 10.5467 11.0448
~~~
{: .output}


Check the expression data and make it a matrix.

~~~
geneIDs = data[,1]
data=data[,-1]
rownames(data)=geneIDs
data[1:5,1:5]
~~~
{: .r}

Let us create data frame containing key experimental design factors for the experiment.

~~~
exp_design = data.frame(mouseIDs=colnames(exp.all),
                        diet=covariates.rna.192$Diet,
                        sex=covariates.rna.192$Sex,
                        coat_color=covariates.rna.192$Coat.Color)
~~~
{: .r}



~~~
Error in is.data.frame(x): object 'exp.all' not found
~~~
{: .error}


~~~
all(colnames(exp.all)==exp_design$mouseIDs)
~~~
{: .r}



~~~
Error in is.data.frame(x): object 'exp.all' not found
~~~
{: .error}

A quick check for sample mixup
------------------------------
Let us do a quick sample mixup check using **Xist** gene expression. Xist is non-coding RNA 

~~~
geneID="ENSMUSG00000086503"
geneName="Xist"
gIndex = which(rownames(exp.all)==geneID)
~~~
{: .r}



~~~
Error in rownames(exp.all): object 'exp.all' not found
~~~
{: .error}



~~~
data= data.frame(exp_design, 
                 exp=as.numeric(exp.all[gIndex,]))
~~~
{: .r}



~~~
Error in data.frame(exp_design, exp = as.numeric(exp.all[gIndex, ])): object 'exp_design' not found
~~~
{: .error}

~~~
head(data)
~~~
{: .r}



~~~
                                                                     
1 function (..., list = character(), package = NULL, lib.loc = NULL, 
2     verbose = getOption("verbose"), envir = .GlobalEnv)            
3 {                                                                  
4     fileExt <- function(x) {                                       
5         db <- grepl("\\\\.[^.]+\\\\.(gz|bz2|xz)$", x)              
6         ans <- sub(".*\\\\.", "", x)                               
~~~
{: .output}

Let us plot **Xist** expression in all samples against sex.

~~~
p <- ggplot(data,aes(x=sex,y=exp)) 
~~~
{: .r}



~~~
Error in ggplot(data, aes(x = sex, y = exp)): could not find function "ggplot"
~~~
{: .error}



~~~
p <- p + geom_point(position = position_jitter(width = 0.1),size=3,
                    aes(colour = factor(sex)))
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}



~~~
p <- p + stat_summary(fun.y=mean, geom="point", shape=5, size=4)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}



~~~
p <- p + ylab("Gene Expression (Read Counts)")
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}



~~~
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=20,face="bold", colour = "blue"),
               plot.title = element_text(size = rel(2)))
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}



~~~
p <- p + ggtitle("Xist: ENSMUSG00000086503")
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}



~~~
p
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}


Let us start with an example identifyingthe genes that are differentially expressed in the samples between two diets **Chow** and **High fat**.


~~~
male_index = which(exp_design$sex=="M")
~~~
{: .r}



~~~
Error in which(exp_design$sex == "M"): object 'exp_design' not found
~~~
{: .error}



~~~
female_index = which(exp_design$sex=="F")
~~~
{: .r}



~~~
Error in which(exp_design$sex == "F"): object 'exp_design' not found
~~~
{: .error}



~~~
chow_index= which(exp_design$diet=="chow")
~~~
{: .r}



~~~
Error in which(exp_design$diet == "chow"): object 'exp_design' not found
~~~
{: .error}



~~~
hf_index= which(exp_design$diet=="HF")
~~~
{: .r}



~~~
Error in which(exp_design$diet == "HF"): object 'exp_design' not found
~~~
{: .error}



~~~
male_chow = intersect(male_index,chow_index)
~~~
{: .r}



~~~
Error in as.vector(y): object 'chow_index' not found
~~~
{: .error}



~~~
male_hf = intersect(male_index,hf_index)
~~~
{: .r}



~~~
Error in as.vector(y): object 'hf_index' not found
~~~
{: .error}

Differential Expression Analysis with **three** samples in each group
------------------------------------------------------------------
To make the example simple, let us subset our expression data such that we have **3 DO mice** under **Chow diet** and 3 DO mice under **High Fat diet**.

~~~
sampleSize = 3
~~~
{: .r}
Later on we will see the effect of sample size by varying it.

~~~
diet_DE = c(male_chow[1:sampleSize],male_hf[1:sampleSize])
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'male_chow' not found
~~~
{: .error}



~~~
exp_design_diet_DE= exp_design[diet_DE,]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_design' not found
~~~
{: .error}



~~~
exp_design_diet_DE
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_design_diet_DE' not found
~~~
{: .error}



~~~
exp_diet_DE=exp.all[,diet_DE]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp.all' not found
~~~
{: .error}



~~~
all(colnames(exp_diet_DE)==as.vector(exp_design_diet_DE$mouseIDs))
~~~
{: .r}



~~~
Error in is.data.frame(x): object 'exp_diet_DE' not found
~~~
{: .error}

~~~
head(exp_diet_DE)
~~~
{: .r}



~~~
Error in head(exp_diet_DE): object 'exp_diet_DE' not found
~~~
{: .error}

Let us filter out genes with zero and low expression (less than 5 read counts) in 50% of the samples.

~~~
thres= 5
nzIndex= as.vector(which(apply(exp_diet_DE,1,function(x){sum(x>thres)/length(x)})>=0.5))
~~~
{: .r}



~~~
Error in apply(exp_diet_DE, 1, function(x) {: object 'exp_diet_DE' not found
~~~
{: .error}



~~~
head(nzIndex)
~~~
{: .r}



~~~
Error in head(nzIndex): object 'nzIndex' not found
~~~
{: .error}



~~~
exp.dietDE = exp_diet_DE[nzIndex,]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_diet_DE' not found
~~~
{: .error}



~~~
dim(exp.dietDE)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp.dietDE' not found
~~~
{: .error}
Let us create data frames for **DESeq2** object

~~~
### colData contains the condition/group information for Differenetial expression analysis
colData <- DataFrame(group = factor(exp_design_diet_DE$diet))
~~~
{: .r}



~~~
Error in DataFrame(group = factor(exp_design_diet_DE$diet)): could not find function "DataFrame"
~~~
{: .error}

~~~
### Create DESeq2 object using expression and colData
dds <- DESeqDataSetFromMatrix(countData = as.data.frame(round(exp.dietDE)),
         colData = colData, design = ~ group)
~~~
{: .r}



~~~
Error in DESeqDataSetFromMatrix(countData = as.data.frame(round(exp.dietDE)), : could not find function "DESeqDataSetFromMatrix"
~~~
{: .error}



~~~
dds <- DESeq(dds)
~~~
{: .r}



~~~
Error in DESeq(dds): could not find function "DESeq"
~~~
{: .error}



~~~
res = results(dds)
~~~
{: .r}



~~~
Error in results(dds): could not find function "results"
~~~
{: .error}

~~~
### summary of Differential Expression analysis
summary(res)
~~~
{: .r}



~~~
Error in summary(res): object 'res' not found
~~~
{: .error}

~~~
plotMA(res, main="M-A Plot: 3 Samples per group", ylim=c(-2,2))
~~~
{: .r}



~~~
Error in plotMA(res, main = "M-A Plot: 3 Samples per group", ylim = c(-2, : could not find function "plotMA"
~~~
{: .error}

~~~
d<-plotCounts(dds, gene=which.min(res$padj), intgroup="group",
              returnData=TRUE)
~~~
{: .r}



~~~
Error in plotCounts(dds, gene = which.min(res$padj), intgroup = "group", : could not find function "plotCounts"
~~~
{: .error}



~~~
p <- ggplot(d, aes(x=group, y=count)) +
  geom_point(position=position_jitter(w=0.2,h=0),size=3)
~~~
{: .r}



~~~
Error in ggplot(d, aes(x = group, y = count)): could not find function "ggplot"
~~~
{: .error}



~~~
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=20,face="bold", colour = "blue"),
               plot.title = element_text(size = rel(2)))
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}



~~~
p
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}
Let us plot the histogram of p-values. The p-value histogram is a good diagnostic test for the differential expression analysis.


~~~
hist(res$pvalue,breaks=100,col="grey", xlab="p-value",main="p-value histogram: 3 Samples per group")
~~~
{: .r}



~~~
Error in hist(res$pvalue, breaks = 100, col = "grey", xlab = "p-value", : object 'res' not found
~~~
{: .error}

Differential Expression Analysis with **ten** samples in each **diet** group
------------------------------------------------------------------


~~~
sampleSize = 10
~~~
{: .r}
Later on we will see the effect of sample size by varying it.

~~~
diet_DE = c(male_chow[1:sampleSize],male_hf[1:sampleSize])
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'male_chow' not found
~~~
{: .error}



~~~
exp_design_diet_DE= exp_design[diet_DE,]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_design' not found
~~~
{: .error}



~~~
exp_design_diet_DE
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_design_diet_DE' not found
~~~
{: .error}



~~~
exp_diet_DE=exp.all[,diet_DE]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp.all' not found
~~~
{: .error}



~~~
all(colnames(exp_diet_DE)==as.vector(exp_design_diet_DE$mouseIDs))
~~~
{: .r}



~~~
Error in is.data.frame(x): object 'exp_diet_DE' not found
~~~
{: .error}



~~~
head(exp_diet_DE)
~~~
{: .r}



~~~
Error in head(exp_diet_DE): object 'exp_diet_DE' not found
~~~
{: .error}

~~~
head(exp_diet_DE)
~~~
{: .r}



~~~
Error in head(exp_diet_DE): object 'exp_diet_DE' not found
~~~
{: .error}

Let us filter out genes with zero and low expression (less than 5 read counts) in 50% of the samples.

~~~
thres= 5
nzIndex= as.vector(which(apply(exp_diet_DE,1,function(x){sum(x>thres)/length(x)})>=0.5))
~~~
{: .r}



~~~
Error in apply(exp_diet_DE, 1, function(x) {: object 'exp_diet_DE' not found
~~~
{: .error}



~~~
head(nzIndex)
~~~
{: .r}



~~~
Error in head(nzIndex): object 'nzIndex' not found
~~~
{: .error}



~~~
exp.dietDE = exp_diet_DE[nzIndex,]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_diet_DE' not found
~~~
{: .error}



~~~
dim(exp.dietDE)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp.dietDE' not found
~~~
{: .error}
Let us create data frames for DESeq2 object

~~~
### colData contains the condition/group information for Differenetial expression analysis
colData <- DataFrame(group = factor(exp_design_diet_DE$diet))
~~~
{: .r}



~~~
Error in DataFrame(group = factor(exp_design_diet_DE$diet)): could not find function "DataFrame"
~~~
{: .error}

~~~
### Create DESeq2 object using expression and colData
dds <- DESeqDataSetFromMatrix(countData = as.data.frame(round(exp.dietDE)),
         colData = colData, design = ~ group)
~~~
{: .r}



~~~
Error in DESeqDataSetFromMatrix(countData = as.data.frame(round(exp.dietDE)), : could not find function "DESeqDataSetFromMatrix"
~~~
{: .error}



~~~
dds <- DESeq(dds)
~~~
{: .r}



~~~
Error in DESeq(dds): could not find function "DESeq"
~~~
{: .error}



~~~
res = results(dds)
~~~
{: .r}



~~~
Error in results(dds): could not find function "results"
~~~
{: .error}

~~~
### summary of Differential Expression analysis
summary(res)
~~~
{: .r}



~~~
Error in summary(res): object 'res' not found
~~~
{: .error}

~~~
plotMA(res, main="M-A Plot: 10 Samples per group", ylim=c(-2,2))
~~~
{: .r}



~~~
Error in plotMA(res, main = "M-A Plot: 10 Samples per group", ylim = c(-2, : could not find function "plotMA"
~~~
{: .error}

~~~
d<-plotCounts(dds, gene=which.min(res$padj), intgroup="group",
              returnData=TRUE)
~~~
{: .r}



~~~
Error in plotCounts(dds, gene = which.min(res$padj), intgroup = "group", : could not find function "plotCounts"
~~~
{: .error}



~~~
p <- ggplot(d, aes(x=group, y=count)) +
  geom_point(position=position_jitter(w=0.2,h=0),size=3)
~~~
{: .r}



~~~
Error in ggplot(d, aes(x = group, y = count)): could not find function "ggplot"
~~~
{: .error}



~~~
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=20,face="bold", colour = "blue"),
               plot.title = element_text(size = rel(2)))
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}



~~~
p
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'p' not found
~~~
{: .error}

~~~
hist(res$pvalue,breaks=100,col="grey", xlab="p-value",main="p-value Histogram: 10 Samples per group")
~~~
{: .r}



~~~
Error in hist(res$pvalue, breaks = 100, col = "grey", xlab = "p-value", : object 'res' not found
~~~
{: .error}


~~~
svd.obj = svd(apply(exp.dietDE,1,function(x){x-mean(x)}))
~~~
{: .r}



~~~
Error in apply(exp.dietDE, 1, function(x) {: object 'exp.dietDE' not found
~~~
{: .error}



~~~
plot(svd.obj$d^2/sum(svd.obj$d^2),ylab="Percent Variance Explained", main="PC of expression data")
~~~
{: .r}



~~~
Error in plot(svd.obj$d^2/sum(svd.obj$d^2), ylab = "Percent Variance Explained", : object 'svd.obj' not found
~~~
{: .error}

~~~
print(cor(svd.obj$u[,1],as.numeric(as.factor(exp_design_diet_DE$diet))))
~~~
{: .r}



~~~
Error in is.factor(x): object 'exp_design_diet_DE' not found
~~~
{: .error}



~~~
print(cor(svd.obj$u[,2],as.numeric(as.factor(exp_design_diet_DE$diet))))
~~~
{: .r}



~~~
Error in is.factor(x): object 'exp_design_diet_DE' not found
~~~
{: .error}



~~~
print(cor(svd.obj$u[,5],as.numeric(as.factor(exp_design_diet_DE$diet))))
~~~
{: .r}



~~~
Error in is.factor(x): object 'exp_design_diet_DE' not found
~~~
{: .error}



~~~
print(cor(svd.obj$u[,1],
    as.numeric(as.factor(covariates.rna.192$Coat.Color[diet_DE]))))
~~~
{: .r}



~~~
Error in is.factor(x): object 'covariates.rna.192' not found
~~~
{: .error}



~~~
print(cor(as.numeric(as.factor(exp_design_diet_DE$diet)),
                         as.numeric(as.factor(covariates.rna.192$Coat.Color[diet_DE]))))
~~~
{: .r}



~~~
Error in is.factor(x): object 'covariates.rna.192' not found
~~~
{: .error}
