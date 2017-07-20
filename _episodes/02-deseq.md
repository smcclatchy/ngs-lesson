---
title: "Differential Expression Analysis with DESeq2"
teaching: 0
exercises: 0
questions:
- "Key question"
objectives:
- "First objective."
keypoints:
- "First key point."
source: Rmd
---



One of the most common applications of RNA-seq technology is using it for identifying genes that are differentially expressed between sample groups, for example, wild type vs mutant, or cancer vs normal samples. 

In the last section, We saw an example of how we can use **EMASE** powereed by **Kallisto** pseudo-alignment (instead of bowtie alignment) and quantify expression abundances at isoform/gene level for a single animal from Diversity Outbred mouse population. 

Let us assume that we have used the same pipeline and quantified expression abundances for all 192 DO samples.

We will be using read counts at gene level and the software tool **DESeq2** for doing differential expression analysis on a subset of the DO mice. 

R Libraries and Data Import
------------------------------------
Let us load the R packages and the data needed for the differential expression analysis.
```
{r R_package, results="hide"}
library("DESeq2")
library(ggplot2)
library(dplyr)
```

Let us load the R object file containing all the data we need.

~~~
load("/data/RData/DO192_DataforSysGenCourse.rdata")
~~~
{: .r}



~~~
Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
'/data/RData/DO192_DataforSysGenCourse.rdata', probable reason 'No such
file or directory'
~~~
{: .error}



~~~
Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
~~~
{: .error}



~~~
exp.all = read.table("/data/emase/expected_read_counts_gene_level.txt", header=T)
~~~
{: .r}



~~~
Warning in file(file, "rt"): cannot open file '/data/emase/
expected_read_counts_gene_level.txt': No such file or directory
~~~
{: .error}



~~~
Error in file(file, "rt"): cannot open the connection
~~~
{: .error}

Let us check the expression data and make it a matrix.

~~~
geneIDs = exp.all[,1]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp.all' not found
~~~
{: .error}



~~~
exp.all=exp.all[,-1]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp.all' not found
~~~
{: .error}



~~~
rownames(exp.all)=geneIDs
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'geneIDs' not found
~~~
{: .error}



~~~
exp.all[1:5,1:5]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp.all' not found
~~~
{: .error}

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
