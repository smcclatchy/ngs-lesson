---
title: "Narayanan's Differential Gene Expression Analysis"
teaching: 0
exercises: 0
questions:
- "How do I find genes differentially expressed between treatment groups?"
objectives:
- "To check for sample mixup in data."
- "To identify genes differentially expressed between treatment groups."
- "To describe the effect of sample size on differential gene expression analysis."
keypoints:
- "Non-coding RNA can be used to check for sample mixups in RNA sequencing data."
source: Rmd
---



One of the most common applications of RNA sequencing technology is to identify genes that are differentially expressed between sample groups, for example, between wild type and mutant, or between tumor and normal samples. Count data report the number of sequence reads (fragments) assigned to each gene, which describes the expression abundance of a gene. Similar data can be found in ChIP-Seq, HiC, shRNA screening, or mass spectrometry experiments.

![](../fig/RNAseq-workflow.png)

Once we have aligned sequence reads and have quantified expression, we can continue the pipeline with differential expression analysis. We will use read counts at the gene level and the R package [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html), among other packages. 
  
We will use quantified total liver gene expression data from 192 male and female Diversity Outbred (DO) mice [Chick, J.M., et al. (2016) *Nature* 534(7608):500-505.](https://www.nature.com/nature/journal/v534/n7608/abs/nature18270.html) Half of the animals were fed a standard rodent chow diet, and the other half fed a high-fat diet.

R Libraries and Data Import
------------------------------------
#### Load packages
Let us load the R packages and the data needed for the differential expression analysis.


~~~
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(limma)
~~~
{: .r}



~~~
Warning: package 'limma' was built under R version 3.4.2
~~~
{: .error}



~~~

Attaching package: 'limma'
~~~
{: .output}



~~~
The following object is masked from 'package:DESeq2':

    plotMA
~~~
{: .output}



~~~
The following object is masked from 'package:BiocGenerics':

    plotMA
~~~
{: .output}

#### Load gene information

Load the data file containing basic gene information used in the analysis.


~~~
gene_info=read.csv(url("ftp://ftp.jax.org/dgatti/ShortCourse2015/ENSMUSG-gene-info-R84.tsv"),header=FALSE,sep="\t")
colnames(gene_info)=c("gene_id","gene_name","chr","strand","start","end")
head(gene_info)
~~~
{: .r}


#### Load R Data files

Load the R data files containing expression data and experimental design information needed for doing  differential expression analysis.


~~~
load("data/ChickMungeretal2016_DiversityOutbred.Rdata")
~~~
{: .r}



~~~
Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
'data/ChickMungeretal2016_DiversityOutbred.Rdata', probable reason 'No such
file or directory'
~~~
{: .error}



~~~
Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
~~~
{: .error}



~~~
load("data/DO192_RNAseq_EMASE_RawCounts.Rdata")
~~~
{: .r}



~~~
Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
'data/DO192_RNAseq_EMASE_RawCounts.Rdata', probable reason 'No such file or
directory'
~~~
{: .error}



~~~
Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
~~~
{: .error}

#### Data munging

Let us munge the data :-)


~~~
exp_all= t(expr.rna.192.rawcounts)
~~~
{: .r}



~~~
Error in t(expr.rna.192.rawcounts): object 'expr.rna.192.rawcounts' not found
~~~
{: .error}



~~~
exp_design=data.frame(Sample_ID=rownames(covariates.rna.192),covariates.rna.192,stringsAsFactors = FALSE)
~~~
{: .r}



~~~
Error in rownames(covariates.rna.192): object 'covariates.rna.192' not found
~~~
{: .error}



~~~
head(exp_design)
~~~
{: .r}



~~~
Error in head(exp_design): object 'exp_design' not found
~~~
{: .error}



~~~
exp_all[1:5,1:5]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_all' not found
~~~
{: .error}



~~~
exp_design[1:5,]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_design' not found
~~~
{: .error}


~~~
exp_design=covariates.rna.192  %>% 
        select(Diet,Sex) %>%   rownames_to_column('Sample_ID')
~~~
{: .r}



~~~
Error in eval(lhs, parent, parent): object 'covariates.rna.192' not found
~~~
{: .error}



~~~
### print few rows of the exp_design 
head(exp_design)
~~~
{: .r}



~~~
Error in head(exp_design): object 'exp_design' not found
~~~
{: .error}


Let us check the samples in the expression data and design data are in the same order.


~~~
all(colnames(exp_all)==exp_design$sample_ID)
~~~
{: .r}



~~~
Error in colnames(exp_all): object 'exp_all' not found
~~~
{: .error}

> ## Challenge 1 Familiarize yourself with the data
> 1). Find the number of samples in the data set.  
> 2). Find the number of genes in the study.  
>
> > ## Solution to Challenge 1
> > 1). 
> > 2).   
> {: .solution}
{: .challenge}


A quick check for sample mixup
------------------------------
Let us do a quick sample mixup check using **Xist** gene expression. Xist is non-coding RNA expressed in females.

Let us plot **Xist** expression in all samples against sex.


~~~
plot_exp <- function(exp_design, gexp, g_id, g_info, variable="Sex"){
      # plots gene expression (raw) counts by Sex variable
      # Arguments:
      #    exp_design: experimental design data frame containing
      #                sample IDs, Diet and Sex information           
      #    gexp:       gene expression data  
      #    g_id:       ensembl gene ID
      #    variable:   variable to plot
      # Output:
      #     gene expression plot 
      #
      if (g_id %in% rownames(gexp)){
        g_ind = which(as.vector(g_info$gene_id)==g_id)
        g_name = as.vector(g_info$gene_name)[g_ind]
        chro = as.vector(g_info$chr)[g_ind]
        g_index = which(rownames(gexp)==g_id)
        exp_data= data.frame(exp_design, 
                     exp=as.numeric(gexp[g_index,]))
        if (variable=="Sex"){
            p <- ggplot(exp_data,aes(x=Sex,y=exp)) 
            p <- p + geom_point(position = position_jitter(width = 0.2),size=3,
                    aes(colour = Sex))
        }else{
            p <- ggplot(exp_data,aes(x=Diet,y=exp)) 
            p <- p + geom_point(position = position_jitter(width = 0.2),size=3,
                    aes(colour = Diet))
        }
        p <- p + stat_summary(fun.y=mean, geom="point", shape=5, size=4)
        p <- p + ylab("Gene Expression (Read Counts)")
        p <- p + theme(axis.text=element_text(size=12),
                axis.title=element_text(size=12,face="bold", 
                                        colour = "blue"),
                plot.title = element_text(size = rel(1.5)))
        p <- p + ggtitle(paste0(g_id,": ",g_name," Chr",chro))
        print(p)
      }else{
        print(paste0(gene_id, " not expressed"))
      }
}
~~~
{: .r}



~~~
### Xist ensembl ID
gene_id='ENSMUSG00000086503'
### plot Xist expression by Sex using the function
plot_exp(exp_design, exp_all, gene_id,gene_info)
~~~
{: .r}



~~~
Error in rownames(gexp): object 'exp_all' not found
~~~
{: .error}



~~~
plot_exp(exp_design, exp_all, gene_id,gene_info,variable="Diet")
~~~
{: .r}



~~~
Error in rownames(gexp): object 'exp_all' not found
~~~
{: .error}



~~~
gene_id='ENSMUSG00000025907'
plot_exp(exp_design, exp_all, gene_id, gene_info)
~~~
{: .r}



~~~
Error in rownames(gexp): object 'exp_all' not found
~~~
{: .error}

> ## Challenge 2 Plot your favorite gene
> Pick your favorite gene (by ensembl ID) and plot its expression by:  
> 1). sex.    
> 2). diet.    
>
> > ## Solution to Challenge 2
> > 1). 
> > 2).   
> {: .solution}
{: .challenge}
  
Differential Expression Analysis with **three** samples in each group
------------------------------------------------------------------

Let us start with an example identifying the genes that are differentially expressed between two diets **Chow** and **High fat**.

Let us first get the sample IDs (mouse IDs).


~~~
head(exp_design)
~~~
{: .r}



~~~
Error in head(exp_design): object 'exp_design' not found
~~~
{: .error}



~~~
exp_design
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_design' not found
~~~
{: .error}



~~~
male_chow_ids= exp_design %>% filter(Sex=='M' & Diet=='chow') %>% pull(Sample_ID) 
~~~
{: .r}



~~~
Error in eval(lhs, parent, parent): object 'exp_design' not found
~~~
{: .error}



~~~
male_chow_ids
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'male_chow_ids' not found
~~~
{: .error}



~~~
male_hf_ids = exp_design %>% filter(Sex=='M' & Diet=='HF') %>% pull(Sample_ID)
~~~
{: .r}



~~~
Error in eval(lhs, parent, parent): object 'exp_design' not found
~~~
{: .error}

To make the example simple, let us subset our expression data such that we have **3 Male DO mice** under **Chow diet** and 3 DO mice under **High Fat diet**.

~~~
sampleSize = 3
~~~
{: .r}
Later on we will see the effect of sample size by changing it to **10**.

~~~
diet_3 = c(male_chow_ids[1:sampleSize],male_hf_ids[1:sampleSize])
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'male_chow_ids' not found
~~~
{: .error}



~~~
exp_design_diet_3 = exp_design %>% filter(Sample_ID %in% diet_3)
~~~
{: .r}



~~~
Error in eval(lhs, parent, parent): object 'exp_design' not found
~~~
{: .error}



~~~
exp_diet_3=exp_all[,diet_3]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_all' not found
~~~
{: .error}



~~~
all(colnames(exp_diet_3)==as.vector(exp_design_diet_3$Sample_ID))
~~~
{: .r}



~~~
Error in colnames(exp_diet_3): object 'exp_diet_3' not found
~~~
{: .error}

~~~
as.data.frame(head(exp_diet_3))
~~~
{: .r}



~~~
Error in head(exp_diet_3): object 'exp_diet_3' not found
~~~
{: .error}



~~~
dim(exp_diet_3)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_diet_3' not found
~~~
{: .error}

Let us filter out genes with zero and low expression.


~~~
threshold = 200
exp_mat_diet_3= as.data.frame(exp_diet_3) %>%
    rownames_to_column('gene_id') %>%
    filter(rowSums(.[,2:7], na.rm=TRUE)>threshold) %>%
    column_to_rownames('gene_id')
~~~
{: .r}



~~~
Error in as.data.frame(exp_diet_3): object 'exp_diet_3' not found
~~~
{: .error}



~~~
dim(exp_mat_diet_3)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_mat_diet_3' not found
~~~
{: .error}



~~~
head(exp_mat_diet_3)
~~~
{: .r}



~~~
Error in head(exp_mat_diet_3): object 'exp_mat_diet_3' not found
~~~
{: .error}
## Differential expression analysis with DESeq2

Let us create data frames for **DESeq2** object 


~~~
### colData contains the condition/group information for Differenetial expression analysis
colData <- DataFrame(group = factor(exp_design_diet_3$Diet))
~~~
{: .r}



~~~
Error in factor(exp_design_diet_3$Diet): object 'exp_design_diet_3' not found
~~~
{: .error}

~~~
### Create DESeq2 object using expression and colData
dds_3reps <- DESeqDataSetFromMatrix(countData = as.data.frame(round(exp_mat_diet_3)),
         colData = colData, design = ~ group)
~~~
{: .r}



~~~
Error in as.data.frame(round(exp_mat_diet_3)): object 'exp_mat_diet_3' not found
~~~
{: .error}



~~~
dds_3reps <- DESeq(dds_3reps)
~~~
{: .r}



~~~
Error in is(object, "DESeqDataSet"): object 'dds_3reps' not found
~~~
{: .error}



~~~
res_3reps = results(dds_3reps)
~~~
{: .r}



~~~
Error in mcols(object): object 'dds_3reps' not found
~~~
{: .error}



~~~
resOrdered_3reps <- res_3reps[order(res_3reps$padj),]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'res_3reps' not found
~~~
{: .error}



~~~
head(resOrdered_3reps)
~~~
{: .r}



~~~
Error in head(resOrdered_3reps): object 'resOrdered_3reps' not found
~~~
{: .error}

DE Analysis summary with **3** samples per group 
-------------------------------------------------


~~~
### summary of Differential Expression analysis
summary(res_3reps)
~~~
{: .r}



~~~
Error in summary(res_3reps): object 'res_3reps' not found
~~~
{: .error}



~~~
sig_genes_3reps = as.data.frame(res_3reps) %>% 
                  rownames_to_column('gene_id') %>%
                  filter(padj<0.1) %>% pull(gene_id)
~~~
{: .r}



~~~
Error in as.data.frame(res_3reps): object 'res_3reps' not found
~~~
{: .error}



~~~
length(sig_genes_3reps)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'sig_genes_3reps' not found
~~~
{: .error}

### P-value histogram

Let us plot the histogram of p-values. The p-value histogram is a good diagnostic test for the differential expression analysis.


~~~
hist(res_3reps$pvalue,breaks=100,col="grey",ylim=c(0,800), xlab="p-value",main="p-value histogram: 3 Samples per group")
~~~
{: .r}



~~~
Error in hist(res_3reps$pvalue, breaks = 100, col = "grey", ylim = c(0, : object 'res_3reps' not found
~~~
{: .error}

### **Top** differentially expressed genes

~~~
### helper functions to get gene information for a gene
### plot DESEQ2 object
deseq_gene_exp_plot <- function(deseq_obj, g_id, g_info){
      # plots gene expression counts by group variable
      # used in DESEQ2 object
      #
      # Arguments:
      #    deseq_obj: DESEQ2 object
      #    g_id:      ensembl gene ID
      #    g_info:    gene information dataframe
      #
      # Output:
      #     gene expression plot 
      #
      g_ind = which(as.vector(g_info$gene_id)==g_id)
      g_name = as.vector(g_info$gene_name)[g_ind]
      chro = as.vector(g_info$chr)[g_ind]
      data <- plotCounts(deseq_obj, gene=g_id, intgroup=c("group"), returnData=TRUE)
      p <- ggplot(data, aes(x=group, y=count, color=group))
      p <- p+ ggtitle(paste0(g_id,": ",g_name," Chr",chro))
      p <- p+ geom_point(position=position_jitter(width=.1,height=0), size=3)
      p <- p + theme(axis.text=element_text(size=12),              axis.title=element_text(size=20,face="bold", colour = "blue"), 
            plot.title = element_text(size=rel(1.5)))
      print(p)
}
~~~
{: .r}


~~~
#par(mfrow=c(2,3),las=1)
n=3
top_genes= rownames(resOrdered_3reps[1:n,])
~~~
{: .r}



~~~
Error in rownames(resOrdered_3reps[1:n, ]): object 'resOrdered_3reps' not found
~~~
{: .error}



~~~
for (i in 1:length(top_genes)){
  g_id = top_genes[i]
  deseq_gene_exp_plot(dds_3reps, g_id, gene_info)
}
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'top_genes' not found
~~~
{: .error}

Differential Expression Analysis with **ten** samples in each **diet** group
------------------------------------------------------------------


~~~
sampleSize=10
diet_10 = c(male_chow_ids[1:sampleSize],male_hf_ids[1:sampleSize])
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'male_chow_ids' not found
~~~
{: .error}



~~~
exp_design_diet_10 = exp_design %>% filter(Sample_ID %in% diet_10)
~~~
{: .r}



~~~
Error in eval(lhs, parent, parent): object 'exp_design' not found
~~~
{: .error}



~~~
head(exp_design_diet_10)
~~~
{: .r}



~~~
Error in head(exp_design_diet_10): object 'exp_design_diet_10' not found
~~~
{: .error}



~~~
exp_diet_10=exp_all[,diet_10]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_all' not found
~~~
{: .error}



~~~
dim(exp_diet_10)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_diet_10' not found
~~~
{: .error}



~~~
head(exp_diet_10)
~~~
{: .r}



~~~
Error in head(exp_diet_10): object 'exp_diet_10' not found
~~~
{: .error}



~~~
all(colnames(exp_diet_10)==as.vector(exp_design_diet_10$Sample_ID))
~~~
{: .r}



~~~
Error in colnames(exp_diet_10): object 'exp_diet_10' not found
~~~
{: .error}

Let us filter out genes with zero and low expression (less than 5 read counts) in 50% of the samples.


~~~
threshold = 2000
head(exp_diet_10)
~~~
{: .r}



~~~
Error in head(exp_diet_10): object 'exp_diet_10' not found
~~~
{: .error}



~~~
dim(exp_diet_10)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_diet_10' not found
~~~
{: .error}



~~~
exp_mat_diet_10= as.data.frame(exp_diet_10) %>%
    rownames_to_column('gene_id') %>%
    filter(rowSums(.[,2:ncol(exp_diet_10)+1], na.rm=TRUE)>threshold) %>%
    column_to_rownames('gene_id')
~~~
{: .r}



~~~
Error in as.data.frame(exp_diet_10): object 'exp_diet_10' not found
~~~
{: .error}



~~~
dim(exp_mat_diet_10)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'exp_mat_diet_10' not found
~~~
{: .error}



~~~
head(exp_mat_diet_10)
~~~
{: .r}



~~~
Error in head(exp_mat_diet_10): object 'exp_mat_diet_10' not found
~~~
{: .error}
Let us create data frames for DESeq2 object

~~~
### colData contains the condition/group information for Differenetial expression analysis
colData <- DataFrame(group = factor(exp_design_diet_10$Diet))
~~~
{: .r}



~~~
Error in factor(exp_design_diet_10$Diet): object 'exp_design_diet_10' not found
~~~
{: .error}

~~~
### Create DESeq2 object using expression and colData
dds_10reps <- DESeqDataSetFromMatrix(countData = as.data.frame(round(exp_mat_diet_10)),
         colData = colData, design = ~ group)
~~~
{: .r}



~~~
Error in as.data.frame(round(exp_mat_diet_10)): object 'exp_mat_diet_10' not found
~~~
{: .error}



~~~
dds_10reps <- estimateSizeFactors(dds_10reps)
~~~
{: .r}



~~~
Error in estimateSizeFactors(dds_10reps): object 'dds_10reps' not found
~~~
{: .error}



~~~
counts_10reps= counts(dds_10reps, normalized=TRUE)
~~~
{: .r}



~~~
Error in counts(dds_10reps, normalized = TRUE): object 'dds_10reps' not found
~~~
{: .error}



~~~
dds_10reps <- DESeq(dds_10reps)
~~~
{: .r}



~~~
Error in is(object, "DESeqDataSet"): object 'dds_10reps' not found
~~~
{: .error}



~~~
res_10reps = results(dds_10reps)
~~~
{: .r}



~~~
Error in mcols(object): object 'dds_10reps' not found
~~~
{: .error}



~~~
resOrdered_10reps <- res_10reps[order(res_10reps$padj),]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'res_10reps' not found
~~~
{: .error}



~~~
head(resOrdered_10reps)
~~~
{: .r}



~~~
Error in head(resOrdered_10reps): object 'resOrdered_10reps' not found
~~~
{: .error}


~~~
n=3
top_genes= rownames(resOrdered_10reps[1:n,])
~~~
{: .r}



~~~
Error in rownames(resOrdered_10reps[1:n, ]): object 'resOrdered_10reps' not found
~~~
{: .error}



~~~
par(mfrow=c(2,3),las=1)
for (i in 1:length(top_genes)){
  g_id = top_genes[i]
  deseq_gene_exp_plot(dds_10reps, g_id, gene_info)
}
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'top_genes' not found
~~~
{: .error}

Differential Gene Expression Analysis Summary  
---------------------------------------------

### P-value histogram comparison

~~~
par(mfrow=c(1,2))
hist(res_10reps$pvalue,breaks=100,col="grey", xlab="p-value",ylim=c(0,1000), main="p-value histogram: 10 Samples per group")
~~~
{: .r}



~~~
Error in hist(res_10reps$pvalue, breaks = 100, col = "grey", xlab = "p-value", : object 'res_10reps' not found
~~~
{: .error}



~~~
hist(res_3reps$pvalue,breaks=100,ylim=c(0,1000),col="grey", xlab="p-value",main="p-value histogram: 3 Samples per group")
~~~
{: .r}



~~~
Error in hist(res_3reps$pvalue, breaks = 100, ylim = c(0, 1000), col = "grey", : object 'res_3reps' not found
~~~
{: .error}



~~~
#rld_10reps <- rlog(dds_10reps, blind = FALSE)
#plotPCA(rld_10reps, intgroup = c("Diet"))
~~~
{: .r}

### DESEQ2 Summary: 3 samples per group

~~~
### summary of Differential Expression analysis
summary(res_3reps)
~~~
{: .r}



~~~
Error in summary(res_3reps): object 'res_3reps' not found
~~~
{: .error}

### DESEQ2 Summary: 10 samples per group

~~~
### summary of Differential Expression analysis
summary(res_10reps)
~~~
{: .r}



~~~
Error in summary(res_10reps): object 'res_10reps' not found
~~~
{: .error}



~~~
sig_genes_10reps = as.data.frame(res_10reps) %>% 
                  rownames_to_column('gene_id') %>%
                  filter(padj<0.1) %>% pull(gene_id)
~~~
{: .r}



~~~
Error in as.data.frame(res_10reps): object 'res_10reps' not found
~~~
{: .error}



~~~
length(sig_genes_10reps)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'sig_genes_10reps' not found
~~~
{: .error}



~~~
length(union(sig_genes_10reps, sig_genes_3reps))
~~~
{: .r}



~~~
Error in union(sig_genes_10reps, sig_genes_3reps): object 'sig_genes_10reps' not found
~~~
{: .error}



~~~
length(intersect(sig_genes_10reps, sig_genes_3reps))
~~~
{: .r}



~~~
Error in intersect(sig_genes_10reps, sig_genes_3reps): object 'sig_genes_10reps' not found
~~~
{: .error}



~~~
# Combining the two above..
comb <- unique(c(sig_genes_10reps, sig_genes_3reps))
~~~
{: .r}



~~~
Error in unique(c(sig_genes_10reps, sig_genes_3reps)): object 'sig_genes_10reps' not found
~~~
{: .error}



~~~
length(comb)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'comb' not found
~~~
{: .error}



~~~
# Comparing comb with the above two
sig_genes_10reps_2 <- comb %in% sig_genes_10reps
~~~
{: .r}



~~~
Error in comb %in% sig_genes_10reps: object 'comb' not found
~~~
{: .error}



~~~
sig_genes_3reps_2 <- comb %in% sig_genes_3reps
~~~
{: .r}



~~~
Error in comb %in% sig_genes_3reps: object 'comb' not found
~~~
{: .error}



~~~
# Generating venn counts to plot venn diagram
sig_genes <- cbind(sig_genes_10reps_2, sig_genes_3reps_2)
~~~
{: .r}



~~~
Error in eval(quote(list(...)), env): object 'sig_genes_10reps_2' not found
~~~
{: .error}



~~~
head(sig_genes)
~~~
{: .r}



~~~
Error in head(sig_genes): object 'sig_genes' not found
~~~
{: .error}



~~~
sig_genes_counts <- vennCounts(sig_genes)
~~~
{: .r}



~~~
Error in as.matrix(x): object 'sig_genes' not found
~~~
{: .error}



~~~
vennDiagram(sig_genes_counts, cex = 1,names = c("10 reps","3 reps"), circle.col = c("red", "blue"))
~~~
{: .r}



~~~
Error in is(object, "VennCounts"): object 'sig_genes_counts' not found
~~~
{: .error}
