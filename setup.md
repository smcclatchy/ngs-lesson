---
layout: page
title: Setup
permalink: /setup/
---
## Installation

This lesson assumes you have the R, RStudio software installed on your computer.

R can be downloaded [here](https://cran.r-project.org/mirrors.html).

RStudio is an environment for developing using R.
It can be downloaded [here](https://www.rstudio.com/products/rstudio/download/).
You will need the Desktop version for your computer.

This workshop also requires package installation from [CRAN](https://cran.r-project.org) 
and from [Bioconductor](http://www.bioconductor.org/).

### CRAN

To install required packages from CRAN, run this code in the Console.

~~~
install.packages(pkgs = c("ggplot2", "plyr", "dplyr", "tibble", "ape", "flashClust", "WGCNA"))
~~~
{: .r}

Alternatively, you can use the Packages tab in RStudio.
Select the Packages tab and click the Install button.
Type this comma-separated list of package names into 
the pop-up window.

ggplot2, plyr, dplyr, tibble, ape, flashClust, WGCNA

### Bioconductor

The `biocLite.R` script installs Bioconductor packages.
Run the following in the Console to source this script 
and run it to install required Bioconductor packages.

~~~
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(pkgs = c("DESeq2", "org.Mm.eg.db", "genefilter", "topGO", "DT", "biomaRt", "limma"))
~~~
{: .r}

## Data files and project organization

1. Make a new folder in your Desktop called `rna-seq`. Move into this new folder.

2. Create  a `data` folder to hold the data, a `scripts` folder to house your scripts, and a `results` folder to hold results. 

Alternatively, you can use the R console to run the following commands for steps 1 and 2.

~~~
setwd("~/Desktop")
dir.create("./rna-seq")
setwd("~/Desktop/rna-seq")
dir.create("./data")
dir.create("./scripts")
dir.create("./results")
~~~
{: .r}

3. Please download the following large files **before the workshop**, and place them in your `data` folder.

[liver data for eQTL, pQTL, and mediation analyses](ftp://ftp.jax.org/scm/ChickMungeretal2016_DiversityOutbred.Rdata)
[raw count data](ftp://ftp.jax.org/scm/DO192_RNAseq_EMASE_RawCounts.Rdata)
[addiction data](ftp://ftp.jax.org/dgatti/AddictionCourse2017/DO_striatum_addiction2017_0912.Rdata)

