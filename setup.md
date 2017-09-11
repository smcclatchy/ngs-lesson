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
install.packages(c("ggplot2", "plyr", "dplyr"))
~~~
{: .r}

Alternatively, you can use the Packages tab in RStudio.
Select the Packages tab and click the Install button.
Type this comma-separated list of package names into 
the pop-up window.

ggplot2, plyr, dplyr

### Bioconductor

The `biocLite.R` script installs Bioconductor packages.
Run the following in the Console to source this script 
and run it to install required Bioconductor packages.

~~~
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("DESeq2", "org.Mm.eg.db", "genefilter", "topGO", "DT"))
~~~
{: .r}


## Data files and project organization
Please download the following large files before the workshop. 

1. Create a project in RStudio. Refer to 
[this lesson](http://swcarpentry.github.io/r-novice-gapminder/02-project-intro/) 
to do this.

2. Create a data directory (folder) to hold the data, a script 
directory to house your scripts, and a results directory to hold results.
You can do this in the  RStudio Files tab, or use Finder on a Mac, 
or go to the Start menu and select (My) Computer on a Windows machine.

3. Download the following data files into your data directory. 

[Data](ftp://ftp.jax.org/dgatti/AddictionCourse2017/DO_striatum_addiction2017_0912.Rdata)

