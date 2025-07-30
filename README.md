# ordinalGWAS: genome-wide association analysis of ordinal phenotypes
The `ordinalGWAS` package is an R package for genetic association analysis of ordinal phenotypes. It is designed to analyze a batch of SNPs, and can be applied to genome-wide analysis by splitting the data into smaller datasets (e.g. by chromosomes).
Current version is v0.1 published on December 23, 2019.

## Note from author (July 30, 2025)
This R package is not under active development or maintenance and may be outdated. It was NOT optimized for large-scale data processing. In particular, it is NOT recommended for GWAS or other computationally intensive applications. Please consider alternative tools that are actively maintained and better suited for such tasks. For additional questions, please contact the author at feifeitao99@gmail.com

## Install the package
You can use `devtools` to install the package in R.
```
library(devtools)
install_github("feifeitao/ordinalGWAS")
```

## User guide
User guide for `ordinalGWAS` is available [here](https://htmlpreview.github.io/?https://github.com/feifeitao/ordinalGWAS/blob/master/vignettes/ordinalGWAS_vignette.html). The user guide describes the input file format, analysis workflow, output format, and provides a set of mock data and analysis examples.
