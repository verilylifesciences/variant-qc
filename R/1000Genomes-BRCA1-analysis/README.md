Data Analysis using Google Genomics
===================================

In this codelab, you will use [Google Genomics](https://cloud.google.com/genomics/), [BigQuery](https://cloud.google.com/bigquery/what-is-bigquery), [Spark](http://spark.apache.org/), and [R](http://www.r-project.org/) to explore [1000 Genomes dataset](https://cloud.google.com/genomics/data/1000-genomes). You will:
* run a principal component analysis (either from scratch or using pre-computed results)
* use BigQuery to explore population variation
* zoom in to specific genome regions, including using the Genomics API to look all the way down to raw reads
* run a GWAS over the BRCA1 variants
* visualize and annotate results using various R packages, including [BioConductor](http://www.bioconductor.org)

### Instructions
1. Read the rendered version of the RMarkdown [AllModalitiesDemo.md](./AllModalitiesDemo.md)

2. Execute the code chunk by chunk in RStudio [AllModalitiesDemo.Rmd](./AllModalitiesDemo.Rmd) or line by line in R.

### Notes

If you prefer, `AllModalitiesDemo.R` can be created from [AllModalitiesDemo.Rmd](./AllModalitiesDemo.Rmd) via
```
require(knitr)
purl("./AllModalitiesDemo.Rmd", documentation=1)
```
