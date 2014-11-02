Data Analysis using Google Genomics
===================================

In this codelab, you will use [Google Genomics](https://cloud.google.com/genomics/), [BigQuery](https://cloud.google.com/bigquery/what-is-bigquery), [Spark](http://spark.apache.org/), and [R](http://www.r-project.org/) to explore the [1000 Genomes dataset](https://cloud.google.com/genomics/data/1000-genomes). Specifically, you will:
* run a principal component analysis (either from scratch or using pre-computed results)
* use BigQuery to explore population variation
* zoom in to specific genome regions, including using the Genomics API to look all the way down to raw reads
* run a GWAS over the variants within BRCA1
* visualize and annotate results using various R packages, including [BioConductor](http://www.bioconductor.org)

### Instructions
1. Make sure you have completed the necessary [prerequisites](../README.md).

2. Read the rendered version of the RMarkdown [AllModalitiesDemo.md](./AllModalitiesDemo.md)

3. Execute the code *chunk by chunk* in RStudio [AllModalitiesDemo.Rmd](./AllModalitiesDemo.Rmd) or *line by line* in R.  _Note:_ there are a few additional R packages used by this demo but you'll see the necessary lines to install them as you execute this step by step from the [unrendered code](./AllModalitiesDemo.Rmd).

### Notes

If you prefer, `AllModalitiesDemo.R` can be created from [AllModalitiesDemo.Rmd](./AllModalitiesDemo.Rmd) via
```
require(knitr)
purl("./AllModalitiesDemo.Rmd", documentation=1)
```
### Troubleshooting

* `cannot open file './data/1kg-pca.tsv': No such file or directory`
 * The codelabs assume that the current working is the directory in which Rmd file resides.

* `Daily Limit for Unauthenticated Use Exceeded. Continued use requires signup.`
 * Make sure you can retrieve variants and reads via the GoogleGenomics R package. Please see [Getting Started with the Google Genomics R Client](https://github.com/googlegenomics/api-client-r).

* `oauth_listener() needs an interactive environment`
 * Run this demo *interactively* the first time around so that httr can invoke the OAuth flow and store credentials in the current working directory.

* `Access Denied: Job genomics-public-data ...`
 * There's a spot in the codelab that says `# put your projectID here` and you'll want to use your Google Cloud Platform project id as the value of that variable.
