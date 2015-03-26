Data Analysis using Google Genomics
===================================

In this codelab, you will use [Google Genomics](https://cloud.google.com/genomics/), [BigQuery](https://cloud.google.com/bigquery/what-is-bigquery), [Spark](http://spark.apache.org/), and [R](http://www.r-project.org/) to explore the [1000 Genomes dataset](https://cloud.google.com/genomics/data/1000-genomes). Specifically, you will:
* run a principal component analysis (either from scratch or using pre-computed results)
* use BigQuery to explore population variation
* zoom in to specific genome regions, including using the Genomics API to look all the way down to raw reads
* run a GWAS over the variants within BRCA1
* visualize and annotate results using various R packages, including [BioConductor](http://www.bioconductor.org)

### Instructions
You can watch, read or execute the analysis:

* **Watch** the brief video [Google Genomics: Data Analysis Overview](https://www.youtube.com/watch?v=vINpqxhcTt0) or watch the extended video [Google Genomics Codelab: Data Analysis in R](https://www.youtube.com/watch?v=tPH5PwjzhBM)
* **Read** the rendered version of the RMarkdown [AllModalitiesDemo.md](./AllModalitiesDemo.md)
* **Execute** the code *chunk by chunk* in RStudio or *line by line* in R.
 1. Clone this github repository.
 2. Make sure you have completed the necessary [prerequisites](../README.md#prerequisites).
 3. Run the RMarkdown [AllModalitiesDemo.Rmd](./AllModalitiesDemo.Rmd) file.
 4. If you prefer, `AllModalitiesDemo.R` can be created from [AllModalitiesDemo.Rmd](./AllModalitiesDemo.Rmd) via
```
require(knitr)
purl("./AllModalitiesDemo.Rmd", documentation=2)
```
### Troubleshooting

* `cannot open file './data/1kg-pca.tsv': No such file or directory`
 * The codelabs assume that the current working is the directory in which Rmd file resides.
 * `setwd("path/to/codelabs/R/1000Genomes-BRCA1-analysis")`

* `Daily Limit for Unauthenticated Use Exceeded. Continued use requires signup.`
 * This error occurs when the GoogleGenomics package has not been authenticated properly. Make sure you have followed the [prerequisites](../README.md#required) and have called the `GoogleGenomics::authenticate` method.

* `oauth_listener() needs an interactive environment`
 * Run this demo *interactively* the first time around so that httr can invoke the OAuth flow and store credentials in the current working directory.

* `Access Denied: Job genomics-public-data ...`
 * There's a spot in the codelab that says `# put your projectID here` and you'll want to use your Google Cloud Platform project id as the value of that variable.
 
* `Error: Quota exceeded: Your project exceeded quota for concurrent query bytes scanned.`
 * This codelab requires billing to be enabled on the BigQuery project used. Skip all calls to `DisplayAndDispatchQuery` if turning on billing isn't an option.
