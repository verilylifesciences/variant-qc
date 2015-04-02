# R Codelabs

* [Data Analysis using Google Genomics](./1000Genomes-BRCA1-analysis)
* [Quality Control using Google Genomics](./PlatinumGenomes-QC)

## Prerequisites

### Required
1. Ensure that [BigQuery is enabled](https://console.developers.google.com/flows/enableapi?apiid=bigquery)
   for your Google Cloud Platform project in the Google Developers Console.

1. Ensure that [GoogleGenomics is enabled](https://console.developers.google.com/flows/enableapi?apiid=genomics)
   for your Google Cloud Platform project in the Google Developers Console.

1. Install the `bigrquery` package
  ```
  install.packages("bigrquery")
  ```
  
1. Install the GoogleGenomics R package. 
  ```
  source("http://bioconductor.org/biocLite.R") 
  biocLite(suppressUpdates=TRUE) 
  options(repos=biocinstallRepos())
  devtools::install_github("googlegenomics/api-client-r")
  library(GoogleGenomics)
  ```
  
1. If you have not already done so, follow the Google Genomics [sign up instructions](https://cloud-dot-devsite.googleplex.com/genomics/install-genomics-tools#authenticate) to generate and download a valid ``client_secrets.json`` file.

1. Now use that file to authenticate the GoogleGenomics library:
  ```
  GoogleGenomics::authenticate("/path/to/client_secrets.json")
  ```


### Optional
1. The codelabs load results of previously run Google Cloud Dataflow and Apache Spark jobs.
  * If you would like to run the jobs yourself, please see [Getting Started with Dataflow](https://github.com/googlegenomics/dataflow-examples) and/or [Getting Started with Spark](https://github.com/googlegenomics/spark-examples).

  * To grab job results from Google Cloud Storage and turn them into a TSV file.  For example:
    ```
    gsutil cat gs://<bucket-name>/<output-path>-pca.tsv/part* > pca-results.tsv
    ```

  * To then load the TSV into R.  For example:
    ```
    pcaResults <- read.delim('pca-results.tsv', col.names=c("Sample", "PC1", "PC2"))
    ```

1. pre-install additional R packages referenced by the codelabs:

  ```
  install.packages("dplyr")
  install.packages("ggplot2")
  biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
  biocLite("ggbio")
  ```


## Tips

The [httr](https://github.com/hadley/httr) package is used for the OAuth flow for R packages [bigrquery](https://github.com/hadley/bigrquery) and [GoogleGenomics](https://github.com/googlegenomics/api-client-r).  By default it will run the OAuth flow and store the credentials in the _current working directory_.  You can edit your [.Rprofile](http://www.statmethods.net/interface/customizing.html) to specify an explicit location for it to store and locate cached credentials.
```
options("httr_oauth_cache"="~/.httr-oauth")
```

## Further reading
* [Getting Started with BigQuery](https://github.com/googlegenomics/getting-started-bigquery)
* [Getting Started with RMarkdown and BigQuery](https://github.com/googlegenomics/getting-started-bigquery/tree/master/RMarkdown)
* [Getting Started with the Google Genomics R Client](https://github.com/googlegenomics/api-client-r)
