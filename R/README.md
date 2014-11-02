# R Codelabs

## Prerequisites

### Required
1. Make sure you can execute a BigQuery query using R package `bigrquery`.  Please see (1) [Getting Started with BigQuery](https://github.com/googlegenomics/getting-started-bigquery) and (2) [Getting Started with RMarkdown and BigQuery](https://github.com/googlegenomics/getting-started-bigquery/tree/master/RMarkdown).
1. Make sure you can retrieve variants and reads via the GoogleGenomics R package. Please see [Getting Started with the Google Genomics R Client](https://github.com/googlegenomics/api-client-r).

### Optional
1. The codelabs load results of previously run Spark jobs.  If you would like to run the Spark jobs yourself, please see [Getting Started with Spark](https://github.com/googlegenomics/spark-examples).

## Tips

The [httr](https://github.com/hadley/httr) package is used for the OAuth flow for R packages [bigrquery](https://github.com/hadley/bigrquery) and [GoogleGenomics](https://github.com/googlegenomics/api-client-r).  By default it will run the OAuth flow and store the credentials in the _current working directory_.  You can edit your [.Rprofile](http://www.statmethods.net/interface/customizing.html) to specify an explicit location for it to store and locate cached credentials.
```
options("httr_oauth_cache"="~/.httr-oauth")
```
