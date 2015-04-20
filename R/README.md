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
  biocLite("GoogleGenomics")
  ```
  
1. If you have not already done so, follow the Google Genomics [sign up instructions](https://cloud.google.com/genomics/install-genomics-tools#authenticate) to generate and download a valid ``client_secrets.json`` file.

1. Now use that file to authenticate the GoogleGenomics library:
  ```
  require(GoogleGenomics)
  GoogleGenomics::authenticate("/PATH/TO/YOUR/client_secrets.json")
  ```

1. If you do not already have them, install additional R packages referenced by the codelabs.
  ```
  install.packages("ggplot2")
  install.packages("dplyr")
  install.packages("scales")
  ```

1. If you do not already have them, install additional Bioconductor R packages referenced by the codelabs _(Note that these are only used in the [Data Analysis using Google Genomics](./1000Genomes-BRCA1-analysis) codelab.)_
```
  biocLite("ggbio", suppressUpdates=TRUE)
  biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene", suppressUpdates=TRUE)
  # This is a large download and takes a little while to install.
  biocLite("BSgenome.Hsapiens.UCSC.hg19", suppressUpdates=TRUE)
```

### Optional

The codelabs load previously computed Principal Coordinate Analysis and Identity-By-State results.  You can run these jobs yourself, giving an opportunity to use Google Cloud Dataflow or Apache Spark.

Instructions are provided to run them over a small portion of the genome, only taking a few minutes, and also how to run them over the whole genome, which can take a few hours depending upon how many machines are running concurrently. For more detail, please see:
   * [Principal Coordinate Analysis](http://googlegenomics.readthedocs.org/en/latest/use_cases/compute_principal_coordinate_analysis/index.html)
   * [Identity-By-State](http://googlegenomics.readthedocs.org/en/latest/use_cases/compute_identity_by_state/index.html)


## Tips

The [httr](https://github.com/hadley/httr) package is used for the OAuth flow for R packages [bigrquery](https://github.com/hadley/bigrquery) and [GoogleGenomics](https://github.com/googlegenomics/api-client-r).  By default it will run the OAuth flow and store the credentials in the _current working directory_.  You can edit your [.Rprofile](http://www.statmethods.net/interface/customizing.html) to specify an explicit location for it to store and locate cached credentials.
```
options("httr_oauth_cache"="~/.httr-oauth")
```

## Further reading
* [Getting Started with BigQuery](https://github.com/googlegenomics/getting-started-bigquery)
* [Getting Started with RMarkdown and BigQuery](https://github.com/googlegenomics/getting-started-bigquery/tree/master/RMarkdown)
* [Getting Started with the Google Genomics R Client](https://github.com/googlegenomics/api-client-r)
