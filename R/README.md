# R Codelabs

* [Data Analysis using Google Genomics](./1000Genomes-BRCA1-analysis)
* [Quality Control using Google Genomics](./PlatinumGenomes-QC)

## Prerequisites

### Required
(1) Ensure that [BigQuery is enabled](https://console.developers.google.com/flows/enableapi?apiid=bigquery)
   for your Google Cloud Platform project in the Google Developers Console.

(2) Ensure that [GoogleGenomics is enabled](https://console.developers.google.com/flows/enableapi?apiid=genomics)
   for your Google Cloud Platform project in the Google Developers Console.

(3) Install the [bigrquery](https://github.com/rstats-db/bigrquery) package and other R packages referenced by the codelabs.
```
install.packages(c("bigrquery", "devtools", "dplyr", "reshape2", "ggplot2", "scales"))
```

(4) _(Only needed for the [Data Analysis using Google Genomics](./1000Genomes-BRCA1-analysis) codelab)_ Install the [GoogleGenomics](https://github.com/Bioconductor/GoogleGenomics) Bioconductor and other Bioconductor packages referenced by the codelabs.
```
source("http://bioconductor.org/biocLite.R") 
useDevel(TRUE) # The latest version of GoogleGenomics is not in the 3.5 release of Bioconductor.
# If gRPC is not installed.
biocLite("GoogleGenomics", suppressUpdates=TRUE)
# If gRPC is installed.
#biocLite("GoogleGenomics", type="source", suppressUpdates=TRUE)
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


## Tips and troubleshooting.

The [httr](https://github.com/hadley/httr) package is used for the OAuth flow for R packages [bigrquery](https://github.com/hadley/bigrquery) and [GoogleGenomics](https://github.com/Bioconductor/GoogleGenomics).  By default it will run the OAuth flow and store the credentials in the _current working directory_.  You can edit your [.Rprofile](http://www.statmethods.net/interface/customizing.html) to specify an explicit location for it to store and locate cached credentials.
```
options("httr_oauth_cache"="~/.httr-oauth")
```

If you have any trouble with the OAuth flow:

1. At the R prompt: `options(httr_oob_default = TRUE)`
2. Delete file .httr-oauth.
3. Retry the query via bigrquery.

For more help, see https://github.com/rstats-db/bigrquery/issues.

## Further reading
* [Getting Started with BigQuery](https://github.com/googlegenomics/getting-started-bigquery)
* [Getting Started with RMarkdown and BigQuery](https://github.com/googlegenomics/getting-started-bigquery/tree/master/RMarkdown)
* [Getting Started with the Google Genomics R Client](https://github.com/googlegenomics/api-client-r)
