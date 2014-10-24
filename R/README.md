# R Codelabs

## Prerequisites

* [Getting Started with BigQuery](https://github.com/googlegenomics/getting-started-bigquery)
 * [Getting Started with RMarkdown and BigQuery](https://github.com/googlegenomics/getting-started-bigquery/tree/master/RMarkdown)
* [Getting Started with the Google Genomics R Client](https://github.com/googlegenomics/api-client-r)
  * TEMPORARY CAVEAT `devtools::install_github("deflaux/api-client-r")` instead of `devtools::install_github("googlegenomics/api-client-r")`
* [Getting Started with Spark](https://github.com/googlegenomics/spark-examples)

## Tips
* These codelabs don't enumerate all the R packages needed.  Just install them on a case-by-case basis as needed.
  * If you are new to BioConductor, see [Installing BioConductor packages](http://www.bioconductor.org/install/)
* For convenient authentication for the R Client, in your .Rprofile
```
setHook(packageEvent("GoogleGenomics", "attach"), function(...) {
  GoogleGenomics::authenticate(file="YOUR/PATH/TO/client_secrets.json")
})
```
