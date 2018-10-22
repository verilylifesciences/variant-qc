# R Reports

RMarkdown code for example quality control reports. The corresponding SQL and tests are in the [sql](../sql) subdirectory.

## To run the reports

(1) Ensure that [BigQuery is
enabled](https://console.cloud.google.com/flows/enableapi?apiid=bigquery) for
your Google Cloud Platform project.

(1) Install the [bigrquery](https://github.com/rstats-db/bigrquery) package and other R packages referenced by the reports.
```
install.packages(c("tidyverse", "bigrquery", "reticulate", "scales"))
```

(1) Install the [Jinja2](http://jinja.pocoo.org/docs/2.10/) Python package. This is used to interpolate the parameters in the [SQL templates](./sql/).
```
pip install jinja2
```

(1) In R, set your working directory to be where the `.Rmd` file you want to run resides.

(1) In R, source any of the dataset-specific `.R` files to run the  [RMarkdown parameterized report](https://bookdown.org/yihui/rmarkdown/parameterized-reports.html) on that particular dataset.

## Tips and troubleshooting.

The [httr](https://github.com/r-lib/httr) package is used for the OAuth flow for R packages [bigrquery](https://github.com/rstats-db/bigrquery). By default it will run the OAuth flow and store the credentials in the _current working directory_.  You can edit your [.Rprofile](http://www.statmethods.net/interface/customizing.html) to specify an explicit location for it to store and locate cached credentials.
```
options("httr_oauth_cache"="~/.httr-oauth")
```

If you have any trouble with the OAuth flow:

1. At the R prompt: `options(httr_oob_default = TRUE)`
2. Delete file `.httr-oauth`.
3. Retry the query via bigrquery.

For more help, see https://github.com/rstats-db/bigrquery/issues.

## Further reading
* [Getting Started with BigQuery](https://github.com/googlegenomics/getting-started-bigquery)
* [Getting Started with RMarkdown and BigQuery](https://github.com/googlegenomics/getting-started-bigquery/tree/master/RMarkdown)
