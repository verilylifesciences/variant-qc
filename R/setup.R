# Copyright 2015 Google Inc., Verily Life Sciences LLC. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

library(tidyverse)
library(scales)
library(bigrquery)
library(reticulate)

jinja <- import("jinja2")
# If you get an error, in the shell run:
#    pip install jinja2
py <- import_builtins()

#' This is a helper method to perform common tasks when dispatching SQL from external files
#' to BigQuery.
#'
#' This method will obtain SQL from a file, optionally perform text replacements,
#' display it, and send it to BigQuery for execution.
#'
#' There are many ways to facilitate templated queries.  Here we use Python via
#' [reticulate](https://github.com/rstudio/reticulate) and
#' [Jinja2](http://jinja.pocoo.org/docs/2.9/).

#' @param sql A character string containing the SQL to run. Use this or the 'sql_path' parameter.
#' @param sql_path The path to the file containing the SQL. Use this or the 'sql' parameter.
#' @param params A list of key/value pairs to be search/replace on an exact match basis
#'    within the text of the query.
#' @param ... Any additional arguments will be passed to the bigrquery::query_exec method.  This is useful
#'    if you want to pass `destination_table=YOUR-DATASET.YOUR-NEW-TABLE`, `warn=FALSE` or
#'    `max_pages=Inf` to query_exec.
#' @return the dataframe holding the query result
perform_bqquery <- function(sql, sql_path, params, ...) {
  if (missing(sql)) {
    sql <- py$open(sql_path, "r")$read()
  }
  sql <- jinja$Template(sql)$render(params)
  cat(sql)
  result <- try(query_exec(sql, use_legacy_sql = FALSE, project = params$PROJECT_ID, ...))
  if ("try-error" == class(result)) {
    warning(result)
    return(data.frame())  # Empty dataframe.
  }
  return(result)
}

#' This is a helper method to display the first few rows of the query results.
#'
#' @param result A dataframe (e.g., results from bigrquery)
#' @param n The number of results to display.
#' @return an HTML table of results.
#'
DisplayQueryResults <- function(result, n = 6) {
  if (nrow(result) == 0) {
    cat("**None**")
  } else {
    print(knitr::kable(head(result, n = n)))
  }
}

#' Returns a data frame with a y position and a label, for use annotating ggplot boxplots.
#'
#' @param d A data frame.
#' @return A data frame with column y as max and column label as length.
#'
get_boxplot_fun_data <- function(d) {
  return(data.frame(y = max(d), label = stringr::str_c("N = ", length(d))))
}
