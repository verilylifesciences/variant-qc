# Copyright 2015 Google Inc. All rights reserved.
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

library(xtable)
library(RCurl)
library(dplyr)
library(scales)
library(bigrquery)
library(ggplot2)

#' This is a helper method to perform common tasks when dispatching SQL from external files
#' to BigQuery.
#'
#' This method will obtain SQL from a local or remote file, optionally perform text replacements,
#' display it, and send it to BigQuery for execution.
#'
#' @param queryUri File path or URI to the file containing the query.
#' @param project The Google Cloud Platform project within which this query should run.
#' @param replacements A list of key/value pairs to be search/replace on an exact match basis
#'    within the text of the query.
#' @param ... Any additional arguments will be passed to the bigrquery::query_exec method.  This is useful
#'    if you want to pass `destination_table=YOUR-DATASET.YOUR-NEW-TABLE`, `warn=FALSE` or
#'    `max_pages=Inf` to query_exec.
#' @return the dataframe holding the query result or null
#'
#' @examples
#' DisplayAndDispatchQuery('https://raw.githubusercontent.com/googlegenomics/codelabs/master/R/PlatinumGenomes-QC/sql/private-variants.sql',
#'                         project=YOUR-PROJECT_ID,
#'                         replacements=list('_GENOME_CALL_TABLE_'='genomics-public-data:platinum_genomes.variants'))
#'
#' DisplayAndDispatchQuery('https://raw.githubusercontent.com/googlegenomics/codelabs/master/R/PlatinumGenomes-QC/sql/private-variants.sql',
#'                         project=YOUR-PROJECT_ID,
#'                         replacements=list('_GENOME_CALL_TABLE_'='genomics-public-data:platinum_genomes.variants'),
#'                         destination_table='YOUR-DATASET.YOUR-NEW-TABLE')
#'
DisplayAndDispatchQuery <- function(queryUri, project, replacements=list(), ...) {
  if (missing(queryUri)) {
    stop("Pass the file path or url to the file containing the query.")
  }
  if(missing(project)) {
    stop("Pass the project id of your Google Cloud Platform project.")
  }

  if (grepl("^https.*", queryUri)) {
    # Read the query from a remote location.
    querySql <- getURL(queryUri, ssl.verifypeer=FALSE)
  } else {
    # Read the query from the local filesystem.
    querySql <- readChar(queryUri, nchars=1e6)
  }

  # If applicable, substitute values in the query template.
  for(replacement in names(replacements)) {
    querySql <- gsub(replacement, replacements[[replacement]], querySql, fixed=TRUE)
  }

  # Display the query to the terminal.
  cat(querySql)

  # Dispatch the query to BigQuery.
  query_exec(querySql, project, useLegacySql=FALSE, ...)
}

#' This is a helper method to display the first few rows of the query results.
#'
#' It prints the results in a format suitable for HTML display such as on github.  Additionally,
#' if the results contain commonly named columns for genomic data, it attempts to suppress the display
#' of callset identifier with the allele value.
#'
#' @param result A dataframe (e.g., results from bigrquery)
#' @param n The number of results to display.
#' @return an HTML table of results.
#'
DisplayQueryResults <- function(result, n=6) {
  # This accommodates the empty result handling of both recent and old versions of bigrquery.
  if (is.null(result) || nrow(result) == 0) {
    cat("**None**")
  } else {
    # Suppress the printing of sample identifiers when the data contains allele
    # values.  This is brittle since the query could use a different name for the columns
    # containing the identifier or allele.
    if(2 <= length(intersect(c("call_call_set_name", "call_set_name",
                               "alternate_bases", "alt.alternate_bases", "alt_concat"),
                             colnames(result)))) {
      if("call_call_set_name" %in% colnames(result)) {
        result <- mutate(result,
                         call_call_set_name = "not displayed")
      } else {
        result <- mutate(result,
                         call_set_name = "not displayed")
      }
    }
    print(xtable(head(result, n=n)), type="html", include.rownames=F)
  }
}

#' Returns a data frame with a y position and a label, for use annotating ggplot boxplots.
#'
#' @param d A data frame.
#' @return A data frame with column y as max and column label as length.
#'
get_boxplot_fun_data <- function(d) {
  return(data.frame(y=max(d), label=paste0("N = ", length(d))))
}
