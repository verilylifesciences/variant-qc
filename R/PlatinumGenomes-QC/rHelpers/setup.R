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

require(xtable)
require(RCurl)
require(dplyr)
require(scales)
require(bigrquery)
require(ggplot2)

DisplayAndDispatchQuery <- function(queryUri, project, replacements=list()) {
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
  query_exec(querySql, project)
}

DisplayQueryResults <- function(result, n=6) {
  if(is.null(result)) {
    cat("**None**")
  } else {
    # Suppress the printing of sample identifiers when the data contains allele
    # values.  This is brittle since the query could use a different name for the columns
    # containing the identifier or allele.
    if(2 == length(intersect(c("call_call_set_name", "alternate_bases"), colnames(result)))) {
      result <- mutate(result, call_call_set_name = "not displayed")
    }
    print(xtable(head(result, n=n)), type="html", include.rownames=F)
  }
}
