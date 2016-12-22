# Copyright 2016 Google Inc. All rights reserved.
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

kResultsDir = "./platinum_genomes"
kGenomeCallTableName = "genomics-public-data.platinum_genomes.variants"

# For illustrative purposes, we generated the multisample variants table using two
# different schemas.  For more detail, see:
# https://github.com/googlegenomics/codelabs/blob/summarizeRefMatches/Java/PlatinumGenomes-variant-transformation/README.md
kDenseMultiSampleTableName = "google.com:biggene.platinum_genomes.multisample_variants_dense_schema"
kOptimizedMultiSampleTableName = "google.com:biggene.platinum_genomes.multisample_variants_optimized_schema"

# This boolean controls which table we will use for the QC.
kMultiSampleTableSchemaIsOptimized = FALSE

if (kMultiSampleTableSchemaIsOptimized) {
  kMultiSampleTableName = kOptimizedMultiSampleTableName
} else {
  kMultiSampleTableName = kDenseMultiSampleTableName
}

queryReplacements <- list(
  "@GENOME_CALL_TABLE" = kGenomeCallTableName,
  "@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE" = kGenomeCallTableName,
  "@MULTISAMPLE_VARIANT_TABLE" = kMultiSampleTableName
)

sampleData <- read.csv("http://storage.googleapis.com/genomics-public-data/platinum-genomes/other/platinum_genomes_sample_info.csv")
sampleInfo <- dplyr::select(sampleData,
                            call_set_name=Catalog.ID,
                            sex=Gender,
                            ethnicity=Description)

ibs <- read.table("./data/platinum-genomes-ibs.tsv",
                  col.names=c("sample1", "sample2", "ibsScore", "similar", "observed"))

# Read in the results of the 2-way PCA over BRCA1.
pca <- read.table("./data/platinum-genomes-X-1kg-brca1-pca.tsv",
                  col.names=c("call_call_set_name", "PC1", "PC2", "count"))

# Use a specific directory for the figures from this dataset.
knitr::opts_chunk$set(fig.path=file.path(kResultsDir, "figure/"))
# When knitting, stop if any failures occur.
knitr::opts_chunk$set(error=FALSE)
