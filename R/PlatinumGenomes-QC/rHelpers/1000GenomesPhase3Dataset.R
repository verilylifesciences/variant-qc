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

kResultsDir = "./1000genomes"
kGenomeCallTableName = "genomics-public-data.1000_genomes_phase_3.variants_20150220_release"
kMultiSampleTableName = kGenomeCallTableName

queryReplacements <- list(
  "@GENOME_CALL_TABLE" = kGenomeCallTableName,
  "@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE" = kGenomeCallTableName,
  "@MULTISAMPLE_VARIANT_TABLE" = kMultiSampleTableName,
  # Work around numeric overflow in R since 1000 genomes has so many calls.
  "COUNT(genotype) AS genotype_count" = "CAST(COUNT(genotype) AS FLOAT64) AS genotype_count",
  #  1000 Genomes does not have depth.  Just hardcode a value.
  "call.DP" = "'not applicable' AS DP",
  # There is no field identifying low quality variant calls because all calls are high quality.
  "EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))" = "FALSE",
  "EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft IN ('PASS', '.'))" = "TRUE",
  # Use precomputed fields to speed up homozygosity query.
  "(SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt WHERE gt >= 0)) FROM v.call) AS AN" = "AN",
  "(SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt WHERE gt = 1)) FROM v.call) AS AC" = "v.AC[ORDINAL(1)] AS AC"
)

sampleData <- read.csv("http://storage.googleapis.com/genomics-public-data/1000-genomes/other/sample_info/sample_info.csv")
sampleInfo <- dplyr::select(sampleData,
                            call_set_name=Sample,
                            sex=Gender,
                            ethnicity=Super_Population)

# Use a specific directory for the figures from this dataset.
knitr::opts_chunk$set(fig.path=file.path(kResultsDir, "figure/"))
# When knitting, DON'T stop if any failures occur.
knitr::opts_chunk$set(error=TRUE)
