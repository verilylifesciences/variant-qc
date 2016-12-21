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

kResultsDir = "./pgp"
kGenomeCallTableName = "google.com:biggene.pgp_20150205.genome_calls"
kMultiSampleTableName = "google.com:biggene.pgp_20150205.multisample_variants"
kMultiSampleTableSchemaIsOptimized = FALSE

queryReplacements <- list(
  "@GENOME_CALL_TABLE" = kGenomeCallTableName,
  "@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE" = kGenomeCallTableName,
  "@MULTISAMPLE_VARIANT_TABLE" = kMultiSampleTableName,
  #  Use correct Complete Genomics field for Ti/Tv by depth.
  "DP" = "totalReadCount",
  # Complete Genomics male samples have a second genotype of -1 on the X chromosome.
  "CAST((SELECT LOGICAL_AND(gt > 0) FROM UNNEST(call.genotype) gt) AS INT64) AS hom_AA" =
    "CAST(call.genotype[ORDINAL(1)] > 0 AND (call.genotype[ORDINAL(2)] = -1 OR call.genotype[ORDINAL(2)] = call.genotype[ORDINAL(1)]) AS INT64) AS hom_AA",
  # Use the correct field to identify low quality variant calls from Complete Genomics data.
  "EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))" =
    "(call.allele1VariantQuality != 'VQHIGH' OR IFNULL(call.allele2VariantQuality != 'VQHIGH', FALSE))",
  "EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft IN ('PASS', '.'))" =
    "(call.allele1VariantQuality = 'VQHIGH' AND IFNULL(call.allele2VariantQuality = 'VQHIGH', TRUE))"
)

sampleData <- read.csv(textConnection(getURL("https://my.pgp-hms.org/google_surveys/1/download")))

# Add newlines to long levels.
levels(sampleData$Race.ethnicity) <- gsub("American Indian / Alaska Native,",
                                          "American Indian / Alaska Native,\n",
                                          levels(sampleData$Race.ethnicity))

# Combine the empty string and "No response" into a single level.
levels(sampleData$Race.ethnicity) <- gsub("^$",
                                          "No response",
                                          levels(sampleData$Race.ethnicity))

sampleInfo <- dplyr::select(sampleData,
                            call_set_name=Participant,
                            sex=Sex.Gender,
                            ethnicity=Race.ethnicity)

# Use a specific directory for the figures from this dataset.
knitr::opts_chunk$set(fig.path=file.path(kResultsDir, "figure/"))
# When knitting, DON'T stop if any failures occur.
knitr::opts_chunk$set(error=TRUE)
