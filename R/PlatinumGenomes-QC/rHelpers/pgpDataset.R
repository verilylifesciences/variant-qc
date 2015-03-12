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

require(RCurl)
require(dplyr)

# These tables contain the Complete Genomics PGP genomes and one variants-only Illumina genome.
queryReplacements <- list("_THE_TABLE_"="google.com:biggene:pgp_20150205.variants",
                          "_THE_EXPANDED_TABLE_"="google.com:biggene:pgp_20150205.expanded_variants")

sampleData <- read.csv(textConnection(getURL("https://my.pgp-hms.org/google_surveys/1/download")))
sampleInfo <- select(sampleData, call_call_set_name=Participant, gender=Sex.Gender)

# Show only a subset of PGP ids in the Identity By State results to keep the
# heatmap to a reasonable size.
sampleIds <- c("hu03E3D2",
               "hu03E3D2-lfr",
               "hu040C0A",
               "hu0486D6",
               "hu0486D6-ilm",
               "hu0486D6-lfr",
               "hu04DF3C",
               "hu04F220",
               "hu050E9C",
               "hu05FD49",
               "hu085B6D",
               "hu085B6D-lfr",
               "hu24A473",
               "hu24A473-lfr",
               "hu34D5B9-1",
               "hu34D5B9-2",
               "hu6A01AF",
               "hu6A01AF-lfr",
               "huC1F1FB",
               "huC1F1FB-lfr",
               "huD2B804",
               "huD2B804-lfr")
ibs <- read.table("http://storage.googleapis.com/genomics-public-data/personal-genome-project/other/personal-genome-project-ibs.tsv",
                  col.names=c("sample1", "sample2", "ibsScore", "similar", "observed"))
ibs <- filter(ibs, sample1 %in% sampleIds & sample2 %in% sampleIds)
