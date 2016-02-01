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

# These tables contain the Complete Genomics PGP genomes.
queryReplacements <- list("_GENOME_CALL_TABLE_"="google.com:biggene:pgp_20150205.genome_calls",
                          "_MULTISAMPLE_VARIANT_TABLE_"="google.com:biggene:pgp_20150205.multisample_variants",
                          "DP"="totalReadCount" # for Ti/Tv by depth
                          )

sampleData <- read.csv(textConnection(getURL("https://my.pgp-hms.org/google_surveys/1/download")))
sampleInfo <- select(sampleData, call_call_set_name=Participant, sex=Sex.Gender)

# Google Genomics variant set id for dataset 'pgp_20150205'
variantSetId <- 9170389916365079788

# Show only a subset of PGP ids in the Identity By State results to keep the
# heatmap to a reasonable size.
sampleIds <- c("hu040C0A",
               "hu04DF3C",
               "hu04F220",
               "hu050E9C",
               "hu05FD49",
               "hu34D5B9-1",
               "hu34D5B9-2",
               "hu089792",
               "huA49E22",
               "huFFAD87",
               "huC14AE1",
               "hu627574",
               "hu5FA322",
               "huE58004",
               "hu2843C9",
               "hu7B594C",
               "hu67EBB3",
               "hu7852C5",
               "hu2D6140",
               "huA02824",
               "hu82E689",
               "huD52556")
# TODO: IBS
ibs <- data.frame(col.names=c("sample1", "sample2", "ibsScore", "similar", "observed"))
#ibs <- filter(ibs, sample1 %in% sampleIds & sample2 %in% sampleIds)

# TODO: 2-way PCA against 1,000 Genomes
pca <- data.frame(col.names=c("call_call_set_name", "PC1", "PC2", "count"))

