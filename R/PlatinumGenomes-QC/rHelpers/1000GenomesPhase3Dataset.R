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

# Notice that both tables are the same since 1,000 Genomes already has non-variant
# calls within the variant records.
queryReplacements <- list("_THE_TABLE_"="genomics-public-data:1000_genomes_phase_3.variants",
                          "_THE_EXPANDED_TABLE_"="genomics-public-data:1000_genomes_phase_3.variants",
                          # Chromosomes do not start with the 'chr' prefix in this dataset
                          "chr"="")

sampleData <- read.csv("http://storage.googleapis.com/genomics-public-data/1000-genomes/other/sample_info/sample_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Sample, gender=Gender)

# TODO: Identity By State results for 1,000 Genomes.
ibs <- data.frame(col.names=c("sample1", "sample2", "ibsScore", "similar", "observed"))
