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

rmarkdown::render(
  "Data-Overview.Rmd",
  params = list(
    PROJECT_ID = "YOUR-PROJECT-ID",
    DATASET_NAME = "1000 Genomes phase 3",
    DATASET_DESCRIPTION = "1000 Genomes phase 3 https://cloud.google.com/genomics/docs/public-datasets/1000-genomes",
    # Zero-based b37 coordinates per https://www.ncbi.nlm.nih.gov/gene/672
    BRCA1_START = 41196311,
    BRCA1_END = 41277500,
    GENOME_CALL_TABLE = "bigquery-public-data.human_genome_variants.1000_genomes_phase_3_variants_20150220",
    GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE = "bigquery-public-data.human_genome_variants.1000_genomes_phase_3_optimized_schema_variants_20150220"
    ),
  output_file = "Data-Overview-1000-Genomes.html")


rmarkdown::render(
  "Sample-Level-QC.Rmd",
  params = list(
    PROJECT_ID = "YOUR-PROJECT-ID",
    DATASET_NAME = "1000 Genomes phase 3",
    DATASET_DESCRIPTION = "1000 Genomes phase 3 https://cloud.google.com/genomics/docs/public-datasets/1000-genomes",
    # Zero-based b37 coordinates per https://www.ncbi.nlm.nih.gov/grc/human
    PAR1_START = 60000,
    PAR1_END = 2699520,
    PAR2_START = 154931043,
    PAR2_END = 155260560,
    GENOME_CALL_TABLE = "bigquery-public-data.human_genome_variants.1000_genomes_phase_3_variants_20150220",
    GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE = "bigquery-public-data.human_genome_variants.1000_genomes_phase_3_optimized_schema_variants_20150220",
    MULTISAMPLE_VARIANT_TABLE = "bigquery-public-data.human_genome_variants.1000_genomes_phase_3_optimized_schema_variants_20150220",
    MULTISAMPLE_IS_OPTIMIZED = TRUE,
    HIGH_QUALITY_CALLS_FILTER = "TRUE",
    LOW_QUALITY_CALLS_FILTER = "FALSE",
    SAMPLE_INFORMATION_QUERY = "
SELECT
  Sample AS name,
  Gender AS sex,
  Super_Population AS ancestry
FROM
  `bigquery-public-data.human_genome_variants.1000_genomes_sample_info`
"
    ),
  output_file = "Sample-Level-QC-1000-Genomes.html")

rmarkdown::render(
  "Variant-Level-QC.Rmd",
  params = list(
    PROJECT_ID = "YOUR-PROJECT-ID",
    DATASET_NAME = "1000 Genomes phase 3",
    DATASET_DESCRIPTION = "1000 Genomes phase 3 https://cloud.google.com/genomics/docs/public-datasets/1000-genomes",
    GENOME_CALL_TABLE = "bigquery-public-data.human_genome_variants.1000_genomes_phase_3_variants_20150220",
    GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE = "bigquery-public-data.human_genome_variants.1000_genomes_phase_3_optimized_schema_variants_20150220",
    MULTISAMPLE_VARIANT_TABLE = "bigquery-public-data.human_genome_variants.1000_genomes_phase_3_optimized_schema_variants_20150220",
    MULTISAMPLE_IS_OPTIMIZED = TRUE,
    HIGH_QUALITY_VARIANTS_FILTER = "TRUE",
    HIGH_QUALITY_CALLS_FILTER = "TRUE",
    MALE_SAMPLES_QUERY = "
SELECT
  Sample AS name
FROM
  `bigquery-public-data.human_genome_variants.1000_genomes_sample_info`
WHERE
  Gender = 'male'
"
    ),
  output_file = "Variant-Level-QC-1000-Genomes.html")
