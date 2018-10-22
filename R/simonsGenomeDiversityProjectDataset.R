# Copyright 2018 Verily Life Sciences LLC.
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
    DATASET_NAME = "Simons Genome Diversity Project",
    DATASET_DESCRIPTION = "Simons Genome Diversity Project https://cloud.google.com/genomics/docs/public-datasets/simons",
    # Zero-based b37 coordinates per https://www.ncbi.nlm.nih.gov/gene/672
    BRCA1_START = 41196311,
    BRCA1_END = 41277500,
    GENOME_CALL_TABLE = "bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_variants",
    GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE = "bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_variants"
    ),
  output_file = "Data-Overview-Simons-Genome-Diversity-Project.html")

rmarkdown::render(
  "Sample-Level-QC.Rmd",
  params = list(
    PROJECT_ID = "YOUR-PROJECT-ID",
    DATASET_NAME = "Simons Genome Diversity Project",
    DATASET_DESCRIPTION = "Simons Genome Diversity Project https://cloud.google.com/genomics/docs/public-datasets/simons",
    # Zero-based b37 coordinates per https://www.ncbi.nlm.nih.gov/grc/human
    PAR1_START = 60000,
    PAR1_END = 2699520,
    PAR2_START = 154931043,
    PAR2_END = 155260560,
    GENOME_CALL_TABLE = "bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_variants",
    GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE = "bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_variants",
    HIGH_QUALITY_CALLS_FILTER = "c.FL IN ('','2','3','4','5','6','7','8','9')",
    LOW_QUALITY_CALLS_FILTER = "c.FL IN ('0','1')",
    SAMPLE_INFORMATION_QUERY = "
SELECT
  id_from_vcf AS name,
  sex,
  region AS ancestry
FROM
  `bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_attributes`
"
    ),
  output_file = "Sample-Level-QC-Simons-Genome-Diversity-Project.html")

rmarkdown::render(
  "Variant-Level-QC.Rmd",
  params = list(
    PROJECT_ID = "YOUR-PROJECT-ID",
    DATASET_NAME = "Simons Genome Diversity Project",
    DATASET_DESCRIPTION = "Simons Genome Diversity Project https://cloud.google.com/genomics/docs/public-datasets/simons",
    GENOME_CALL_TABLE = "bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_variants",
    GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE = "bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_variants",
    HIGH_QUALITY_VARIANTS_FILTER = "TRUE",
    HIGH_QUALITY_CALLS_FILTER = "c.FL IN ('','2','3','4','5','6','7','8','9')",
    MALE_SAMPLES_QUERY = "
SELECT
  id_from_vcf AS name
FROM
  `bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_attributes`
WHERE
sex = 'male'
"
    ),
  output_file = "Variant-Level-QC-Simons-Genome-Diversity-Project.html")
