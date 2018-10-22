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
    PROJECT_ID = "YOUR-PROJECT-ID"
    ),
  output_file = "Data-Overview-DeepVariant-Platinum-Genomes.html")

rmarkdown::render(
  "Sample-Level-QC.Rmd",
  params = list(
    PROJECT_ID = "YOUR-PROJECT-ID"
    ),
  output_file = "Sample-Level-QC-DeepVariant-Platinum-Genomes.html")

rmarkdown::render(
  "Variant-Level-QC.Rmd",
  params = list(
    PROJECT_ID = "YOUR-PROJECT-ID"
    ),
  output_file = "Variant-Level-QC-DeepVariant-Platinum-Genomes.html")
