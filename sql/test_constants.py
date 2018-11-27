# Copyright 2018 Verily Life Sciences LLC.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Test data to be used commonly across all tests."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

class TestConstants(object):
  """A collection of constants for use in query tests.

  VCFs can encode the same data in several ways. Below we try to succinctly
  capture sufficient variety. Each query should be tested against at least one
  of the tables below which each reflect the same set of variants, samples, and
  the genotypes of those samples for each variant.
  """

  # Basic input for a 'genome call table'.
  # This includes a variety of variant types, genotypes, and non-variant
  # segments.
  GENOME_CALL_TABLE = "basic_genome_call_table"
  GENOME_CALL_TABLE_INPUT = """
SELECT * FROM UNNEST([
STRUCT<reference_name STRING,
       start_position INT64,
       end_position INT64,
       reference_bases STRING,
       alternate_bases ARRAY<STRUCT<alt STRING>>,
       call ARRAY<STRUCT<name STRING,
                         genotype ARRAY<INT64>>> >
    -- SNP with ref-homozygous, heterozygous, and alt-homozygous genotypes.
    ('chr1', 123, 124, 'A', ARRAY[STRUCT('T'),
                                  STRUCT('<*>')], ARRAY[STRUCT('sample1', ARRAY[0, 0]),
                                                        STRUCT('sample2', ARRAY[0, 1]),
                                                        STRUCT('sample3', ARRAY[1, 1])]),
    -- SNP with no-calls
    ('chr1', 456, 457, 'C', ARRAY[STRUCT('T'),
                                  STRUCT('<*>')], ARRAY[STRUCT('sample1', ARRAY[-1, 1]),
                                                        STRUCT('sample2', ARRAY[-1, -1])]),
    -- Multiallelic SNP
    ('chr1', 456, 457, 'C', ARRAY[STRUCT('T'),
                                  STRUCT('G'),
                                  STRUCT('<*>')], ARRAY[STRUCT('sample3', ARRAY[1, 2])]),
    -- Insertion
    ('chr1', 456, 457, 'C', ARRAY[STRUCT('CC'),
                                  STRUCT('<*>')], ARRAY[STRUCT('sample4', ARRAY[1, 1])]),
    -- Deletion
    ('chr1', 456, 458, 'CA', ARRAY[STRUCT('C'),
                                   STRUCT('<*>')], ARRAY[STRUCT('sample5', ARRAY[0, 1])]),
    -- Overlapping non-variant segments
    ('chr1', 100, 800, 'C', NULL,                 ARRAY[STRUCT('sample6', ARRAY[0, 0])]),
    ('chr1', 100, 800, 'C', ARRAY[STRUCT('<*>')], ARRAY[STRUCT('sample7', ARRAY[0, 0])])
])"""

  # Basic input for a 'multisample variants table'.
  # This includes a variety of variant types and genotypes.
  MULTISAMPLE_VARIANTS_TABLE = "basic_multisample_variants_table"
  MULTISAMPLE_VARIANTS_TABLE_INPUT = """
SELECT * FROM UNNEST([
STRUCT<reference_name STRING,
       start_position INT64,
       end_position INT64,
       reference_bases STRING,
       alternate_bases ARRAY<STRUCT<alt STRING>>,
       call ARRAY<STRUCT<name STRING,
                         genotype ARRAY<INT64>>> >
    -- SNP with ref-homozygous, heterozygous, and alt-homozygous genotypes.
    ('chr1', 123, 124, 'A', ARRAY[STRUCT('T'),
                                  STRUCT('<*>')], ARRAY[STRUCT('sample1', ARRAY[0, 0]),
                                                        STRUCT('sample2', ARRAY[0, 1]),
                                                        STRUCT('sample3', ARRAY[1, 1]),
                                                        STRUCT('sample4', ARRAY[-1, -1]),
                                                        STRUCT('sample5', ARRAY[-1, -1])]),
    -- SNP with no-calls
    ('chr1', 456, 457, 'C', ARRAY[STRUCT('T'),
                                  STRUCT('<*>')], ARRAY[STRUCT('sample1', ARRAY[-1, 1]),
                                                        STRUCT('sample2', ARRAY[-1, -1]),
                                                        STRUCT('sample3', ARRAY[-1, -1]),
                                                        STRUCT('sample4', ARRAY[-1, -1]),
                                                        STRUCT('sample5', ARRAY[-1, -1])]),
    -- Multiallelic SNP
    ('chr1', 456, 457, 'C', ARRAY[STRUCT('T'),
                                  STRUCT('G'),
                                  STRUCT('<*>')], ARRAY[STRUCT('sample1', ARRAY[-1, -1]),
                                                        STRUCT('sample2', ARRAY[-1, -1]),
                                                        STRUCT('sample3', ARRAY[1, 2]),
                                                        STRUCT('sample4', ARRAY[-1, -1]),
                                                        STRUCT('sample5', ARRAY[-1, -1])]),
    -- Insertion
    ('chr1', 456, 457, 'C', ARRAY[STRUCT('CC'),
                                  STRUCT('<*>')], ARRAY[STRUCT('sample1', ARRAY[-1, -1]),
                                                        STRUCT('sample2', ARRAY[-1, -1]),
                                                        STRUCT('sample3', ARRAY[-1, -1]),
                                                        STRUCT('sample4', ARRAY[1, 1]),
                                                        STRUCT('sample5', ARRAY[-1, -1])]),
    -- Deletion
    ('chr1', 456, 458, 'CA', ARRAY[STRUCT('C'),
                                   STRUCT('<*>')], ARRAY[STRUCT('sample1', ARRAY[-1, -1]),
                                                         STRUCT('sample2', ARRAY[-1, -1]),
                                                         STRUCT('sample3', ARRAY[-1, -1]),
                                                         STRUCT('sample4', ARRAY[-1, -1]),
                                                         STRUCT('sample5', ARRAY[0, 1])])
])"""

  # Basic input for a 'multisample variants table' with a schema optimized for
  # large cohorts. This includes a variety of variant types and genotypes.
  OPT_MULTISAMPLE_VARIANTS_TABLE = "basic_optimized_multisample_variants_table"
  OPT_MULTISAMPLE_VARIANTS_TABLE_INPUT = """
SELECT * FROM UNNEST([
STRUCT<reference_name STRING,
       start_position INT64,
       end_position INT64,
       reference_bases STRING,
       AN INT64,
       alternate_bases ARRAY<STRUCT<alt STRING, AC INT64, AF FLOAT64>>,
       call ARRAY<STRUCT<name STRING,
                         genotype ARRAY<INT64>>>,
       hom_ref_call ARRAY<STRING> >
    -- SNP with ref-homozygous, heterozygous, and alt-homozygous genotypes.
    ('chr1', 123, 124, 'A', 6, ARRAY[STRUCT('T', 3, .5),
                                     STRUCT('<*>', 0, 0.0)], ARRAY[STRUCT('sample2', ARRAY[0, 1]),
                                                                   STRUCT('sample3', ARRAY[1, 1]),
                                                                   STRUCT('sample4', ARRAY[-1, -1]),
                                                                   STRUCT('sample5', ARRAY[-1, -1])],
                                                             ARRAY['sample1']),
    -- SNP with no-calls
    ('chr1', 456, 457, 'C', 1, ARRAY[STRUCT('T', 1, 1),
                                     STRUCT('<*>', 0, 0.0)], ARRAY[STRUCT('sample1', ARRAY[-1, 1]),
                                                                   STRUCT('sample2', ARRAY[-1, -1]),
                                                                   STRUCT('sample3', ARRAY[-1, -1]),
                                                                   STRUCT('sample4', ARRAY[-1, -1]),
                                                                   STRUCT('sample5', ARRAY[-1, -1])],
                                                             NULL),
    -- Multiallelic SNP
    ('chr1', 456, 457, 'C', 2, ARRAY[STRUCT('T', 1, .5),
                                     STRUCT('G', 1, .5),
                                     STRUCT('<*>', 0, 0.0)], ARRAY[STRUCT('sample1', ARRAY[-1, -1]),
                                                                   STRUCT('sample2', ARRAY[-1, -1]),
                                                                   STRUCT('sample3', ARRAY[1, 2]),
                                                                   STRUCT('sample4', ARRAY[-1, -1]),
                                                                   STRUCT('sample5', ARRAY[-1, -1])],
                                                             NULL),
    -- Insertion
    ('chr1', 456, 457, 'C', 2, ARRAY[STRUCT('CC', 2, 1),
                                     STRUCT('<*>', 0, 0.0)], ARRAY[STRUCT('sample1', ARRAY[-1, -1]),
                                                                   STRUCT('sample2', ARRAY[-1, -1]),
                                                                   STRUCT('sample3', ARRAY[-1, -1]),
                                                                   STRUCT('sample4', ARRAY[1, 1]),
                                                                   STRUCT('sample5', ARRAY[-1, -1])],
                                                             NULL),
    -- Deletion
    ('chr1', 456, 458, 'CA', 2, ARRAY[STRUCT('C', 1, .5),
                                      STRUCT('<*>', 0, 0.0)], ARRAY[STRUCT('sample1', ARRAY[-1, -1]),
                                                                    STRUCT('sample2', ARRAY[-1, -1]),
                                                                    STRUCT('sample3', ARRAY[-1, -1]),
                                                                    STRUCT('sample4', ARRAY[-1, -1]),
                                                                    STRUCT('sample5', ARRAY[0, 1])],
                                                              NULL)
])"""

  # Test against all three tables when applicable (e.g., when the query ignores
  # non-variant segments and does not need to know which samples match the
  # reference).
  GENOME_CALL_OR_MULTISAMPLE_TABLES = {
    GENOME_CALL_TABLE: GENOME_CALL_TABLE_INPUT,
    MULTISAMPLE_VARIANTS_TABLE: MULTISAMPLE_VARIANTS_TABLE_INPUT,
    OPT_MULTISAMPLE_VARIANTS_TABLE: OPT_MULTISAMPLE_VARIANTS_TABLE_INPUT
  }
