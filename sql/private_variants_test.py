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

"""Tests for query private_variants.sql.

See https://github.com/verilylifesciences/analysis-py-utils for more details
about the testing framework.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ddt import data
from ddt import ddt
from jinja2 import Template
import os
import unittest
from verily.bigquery_wrapper import bq_test_case
from test_constants import TestConstants


@ddt
class QueryTest(bq_test_case.BQTestCase):

  @classmethod
  def setUpClass(cls):
    """Set up class."""
    super(QueryTest, cls).setUpClass(use_mocks=False)
    cls.sql_to_test = open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "private_variants.sql"),
        "r").read()

  @classmethod
  def create_mock_tables(cls):
    """Create mock tables."""
    # Basic input tables.
    cls.basic_tables = {}
    for tbl_name, tbl_input in (TestConstants
                                .GENOME_CALL_OR_MULTISAMPLE_TABLES.iteritems()):
      cls.basic_tables[tbl_name] = cls.client.path(tbl_name)
      cls.client.create_table_from_query(tbl_input, tbl_name)

    # Table created specifically for this test.
    cls.src_table_name = cls.client.path("genome_call_table")
    cls.client.create_table_from_query("""
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
  ])
    """, cls.src_table_name)

  def test(self):
    sql = Template(self.sql_to_test).render(
        {"GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE": self.src_table_name,
         "HIGH_QUALITY_CALLS_FILTER": "TRUE"})
    expected = [
        # name, private_variant_count
        ("sample3", 1),
        ("sample4", 1),
        ("sample5", 1),
        ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)

  @data(TestConstants.GENOME_CALL_TABLE,
        TestConstants.MULTISAMPLE_VARIANTS_TABLE,
        TestConstants.OPT_MULTISAMPLE_VARIANTS_TABLE)
  def test_basic_input(self, table):
    sql = Template(self.sql_to_test).render(
        {"GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE": self.basic_tables[table],
         "HIGH_QUALITY_CALLS_FILTER": "TRUE"})
    expected = [
        # name, private_variant_count
        ("sample3", 1),
        ("sample4", 1),
        ("sample5", 1),
        ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)

if __name__ == "__main__":
  unittest.main()
