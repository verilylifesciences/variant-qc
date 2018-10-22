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

"""Tests for query sample_level_missingness.sql.

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
            "sample_level_missingness.sql"),
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
    -- SNP with no-calls
    ('chr1', 123, 124, 'A', ARRAY[STRUCT('T')],  ARRAY[STRUCT('sample1', ARRAY[1, -1]),
                                                       STRUCT('sample2', ARRAY[0, 1])]),
    -- Deletion with no-calls
    ('chr1', 456, 458, 'AA', ARRAY[STRUCT('A')], ARRAY[STRUCT('sample1', ARRAY[1, 1]),
                                                       STRUCT('sample2', ARRAY[-1, -1])]),
    -- Non-variant segments
    ('chr1', 500, 800, 'C', NULL,                ARRAY[STRUCT('sample1', ARRAY[0, 0])]),
    ('chr1', 750, 900, 'C', NULL,                ARRAY[STRUCT('sample2', ARRAY[0, 0])])
])
    """, cls.src_table_name)

  def test(self):
    sql = Template(self.sql_to_test).render(
        {"GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE": self.src_table_name,
         "LOW_QUALITY_CALLS_FILTER": "FALSE"})
    expected = [
        # name, no_calls, all_calls, missingness_rate
        ("sample1", 1, 303, 1/303),
        ("sample2", 2, 153, 2/153)
        ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)

  @data(TestConstants.GENOME_CALL_TABLE)
  def test_basic_input_genome_call_table(self, table):
    sql = Template(self.sql_to_test).render(
        {"GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE": self.basic_tables[table],
         "LOW_QUALITY_CALLS_FILTER": "FALSE"})
    expected = [
        # name, no_calls, all_calls, missingness_rate
        ("sample1", 1, 2, 0.5),
        ("sample2", 1, 2, 0.5),
        ("sample3", 0, 2, 0),
        ("sample4", 0, 1, 0),
        ("sample5", 0, 2, 0),
        ("sample6", 0, 700, 0),
        ("sample7", 0, 700, 0),
    ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)

  @data(TestConstants.MULTISAMPLE_VARIANTS_TABLE)
  def test_basic_input_multisample_variants_table(self, table):
    sql = Template(self.sql_to_test).render(
        {"GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE": self.basic_tables[table],
         "LOW_QUALITY_CALLS_FILTER": "FALSE"})
    expected = [
        # name, no_calls, all_calls, missingness_rate
        ("sample1", 5, 6, 5/6),
        ("sample2", 5, 6, 5/6),
        ("sample3", 4, 6, 4/6),
        ("sample4", 5, 6, 5/6),
        ("sample5", 4, 6, 4/6)
    ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)


if __name__ == "__main__":
  unittest.main()
