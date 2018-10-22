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

"""Tests for query homozygosity_coefficient.sql.

See https://github.com/verilylifesciences/analysis-py-utils for more details
about the testing framework.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from jinja2 import Template
import os
import unittest
from verily.bigquery_wrapper import bq_test_case
from test_constants import TestConstants


class QueryTest(bq_test_case.BQTestCase):

  @classmethod
  def setUpClass(cls):
    """Set up class."""
    super(QueryTest, cls).setUpClass(use_mocks=False)
    cls.sql_to_test = open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "homozygosity_coefficient.sql"),
        "r").read()

  @classmethod
  def create_mock_tables(cls):
    """Create mock tables."""
    # Basic input table.
    cls.basic_table_name = cls.client.path(
        TestConstants.MULTISAMPLE_VARIANTS_TABLE)
    cls.client.create_table_from_query(
        TestConstants.MULTISAMPLE_VARIANTS_TABLE_INPUT, cls.basic_table_name)

    # Table created specifically for this test.
    cls.src_table_name = cls.client.path("multisample_variant_table")
    cls.client.create_table_from_query("""
SELECT * FROM UNNEST([
STRUCT<reference_name STRING,
       start_position INT64,
       end_position INT64,
       reference_bases STRING,
       alternate_bases ARRAY<STRUCT<alt STRING>>,
       call ARRAY<STRUCT<name STRING,
                         genotype ARRAY<INT64>>> >
    ('chr1', 123, 124, 'A', ARRAY[STRUCT('T')], ARRAY[STRUCT('sample1', ARRAY[0, -1]),
                                                      STRUCT('sample2', ARRAY[0, 0]),
                                                      STRUCT('sample3', ARRAY[0, 1]),
                                                      STRUCT('sample4', ARRAY[0, 0]),
                                                      STRUCT('sample5', ARRAY[-1, -1]),
                                                      STRUCT('sample6', ARRAY[0, 0])]),
    ('chr1', 456, 457, 'G', ARRAY[STRUCT('C')], ARRAY[STRUCT('sample1', ARRAY[0, 0]),
                                                      STRUCT('sample2', ARRAY[-1, -1]),
                                                      STRUCT('sample3', ARRAY[0, 0]),
                                                      STRUCT('sample4', ARRAY[1, -1]),
                                                      STRUCT('sample5', ARRAY[1, 0]),
                                                      STRUCT('sample6', ARRAY[1, 1])]),
    ('chr1', 789, 790, 'C', ARRAY[STRUCT('G')], ARRAY[STRUCT('sample1', ARRAY[0, 0]),
                                                      STRUCT('sample2', ARRAY[-1, -1]),
                                                      STRUCT('sample3', ARRAY[0, 0]),
                                                      STRUCT('sample4', ARRAY[1, 1]),
                                                      STRUCT('sample5', ARRAY[1, 1]),
                                                      STRUCT('sample6', ARRAY[1, 1])])
])
    """, cls.src_table_name)

  def test(self):
    sql = Template(self.sql_to_test).render(
        {"MULTISAMPLE_VARIANT_TABLE": self.src_table_name,
         "HIGH_QUALITY_CALLS_FILTER": "TRUE"})
    expected = [
        # name, O_HOM, E_HOM, N_SITES, F
        ("sample1", 2, 0.93, 2, 1.0),
        ("sample2", 1, 0.75, 1, 1.0),
        ("sample3", 2, 1.68, 3, 0.24188),
        ("sample4", 2, 1.22, 2, 1.0),
        ("sample5", 1, 0.93, 2, 0.06459),
        ("sample6", 3, 1.68, 3, 1.0)
    ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)

  def test_basic_input(self):
    sql = Template(self.sql_to_test).render(
        {"MULTISAMPLE_VARIANT_TABLE": self.basic_table_name,
         "HIGH_QUALITY_CALLS_FILTER": "TRUE"})
    expected = [
        # name, O_HOM, E_HOM, N_SITES, F
        ("sample1", 1, 0.4, 1, 1),
        ("sample2", 0, 0.4, 1, -0.66667),
        ("sample3", 1, 0.4, 1, 1)
    ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)

if __name__ == "__main__":
  unittest.main()
