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

"""Tests for query characterize_alts.sql.

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
            "characterize_alts.sql"),
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
       alternate_bases ARRAY<STRUCT<alt STRING>> >
    -- SNP and insertion
    ('chr1', 123, 124, 'A', ARRAY[STRUCT('T'), STRUCT('TTTTTT')]),
    -- Deletion
    ('chr1', 123, 127, 'ATAG', ARRAY[STRUCT('A')]),
    -- An alternate with special characters.
    ('chr1', 123, 200, 'A', ARRAY[STRUCT('<DEL>')])
  ])
  """, cls.src_table_name)

  def test(self):
    sql = Template(self.sql_to_test).render(
        {"GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE": self.src_table_name})
    expected = [
        # number_of_variant_records, alt_contains_no_special_characters,
        # max_ref_len, max_alt_len
        (3, True, 4, 6),
        (1, False, 1, 5),
        ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)

  @data(TestConstants.GENOME_CALL_TABLE,
        TestConstants.MULTISAMPLE_VARIANTS_TABLE,
        TestConstants.OPT_MULTISAMPLE_VARIANTS_TABLE)
  def test_basic_input(self, table):
    sql = Template(self.sql_to_test).render(
        {"GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE": self.basic_tables[table]})
    expected = [
        # number_of_variant_records, alt_contains_no_special_characters,
        # max_ref_len, max_alt_len
        (6, True, 2, 2),
        ]
    self.expect_query_result(
        query=sql, expected=expected, enforce_ordering=False)

if __name__ == "__main__":
  unittest.main()
