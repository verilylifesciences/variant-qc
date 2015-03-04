#!/usr/bin/env python

# Copyright 2014 Google Inc. All Rights Reserved.
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

"""A reducer for expansion of gVCF or CGI data.

To test locally, set both BIG_QUERY_SOURCE and BIG_QUERY_SINK to False and run:
  cat ./data/platinum-genomes-brca1.json | ./gvcf-expand-mapper.py | sort \
  | ./gvcf-expand-reducer.py > out.json
"""

import json
import sys

from gvcf_expander import GvcfExpander
from gvcf_expander import Pair

# Whether the source data from this job is coming from the BigQuery connector
# for Hadoop Streaming
BIG_QUERY_SINK = True


def main():
  """Entry point to the script."""

  # Basic parsing of command line arguments to allow a filename
  # to be passed when running this code in the debugger.
  file_handle = sys.stdin
  if 2 <= len(sys.argv):
    file_handle = open(sys.argv[1], "r")

  expander = GvcfExpander(bin_size=1000,
                          filter_ref_matches=False,
                          emit_ref_blocks=False)

  line = file_handle.readline()
  while line:
    line = line.strip()
    if not line:
      line = file_handle.readline()
      continue

    (key, value) = line.split("\t")
    fields = json.loads(value)
    results = expander.reduce(pair=Pair(key, fields))

    for result in results:
      emit(result)

    line = file_handle.readline()

  results = expander.finalize()

  for result in results:
    emit(result)


def emit(fields):
  """Emits a reduced value to stdout.

  Args:
    fields: (dict)

  Returns: n/a

  Side Effects:
    a value is written to stdout
  """
  if BIG_QUERY_SINK:
    print "0\t%s" % (json.dumps(fields))
  else:
    print "%s" % (json.dumps(fields))


if __name__ == "__main__":
  main()
