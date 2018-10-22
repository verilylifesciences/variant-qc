#!/bin/bash

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

set -o nounset
set -o errexit
set -o xtrace

virtualenv --system-site-packages virtualTestEnv
# Work around virtual env error 'PS1: unbound variable'
set +o nounset
source virtualTestEnv/bin/activate
set -o nounset

echo Setting up
pip install --upgrade pip
pip install --upgrade setuptools
# The test scripts in this package make use of the BigQuery test framework
# in this repository. Its not currently available on https://pypi.org/ so we
# install a particular release of it from GitHub instead.
pip install git+https://github.com/verilylifesciences/analysis-py-utils.git@v0.3.0

tests=`find ./sql -name '*_test.py'`

for test in $tests; do
  echo Running $test
  python $test
done
