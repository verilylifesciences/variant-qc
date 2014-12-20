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

"""A library for expansion of gVCF data.

Given gVCF data as input, this library will find the reference
matching blocks which overlap a variant and add the reference-matching
sample genotypes to the variant record.

At this current time, only SNPs are expanded.  Other variant types are
left as-is.

"""

from collections import namedtuple
import json
import math
import unittest


Pair = namedtuple('Pair', 'k v')


class GvcfExpander(object):
  """Common logic for gVCF expansion."""

  def __init__(self, bin_size=100,
               filter_ref_matches=False,
               emit_ref_blocks=True):
    self.bin_size = bin_size
    self.filter_ref_matches = filter_ref_matches
    self.emit_ref_blocks = emit_ref_blocks
    self.current_key = None
    self.binned_calls = []
    self.sample_refs = {}

  def is_variant(self, fields):
    """Determines whether or not the VCF fields constitute a variant."""
    if 'alternate_bases' in fields and fields['alternate_bases']:
      return True
    return False

  def is_snp(self, fields):
    """Determines whether or not the VCF fields constitute a SNP."""
    if (self.is_variant(fields)
        and fields['reference_bases'] in ['A', 'C', 'G', 'T']
        and all(alt in ['A', 'C', 'G', 'T']
                for alt in fields['alternate_bases'])):
      return True
    return False

  def get_start(self, fields):
    return int(fields['start'])

  def get_end(self, fields):
    if 'END' in fields:
      return int(fields['END'])
    else:
      return int(fields['end'])

  def compute_start_bin(self, fields):
    return int(math.floor(self.get_start(fields)/self.bin_size))

  def compute_end_bin(self, fields):
    return int(math.floor(self.get_end(fields)/self.bin_size))

  def map(self, fields):
    results = []

    start_bin = self.compute_start_bin(fields)
    if self.is_variant(fields):
      end_bin = start_bin
    else:
      end_bin = self.compute_end_bin(fields)

    for current_bin in range(start_bin, end_bin+1):
      key = '%s:%s' % (fields['reference_name'], str(current_bin))
      results.append(Pair(key, fields))

    return results

  def reduce(self, pair):
    expanded_calls = []

    if None is self.current_key:
      self.current_key = pair.k

    if pair.k != self.current_key:
      expanded_calls = self.expand_binned_calls()
      self.current_key = pair.k
      self.binned_calls = []
      self.sample_refs = {}

    self.binned_calls.append(pair.v)

    return expanded_calls

  def finalize(self):
    return self.expand_binned_calls()

  def expand_binned_calls(self):
    expanded_calls = []
    current_bin = int(self.current_key.split(':')[1])

    # Sort by start position descending and ensure that if a variant
    # and a ref-matching block are at the same position, the
    # ref-matching block comes first.
    calls = sorted(self.binned_calls,
                   key=lambda k: (int(k['start']),
                                  len(k.get('alternate_bases', []))))

    for call in calls:
      if self.is_variant(call):
        expanded_calls.append(self.expand_variant(call))
      else:
        self.accumulate_block(call)
        if self.emit_ref_blocks:
          # Don't output ref-matches that we previously emitted
          if self.compute_start_bin(call) == current_bin:
            expanded_calls.append(call)

    return expanded_calls

  def expand_variant(self, variant_call):
    # Only expand SNPs.  Leave indels as-is.
    if not self.is_snp(variant_call):
      return variant_call

    expansion_calls = []
    for sample_id in self.sample_refs.keys():
      ref_call = self.sample_refs[sample_id]
      # This criteria is correct for SNPs.  It would need to be modified
      # if we wish to support insertions or deletions
      if (self.get_start(ref_call) <= self.get_start(variant_call) and
          self.get_end(ref_call) >= self.get_start(variant_call) + 1):
        expansion_calls.extend(ref_call['call'])
      else:
        # This ref_block is now outside of the range we are
        # considering and therefore obsolete, nuke it
        del self.sample_refs[sample_id]

    if self.filter_ref_matches:
      # Get the sample_ids already called in this variant
      variant_sample_names = [call['call_set_name'] for call in
                              variant_call['call']]

      # Filter out ref_calls for samples called in this variant; this
      # is naive; a better approach might be to compare the quality
      # score of the referenece call to that of the variant call for
      # the sample in question
      variant_call['call'].extend(
          [call for call in expansion_calls
           if call['call_set_name'] not in variant_sample_names]
          )
    else:
      variant_call['call'].extend(expansion_calls)

    return variant_call

  def accumulate_block(self, ref_call):
    # Since everything is sorted by start, we only need to stash
    # at most one ref block per sample in the ref block at a time
    self.sample_refs[ref_call['call'][0]['call_set_name']] = ref_call


class GvcfExpanderTest(unittest.TestCase):
  """Unit tests for common logic for gVCF expansion."""

  def test_is_variant(self):
    expander = GvcfExpander()

    self.assertTrue(expander.is_variant(json.loads(self.snp_1)))
    self.assertTrue(expander.is_variant(json.loads(self.snp_2)))
    self.assertTrue(expander.is_variant(json.loads(self.insertion_1)))
    self.assertTrue(expander.is_variant(json.loads(self.deletion_1)))
    self.assertFalse(expander.is_variant(json.loads(self.ref_a)))
    self.assertFalse(expander.is_variant(json.loads(self.ref_b)))
    self.assertFalse(expander.is_variant(json.loads(self.ref_c)))
    self.assertFalse(expander.is_variant(json.loads(self.ref_d)))
    self.assertFalse(expander.is_variant(json.loads(self.ref_ambiguous)))
    self.assertFalse(expander.is_variant(json.loads(self.no_call_1)))

  def test_is_snp(self):
    expander = GvcfExpander()

    self.assertTrue(expander.is_snp(json.loads(self.snp_1)))
    self.assertTrue(expander.is_snp(json.loads(self.snp_2)))
    self.assertFalse(expander.is_snp(json.loads(self.insertion_1)))
    self.assertFalse(expander.is_snp(json.loads(self.deletion_1)))
    self.assertFalse(expander.is_snp(json.loads(self.ref_a)))
    self.assertFalse(expander.is_snp(json.loads(self.ref_b)))
    self.assertFalse(expander.is_snp(json.loads(self.ref_c)))
    self.assertFalse(expander.is_snp(json.loads(self.ref_d)))
    self.assertFalse(expander.is_snp(json.loads(self.ref_ambiguous)))
    self.assertFalse(expander.is_snp(json.loads(self.no_call_1)))

  def test_mapper_variant(self):
    test_input = json.loads(self.snp_1)

    expected_output = json.loads(self.snp_1)  # no changes

    expander = GvcfExpander()
    result = expander.map(fields=test_input)
    self.assertEqual(1, len(result))
    (key, value) = result[0]
    self.assertEqual('13:1022656', key)
    self.assertDictEqual(expected_output, value)

  def test_mapper_ref_block(self):
    test_input = json.loads(self.ref_a)

    expected_output = json.loads(self.ref_a)  # no changes

    expander = GvcfExpander()
    result = expander.map(fields=test_input)
    self.assertEqual(3, len(result))
    self.assertEqual('13:1022656', result[0].k)
    self.assertEqual('13:1022657', result[1].k)
    self.assertEqual('13:1022658', result[2].k)
    for i in range(0, 3):
      self.assertDictEqual(expected_output, result[i].v)

  def test_mapper_no_call(self):
    test_input = json.loads(self.no_call_1)

    expected_output = json.loads(self.no_call_1)  # no changes

    expander = GvcfExpander()
    result = expander.map(fields=test_input)
    self.assertEqual(1, len(result))
    (key, value) = result[0]
    self.assertEqual('13:1022656', key)
    self.assertDictEqual(expected_output, value)

  def test_mr_variant(self):
    test_input = json.loads(self.snp_1)

    expected_output = json.loads(self.snp_1)  # no changes

    expander = GvcfExpander()
    result = expander.map(fields=test_input)
    self.assertEqual(1, len(result))
    (key, value) = result[0]
    self.assertEqual('13:1022656', key)
    self.assertDictEqual(expected_output, value)

    result = expander.reduce(pair=result[0])
    self.assertEqual(0, len(result))
    result = expander.finalize()
    self.assertEqual(1, len(result))
    value = result[0]
    self.assertDictEqual(expected_output, value)

  def test_mr_ref(self):
    test_input = json.loads(self.ref_a)

    expected_output = json.loads(self.ref_a)  # no changes

    expander = GvcfExpander()
    pairs = expander.map(fields=test_input)
    self.assertEqual(3, len(pairs))

    result = expander.reduce(pairs[0])
    self.assertEqual(0, len(result))
    result = expander.reduce(pairs[1])
    self.assertEqual(1, len(result))
    value = result[0]
    self.assertDictEqual(expected_output, value)

    for i in range(2, 3):
      result = expander.reduce(pairs[i])
      self.assertEqual(0, len(result))

    result = expander.finalize()
    self.assertEqual(0, len(result))

  def test_mr_no_call(self):
    test_input = json.loads(self.no_call_1)

    expected_output = json.loads(self.no_call_1)  # no changes

    expander = GvcfExpander()
    result = expander.map(fields=test_input)
    self.assertEqual(1, len(result))
    (key, value) = result[0]
    self.assertEqual('13:1022656', key)
    self.assertDictEqual(expected_output, value)

    result = expander.reduce(pair=result[0])
    self.assertEqual(0, len(result))
    result = expander.finalize()
    self.assertEqual(1, len(result))
    value = result[0]
    self.assertDictEqual(expected_output, value)

  def test_mr(self):

    for filter_ref_matches in [True, False]:
      expander = GvcfExpander(filter_ref_matches=filter_ref_matches)
      pairs = []

      pairs.extend(expander.map(fields=json.loads(self.ref_a)))
      pairs.extend(expander.map(fields=json.loads(self.ref_b)))
      pairs.extend(expander.map(fields=json.loads(self.ref_c)))
      pairs.extend(expander.map(fields=json.loads(self.ref_ambiguous)))
      pairs.extend(expander.map(fields=json.loads(self.snp_1)))
      pairs.extend(expander.map(fields=json.loads(self.snp_2)))
      pairs.extend(expander.map(fields=json.loads(self.no_call_1)))
      self.assertEqual(11, len(pairs))

      # Sort these by key so that all pairs in the same bin are in a
      # row, but not necessarily in any order within that bin
      pairs = sorted(pairs)

      self.assertEqual(0, len(expander.reduce(pairs[0])))
      self.assertEqual(0, len(expander.reduce(pairs[1])))
      self.assertEqual(0, len(expander.reduce(pairs[2])))
      self.assertEqual(0, len(expander.reduce(pairs[3])))
      self.assertEqual(0, len(expander.reduce(pairs[4])))
      self.assertEqual(0, len(expander.reduce(pairs[5])))
      self.assertEqual(0, len(expander.reduce(pairs[6])))
      result = expander.reduce(pairs[7])
      self.assertEqual(7, len(result))
      self.assertIn(json.loads(self.ref_b), result)
      self.assertIn(json.loads(self.ref_c), result)
      self.assertIn(json.loads(self.ref_a), result)
      self.assertIn(json.loads(self.ref_ambiguous), result)
      self.assertIn(json.loads(self.no_call_1), result)

      if filter_ref_matches:
        self.assertIn(json.loads(self.expanded_snp_1_filtered), result)
        self.assertIn(json.loads(self.expanded_snp_2), result)
      else:
        self.assertIn(json.loads(self.expanded_snp_1), result)
        self.assertIn(json.loads(self.expanded_snp_2), result)

      for i in range(7, 11):
        result = expander.reduce(pairs[i])
        self.assertEqual(0, len(result))

  def test_same_start(self):
    expander = GvcfExpander(emit_ref_blocks=False)
    expander.reduce(expander.map(fields=json.loads('{"reference_name":"chr17","start":"41196822","end":"41196841","reference_bases":"T","alternate_bases":[],"quality":0,"filter":["PASS"],"names":[],"call":[{"call_set_id":"3049512673186936334-5","call_set_name":"NA12883","genotype":[0,0],"genotype_likelihood":[],"AD":[],"DP":"39","FILTER":["PASS"],"GQX":"78","MQ":"48","PL":[],"QUAL":0}],"AC":[],"AF":[],"BLOCKAVG_min30p3a":true}'))[0])
    expander.reduce(expander.map(fields=json.loads('{"reference_name":"chr17","start":"41196839","end":"41196841","reference_bases":"T","alternate_bases":[],"quality":0,"filter":["PASS"],"names":[],"call":[{"call_set_id":"3049512673186936334-15","call_set_name":"NA12892","genotype":[0,0],"genotype_likelihood":[],"AD":[],"DP":"32","FILTER":["PASS"],"GQX":"63","MQ":"47","PL":[],"QUAL":0}],"AC":[],"AF":[],"BLOCKAVG_min30p3a":true}'))[0])
    expander.reduce(expander.map(fields=json.loads('{"reference_name":"chr17","start":"41196823","end":"41196841","reference_bases":"T","alternate_bases":[],"quality":0,"filter":["PASS"],"names":[],"call":[{"call_set_id":"3049512673186936334-7","call_set_name":"NA12880","genotype":[0,0],"genotype_likelihood":[],"AD":[],"DP":"56","FILTER":["PASS"],"GQX":"99","MQ":"50","PL":[],"QUAL":0}],"AC":[],"AF":[],"BLOCKAVG_min30p3a":true}'))[0])
    expander.reduce(expander.map(fields=json.loads('{"reference_name":"chr17","start":"41196840","end":"41196843","reference_bases":"G","alternate_bases":[],"quality":0,"filter":["PASS"],"names":[],"call":[{"call_set_id":"3049512673186936334-2","call_set_name":"NA12889","genotype":[0,0],"genotype_likelihood":[],"AD":[],"DP":"27","FILTER":["PASS"],"GQX":"63","MQ":"54","PL":[],"QUAL":0}],"AC":[],"AF":[],"BLOCKAVG_min30p3a":true}'))[0])
    expander.reduce(expander.map(fields=json.loads('{"reference_name":"chr17","start":"41196840","end":"41196841","reference_bases":"G","alternate_bases":[],"quality":91.489999999999995,"filter":["PASS"],"names":[],"call":[{"call_set_id":"3049512673186936334-0","call_set_name":"NA12882","genotype":[0,0],"genotype_likelihood":[],"AD":[],"DP":"40","FILTER":["LowGQX"],"GQX":"1","MQ":"53","PL":[],"QUAL":0},{"call_set_id":"3049512673186936334-1","call_set_name":"NA12877","genotype":[0,0],"genotype_likelihood":[],"AD":["40"],"DP":"48","FILTER":["PASS"],"GQ":61.490000000000002,"GQX":"61","MQ":"52","PL":[],"QUAL":91.489999999999995},{"call_set_id":"3049512673186936334-3","call_set_name":"NA12885","genotype":[0,0],"genotype_likelihood":[],"AD":["34"],"DP":"43","FILTER":["PASS"],"GQ":62.57,"GQX":"63","MQ":"48","PL":[],"QUAL":92.560000000000002},{"call_set_id":"3049512673186936334-8","call_set_name":"NA12891","genotype":[0,0],"genotype_likelihood":[],"AD":["31"],"DP":"39","FILTER":["LowGQX"],"GQ":0.01,"GQX":"0","MQ":"51","PL":[],"QUAL":6.4699999999999998},{"call_set_id":"3049512673186936334-10","call_set_name":"NA12886","genotype":[0,0],"genotype_likelihood":[],"AD":["34"],"DP":"42","FILTER":["PASS"],"GQ":93.260000000000005,"GQX":"93","MQ":"48","PL":[],"QUAL":123.26000000000001},{"call_set_id":"3049512673186936334-11","call_set_name":"NA12884","genotype":[0,0],"genotype_likelihood":[],"AD":["23"],"DP":"29","FILTER":["LowGQX"],"GQ":0.22,"GQX":"0","MQ":"49","PL":[],"QUAL":17.329999999999998},{"call_set_id":"3049512673186936334-12","call_set_name":"NA12890","genotype":[0,0],"genotype_likelihood":[],"AD":[],"DP":"28","FILTER":["PASS"],"GQX":"39","MQ":"55","PL":[],"QUAL":0},{"call_set_id":"3049512673186936334-14","call_set_name":"NA12878","genotype":[0,0],"genotype_likelihood":[],"AD":["33"],"DP":"41","FILTER":["LowGQX"],"GQ":0.34999999999999998,"GQX":"0","MQ":"49","PL":[],"QUAL":19.309999999999999},{"call_set_id":"3049512673186936334-16","call_set_name":"NA12881","genotype":[0,0],"genotype_likelihood":[],"AD":["33"],"DP":"43","FILTER":["PASS"],"GQ":93.260000000000005,"GQX":"93","MQ":"50","PL":[],"QUAL":123.26000000000001}],"AC":[],"AF":[],"DP":"50","MQ":"52","MQ0":"0"}'))[0])
    expander.reduce(expander.map(fields=json.loads('{"reference_name":"chr17","start":"41196840","end":"41196844","reference_bases":"G","alternate_bases":[],"quality":0,"filter":["PASS"],"names":[],"call":[{"call_set_id":"3049512673186936334-9","call_set_name":"NA12888","genotype":[0,0],"genotype_likelihood":[],"AD":[],"DP":"24","FILTER":["PASS"],"GQX":"66","MQ":"50","PL":[],"QUAL":0}],"AC":[],"AF":[],"BLOCKAVG_min30p3a":true}'))[0])
    expander.reduce(expander.map(fields=json.loads('{"reference_name":"chr17","start":"41196840","end":"41196841","reference_bases":"G","alternate_bases":["T"],"quality":85.680000000000007,"filter":["TruthSensitivityTranche99.90to100.00","LowQD"],"names":[],"call":[{"call_set_id":"3049512673186936334-6","call_set_name":"NA12879","genotype":[0,1],"genotype_likelihood":[],"AD":["30","13"],"DP":"44","FILTER":["TruthSensitivityTranche99.90to100.00","LowQD"],"GQ":99,"GQX":"86","MQ":"49","PL":["116","0","744"],"QUAL":85.680000000000007,"VF":0.30199999999999999},{"call_set_id":"3049512673186936334-13","call_set_name":"NA12893","genotype":[0,1],"genotype_likelihood":[],"AD":["31","6"],"DP":"37","FILTER":["TruthSensitivityTranche99.90to100.00","LowQD"],"GQ":74.269999999999996,"GQX":"44","MQ":"57","PL":["74","0","1010"],"QUAL":44.280000000000001,"VF":0.16200000000000001}],"AC":["1"],"AF":[0.5],"AN":"2","BaseQRankSum":-2.2509999999999999,"DP":"48","Dels":0.080000000000000002,"FS":36.518000000000001,"HRun":"19","HaplotypeScore":23.5017,"MQ":"49","MQ0":"0","MQRankSum":0.51500000000000001,"QD":1.79,"ReadPosRankSum":-2.8980000000000001,"SB":-0.01,"VQSLOD":-1.5014000000000001,"culprit":"QD","set":"FilteredInAll"}'))[0])
    expander.reduce(expander.map(fields=json.loads('{"reference_name":"chr17","start":"41196840","end":"41196845","reference_bases":"G","alternate_bases":[],"quality":0,"filter":["PASS"],"names":[],"call":[{"call_set_id":"3049512673186936334-4","call_set_name":"NA12887","genotype":[0,0],"genotype_likelihood":[],"AD":[],"DP":"37","FILTER":["PASS"],"GQX":"81","MQ":"50","PL":[],"QUAL":0}],"AC":[],"AF":[],"BLOCKAVG_min30p3a":true}'))[0])

    result = expander.finalize()
    self.assertEqual(1, len(result))
    self.assertEqual(17, len(result[0]['call']))

  def setUp(self):
    self.maxDiff = None

    self.ref_a = """
{
  "reference_name": "13",
  "start": "102265642",
  "reference_bases": "A",
  "END": "102265842",
  "call": [
    {
      "callset_id": "1",
      "call_set_name": "same_start",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id": "2",
      "call_set_name": "same_start_second_sample",
      "genotype":[
        0,
        0
      ]
    }
  ]
}
"""

    self.ref_b = """
{
  "reference_name": "13",
  "start": "102265602",
  "reference_bases": "A",
  "END": "102265842",
  "call": [
    {
      "callset_id": "1",
      "call_set_name": "different_start",
      "genotype":[
        0,
        0
      ]
    }
  ]
}
"""

    self.ref_ambiguous = """
{
  "reference_name": "13",
  "start": "102265642",
  "reference_bases": "A",
  "END": "102265650",
  "call": [
    {
      "callset_id": "3",
      "call_set_name": "ambiguous",
      "genotype":[
        0,
        0
      ]
    }
  ]
}
"""

    self.ref_c = """
{
  "reference_name": "13",
  "start": "102265602",
  "reference_bases": "A",
  "END": "102265642",
  "call": [
    {
      "callset_id": "1",
      "call_set_name": "does_not_overlap_var_1",
      "genotype":[
        0,
        0
      ]
    }
  ]
}
"""

    self.ref_d = """
{
  "reference_name": "chr14",
  "start": "22973137",
  "end": "22973138",
  "reference_bases": "A",
  "alternate_bases": [
  ],
  "END": "22973211",
  "call": [
    {
      "callset_id": "715080930289-10",
      "call_set_name": "foo1",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id": "715080930289-14",
      "call_set_name": "foo2",
      "genotype":[
        0,
        0
      ]
    }
  ]
}
"""

    self.no_call_1 = """
{
  "reference_name": "13",
  "start": "102265642",
  "end": "102265645",
  "reference_bases":"TGA",
  "alternate_bases":[

  ],
  "call":[
    {
      "callset_id":"7122130836277736291-116",
      "phaseset":"7278593",
      "call_set_name":"no_call",
      "genotype":[
        -1,
        -1
      ]
    }
  ]
}
"""

    self.snp_1 = """
{
  "reference_name": "13",
  "start": "102265642",
  "end": "102265643",
  "reference_bases": "A",
  "alternate_bases": [
    "G"
  ],
  "call": [
    {
      "callset_id": "383928317087-12",
      "call_set_name": "hu52B7E5",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-34",
      "call_set_name": "hu1187FF",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-38",
      "call_set_name": "huC434ED",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "3",
      "call_set_name": "ambiguous",
      "genotype":[
        1,
        0
      ]
    }
  ]
}
"""

    self.expanded_snp_1 = """
{
  "reference_name": "13",
  "start": "102265642",
  "end": "102265643",
  "reference_bases": "A",
  "alternate_bases": [
    "G"
  ],
  "call": [
    {
      "callset_id": "383928317087-12",
      "call_set_name": "hu52B7E5",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-34",
      "call_set_name": "hu1187FF",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-38",
      "call_set_name": "huC434ED",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "3",
      "call_set_name": "ambiguous",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "1",
      "call_set_name": "different_start",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id": "1",
      "call_set_name": "same_start",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id": "2",
      "call_set_name": "same_start_second_sample",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id":"7122130836277736291-116",
      "phaseset":"7278593",
      "call_set_name":"no_call",
      "genotype":[
        -1,
        -1
      ]
    },
    {
      "callset_id": "3",
      "call_set_name": "ambiguous",
      "genotype":[
        0,
        0
      ]
    }
  ]
}
"""

    self.expanded_snp_1_filtered = """
{
  "reference_name": "13",
  "start": "102265642",
  "end": "102265643",
  "reference_bases": "A",
  "alternate_bases": [
    "G"
  ],
  "call": [
    {
      "callset_id": "383928317087-12",
      "call_set_name": "hu52B7E5",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-34",
      "call_set_name": "hu1187FF",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-38",
      "call_set_name": "huC434ED",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "3",
      "call_set_name": "ambiguous",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "1",
      "call_set_name": "different_start",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id": "1",
      "call_set_name": "same_start",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id": "2",
      "call_set_name": "same_start_second_sample",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id":"7122130836277736291-116",
      "phaseset":"7278593",
      "call_set_name":"no_call",
      "genotype":[
        -1,
        -1
      ]
    }
  ]
}
"""

    self.snp_2 = """
{
  "reference_name": "13",
  "start": "102265640",
  "end": "102265641",
  "reference_bases": "A",
  "alternate_bases": [
    "T"
  ],
  "call": [
    {
      "callset_id": "383928317087-12",
      "call_set_name": "hu52B7E5",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-34",
      "call_set_name": "hu1187FF",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-38",
      "call_set_name": "huC434ED",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-48",
      "call_set_name": "hu0211D6",
      "genotype":[
        1,
        0
      ]
    }
  ]
}
"""

    self.expanded_snp_2 = """
{
  "reference_name": "13",
  "start": "102265640",
  "end": "102265641",
  "reference_bases": "A",
  "alternate_bases": [
    "T"
  ],
  "call": [
    {
      "callset_id": "383928317087-12",
      "call_set_name": "hu52B7E5",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-34",
      "call_set_name": "hu1187FF",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-38",
      "call_set_name": "huC434ED",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "383928317087-48",
      "call_set_name": "hu0211D6",
      "genotype":[
        1,
        0
      ]
    },
    {
      "callset_id": "1",
      "call_set_name": "different_start",
      "genotype":[
        0,
        0
      ]
    },
    {
      "callset_id": "1",
      "call_set_name": "does_not_overlap_var_1",
      "genotype":[
        0,
        0
      ]
    }
  ]
}
"""

    self.insertion_1 = """
{
  "reference_name": "13",
  "start": "102265642",
  "end": "102265643",
  "reference_bases": "A",
  "alternate_bases": [
    "AGG"
  ],
  "call": [
    {
      "callset_id": "1",
      "call_set_name": "sample_with_an_insertion",
      "genotype":[
        1,
        0
      ]
    }
  ]
}
"""

    self.deletion_1 = """
{
  "reference_name": "13",
  "start": "102265642",
  "end": "102265644",
  "reference_bases": "AT",
  "alternate_bases": [
    "A"
  ],
  "call": [
    {
      "callset_id": "1",
      "call_set_name": "sample_with_a_deletion",
      "genotype":[
        1,
        0
      ]
    }
  ]
}
"""

if __name__ == '__main__':
  unittest.main()
