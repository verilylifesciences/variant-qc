/*
 * Copyright (C) 2015 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */
package com.google.cloud.genomics.examples;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hamcrest.CoreMatchers;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import com.google.api.services.bigquery.model.TableRow;
import com.google.cloud.dataflow.sdk.transforms.DoFnTester;
import com.google.cloud.genomics.examples.TransformNonVariantSegmentData.FilterCallsFn;
import com.google.cloud.genomics.examples.TransformNonVariantSegmentData.FlagVariantsWithAmbiguousCallsFn;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.google.protobuf.ListValue;
import com.google.protobuf.Value;

@RunWith(JUnit4.class)
public class TransformNonVariantSegmentDataTest {

  @Test
  public void testFilterVariantCallsFn() {

    Map<String, ListValue> passingFilter = new HashMap<String, ListValue>();
    passingFilter.put(
        TransformNonVariantSegmentData.FILTER_FIELD,
        ListValue.newBuilder()
            .addValues(Value.newBuilder().setStringValue(TransformNonVariantSegmentData.PASSING_FILTER).build())
            .build());
    VariantCall call1 = VariantCall.newBuilder().putAllInfo(passingFilter).build();
    
    Map<String, ListValue> failingFilter = new HashMap<String, ListValue>();
    failingFilter.put(
        TransformNonVariantSegmentData.FILTER_FIELD,
        ListValue.newBuilder()
            .addValues(Value.newBuilder().setStringValue("VQSRTrancheSNP99.90to100.00").build())
            .build());
    VariantCall call2 = VariantCall.newBuilder().putAllInfo(failingFilter).build();

    Map<String, ListValue> ambiguousFilter = new HashMap<String, ListValue>();
    ambiguousFilter.put(
        TransformNonVariantSegmentData.FILTER_FIELD,
        ListValue.newBuilder()
            .addValues(Value.newBuilder().setStringValue("VQSRTrancheSNP99.90to100.00").build())
            .addValues(Value.newBuilder().setStringValue(TransformNonVariantSegmentData.PASSING_FILTER).build())
            .build());    
    VariantCall call3 = VariantCall.newBuilder().putAllInfo(ambiguousFilter).build();
    
    Variant inputVariant =
        Variant.newBuilder().addAllCalls(Arrays.asList(call1, call2, call3)).build();
    
    Variant expectedVariant = Variant.newBuilder().addAllCalls(Arrays.asList(call1, call3)).build();

    DoFnTester<Variant, Variant> filterCallsFn = DoFnTester.of(new FilterCallsFn());
    Assert.assertThat(filterCallsFn.processBatch(inputVariant),
        CoreMatchers.allOf(CoreMatchers.hasItems(expectedVariant)));
  }

  @Test
  public void testAmbiguousVariantCallsFns() {

    DoFnTester<Variant, Variant> flagVariantsFn =
        DoFnTester.of(new TransformNonVariantSegmentData.FlagVariantsWithAmbiguousCallsFn());

    VariantCall call1 =
        VariantCall.newBuilder().setCallSetName("sample1").addAllGenotype(Arrays.asList(0, 1))
            .build();
    VariantCall call2a =
        VariantCall.newBuilder().setCallSetName("sample2").addAllGenotype(Arrays.asList(0, 1))
            .build();
    VariantCall call2b =
        VariantCall.newBuilder().setCallSetName("sample2").addAllGenotype(Arrays.asList(-1, -1))
            .build();

    Variant inputVariant = Variant.newBuilder().addAllCalls(Arrays.asList(call1, call2a)).build();
    Variant ambiguousInputVariant =
        Variant.newBuilder().addAllCalls(Arrays.asList(call1, call2a, call2b)).build();

    Variant expectedVariant =
        Variant.newBuilder().addAllCalls(Arrays.asList(call1, call2a))
            .putAllInfo(FlagVariantsWithAmbiguousCallsFn.NO_AMBIGUOUS_CALLS_INFO).build();
    Variant ambiguousExpectedVariant =
        Variant.newBuilder().addAllCalls(Arrays.asList(call1, call2a, call2b))
            .putAllInfo(FlagVariantsWithAmbiguousCallsFn.HAS_AMBIGUOUS_CALLS_INFO).build();

    Assert.assertThat(flagVariantsFn.processBatch(inputVariant, ambiguousInputVariant),
        CoreMatchers.allOf(CoreMatchers.hasItems(expectedVariant, ambiguousExpectedVariant)));
        
    DoFnTester<Variant, TableRow> formatVariantsFn = DoFnTester.of(new TransformNonVariantSegmentData.FormatVariantsFn());
    List<TableRow> rows = formatVariantsFn.processBatch(expectedVariant, ambiguousExpectedVariant);
    assertEquals(2, rows.size());
    assertEquals("false", rows.get(0).get(TransformNonVariantSegmentData.HAS_AMBIGUOUS_CALLS_FIELD).toString());
    assertEquals("true", rows.get(1).get(TransformNonVariantSegmentData.HAS_AMBIGUOUS_CALLS_FIELD).toString());
  }

}
