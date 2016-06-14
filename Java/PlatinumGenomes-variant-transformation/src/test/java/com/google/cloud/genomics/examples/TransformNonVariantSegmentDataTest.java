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
import static org.junit.Assert.assertFalse;

import com.google.api.services.bigquery.model.TableRow;
import com.google.cloud.dataflow.sdk.transforms.DoFnTester;
import com.google.cloud.genomics.examples.TransformNonVariantSegmentData.FilterCallsFn;
import com.google.cloud.genomics.examples.TransformNonVariantSegmentData.FlagVariantsWithAmbiguousCallsFn;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.google.protobuf.ListValue;
import com.google.protobuf.Value;

import org.hamcrest.CoreMatchers;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

@RunWith(JUnit4.class)
public class TransformNonVariantSegmentDataTest {

  @Test
  public void testFilterVariantCallsFn() {
    DoFnTester<Variant, Variant> filterCallsFn = DoFnTester.of(new FilterCallsFn());

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

    // Test a variant.
    Variant inputVariant = Variant.newBuilder()
        .setReferenceBases("A")
        .addAlternateBases("T")
        .addAllCalls(Arrays.asList(call1, call2, call3))
        .build();

    Variant expectedVariant = Variant.newBuilder()
        .setReferenceBases("A")
        .addAlternateBases("T")
        .addAllCalls(Arrays.asList(call1, call3))
        .build();

    Iterator<Variant> filtered1 = filterCallsFn.processBatch(inputVariant).iterator();
    assertEquals(filtered1.next(), expectedVariant);
    assertFalse(filtered1.hasNext());

    // Also test a non-variant segment.  These are not filtered.
    Variant inputBlockRecord = Variant.newBuilder()
        .setReferenceBases("A")
        .addAllCalls(Arrays.asList(call1, call2, call3))
        .build();

    Variant expectedBlockRecord = Variant.newBuilder()
        .setReferenceBases("A")
        .addAllCalls(Arrays.asList(call1, call2, call3))
        .build();

    Iterator<Variant> filtered2 = filterCallsFn.processBatch(inputBlockRecord).iterator();
    assertEquals(filtered2.next(), expectedBlockRecord);
    assertFalse(filtered2.hasNext());
  }

  @Test
  public void testFilterVariantCallsFn_AllCallsRemoved() {

    Map<String, ListValue> passingFilter = new HashMap<String, ListValue>();
    passingFilter.put(
        TransformNonVariantSegmentData.FILTER_FIELD,
        ListValue.newBuilder()
            .addValues(Value.newBuilder().setStringValue(TransformNonVariantSegmentData.PASSING_FILTER).build())
            .build());
    VariantCall passingCall = VariantCall.newBuilder().putAllInfo(passingFilter).build();

    Map<String, ListValue> failingFilter = new HashMap<String, ListValue>();
    failingFilter.put(
        TransformNonVariantSegmentData.FILTER_FIELD,
        ListValue.newBuilder()
            .addValues(Value.newBuilder().setStringValue("VQSRTrancheSNP99.90to100.00").build())
            .build());
    VariantCall failingCall = VariantCall.newBuilder().putAllInfo(failingFilter).build();

    DoFnTester<Variant, Variant> filterCallsFn = DoFnTester.of(new FilterCallsFn());
    List filteredVariants = filterCallsFn.processBatch(
        Variant.newBuilder().addAllCalls(Arrays.asList(failingCall)).build(),
        Variant.newBuilder().addAllCalls(Arrays.asList(passingCall)).build());
    assertEquals(1, filteredVariants.size());
  }

  @Test
  public void testAmbiguousVariantCallsFn() {

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

    DoFnTester<Variant, TableRow> formatVariantsFn = DoFnTester.of(new TransformNonVariantSegmentData.FormatVariantsFn(true, false));
    List<TableRow> rows = formatVariantsFn.processBatch(expectedVariant, ambiguousExpectedVariant);
    assertEquals(2, rows.size());
    assertEquals("false", rows.get(0).get(TransformNonVariantSegmentData.HAS_AMBIGUOUS_CALLS_FIELD).toString());
    assertEquals("true", rows.get(1).get(TransformNonVariantSegmentData.HAS_AMBIGUOUS_CALLS_FIELD).toString());
  }

  @Test
  public void testFormatAlt() {

    VariantCall noCall = VariantCall.newBuilder().setCallSetName("noCall").addAllGenotype(Arrays.asList(-1, -1)).build();
    VariantCall noCallRef = VariantCall.newBuilder().setCallSetName("noCallRef").addAllGenotype(Arrays.asList(-1, 0)).build();
    VariantCall noCallAlt = VariantCall.newBuilder().setCallSetName("noCallAlt").addAllGenotype(Arrays.asList(-1, 1)).build();
    VariantCall refMatch = VariantCall.newBuilder().setCallSetName("refMatch").addAllGenotype(Arrays.asList(0, 0)).build();
    VariantCall hetAlt = VariantCall.newBuilder().setCallSetName("hetAlt").addAllGenotype(Arrays.asList(0, 1)).build();
    VariantCall homAlt = VariantCall.newBuilder().setCallSetName("homAlt").addAllGenotype(Arrays.asList(1, 1)).build();
    VariantCall multiAllelic = VariantCall.newBuilder().setCallSetName("multiAllelic").addAllGenotype(Arrays.asList(1, 2)).build();

    Variant snpVariant = Variant.newBuilder()
        .setReferenceBases("T")
        .addAllAlternateBases(Arrays.asList("A", "C"))
        .putAllInfo(FlagVariantsWithAmbiguousCallsFn.NO_AMBIGUOUS_CALLS_INFO)
        .addAllCalls(Arrays.asList(noCall, noCallRef, noCallAlt, refMatch, hetAlt, homAlt, multiAllelic))
        .build();

    // This should have been weeded out earlier in the pipeline, but just checking correctness of calculations here.
    Variant allNoCallsVariant = Variant.newBuilder()
        .setReferenceBases("T")
        .addAllAlternateBases(Arrays.asList("A", "C"))
        .putAllInfo(FlagVariantsWithAmbiguousCallsFn.NO_AMBIGUOUS_CALLS_INFO)
        .addAllCalls(Arrays.asList(noCall, noCall, noCall, noCall, noCall))
        .build();

    Variant noAltGenotypesVariant = Variant.newBuilder()
        .setReferenceBases("T")
        .addAllAlternateBases(Arrays.asList("A", "C"))
        .putAllInfo(FlagVariantsWithAmbiguousCallsFn.NO_AMBIGUOUS_CALLS_INFO)
        .addAllCalls(Arrays.asList(noCall, noCallRef, refMatch))
        .build();

    Variant deletionVariant = Variant.newBuilder()
        .setReferenceBases("AA")
        .addAllAlternateBases(Arrays.asList("A", "C"))
        .putAllInfo(FlagVariantsWithAmbiguousCallsFn.NO_AMBIGUOUS_CALLS_INFO)
        .addAllCalls(Arrays.asList(noCall, noCallRef, noCallAlt, refMatch, hetAlt, homAlt, multiAllelic))
        .build();

    DoFnTester<Variant, TableRow> formatVariantsFn = DoFnTester.of(new TransformNonVariantSegmentData.FormatVariantsFn(true, false));
    List<TableRow> rows = formatVariantsFn.processBatch(snpVariant, allNoCallsVariant, noAltGenotypesVariant, deletionVariant);
    assertEquals(4, rows.size());

    assertEquals(10, rows.get(0).get(TransformNonVariantSegmentData.ALLELE_NUMBER_FIELD));
    assertEquals("[{alternate_bases=A, AC=5, AF=0.5}, {alternate_bases=C, AC=1, AF=0.1}]",
        rows.get(0).get("alt").toString());
    assertEquals(Arrays.asList("refMatch"), rows.get(0).get(TransformNonVariantSegmentData.REF_MATCH_CALLSETS_FIELD));

    assertEquals(0, rows.get(1).get(TransformNonVariantSegmentData.ALLELE_NUMBER_FIELD));
    assertEquals("[{alternate_bases=A, AC=0, AF=0.0}, {alternate_bases=C, AC=0, AF=0.0}]",
        rows.get(1).get("alt").toString());
    assertEquals(Arrays.asList(), rows.get(1).get(TransformNonVariantSegmentData.REF_MATCH_CALLSETS_FIELD));

    assertEquals(3, rows.get(2).get(TransformNonVariantSegmentData.ALLELE_NUMBER_FIELD));
    assertEquals("[{alternate_bases=A, AC=0, AF=0.0}, {alternate_bases=C, AC=0, AF=0.0}]",
        rows.get(2).get("alt").toString());
    assertEquals(Arrays.asList("refMatch"), rows.get(2).get(TransformNonVariantSegmentData.REF_MATCH_CALLSETS_FIELD));

    assertEquals(0, rows.get(3).get(TransformNonVariantSegmentData.ALLELE_NUMBER_FIELD));
    assertEquals("[{alternate_bases=A, AC=5, AF=0.0}, {alternate_bases=C, AC=1, AF=0.0}]",
        rows.get(3).get("alt").toString());
    assertEquals(Arrays.asList("refMatch"), rows.get(3).get(TransformNonVariantSegmentData.REF_MATCH_CALLSETS_FIELD));
  }

  @Test
  public void testFormatCalls() {

    Map depthInfo = new HashMap<String, List<String>>();
    depthInfo.put(
        TransformNonVariantSegmentData.DEPTH_FIELD,
        ListValue.newBuilder()
            .addValues(Value.newBuilder().setStringValue("30").build())
            .build());

    Map dotDepthInfo = new HashMap<String, List<String>>();
    dotDepthInfo.put(
        TransformNonVariantSegmentData.DEPTH_FIELD,
        ListValue.newBuilder()
            .addValues(Value.newBuilder().setStringValue(".").build())
            .build());

    VariantCall callWithValidDepth = VariantCall.newBuilder()
        .setCallSetName("hetAlt").addAllGenotype(Arrays.asList(0, 1))
        .putAllInfo(depthInfo)
        .build();

    VariantCall callWithDotDepth = VariantCall.newBuilder()
        .setCallSetName("homAlt").addAllGenotype(Arrays.asList(1, 1))
        .putAllInfo(dotDepthInfo)
        .build();

    Variant variant = Variant.newBuilder()
        .putAllInfo(FlagVariantsWithAmbiguousCallsFn.NO_AMBIGUOUS_CALLS_INFO)
        .addAllCalls(Arrays.asList(callWithValidDepth, callWithDotDepth))
        .build();

    DoFnTester<Variant, TableRow> formatVariantsFn = DoFnTester.of(new TransformNonVariantSegmentData.FormatVariantsFn(true, false));
    List<TableRow> rows = formatVariantsFn.processBatch(variant);
    assertEquals(1, rows.size());

    assertEquals("[{call_set_name=hetAlt, phaseset=, genotype=[0, 1], genotype_likelihood=[], DP=30},"
        + " {call_set_name=homAlt, phaseset=, genotype=[1, 1], genotype_likelihood=[], DP=null}]",
        rows.get(0).get("call").toString());
  }
}
