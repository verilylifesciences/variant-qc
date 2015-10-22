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

import java.util.Arrays;

import org.hamcrest.CoreMatchers;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import com.google.cloud.dataflow.sdk.transforms.DoFnTester;
import com.google.cloud.genomics.examples.TransformNonVariantSegmentData.FilterCallsFn;
import com.google.cloud.genomics.examples.TransformNonVariantSegmentData.FlagVariantsWithAmbiguousCallsFn;
import com.google.common.collect.ImmutableSet;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;

@RunWith(JUnit4.class)
public class TransformNonVariantSegmentDataTest {

  @Test
  public void testFilterVariantCallsFn() {

    DoFnTester<Variant, Variant> filterCallsFn =
        DoFnTester.of(new FilterCallsFn(ImmutableSet.of("filterMeOut", "alsoFilterMeOut")));

    VariantCall call1 = VariantCall.newBuilder().setCallSetName("filterMeOut").build();
    VariantCall call2 = VariantCall.newBuilder().setCallSetName("keepMe").build();
    VariantCall call3 = VariantCall.newBuilder().setCallSetName("alsoFilterMeOut").build();
    Variant inputVariant =
        Variant.newBuilder().addAllCalls(Arrays.asList(call1, call2, call3)).build();

    Variant expectedVariant = Variant.newBuilder().addAllCalls(Arrays.asList(call2)).build();

    Assert.assertThat(filterCallsFn.processBatch(inputVariant),
        CoreMatchers.allOf(CoreMatchers.hasItems(expectedVariant)));
  }

  @Test
  public void testFlagVariantsWithAmbiguousVariantCallsFn() {

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
  }

}
