package com.google.cloud.genomics.examples;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.hamcrest.CoreMatchers;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.dataflow.sdk.transforms.DoFnTester;
import com.google.common.collect.ImmutableSet;

@RunWith(JUnit4.class)
public class TransformNonVariantSegmentDataTest {

  @Test
  public void testFilterCallsFn() {

    DoFnTester<Variant, Variant> filterCallsFn =
        DoFnTester.of(new TransformNonVariantSegmentData.FilterCallsFn(ImmutableSet.of("filterMeOut", 
            "alsoFilterMeOut")));

    Call call1 = new Call().setCallSetName("filterMeOut");
    Call call2 = new Call().setCallSetName("keepMe");
    Call call3 = new Call().setCallSetName("alsoFilterMeOut");
    Variant inputVariant = new Variant().setCalls(Arrays.asList(call1, call2, call3));
    
    Variant expectedVariant = new Variant().setCalls(Arrays.asList(call2));
    
    Assert.assertThat(filterCallsFn.processBatch(inputVariant),
        CoreMatchers.hasItems(expectedVariant));
  }

  @Test
  public void testFlagVariantsWithAmbiguousCallsFn() {

    DoFnTester<Variant, Variant> flagVariantsFn =
        DoFnTester.of(new TransformNonVariantSegmentData.FlagVariantsWithAmbiguousCallsFn());

    Call call1 = new Call().setCallSetName("sample1").setGenotype(Arrays.asList(0,1));
    Call call2a = new Call().setCallSetName("sample2").setGenotype(Arrays.asList(0,1));
    Call call2b = new Call().setCallSetName("sample2").setGenotype(Arrays.asList(-1,-1));
    
    Variant inputVariant = new Variant().setCalls(Arrays.asList(call1, call2a));
    Variant ambiguousInputVariant = new Variant().setCalls(Arrays.asList(call1, call2a, call2b));
    
    Variant expectedVariant = new Variant().setCalls(Arrays.asList(call1, call2a))
        .setInfo(new HashMap<String, List<String>>());
        expectedVariant.getInfo().put(TransformNonVariantSegmentData.HAS_AMBIGUOUS_CALLS_FIELD,
            Arrays.asList(Boolean.FALSE.toString()));
    Variant ambiguousExpectedVariant = new Variant().setCalls(Arrays.asList(call1, call2a, call2b))
        .setInfo(new HashMap<String, List<String>>());
    ambiguousExpectedVariant.getInfo().put(TransformNonVariantSegmentData.HAS_AMBIGUOUS_CALLS_FIELD,
            Arrays.asList(Boolean.TRUE.toString()));

    Assert.assertThat(flagVariantsFn.processBatch(inputVariant, ambiguousInputVariant),
        CoreMatchers.hasItems(expectedVariant, ambiguousExpectedVariant));
  }

}
