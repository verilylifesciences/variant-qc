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

import com.google.api.client.util.Preconditions;
import com.google.api.services.bigquery.model.TableFieldSchema;
import com.google.api.services.bigquery.model.TableRow;
import com.google.api.services.bigquery.model.TableSchema;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.BigQueryIO;
import com.google.cloud.dataflow.sdk.options.Default;
import com.google.cloud.dataflow.sdk.options.Description;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.options.Validation;
import com.google.cloud.dataflow.sdk.transforms.Aggregator;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.Sum;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.functions.JoinNonVariantSegmentsWithVariants;
import com.google.cloud.genomics.dataflow.readers.VariantStreamer;
import com.google.cloud.genomics.dataflow.utils.CallSetNamesOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.dataflow.utils.ShardOptions;
import com.google.cloud.genomics.utils.OfflineAuth;
import com.google.cloud.genomics.utils.ShardBoundary;
import com.google.cloud.genomics.utils.ShardUtils;
import com.google.cloud.genomics.utils.grpc.MergeNonVariantSegmentsWithSnps;
import com.google.cloud.genomics.utils.grpc.VariantUtils;
import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Iterables;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import com.google.genomics.v1.StreamVariantsRequest;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.google.protobuf.ListValue;
import com.google.protobuf.Value;

import java.io.IOException;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;

/**
 * Sample pipeline that transforms data with non-variant segments (such as data that was in source
 * format Genome VCF (gVCF) or Complete Genomics) to variant-only data with calls from
 * non-variant-segments merged into the variants with which they overlap. The resultant data is
 * emitted to a BigQuery table.
 *
 * The specific details of the merge logic depend upon option --variantMergeStrategy=<Class>
 *
 * This pipeline assumes the call set names in the variant set are unique.
 *
 * The data source is the Google Genomics Variants API. The data sink is BigQuery.
 *
 * <p>
 * The sample could be expanded upon to:
 * <ol>
 * <li>emit additional fields from the variants and calls
 * <li>perform additional data munging
 * <li>write data to a different sink such as Google Cloud Storage or Google Datastore
 * </ol>
 */
public class TransformNonVariantSegmentData {

  private static final Logger LOG = Logger
      .getLogger(TransformNonVariantSegmentData.class.getName());

  public static final String VARIANT_API_FIELDS = "variants(alternateBases,calls,end,filter,id,names,quality,referenceBases,referenceName,start)";
  public static final String HAS_AMBIGUOUS_CALLS_FIELD = "ambiguousCalls";
  public static final String REF_MATCH_CALLSETS_FIELD = "refMatchCallsets";
  public static final String ALT_RECORD_FIELD = "alt";
  // AC : allele count in genotypes, for each ALT allele, in the same order as listed
  public static final String ALLELE_COUNT_FIELD = "AC";
  // AF : allele frequency for each ALT allele in the same order as listed: use this when estimated from primary
  // data, not called genotypes
  public static final String ALLELE_FREQUENCY_FIELD = "AF";
  // AN : total number of alleles in called genotypes
  public static final String ALLELE_NUMBER_FIELD = "AN";
  public static final String DEPTH_FIELD = "DP";
  public static final String FILTER_FIELD = "FILTER";
  public static final String PASSING_FILTER = "PASS";

  /**
   * Options supported by {@link TransformNonVariantSegmentData}.
   * <p>
   * Inherits standard configuration options for pipelines operating on shards of genomic data.
   */
  private static interface Options extends
    // Options for call set names.
    CallSetNamesOptions,
    // Options for calculating over regions, chromosomes, or whole genomes.
    ShardOptions,
    // Options for special handling of data with non-variant segment records.  This
    // is needed since IBS must take into account reference-matches in addition
    // to the variants (unlike other analyses such as PCA).
    JoinNonVariantSegmentsWithVariants.Options {

    @Description("BigQuery table to write to, specified as "
        + "<project_id>:<dataset_id>.<table_id>. The BigQuery dataset must already exist.")
    @Validation.Required
    String getOutputTable();
    void setOutputTable(String value);

    @Description("Whether to append to an existing BigQuery table.")
    @Default.Boolean(false)
    boolean getAppendToTable();
    void setAppendToTable(boolean value);

    @Description("Omit low quality variant calls.  Specifically, exclude any variant calls where call.FILTER != \"PASS\".")
    @Default.Boolean(false)
    boolean getOmitLowQualityCalls();
    void setOmitLowQualityCalls(boolean value);

    @Description("Whether to reduce the size of the final output by only listing the call set names that match the reference"
        + "for each variant.")
    @Default.Boolean(true)
    boolean getSummarizeRefMatchCallSets();
    void setSummarizeRefMatchCallSets(boolean value);
  }

  /**
   * Construct the table schema for the output table.
   *
   * @return The schema for the destination table.
   */
  private static TableSchema getTableSchema(boolean summarizeRefMatches) {
    List<TableFieldSchema> callFields = new ArrayList<>();
    callFields.add(new TableFieldSchema().setName("call_set_name").setType("STRING"));
    callFields.add(new TableFieldSchema().setName("phaseset").setType("STRING"));
    callFields.add(new TableFieldSchema().setName("genotype").setType("INTEGER")
        .setMode("REPEATED"));
    callFields.add(new TableFieldSchema().setName("genotype_likelihood").setType("FLOAT")
        .setMode("REPEATED"));
    callFields.add(new TableFieldSchema().setName(DEPTH_FIELD).setType("INTEGER"));

    List<TableFieldSchema> altFields = new ArrayList<>();
    altFields.add(new TableFieldSchema().setName("alternate_bases").setType("STRING"));
    altFields.add(new TableFieldSchema().setName(ALLELE_COUNT_FIELD).setType("INTEGER"));
    altFields.add(new TableFieldSchema().setName(ALLELE_FREQUENCY_FIELD).setType("FLOAT"));

    List<TableFieldSchema> fields = new ArrayList<>();
    fields.add(new TableFieldSchema().setName("variant_id").setType("STRING"));
    fields.add(new TableFieldSchema().setName("reference_name").setType("STRING"));
    fields.add(new TableFieldSchema().setName("start").setType("INTEGER"));
    fields.add(new TableFieldSchema().setName("end").setType("INTEGER"));
    fields.add(new TableFieldSchema().setName("reference_bases").setType("STRING"));
    fields.add(new TableFieldSchema().setName("names").setType("STRING").setMode("REPEATED"));
    fields.add(new TableFieldSchema().setName("filter").setType("STRING").setMode("REPEATED"));
    fields.add(new TableFieldSchema().setName("quality").setType("FLOAT"));
    fields.add(new TableFieldSchema().setName(ALLELE_NUMBER_FIELD).setType("INTEGER"));
    fields.add(new TableFieldSchema().setName(HAS_AMBIGUOUS_CALLS_FIELD).setType("BOOLEAN"));
    if (summarizeRefMatches) {
      fields.add(new TableFieldSchema().setName(REF_MATCH_CALLSETS_FIELD).setType("STRING").setMode("REPEATED"));
    }
    fields.add(new TableFieldSchema().setName(ALT_RECORD_FIELD).setType("RECORD").setMode("REPEATED")
        .setFields(altFields));
    fields.add(new TableFieldSchema().setName("call").setType("RECORD").setMode("REPEATED")
        .setFields(callFields));

    return new TableSchema().setFields(fields);
  }

  /**
   * Pipeline function to prepare the data for writing to BigQuery, including computing allelic frequency.
   *
   * It builds a TableRow object containing data from the variant mapped onto the schema to be used
   * for the destination table.
   */
  static class FormatVariantsFn extends DoFn<Variant, TableRow> {
    final boolean summarizeRefMatches;
    final boolean computeFrequencyForNonSnps;

    public FormatVariantsFn(boolean summarizeRefMatches, boolean computeFrequencyForNonSnps) {
      this.summarizeRefMatches = summarizeRefMatches;
      this.computeFrequencyForNonSnps = computeFrequencyForNonSnps;
    }

    @Override
    public void processElement(ProcessContext c) {
      Variant v = c.element();

      List<String> refMatchCallsets = new ArrayList<>();
      HashMultiset<Integer> genotypeCount = HashMultiset.create();  // This will typically hold counts for genotypes -1,0,1,2.

      List<TableRow> calls = new ArrayList<>();
      for (VariantCall call : v.getCallsList()) {
        genotypeCount.addAll(call.getGenotypeList());

        if (summarizeRefMatches && Iterables.all(call.getGenotypeList(), Predicates.equalTo(0))) {
          refMatchCallsets.add(call.getCallSetName());
        } else {
          String depth = null;
          if (null != call.getInfo().get(DEPTH_FIELD)) {
            String value = call.getInfo().get(DEPTH_FIELD).getValues(0).getStringValue();
            if (!(null == value || value.equals("."))) {
              depth = value;
            }
          }
          calls.add(new TableRow()
              .set("call_set_name", call.getCallSetName())
              .set("phaseset", call.getPhaseset())
              .set("genotype", call.getGenotypeList())
              .set(
                  "genotype_likelihood",
                  (call.getGenotypeLikelihoodList() == null) ? new ArrayList<Double>() : call
                      .getGenotypeLikelihoodList())
              .set(DEPTH_FIELD, depth));
        }
      }

      // Compute AN/AC/AF for SNPs.  Note that no-calls (genotype -1) are excluded from AN.
      List<TableRow> alts = new ArrayList<>();
      int numAlts = v.getAlternateBasesCount();
      int alleleNumber = 0;
      if (computeFrequencyForNonSnps || VariantUtils.IS_SNP.apply(v)) {
        // This for loop iterates over both 0 (reference allele) and over the 1-based alternate allele values.
        for (int i = 0; i <= numAlts; i++) {
          alleleNumber += genotypeCount.count(i);
        }
      }
      for (int j = 0; j < numAlts; j++) {
        alts.add(new TableRow()
            .set("alternate_bases", v.getAlternateBases(j))
            .set(ALLELE_COUNT_FIELD, genotypeCount.count(j + 1))
            .set(ALLELE_FREQUENCY_FIELD,
                (0 < alleleNumber) ? (genotypeCount.count(j + 1) / (double) alleleNumber) : 0.0));
      }

      TableRow row =
          new TableRow()
              .set("variant_id", v.getId())
              .set("reference_name", v.getReferenceName())
              .set("start", v.getStart())
              .set("end", v.getEnd())
              .set("reference_bases", v.getReferenceBases())
              .set("alt", alts)
              .set("names", (v.getNamesList() == null) ? new ArrayList<String>() : v.getNamesList())
              .set("filter", (v.getFilterList() == null) ? new ArrayList<String>() : v.getFilterList())
              .set("quality", v.getQuality())
              .set(ALLELE_NUMBER_FIELD, alleleNumber)
              .set(HAS_AMBIGUOUS_CALLS_FIELD, v.getInfo().get(HAS_AMBIGUOUS_CALLS_FIELD).getValues(0).getStringValue())
              .set("call", calls);
      if (summarizeRefMatches) {
        row.set(REF_MATCH_CALLSETS_FIELD, refMatchCallsets);
      }

      c.output(row);
    }
  }

  /**
   * Pipeline function to filter out low quality calls.
   */
  public static final class FilterCallsFn extends DoFn<Variant, Variant> {
    final Predicate<Value> isPassing = Predicates.equalTo(Value.newBuilder().setStringValue(PASSING_FILTER).build());

    @Override
    public void processElement(ProcessContext context) {
      Variant variant = context.element();

      // Don't filter non-variant segments.
      if(VariantUtils.IS_NON_VARIANT_SEGMENT.apply(variant)) {
        context.output(variant);
        return;
      }

      // We may have variants without calls if any callsets were deleted from the variant set.
      // Skip those.
      if (null == variant.getCallsList() || variant.getCallsList().isEmpty()) {
        return;
      }

      List<VariantCall> filteredCalls =
          Lists.newArrayList(Iterables.filter(variant.getCallsList(), new Predicate<VariantCall>() {
            @Override
            public boolean apply(VariantCall call) {
              ListValue filters = call.getInfo().get(FILTER_FIELD);
              if(null != filters && Iterables.any(filters.getValuesList(), isPassing)) {
                return true;
              }
              return false;
            }
          }));

      // After filtering, the variant may no longer have any calls.  Skip empty variants.
      if (filteredCalls.isEmpty()) {
        return;
      }

      context.output(Variant.newBuilder(variant).clearCalls().addAllCalls(filteredCalls).build());
    }
  }

  /**
   * Pipeline function to flag any variants with more than one call for a particular callSetName.
   *
   * We don't see this for the tidy test datasets such as Platinum Genomes. But in practice, we have
   * seen datasets with the same individual sequenced twice, mistakes in data conversions prior to
   * this step, etc... So this function flags the particular variants with ambiguous calls for the
   * same individual, and also updates the UI with a count of such variants.
   */
  public static final class FlagVariantsWithAmbiguousCallsFn extends DoFn<Variant, Variant> {

    final Aggregator<Long, Long> variantsWithAmbiguousCallsCount = createAggregator("Number of variants containing ambiguous calls", new Sum.SumLongFn());
    public static final Map HAS_AMBIGUOUS_CALLS_INFO;
    public static final Map NO_AMBIGUOUS_CALLS_INFO;

    static {
      HAS_AMBIGUOUS_CALLS_INFO = new HashMap<String, List<String>>();
      HAS_AMBIGUOUS_CALLS_INFO.put(
          HAS_AMBIGUOUS_CALLS_FIELD,
          ListValue.newBuilder()
              .addValues(Value.newBuilder().setStringValue(Boolean.toString(Boolean.TRUE)).build())
              .build());
      NO_AMBIGUOUS_CALLS_INFO = new HashMap<String, List<String>>();
      NO_AMBIGUOUS_CALLS_INFO
          .put(
              HAS_AMBIGUOUS_CALLS_FIELD,
              ListValue
                  .newBuilder()
                  .addValues(
                      Value.newBuilder().setStringValue(Boolean.toString(Boolean.FALSE)).build())
                  .build());
    }

    @Override
    public void processElement(ProcessContext context) {
      Variant variant = context.element();

      // We may have variants without calls if any callsets were deleted from the variant set and/or
      // a filtering step removed all calls. Omit these variants from our output.
      if(null == variant.getCallsList() || variant.getCallsList().isEmpty()) {
        return;
      }

      // Gather calls together for the same callSetName.
      ListMultimap<String, VariantCall> indexedCalls =
          Multimaps.index(variant.getCallsList(), new Function<VariantCall, String>() {
            @Override
            public String apply(final VariantCall c) {
              return c.getCallSetName();
            }
          });

      // Identify and count variants with multiple calls per callSetName.
      boolean isAmbiguous = false;
      for (Entry<String, Collection<VariantCall>> entry : indexedCalls.asMap().entrySet()) {
        if (1 < entry.getValue().size()) {
          LOG.warning("Variant " + variant.getId() + " contains ambiguous calls for at least one genome: "
              + entry.getValue().iterator().next().getCallSetName());
          isAmbiguous = true;
          variantsWithAmbiguousCallsCount.addValue(1l);
          break;  // No need to look for additional ambiguity; one instance is enough to warrant the flag.
        }
      }

      // Add the flag to the variant.
      context.output(Variant.newBuilder(variant)
          .putAllInfo(isAmbiguous ? HAS_AMBIGUOUS_CALLS_INFO : NO_AMBIGUOUS_CALLS_INFO)
          .build());
    }
  }

  public static void main(String[] args) throws IOException, GeneralSecurityException {
    // Register the options so that they show up via --help.
    PipelineOptionsFactory.register(TransformNonVariantSegmentData.Options.class);
    TransformNonVariantSegmentData.Options options =
        PipelineOptionsFactory.fromArgs(args).withValidation()
            .as(TransformNonVariantSegmentData.Options.class);

    Preconditions.checkState(options.getHasNonVariantSegments(),
        "This job is only valid for data containing non-variant segments. "
            + "Set the --hasNonVariantSegments command line option accordingly.");

    // Set up the prototype request and auth.
    StreamVariantsRequest prototype = CallSetNamesOptions.Methods.getRequestPrototype(options);
    final OfflineAuth auth = GenomicsOptions.Methods.getGenomicsAuth(options);

    List<StreamVariantsRequest> requests = options.isAllReferences() ?
        ShardUtils.getVariantRequests(prototype, ShardUtils.SexChromosomeFilter.INCLUDE_XY,
            options.getBasesPerShard(), auth) :
          ShardUtils.getVariantRequests(prototype, options.getBasesPerShard(), options.getReferences());

    Pipeline p = Pipeline.create(options);

    // Create a collection of data with non-variant segments omitted but calls from overlapping
    // non-variant segments added to SNPs and write them to BigQuery.
    PCollection<Variant> variants = p.begin()
        .apply(Create.of(requests))
        .apply(new VariantStreamer(auth, ShardBoundary.Requirement.STRICT, VARIANT_API_FIELDS));

    PCollection<Variant> filteredVariants =
        options.getOmitLowQualityCalls() ? variants.apply(ParDo.of(new FilterCallsFn())) : variants;

    filteredVariants.apply(new JoinNonVariantSegmentsWithVariants.BinShuffleAndCombineTransform())
        .apply(ParDo.of(new FlagVariantsWithAmbiguousCallsFn()))
        .apply(ParDo.of(new FormatVariantsFn(options.getSummarizeRefMatchCallSets(),
            !options.getVariantMergeStrategy().equals(MergeNonVariantSegmentsWithSnps.class))))
        .apply(
            BigQueryIO.Write.to(options.getOutputTable()).withSchema(getTableSchema(options.getSummarizeRefMatchCallSets()))
                .withCreateDisposition(BigQueryIO.Write.CreateDisposition.CREATE_IF_NEEDED)
                .withWriteDisposition(options.getAppendToTable()
                    ? BigQueryIO.Write.WriteDisposition.WRITE_APPEND : BigQueryIO.Write.WriteDisposition.WRITE_EMPTY));

    p.run();
  }
}
