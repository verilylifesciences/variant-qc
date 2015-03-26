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

import java.io.IOException;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.List;

import com.google.api.client.repackaged.com.google.common.base.Preconditions;
import com.google.api.services.bigquery.model.TableFieldSchema;
import com.google.api.services.bigquery.model.TableRow;
import com.google.api.services.bigquery.model.TableSchema;
import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.BigQueryIO;
import com.google.cloud.dataflow.sdk.options.Description;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.options.Validation;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.functions.JoinNonVariantSegmentsWithVariants;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.cloud.genomics.dataflow.utils.GenomicsDatasetOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig.SexChromosomeFilter;
import com.google.cloud.genomics.utils.GenomicsFactory;

/**
 * Sample pipeline that transforms data with non-variant segments (such as data that was in source
 * format Genome VCF (gVCF) or Complete Genomics) to variant-only data with calls from
 * non-variant-segments merged into the variants with which they overlap. The resultant data is
 * emitted to a BigQuery table.
 *
 * This is currently done only for SNP variants. Indels and structural variants are left as-is.
 *
 * The data source is the Google Genomics Variants API. The data sink is BigQuery
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

  /**
   * Options supported by {@link TransformNonVariantSegmentData}.
   * <p>
   * Inherits standard configuration options for Genomics pipelines and datasets.
   */
  private static interface Options extends GenomicsDatasetOptions {
    @Override
    @Description("BigQuery table to write to, specified as "
        + "<project_id>:<dataset_id>.<table_id>. The dataset must already exist.")
    @Validation.Required
    String getOutput();

    @Override
    void setOutput(String value);
  }

  /**
   * Construct the table schema for the output table.
   *
   * @return The schema for the destination table.
   */
  private static TableSchema getTableSchema() {
    List<TableFieldSchema> callFields = new ArrayList<>();
    callFields.add(new TableFieldSchema().setName("call_set_name").setType("STRING"));
    callFields.add(new TableFieldSchema().setName("phaseset").setType("STRING"));
    callFields.add(new TableFieldSchema().setName("genotype").setType("INTEGER")
        .setMode("REPEATED"));
    callFields.add(new TableFieldSchema().setName("genotype_likelihood").setType("FLOAT")
        .setMode("REPEATED"));

    List<TableFieldSchema> fields = new ArrayList<>();
    fields.add(new TableFieldSchema().setName("variant_id").setType("STRING"));
    fields.add(new TableFieldSchema().setName("reference_name").setType("STRING"));
    fields.add(new TableFieldSchema().setName("start").setType("INTEGER"));
    fields.add(new TableFieldSchema().setName("end").setType("INTEGER"));
    fields.add(new TableFieldSchema().setName("reference_bases").setType("STRING"));
    fields.add(new TableFieldSchema().setName("alternate_bases").setType("STRING")
        .setMode("REPEATED"));
    fields.add(new TableFieldSchema().setName("names").setType("STRING")
        .setMode("REPEATED"));
    fields.add(new TableFieldSchema().setName("filter").setType("STRING")
        .setMode("REPEATED"));
    fields.add(new TableFieldSchema().setName("quality").setType("FLOAT"));
    fields.add(new TableFieldSchema().setName("call").setType("RECORD").setMode("REPEATED")
        .setFields(callFields));

    return new TableSchema().setFields(fields);
  }

  /**
   * Prepare the data for writing to BigQuery by building a TableRow object containing data from
   * the variant mapped onto the schema to be used for the destination table.
   */
  static class FormatVariantsFn extends DoFn<Variant, TableRow> {
    @Override
    public void processElement(ProcessContext c) {
      Variant v = c.element();

      List<TableRow> calls = new ArrayList<>();
      for (Call call : v.getCalls()) {
        calls.add(new TableRow()
            .set("call_set_name", call.getCallSetName())
            .set("phaseset", call.getPhaseset())
            .set("genotype", call.getGenotype())
            .set("genotype_likelihood",
                (call.getGenotypeLikelihood() == null) ? new ArrayList<Double>() :
                  call.getGenotypeLikelihood()
                  ));
      }

      TableRow row = new TableRow()
          .set("variant_id", v.getId())
          .set("reference_name", v.getReferenceName())
          .set("start", v.getStart())
          .set("end", v.getEnd())
          .set("reference_bases", v.getReferenceBases())
          .set("alternate_bases", v.getAlternateBases())
          .set("names", (v.getNames() == null) ? new ArrayList<String>() : v.getNames())
          .set("filter", (v.getFilter() == null) ? new ArrayList<String>() : v.getFilter())
          .set("quality", v.getQuality())
          .set("call", calls);
      c.output(row);
    }
  }

  public static void main(String[] args) throws IOException, GeneralSecurityException {
    Options options = PipelineOptionsFactory.fromArgs(args).withValidation().as(Options.class);
    GenomicsDatasetOptions.Methods.validateOptions(options);

    Preconditions.checkState(options.getHasNonVariantSegments(),
        "This job is only valid for data containing non-variant segments. "
            + "Set the --hasNonVariantSegments command line option accordingly.");

    GenomicsFactory.OfflineAuth auth = GenomicsOptions.Methods.getGenomicsAuth(options);
    List<SearchVariantsRequest> requests =
        GenomicsDatasetOptions.Methods.getVariantRequests(options, auth,
            SexChromosomeFilter.INCLUDE_XY);

    Pipeline p = Pipeline.create(options);
    DataflowWorkarounds.registerGenomicsCoders(p);

    PCollection<SearchVariantsRequest> input =
        DataflowWorkarounds.getPCollection(requests, p);

    PCollection<Variant> variants =
        JoinNonVariantSegmentsWithVariants.joinVariantsTransform(input, auth);

    variants.apply(ParDo.of(new FormatVariantsFn())).apply(
        BigQueryIO.Write.to(options.getOutput()).withSchema(getTableSchema())
            .withCreateDisposition(BigQueryIO.Write.CreateDisposition.CREATE_IF_NEEDED)
            .withWriteDisposition(BigQueryIO.Write.WriteDisposition.WRITE_TRUNCATE));

    p.run();
  }
}
