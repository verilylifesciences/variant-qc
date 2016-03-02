# Part 2: Data Transformation

In Part 2 of the codelab, we perform a data transformation to make it more amenable to our QC analyses in Parts 3 and 4 of this codelab.

* [Overview](#overview)
* [Motivation](#motivation)
* [Running the Cluster Compute Job](#running-the-cluster-compute-job)
* [Results](#results)
* [Optional: modify the Cluster Compute Job](#optional-modify-the-cluster-compute-job)

## Overview

We refer to the original table as the "genome calls" table. It contains *all* reference calls and variant calls. To facilitate variant-centric analysis, we generate a second table, the "multi-sample variants" table.

The multi-sample variants table resembles a multi-sample VCF file. In this table:

* every variant record includes calls for *all* callsets
* variants which contained *only reference calls for all callsets* are omitted

## Motivation

Data from source files in [genome VCF](https://sites.google.com/site/gvcftools/home/about-gvcf/gvcf-conventions) (gVCF) format or in Complete Genomics format can be challenging to query due to the presence of non-variant segment records.

For example to lookup [rs9536314](http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=rs9536314) in the Klotho gene, the `WHERE` clause
```
    WHERE
      reference_name = 'chr13'
      AND start = 33628137
```
becomes
```
    WHERE
      reference_name = 'chr13'
      AND start <= 33628137
      AND end >= 33628138
```
to capture not only that variant, but any other records that overlap that genomic position.

Suppose we want to calculate an aggregate for a particular variant, such as the number of samples with the variant on one or both alleles and of samples that match the reference?  The WHERE clause above will do the trick.  But then suppose we want to do this for all SNPs in our dataset?  There are [a few ways to do this](https://github.com/googlegenomics/bigquery-examples/tree/master/pgp/data-stories/schema-comparisons#motivation). In this codelab we will use a cluster computing job to transform data with non-variant segments to variant-only data with calls from non-variant-segments merged into the variants with which they overlap. This is currently done only for SNP variants. Indels and structural variants are left as-is.  

> If you are working with Platinum Genomes its okay to skip ahead to [Part 3: Sample-Level QC](./Sample-Level-QC.md) which makes use of previously computed results from this cluster compute job.

If you are running this codelab against data that does not contain non-variant segments, this data conversion is unnecessary.  Just use the same BigQuery table for both `_GENOME_CALL_TABLE_` and `_MULTISAMPLE_VARIANT_TABLE_`.

## Running the Cluster Compute Job

We use data from variant set [3049512673186936334](https://cloud.google.com/genomics/data/platinum-genomes) as the data source for this job.

* **Java** [Use Cloud Dataflow to transform data with Non-Variant Segments](../../Java/PlatinumGenomes-variant-transformation)

## Results

We have now created table [google.com:biggene:platinum_genomes.multisample_variants](https://bigquery.cloud.google.com/table/google.com:biggene:platinum_genomes.multisample_variants?pli=1)

In the analyses that follow in this codelab, sometimes we work with the original data and sometimes we work with the transformed variant-only data.

## Optional: modify the Cluster Compute Job

Some ideas:

* Add additional fields from the variant to the BigQuery schema.
* merge calls in records with 1/2 genotypes with the same variant found in records with 0/1 genotypes

Reference Name | Start     | End       | Reference Bases | Alternate Bases
---------------|-----------|-----------|-----------------|-----------------
chr6           | 120458771 | 120458773 |TA               |TAA
chr6           | 120458771 | 120458773 |TA               |TAA,T

--------------------------------------------------------
_Next_: [Part 3: Sample-Level QC](./Sample-Level-QC.md)

