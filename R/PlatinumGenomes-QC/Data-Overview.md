<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->

<!-- Copyright 2015 Google Inc. All rights reserved. -->

<!-- Licensed under the Apache License, Version 2.0 (the "License"); -->
<!-- you may not use this file except in compliance with the License. -->
<!-- You may obtain a copy of the License at -->

<!--     http://www.apache.org/licenses/LICENSE-2.0 -->

<!-- Unless required by applicable law or agreed to in writing, software -->
<!-- distributed under the License is distributed on an "AS IS" BASIS, -->
<!-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. -->
<!-- See the License for the specific language governing permissions and -->
<!-- limitations under the License. -->

# Part 1: Data Overview

In Part 1 of the codelab, we perform some queries to acquaint ourselves with the data and determine whether it has any characteristics requiring any additional consideration in the QC checks that follow.

* [Variants](#variants)
* [Non-Variant Segments](#non-variant-segments)
* [Alternative Allele Field](#alternative-allele-field)
* [Genotype Field](#genotype-field)

The following example makes use of [Illumina Platinum Genomes](http://www.illumina.com/platinumgenomes/).  For more detail about how this data was loaded into the Google Genomics API, please see [Google Genomics Public Data](https://cloud.google.com/genomics/data/platinum-genomes).





By default this codelab runs on the Illumina Platinum Genomes Variants. Update the table and change the source of sample information here if you wish to run the queries against a different dataset.

```r
source("./rHelpers/platinumGenomesDataset.R")

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpDataset.R")
```

## Variants

Let's take a look at a few of the [variants within BRCA1 via BigQuery](https://github.com/googlegenomics/getting-started-bigquery/blob/master/RMarkdown/literate-programming-demo.md#data-visualization):

```r
result <- DisplayAndDispatchQuery("./sql/variant-level-data-for-brca1.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
#standardSQL
--
-- Retrieve variant-level information for BRCA1 variants.
--
SELECT
  reference_name,
  start,
  `end`,
  reference_bases AS ref,
  ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat,
  quality,
  ARRAY_TO_STRING(v.filter, ',') AS filters,
  ARRAY_TO_STRING(v.names, ',') AS names,
  ARRAY_LENGTH(v.call) AS num_samples
FROM
  `genomics-public-data.platinum_genomes.variants` v
WHERE
  reference_name IN ('17', 'chr17')
  AND start BETWEEN 41196311 AND 41277499 # per GRCh37
  # Skip non-variant segments.
  AND EXISTS (SELECT alt FROM UNNEST(v.alternate_bases) alt WHERE alt NOT IN ("<NON_REF>", "<*>"))
ORDER BY
  start,
  alt_concat
```
Number of rows returned by this query: **281**.

Displaying the first few rows of the dataframe of results:

|reference_name |    start|      end|ref |alt_concat | quality|filters                                           |names | num_samples|
|:--------------|--------:|--------:|:---|:----------|-------:|:-------------------------------------------------|:-----|-----------:|
|chr17          | 41196407| 41196408|G   |A          |  733.47|PASS                                              |      |           3|
|chr17          | 41196820| 41196823|CTT |C,CT       |  287.18|PASS                                              |      |           1|
|chr17          | 41197273| 41197274|C   |A          | 1011.08|PASS                                              |      |           3|
|chr17          | 41197957| 41197958|G   |T          |  178.48|TruthSensitivityTranche99.90to100.00              |      |           4|
|chr17          | 41198182| 41198183|A   |C          |   98.02|TruthSensitivityTranche99.00to99.90               |      |           1|
|chr17          | 41198186| 41198187|A   |C          |    7.68|TruthSensitivityTranche99.90to100.00,LowGQX,LowQD |      |           4|

These are the variant-level fields common to all variant sets exported to BigQuery from Google Genomics.  There are often dataset-specific variant-level fields as well.  For more information about additional fields, see the schema for the table being queried.  

> In this case, see the Platinum Genomes [variants table schema](https://bigquery.cloud.google.com/table/genomics-public-data:platinum_genomes.variants).

## Non-Variant Segments

Let's take a look at a few non-variant segments within BRCA1:

```r
result <- DisplayAndDispatchQuery("./sql/non-variant-segments-brca1.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
#standardSQL
--
-- Retrieve non-variant segments for BRCA1.
--
SELECT
  call.call_set_name,
  (SELECT STRING_AGG(CAST(gt AS STRING)) from UNNEST(call.genotype) gt) AS genotype,
  reference_name,
  start,
  `end`,
  reference_bases AS ref,
  ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat
FROM
  `genomics-public-data.platinum_genomes.variants` v, v.call call
WHERE
  reference_name IN ('17', 'chr17')
  AND start BETWEEN 41196311 AND 41277499 # per GRCh37
  # Skip all variant sites.
  AND NOT EXISTS (SELECT alt FROM UNNEST(v.alternate_bases) alt WHERE alt NOT IN ("<NON_REF>", "<*>"))
ORDER BY
  start,
  call.call_set_name
LIMIT
  10000
```
Number of rows returned by this query: **8123**.

Displaying the first few rows of the dataframe of results:

|call_set_name |genotype |reference_name |    start|      end|ref |alt_concat |
|:-------------|:--------|:--------------|--------:|--------:|:---|:----------|
|not displayed |0,0      |chr17          | 41196321| 41196381|T   |           |
|not displayed |0,0      |chr17          | 41196349| 41196417|A   |           |
|not displayed |0,0      |chr17          | 41196369| 41196407|T   |           |
|not displayed |0,0      |chr17          | 41196376| 41196621|T   |           |
|not displayed |0,0      |chr17          | 41196381| 41196407|T   |           |
|not displayed |0,0      |chr17          | 41196408| 41196543|G   |           |

When the data contains non-variant segments, for any analyses that require us to know for example _"how many samples do and do not have a particular SNP?"_, we'll need to make sure that the non-variant segments are considered in addition to the variants.

> The source Platinum Genomes data loaded into the Google Genomics API was in [genome VCF](https://sites.google.com/site/gvcftools/home/about-gvcf/gvcf-conventions) (gVCF) format and therefore has non-variant segments.  

Note that Complete Genomics data also includes non-variant segments and requires the same consideration.

If this query were run on a different dataset and returned no rows, then the data only contains variant records.

## Alternative Allele Field

And then let's take a look at the domain and range of values for alternate_bases:

```r
result <- DisplayAndDispatchQuery("./sql/characterize-alts.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
#standardSQL
--
-- Check whether variants are only SNPs and INDELs, with no special characters.
--
SELECT
  COUNT(1) AS number_of_variant_records,
  REGEXP_CONTAINS(alt,
             r'^[ACGT]+$') AS alt_contains_no_special_characters,
  MAX(LENGTH(reference_bases)) AS max_ref_len,
  MAX(LENGTH(alt)) AS max_alt_len
FROM
  `genomics-public-data.platinum_genomes.variants` v, v.alternate_bases alt
WHERE
  # Skip non-variant segments.
  EXISTS (SELECT alt FROM UNNEST(v.alternate_bases) alt WHERE alt NOT IN ("<NON_REF>", "<*>"))
GROUP BY
  alt_contains_no_special_characters
```
Number of rows returned by this query: **1**.

Displaying the first few rows of the dataframe of results:

| number_of_variant_records|alt_contains_no_special_characters | max_ref_len| max_alt_len|
|-------------------------:|:----------------------------------|-----------:|-----------:|
|                  11162053|TRUE                               |          56|          44|

> In the case of Platinum Genomes we see from the query results that there are no special charaters in alternate_bases and the maximum length is ~50 base pairs, so just SNPs and small INDELs.

If this query was run on a different dataset, you may wish to run additional queries to understand the domain and range of possible values in the alternate_bases field (e.g., large deletions coded as `<DEL>`, complex structural variants, etc...)

## Genotype Field

And finally let's take a look at the domain and range of values for genotype:

```r
result <- DisplayAndDispatchQuery("./sql/characterize-genotypes.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
#standardSQL
--
-- Query to show the variety of genotypes.
--
SELECT
  genotype,
  COUNT(genotype) AS genotype_count
FROM (
  SELECT
  (SELECT STRING_AGG(CAST(gt AS STRING)) from UNNEST(call.genotype) gt) AS genotype
  FROM
  `genomics-public-data.platinum_genomes.variants` v, v.call call)
GROUP BY
  genotype
ORDER BY
  genotype_count DESC,
  genotype
```
Number of rows returned by this query: **8**.

Displaying the first few rows of the dataframe of results:

|genotype | genotype_count|
|:--------|--------------:|
|0,0      |      268837097|
|0,1      |       22872404|
|1,1      |       10796056|
|-1       |        3024578|
|0        |        2750341|
|-1,-1    |        1033653|
|1,2      |         230089|
|1        |           7473|


> In the case of Platinum Genomes we see from the query results the variety of genotypes:
>
> * no-calls (the -1 values)
> * genotypes higher than 1 indicating that the data is not strictly bi-allelic
> * genotypes consisting of just a single allele

# Summary

> To summarize attributes we need to consider when working with Platinum Genomes data:
>
> * It has non-variant segments which adds complexity above and beyond [similar examples for the 1,000 Genomes dataset](https://github.com/googlegenomics/bigquery-examples/blob/master/1000genomes/sql/README.md).
> * It is comprised only of SNPs and INDELs (contains no structural variants).
> * The values for `alternate_bases` are just comprised of the letters A,C,G,T (e.g., contains no `<DEL>` values).
> * It contains some single-allele and 1/2 genotypes.

--------------------------------------------------------
_Next_: [Part 2: Data Transformation](./Data-Transformation.md)
