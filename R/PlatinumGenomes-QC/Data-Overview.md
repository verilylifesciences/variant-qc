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






```r
# By default this codelab runs upon the Illumina Platinum Genomes Variants.  
# Change the table here if you wish to run these queries against your own data.
queryReplacements <- list("_THE_TABLE_"="genomics-public-data:platinum_genomes.variants")

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpCGIOnlyDataset.R")
```

## Variants

Let's take a look at a few of the [variants within BRCA1 via BigQuery](https://github.com/googlegenomics/getting-started-bigquery/blob/master/RMarkdown/literate-programming-demo.md#data-visualization):

```r
result <- DisplayAndDispatchQuery("https://raw.githubusercontent.com/googlegenomics/getting-started-bigquery/master/sql/variant-level-data-for-brca1.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Retrieve variant-level information for BRCA1 variants.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
  GROUP_CONCAT(names) WITHIN RECORD AS names,
  COUNT(call.call_set_name) WITHIN RECORD AS num_samples,
FROM
  [genomics-public-data:platinum_genomes.variants]
WHERE
  reference_name = 'chr17'
  AND start BETWEEN 41196311
  AND 41277499
OMIT RECORD IF EVERY(alternate_bases IS NULL)
ORDER BY
  start,
  alternate_bases
```
Number of rows returned by this query: **335**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.0 by xtable 1.7-4 package -->
<!-- Thu Aug 13 16:11:05 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> quality </th> <th> filter </th> <th> names </th> <th> num_samples </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td align="right"> 733.47 </td> <td> PASS </td> <td>  </td> <td align="right">   7 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td align="right"> 41196822 </td> <td> CT </td> <td> C </td> <td align="right"> 63.74 </td> <td> LowQD </td> <td>  </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td align="right"> 41196823 </td> <td> CTT </td> <td> C,CT </td> <td align="right"> 314.59 </td> <td> PASS </td> <td>  </td> <td align="right">   3 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196840 </td> <td align="right"> 41196841 </td> <td> G </td> <td> T </td> <td align="right"> 85.68 </td> <td> TruthSensitivityTranche99.90to100.00,LowQD </td> <td>  </td> <td align="right">   2 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td align="right"> 41197274 </td> <td> C </td> <td> A </td> <td align="right"> 1011.08 </td> <td> PASS </td> <td>  </td> <td align="right">   7 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197938 </td> <td align="right"> 41197939 </td> <td> A </td> <td> AT </td> <td align="right"> 86.95 </td> <td> LowQD </td> <td>  </td> <td align="right">   3 </td> </tr>
   </table>

These are the variant-level fields common to all variant sets exported to BigQuery from Google Genomics.  There are often dataset-specific variant-level fields as well.  For more information about additional fields, see the schema for the table being queried.  

> In this case, see the Platinum Genomes [variants table schema](https://bigquery.cloud.google.com/table/genomics-public-data:platinum_genomes.variants).

## Non-Variant Segments

Let's take a look at a few non-variant segments within BRCA1:

```r
result <- DisplayAndDispatchQuery("./sql/non-variant-segments.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Retrieve non-variant segments for BRCA1, flattening by sample.
SELECT
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
FROM
  [genomics-public-data:platinum_genomes.variants]
WHERE
  reference_name = 'chr17'
  AND start BETWEEN 41196311
  AND 41277499
OMIT RECORD IF SOME(alternate_bases IS NOT NULL)
ORDER BY
  start,
  call.call_set_name
Retrieving data:  3.1s
```
Number of rows returned by this query: **22777**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.0 by xtable 1.7-4 package -->
<!-- Thu Aug 13 16:11:10 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> genotype </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th>  </tr>
  <tr> <td> NA12886 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196313 </td> <td align="right"> 41196746 </td> <td> G </td> <td>  </td> </tr>
  <tr> <td> NA12878 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196321 </td> <td align="right"> 41196381 </td> <td> T </td> <td>  </td> </tr>
  <tr> <td> NA12888 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196322 </td> <td align="right"> 41196356 </td> <td> G </td> <td>  </td> </tr>
  <tr> <td> NA12881 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196337 </td> <td align="right"> 41196620 </td> <td> T </td> <td>  </td> </tr>
  <tr> <td> NA12893 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196339 </td> <td align="right"> 41196489 </td> <td> C </td> <td>  </td> </tr>
  <tr> <td> NA12877 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196349 </td> <td align="right"> 41196417 </td> <td> A </td> <td>  </td> </tr>
   </table>

When the data contains non-variant segments, for any analyses that require us to know for example _"how many samples do and do not have a particular SNP?"_, we'll need to make sure that the non-variant segments are considered in addition to the variants.

> The source Platinum Genomes data loaded into the Google Genomics API was in [genome VCF](https://sites.google.com/site/gvcftools/home/about-gvcf/gvcf-conventions) (gVCF) format and therefore has non-variant segments.  

Note that Complete Genomics data also includes non-variant segments and requires the same consideration.

If this query was run on a different dataset and returned no rows, then the data only contains variant records.

## Alternative Allele Field

And then let's take a look at the domain and range of values for alternate_bases:

```r
result <- DisplayAndDispatchQuery("./sql/characterize-alts.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Variants are only SNPs and INDELs, with no special characters.
SELECT
  COUNT(1) AS number_of_variant_records,
  REGEXP_MATCH(alternate_bases,
    r'^[A,C,G,T]+$') AS alt_contains_no_special_characters,
  MAX(LENGTH(reference_bases)) AS max_ref_len,
  MAX(LENGTH(alternate_bases)) AS max_alt_len
FROM
  [genomics-public-data:platinum_genomes.variants]
OMIT RECORD IF EVERY(alternate_bases IS NULL)
GROUP BY
  alt_contains_no_special_characters
```
Number of rows returned by this query: **1**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.0 by xtable 1.7-4 package -->
<!-- Thu Aug 13 16:11:11 2015 -->
<table border=1>
<tr> <th> number_of_variant_records </th> <th> alt_contains_no_special_characters </th> <th> max_ref_len </th> <th> max_alt_len </th>  </tr>
  <tr> <td align="right"> 12634588 </td> <td> TRUE </td> <td align="right">  56 </td> <td align="right">  47 </td> </tr>
   </table>

> In the case of Platinum Genomes we see from the query results that there are no special charaters in alternate_bases and the maximum length is ~50 base pairs, so just SNPs and small INDELs.

If this query was run on a different dataset, you may wish to run additional queries to understand the domain and range of possible values in the alternate_bases field (e.g., large deletions coded as `<DEL>`, complex structural variants, etc...)

## Genotype Field

And finally let's take a look at the domain and range of values for genotype:

```r
result <- DisplayAndDispatchQuery("./sql/genotypes-brca1.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Query to show the variety of genotypes within BRCA1, such as
# single allele genotypes.
SELECT
  genotype,
  COUNT(genotype) AS genotype_count
FROM (
  SELECT
    GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  FROM
  [genomics-public-data:platinum_genomes.variants]
  WHERE
    reference_name = 'chr17'
    AND start BETWEEN 41196311
    AND 41277499
    )
GROUP BY
  genotype
ORDER BY
  genotype_count DESC
```
Number of rows returned by this query: **7**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.0 by xtable 1.7-4 package -->
<!-- Thu Aug 13 16:11:12 2015 -->
<table border=1>
<tr> <th> genotype </th> <th> genotype_count </th>  </tr>
  <tr> <td> 0,0 </td> <td align="right"> 22519 </td> </tr>
  <tr> <td> 0,1 </td> <td align="right"> 1677 </td> </tr>
  <tr> <td> 0 </td> <td align="right"> 226 </td> </tr>
  <tr> <td> 1,1 </td> <td align="right">  73 </td> </tr>
  <tr> <td> -1 </td> <td align="right">  50 </td> </tr>
  <tr> <td> -1,-1 </td> <td align="right">   5 </td> </tr>
  <tr> <td> 1,2 </td> <td align="right">   4 </td> </tr>
   </table>


> In the case of Platinum Genomes we see from the query results the variety of genotypes just within BRCA1:
* no-calls (the -1 values)
* genotypes higher than 1 indicating that the data is not strictly bi-allelic
* genotypes consisting of just a single allele on an autosome

Remove the WHERE clause to run this query over the entire dataset.

# Summary

> To summarize attributes we need to consider when working with Platinum Genomes data:
> * It has non-variant segments which adds complexity above and beyond [similar examples for the 1,000 Genomes dataset](https://github.com/googlegenomics/bigquery-examples/blob/master/1000genomes/sql/README.md).
> * It is comprised only of SNPs and INDELs (contains no structural variants).
> * The values for `alternate_bases` are just comprised of the letters A,C,G,T (e.g., contains no `<DEL>` values).
> * It contains some single-allele and 1/2 genotypes.

--------------------------------------------------------
_Next_: [Part 2: Data Transformation](./Data-Transformation.md)
