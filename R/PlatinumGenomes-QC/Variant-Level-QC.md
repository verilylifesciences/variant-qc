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

# Part 4: Variant-Level QC





In Part 4 of the codelab, we perform some quality control analyses that could help to identify any problematic variants which should be excluded from further analysis.  The appropriate cut off thresholds will depend upon the input dataset and/or other factors.

* [Ti/Tv by Genomic Window](#titv-by-genomic-window)
* [Ti/Tv by Alternate Allele Counts](#titv-by-alternate-allele-counts)
* [Ti/Tv by Depth](#titv-by-depth)
* [Missingness Rate](#missingness-rate)
* [Hardy-Weinberg Equilibrium](#hardy-weinberg-equilibrium)
* [Heterozygous Haplotype](#heterozygous-haplotype)

By default this codelab runs upon the Illumina Platinum Genomes Variants. Update the table and change the source of sample information here if you wish to run the queries against a different dataset.

```r
source("./rHelpers/platinumGenomesDataset.R")

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpDataset.R")
```

## Ti/Tv by Genomic Window

Check whether the ratio of transitions vs. transversions in SNPs appears to be reasonable in each window of genomic positions.  This query may help identify problematic regions.


```r
result <- DisplayAndDispatchQuery("./sql/ti-tv-by-genomic-window.sql",
                                  project=project,
                                  replacements=c("@WINDOW_SIZE"="100000",
                                                 queryReplacements))
```

```
# Compute the Ti/Tv ratio for variants within genomic windows.
WITH filtered_snp_calls AS (
  SELECT
    reference_name,
    CAST(FLOOR(start / 100000) AS INT64) AS genomic_window,
    CONCAT(reference_bases, '->', alternate_bases[ORDINAL(1)]) AS mutation
  FROM
    `genomics-public-data.platinum_genomes.variants` v
  WHERE
    # Only include biallelic snps with at least one passing variant call.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    AND EXISTS (SELECT gt FROM UNNEST(v.call) AS call, UNNEST(call.genotype) AS gt WHERE gt > 0
      AND EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft IN ('PASS', '.')))
),

mutation_type_counts AS (
  SELECT
    reference_name,
    genomic_window,
    SUM(CAST(mutation IN ('A->G', 'G->A', 'C->T', 'T->C') AS INT64)) AS transitions,
    SUM(CAST(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                          'A->T', 'T->A', 'C->G', 'G->C') AS INT64)) AS transversions,
    COUNT(mutation) AS num_variants_in_group
  FROM filtered_snp_calls
  GROUP BY
    reference_name,
    genomic_window
)

SELECT
  reference_name,
  genomic_window * 100000 AS window_start,
  transitions,
  transversions,
  transitions/transversions AS titv,
  num_variants_in_group
FROM mutation_type_counts
WHERE
  transversions > 0
ORDER BY
  reference_name,
  window_start

Retrieving data:  9.4s
Retrieving data: 19.5s
```
Number of rows returned by this query: **28354**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:51:43 2016 -->
<table border=1>
<tr> <th> reference_name </th> <th> window_start </th> <th> transitions </th> <th> transversions </th> <th> titv </th> <th> num_variants_in_group </th>  </tr>
  <tr> <td> chr1 </td> <td align="right"> 500000 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right"> 1.00 </td> <td align="right">   2 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 700000 </td> <td align="right">  51 </td> <td align="right">  31 </td> <td align="right"> 1.65 </td> <td align="right">  82 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 800000 </td> <td align="right"> 162 </td> <td align="right">  67 </td> <td align="right"> 2.42 </td> <td align="right"> 229 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 900000 </td> <td align="right"> 212 </td> <td align="right">  73 </td> <td align="right"> 2.90 </td> <td align="right"> 285 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 1000000 </td> <td align="right"> 106 </td> <td align="right">  45 </td> <td align="right"> 2.36 </td> <td align="right"> 151 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 1100000 </td> <td align="right"> 151 </td> <td align="right">  38 </td> <td align="right"> 3.97 </td> <td align="right"> 189 </td> </tr>
   </table>

Visualizing the results:

```r
ggplot(filter(result, reference_name %in% c("1", "chr1")), aes(x=window_start, y=titv)) +
  geom_point() +
  stat_smooth() +
  scale_x_continuous(labels=comma) +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
  ggtitle("Ti/Tv by 100,000 base pair windows")
```

<img src="./platinum_genomes/figure/titvByWindow-1.png" title="plot of chunk titvByWindow" alt="plot of chunk titvByWindow" style="display: block; margin: auto;" />

## Ti/Tv by Alternate Allele Counts

Check whether the ratio of transitions vs. transversions in SNPs appears to be resonable across the range of rare variants to common variants.  This query may help to identify problems with rare or common variants.


```r
result <- DisplayAndDispatchQuery("./sql/ti-tv-by-alternate-allele-count.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Compute the Ti/Tv ratio for variants binned by alternate allele count.
WITH filtered_snp_calls AS (
  SELECT
    ( -- Compute the allele count for this site of variation.
      SELECT COUNTIF(gt = 1)
      FROM UNNEST(v.call) AS call, UNNEST(call.genotype) AS gt
      WHERE
        # Skip homozygous reference calls, no-calls, and non-passing variants.
        EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
        AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
        AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
    ) AS AC,
    CONCAT(reference_bases, '->', alternate_bases[ORDINAL(1)]) AS mutation
  FROM
    `genomics-public-data.platinum_genomes.variants` v
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
),

mutation_type_counts AS (
  SELECT
    AC,
    SUM(CAST(mutation IN ('A->G', 'G->A', 'C->T', 'T->C') AS INT64)) AS transitions,
    SUM(CAST(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                          'A->T', 'T->A', 'C->G', 'G->C') AS INT64)) AS transversions,
    COUNT(mutation) AS num_variants_in_group
  FROM filtered_snp_calls
  GROUP BY
    AC
)

SELECT
  AC,
  transitions,
  transversions,
  transitions/transversions AS titv,
  num_variants_in_group
FROM mutation_type_counts
WHERE
  transversions > 0 AND AC > 0
ORDER BY
  AC DESC
```
Number of rows returned by this query: **34**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:51:49 2016 -->
<table border=1>
<tr> <th> AC </th> <th> transitions </th> <th> transversions </th> <th> titv </th> <th> num_variants_in_group </th>  </tr>
  <tr> <td align="right">  34 </td> <td align="right"> 300985 </td> <td align="right"> 142987 </td> <td align="right"> 2.10 </td> <td align="right"> 443972 </td> </tr>
  <tr> <td align="right">  33 </td> <td align="right"> 96846 </td> <td align="right"> 45222 </td> <td align="right"> 2.14 </td> <td align="right"> 142068 </td> </tr>
  <tr> <td align="right">  32 </td> <td align="right"> 60764 </td> <td align="right"> 29492 </td> <td align="right"> 2.06 </td> <td align="right"> 90256 </td> </tr>
  <tr> <td align="right">  31 </td> <td align="right"> 26838 </td> <td align="right"> 12751 </td> <td align="right"> 2.10 </td> <td align="right"> 39589 </td> </tr>
  <tr> <td align="right">  30 </td> <td align="right"> 15224 </td> <td align="right"> 7442 </td> <td align="right"> 2.05 </td> <td align="right"> 22666 </td> </tr>
  <tr> <td align="right">  29 </td> <td align="right"> 13632 </td> <td align="right"> 6442 </td> <td align="right"> 2.12 </td> <td align="right"> 20074 </td> </tr>
   </table>

Visualizing the results:

```r
ggplot(result, aes(x=AC, y=titv)) +
  geom_point() +
  stat_smooth() +
  scale_x_continuous(labels=comma) +
  xlab("Total Number of Sample Alleles with the Variant") +
  ylab("Ti/Tv") +
  ggtitle("Ti/Tv by Alternate Allele Count")
```

<img src="./platinum_genomes/figure/titvByAlt-1.png" title="plot of chunk titvByAlt" alt="plot of chunk titvByAlt" style="display: block; margin: auto;" />

## Ti/Tv by Depth

Visualize the ratio of transitions vs. transversions by depth of coverage.


```r
query <- "./sql/ti-tv-by-depth.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Transition/Transversion Ratio by Depth of Coverage.
WITH filtered_snp_calls AS (
  SELECT
    call.call_set_name,
    call.DP,
    CONCAT(reference_bases, '->', alternate_bases[ORDINAL(1)]) AS mutation
  FROM
    `genomics-public-data.platinum_genomes.variants` v, v.call call
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    # Skip homozygous reference calls, no-calls, and non-passing variants.
    AND EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
),

mutation_type_counts AS (
  SELECT
    call_set_name,
    DP,
    SUM(CAST(mutation IN ('A->G', 'G->A', 'C->T', 'T->C') AS INT64)) AS transitions,
    SUM(CAST(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                          'A->T', 'T->A', 'C->G', 'G->C') AS INT64)) AS transversions,
    COUNT(mutation) AS num_variants_in_group
  FROM filtered_snp_calls
  WHERE
    DP IS NOT NULL
  GROUP BY
    call_set_name,
    DP
)

SELECT
  call_set_name,
  DP AS average_depth,
  transitions,
  transversions,
  transitions/transversions AS titv,
  num_variants_in_group
FROM mutation_type_counts
WHERE
  transversions > 0
ORDER BY
  call_set_name,
  average_depth
```

<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:51:55 2016 -->
<table border=1>
<tr> <th> call_set_name </th> <th> average_depth </th> <th> transitions </th> <th> transversions </th> <th> titv </th> <th> num_variants_in_group </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right">   3 </td> <td align="right">   4 </td> <td align="right">   8 </td> <td align="right"> 0.50 </td> <td align="right">  12 </td> </tr>
  <tr> <td> NA12877 </td> <td align="right">   4 </td> <td align="right">   7 </td> <td align="right">  10 </td> <td align="right"> 0.70 </td> <td align="right">  17 </td> </tr>
  <tr> <td> NA12877 </td> <td align="right">   5 </td> <td align="right">  12 </td> <td align="right">  10 </td> <td align="right"> 1.20 </td> <td align="right">  22 </td> </tr>
  <tr> <td> NA12877 </td> <td align="right">   6 </td> <td align="right">   8 </td> <td align="right">  11 </td> <td align="right"> 0.73 </td> <td align="right">  19 </td> </tr>
  <tr> <td> NA12877 </td> <td align="right">   7 </td> <td align="right">   8 </td> <td align="right">  11 </td> <td align="right"> 0.73 </td> <td align="right">  19 </td> </tr>
  <tr> <td> NA12877 </td> <td align="right">   8 </td> <td align="right">  11 </td> <td align="right">   6 </td> <td align="right"> 1.83 </td> <td align="right">  17 </td> </tr>
   </table>


```r
ggplot(result, aes(x=average_depth, y=titv, color=call_set_name)) + 
  geom_point() +
  ggtitle("Ti/Tv Ratio By Depth") +
  xlab("Coverage Depth") + 
  ylab("Ti/Tv")
```

<img src="./platinum_genomes/figure/titv-by-depth-1.png" title="plot of chunk titv-by-depth" alt="plot of chunk titv-by-depth" style="display: block; margin: auto;" />

## Missingness Rate

For each variant, compute the missingness rate.  This query can be used to identify variants with a poor call rate.


```r
sortAndLimit <- "ORDER BY missingness_rate DESC, reference_name, start, reference_bases, alt_concat LIMIT 1000"
result <- DisplayAndDispatchQuery("./sql/variant-level-missingness.sql",
                                  project=project,
                                  replacements=c(
                                    "-- Optionally add a clause here to constrain the results."=sortAndLimit,
                                                 queryReplacements))
```

```
# Compute the ratio of no-calls for each variant.
WITH variant_missingness AS (
  SELECT
    reference_name,
    start,
    `end`,
    reference_bases,
    ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat,
    (SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt)) FROM v.call) AS all_calls,
    (SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt WHERE gt < 0)) FROM v.call) AS no_calls
  FROM
    `google.com:biggene.platinum_genomes.expanded_variants` v
)

SELECT
  reference_name,
  start,
  `end`,
  reference_bases,
  alt_concat,
  no_calls,
  all_calls,
  no_calls/all_calls AS missingness_rate
FROM variant_missingness
WHERE
  all_calls > 0
ORDER BY missingness_rate DESC, reference_name, start, reference_bases, alt_concat LIMIT 1000
```
Number of rows returned by this query: **1000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:52:02 2016 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alt_concat </th> <th> no_calls </th> <th> all_calls </th> <th> missingness_rate </th>  </tr>
  <tr> <td> chr1 </td> <td align="right"> 723799 </td> <td align="right"> 723800 </td> <td> G </td> <td> C </td> <td align="right">  17 </td> <td align="right">  17 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 823721 </td> <td align="right"> 823722 </td> <td> C </td> <td> A </td> <td align="right">  17 </td> <td align="right">  17 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 825391 </td> <td align="right"> 825392 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">  17 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 825393 </td> <td align="right"> 825394 </td> <td> A </td> <td> T </td> <td align="right">  17 </td> <td align="right">  17 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 867994 </td> <td align="right"> 867995 </td> <td> T </td> <td> G </td> <td align="right">  17 </td> <td align="right">  17 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 867996 </td> <td align="right"> 867997 </td> <td> C </td> <td> G </td> <td align="right">  17 </td> <td align="right">  17 </td> <td align="right"> 1.00 </td> </tr>
   </table>

## Hardy-Weinberg Equilibrium

For each variant, compute the expected versus observed relationship between allele frequencies and genotype frequencies per the Hardy-Weinberg Equilibrium.


```r
sortAndLimit <- "ORDER BY ChiSq DESC, reference_name, start, alt LIMIT 1000"
result <- DisplayAndDispatchQuery("./sql/hardy-weinberg.sql",
                                  project=project,
                                  replacements=c("-- Optionally add a clause here to constrain the results."=sortAndLimit,
                                                 queryReplacements))
```

```
# The following query computes the Hardy-Weinberg equilibrium for variants.
WITH variants AS (
  SELECT
    reference_name,
    start,
    `end`,
    reference_bases,
    v.alternate_bases[ORDINAL(1)] AS alt,
    -- Within each variant count the number of HOM_REF/HOM_ALT/HET samples.
    (SELECT SUM(CAST((SELECT LOGICAL_AND(gt = 0)
      FROM UNNEST(call.genotype) gt) AS INT64)) FROM v.call) AS HOM_REF,
    (SELECT SUM(CAST((SELECT LOGICAL_AND(gt = 1)
      FROM UNNEST(call.genotype) gt) AS INT64)) FROM v.call) AS HOM_ALT,
    (SELECT SUM(CAST((SELECT LOGICAL_OR(gt = 0) AND LOGICAL_OR(gt = 1)
      FROM UNNEST(call.genotype) gt) AS INT64)) FROM v.call) AS HET
  FROM
    `google.com:biggene.platinum_genomes.expanded_variants` v
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
),

observations AS (
  SELECT
    reference_name,
    start,
    reference_bases,
    alt,
    HOM_REF AS OBS_HOM1,
    HET AS OBS_HET,
    HOM_ALT AS OBS_HOM2,
    HOM_REF + HET + HOM_ALT AS SAMPLE_COUNT
  FROM variants
),

expectations AS (
  SELECT
    reference_name,
    start,
    reference_bases,
    alt,
    OBS_HOM1,
    OBS_HET,
    OBS_HOM2,

    # Expected AA
    # p^2
    # ((COUNT(AA) + (COUNT(Aa)/2) /
    #  SAMPLE_COUNT) ^ 2) * SAMPLE_COUNT
    POW((OBS_HOM1 + (OBS_HET/2)) /
      SAMPLE_COUNT, 2) * SAMPLE_COUNT
      AS E_HOM1,

    # Expected Aa
    # 2pq
    # 2 * (COUNT(AA) + (COUNT(Aa)/2) / SAMPLE_COUNT) *
    # (COUNT(aa) + (COUNT(Aa)/2) / SAMPLE_COUNT)
    # * SAMPLE_COUNT
    2 * ((OBS_HOM1 + (OBS_HET/2)) / SAMPLE_COUNT) *
      ((OBS_HOM2 + (OBS_HET/2)) / SAMPLE_COUNT)
      * SAMPLE_COUNT
      AS E_HET,

    # Expected aa
    # q^2
    # (COUNT(aa) + (COUNT(Aa)/2) /
    #  SAMPLE_COUNT) ^ 2 * SAMPLE_COUNT
    POW((OBS_HOM2 + (OBS_HET/2)) /
      SAMPLE_COUNT, 2) * SAMPLE_COUNT
      AS E_HOM2

  FROM observations
  WHERE SAMPLE_COUNT > 0
)

SELECT
  reference_name,
  start,
  reference_bases,
  alt,
  OBS_HOM1,
  OBS_HET,
  OBS_HOM2,
  E_HOM1,
  E_HET,
  E_HOM2,

  # Chi Squared Calculation
  # SUM(((Observed - Expected)^2) / Expected )
  ROUND((POW(OBS_HOM1 - E_HOM1, 2) / E_HOM1)
  + (POW(OBS_HET - E_HET, 2) / E_HET)
  + (POW(OBS_HOM2 - E_HOM2, 2) / E_HOM2), 6)
  AS ChiSq,

  # Determine if Chi Sq value is significant
  # The chi-squared score corresponding to a nominal P-value of 0.05
  # for a table with 2 degrees of freedom is 5.991.
  IF((POW(OBS_HOM1 - E_HOM1, 2) / E_HOM1)
  + (POW(OBS_HET - E_HET, 2) / E_HET)
  + (POW(OBS_HOM2 - E_HOM2, 2) / E_HOM2)
  > 5.991, "TRUE", "FALSE") AS PVALUE_SIG

FROM expectations
WHERE
  E_HOM1 > 0 AND E_HET > 0 AND E_HOM2 > 0
ORDER BY ChiSq DESC, reference_name, start, alt LIMIT 1000
```
Number of rows returned by this query: **1000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:52:12 2016 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> reference_bases </th> <th> alt </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th> <th> ChiSq </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr1 </td> <td align="right"> 10176 </td> <td> A </td> <td> C </td> <td align="right">   0 </td> <td align="right">  17 </td> <td align="right">   0 </td> <td align="right"> 4.25 </td> <td align="right"> 8.50 </td> <td align="right"> 4.25 </td> <td align="right"> 17.00 </td> <td> TRUE </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 10688 </td> <td> C </td> <td> G </td> <td align="right">  16 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right"> 15.06 </td> <td align="right"> 1.88 </td> <td align="right"> 0.06 </td> <td align="right"> 17.00 </td> <td> TRUE </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 14698 </td> <td> C </td> <td> G </td> <td align="right">   0 </td> <td align="right">  17 </td> <td align="right">   0 </td> <td align="right"> 4.25 </td> <td align="right"> 8.50 </td> <td align="right"> 4.25 </td> <td align="right"> 17.00 </td> <td> TRUE </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 15117 </td> <td> A </td> <td> G </td> <td align="right">   0 </td> <td align="right">  17 </td> <td align="right">   0 </td> <td align="right"> 4.25 </td> <td align="right"> 8.50 </td> <td align="right"> 4.25 </td> <td align="right"> 17.00 </td> <td> TRUE </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 16377 </td> <td> T </td> <td> C </td> <td align="right">   0 </td> <td align="right">  17 </td> <td align="right">   0 </td> <td align="right"> 4.25 </td> <td align="right"> 8.50 </td> <td align="right"> 4.25 </td> <td align="right"> 17.00 </td> <td> TRUE </td> </tr>
  <tr> <td> chr1 </td> <td align="right"> 54365 </td> <td> A </td> <td> G </td> <td align="right">  14 </td> <td align="right">   0 </td> <td align="right">   3 </td> <td align="right"> 11.53 </td> <td align="right"> 4.94 </td> <td align="right"> 0.53 </td> <td align="right"> 17.00 </td> <td> TRUE </td> </tr>
   </table>

See also [a version of this query](./sql/hardy-weinberg-udf.sql) that uses BigQuery [user-defined javascript functions](https://cloud.google.com/bigquery/user-defined-functions).

## Heterozygous Haplotype
For each variant within the X and Y chromosome, identify heterozygous variants in male genomes.

First we use our sample information to determine which genomes are male.  

```r
maleSampleIds <- paste("'",
                       filter(sampleInfo,
                              sex %in% c("MALE", "Male", "male", "M", "m"))$call_set_name,
                       "'", sep="", collapse=",")
```


```r
sortAndLimit <- "ORDER BY reference_name, start, alt_concat, call_set_name LIMIT 1000"
result <- DisplayAndDispatchQuery("./sql/sex-chromosome-heterozygous-haplotypes.sql",
                                  project=project,
                                  replacements=c("@MALE_SAMPLE_IDS"=maleSampleIds,
                                                 "-- Optionally add a clause here to constrain the results."=sortAndLimit,
                                                 queryReplacements))
```

```
# Retrieve heterozygous haplotype calls on chromosomes X and Y.
SELECT
  reference_name,
  start,
  `end`,
  reference_bases,
  ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat,
  call.call_set_name,
  (SELECT STRING_AGG(CAST(gt AS STRING)) from UNNEST(call.genotype) gt) AS genotype
FROM
  `genomics-public-data.platinum_genomes.variants` v, v.call call
WHERE
  reference_name IN ('chrX', 'X', 'chrY', 'Y')
  AND call_set_name IN ('NA12877','NA12882','NA12883','NA12884','NA12886','NA12888','NA12889','NA12891','NA12893')
  AND (SELECT LOGICAL_OR(gt = 0) AND LOGICAL_OR(gt = 1) FROM UNNEST(call.genotype) gt)
ORDER BY reference_name, start, alt_concat, call_set_name LIMIT 1000
```
Number of rows returned by this query: **1000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:52:16 2016 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alt_concat </th> <th> call_set_name </th> <th> genotype </th>  </tr>
  <tr> <td> chrX </td> <td align="right"> 2701389 </td> <td align="right"> 2701390 </td> <td> T </td> <td> G </td> <td> not displayed </td> <td> 0,1 </td> </tr>
  <tr> <td> chrX </td> <td align="right"> 2701401 </td> <td align="right"> 2701402 </td> <td> T </td> <td> G </td> <td> not displayed </td> <td> 0,1 </td> </tr>
  <tr> <td> chrX </td> <td align="right"> 2702609 </td> <td align="right"> 2702610 </td> <td> C </td> <td> A </td> <td> not displayed </td> <td> 0,1 </td> </tr>
  <tr> <td> chrX </td> <td align="right"> 2702609 </td> <td align="right"> 2702610 </td> <td> C </td> <td> A </td> <td> not displayed </td> <td> 0,1 </td> </tr>
  <tr> <td> chrX </td> <td align="right"> 2703499 </td> <td align="right"> 2703500 </td> <td> A </td> <td> T </td> <td> not displayed </td> <td> 0,1 </td> </tr>
  <tr> <td> chrX </td> <td align="right"> 2703499 </td> <td align="right"> 2703500 </td> <td> A </td> <td> T </td> <td> not displayed </td> <td> 0,1 </td> </tr>
   </table>

# Removing variants from the Cohort

To mark a variant as problematic so that downstream analyses can filter it out:

* See the [variants patch](https://cloud.google.com/genomics/reference/rest/v1/variants/patch) method.

To remove variants from BigQuery only:

* Materialize the results of queries that include the non-problematic variants to a new table.
* Alternatively, write a custom filtering job similar to what we explored in [Part 2: Data Conversion](./Data-Conversion.md) of this codelab.

To entirely remove a variant from a variant set in the Genomics API:

* See the [variants delete](https://cloud.google.com/genomics/reference/rest/v1/variants/delete) method.
* *Note:* deletion cannot be undone.
