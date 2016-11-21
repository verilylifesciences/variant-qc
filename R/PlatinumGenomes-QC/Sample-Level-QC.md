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

# Part 3: Sample-Level QC





In Part 3 of the codelab, we perform some quality control analyses that could help to identify any problematic genomes that should be removed from the cohort before proceeding with further analysis.  The appropriate cut off thresholds will depend upon the input dataset and/or other factors.

* [Genome Variant Call Rate](#genome-variant-call-rate)
* [Missingness Rate](#missingness-rate)
* [Singleton Rate](#singleton-rate)
* [Heterozygosity Rate](#heterozygosity-rate)
* [Homozygosity Rate](#homozygosity-rate)
* [Inbreeding Coefficient](#inbreeding-coefficient)
* [Ti/Tv Ratio per Chromosome](#titv-ratio-per-chromosome)
* [Sex Inference](#sex-inference)
* [Ethnicity Inference](#ethnicity-inference)
* [Genome Similarity](#genome-similarity)

By default this codelab runs upon the Illumina Platinum Genomes Variants. Update the table and change the source of sample information here if you wish to run the queries against a different dataset.

```r
source("./rHelpers/platinumGenomesDataset.R")

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpDataset.R")
```

## Genome Variant Call Rate

For each genome, count the number of variant calls.  Any genomes whose count is far away from the mean may indicate a problem such as sample quality or identical data loaded multiple times.


```r
result <- DisplayAndDispatchQuery("./sql/genome-variant-calls.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Count the number of variant calls per genome.
SELECT
  call.call_set_name,
  COUNT(call.call_set_name) AS number_of_calls
FROM
  `genomics-public-data.platinum_genomes.variants` v, v.call call
WHERE
  # Skip homozygous reference calls, no-calls, and non-passing variants.
  EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
  AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
  AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
GROUP BY
  call_set_name
ORDER BY
  call_set_name
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:07:53 2016 -->
<table border=1>
<tr> <th> call_set_name </th> <th> number_of_calls </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 3859838 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 3874126 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 3875322 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 3877931 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 3903770 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 3888561 </td> </tr>
   </table>

Let's join this with the sample information and visualize the results:

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
ggplot(joinedResult, aes(y=number_of_calls, x=ethnicity, fill=ethnicity)) +
  geom_boxplot() +
  scale_y_continuous(expand = c(0.3, 0)) +
  stat_summary(fun.data=get_boxplot_fun_data, geom="text", position=position_dodge(width=0.9), vjust=-0.5) +
  ylab("Number of Variant Calls") +
  xlab("Ethnicity") +
  ggtitle("Box plot: Count of variant calls per genome by ethnicity") +
  theme(axis.text.x=element_text(angle=50, hjust=1))
```

<img src="./platinum_genomes/figure/genomeCallsSummary-1.png" title="plot of chunk genomeCallsSummary" alt="plot of chunk genomeCallsSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(joinedResult) +
  geom_point(aes(x=call_set_name, y=number_of_calls, color=sex)) +
  scale_x_discrete(expand=c(0.05, 1)) +
  scale_y_continuous(labels=comma) +
  xlab("Sample") +
  ylab("Number of Calls") +
  ggtitle("Scatter Plot: Count of Calls Per Genome")
if(nrow(result) <= 20) {
  p + theme(axis.text.x=element_text(angle=50, hjust=1))
} else {
  p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())
}
```

<img src="./platinum_genomes/figure/genomeCalls-1.png" title="plot of chunk genomeCalls" alt="plot of chunk genomeCalls" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
sampleResults <- result
```

## Missingness Rate

For each genome, determine the percentage of sites explicitly called as a no-call.  If this percentage is too high, the genome may be problematic.


```r
result <- DisplayAndDispatchQuery("./sql/sample-level-missingness.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Compute the ratio of positions corresponding to no-calls versus all positions
# called (reference, variant, and no-calls).
WITH deltas AS (
  SELECT
    `end` - start AS delta,
    call.call_set_name,
    EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
      OR
      EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.')) AS has_no_calls
  FROM
    `genomics-public-data.platinum_genomes.variants` v, v.call call
),

positions_called AS (
  SELECT
    call_set_name,
    SUM(IF(has_no_calls, delta, 0)) AS no_calls,
    SUM(delta) AS all_calls
  FROM deltas
  GROUP BY
    call_set_name
)

SELECT
  call_set_name,
  no_calls,
  all_calls,
  (no_calls/all_calls) AS missingness_rate
FROM positions_called
ORDER BY
  call_set_name
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:07:59 2016 -->
<table border=1>
<tr> <th> call_set_name </th> <th> no_calls </th> <th> all_calls </th> <th> missingness_rate </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 85025808 </td> <td align="right">  </td> <td align="right"> 0.03 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 100372181 </td> <td align="right">  </td> <td align="right"> 0.04 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 101872729 </td> <td align="right">  </td> <td align="right"> 0.04 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 101488526 </td> <td align="right">  </td> <td align="right"> 0.04 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 102372220 </td> <td align="right">  </td> <td align="right"> 0.04 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 84981715 </td> <td align="right">  </td> <td align="right"> 0.03 </td> </tr>
   </table>

Note that for some datasets, we see message "NAs introduced by coercion to integer range" when [bigrquery](https://github.com/rstats-db/bigrquery) converts 64-bit integer results from BigQuery to 32-bit R integers in the dataframe. For this query, the particular column with the issue is not used in our downstream analysis in R, so we can omit it.

```r
.Machine$integer.max
```

```
## [1] 2147483647
```

```r
result <- dplyr::select(result, -all_calls)
```

Let's join this with the sample information and visualize the results:

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
ggplot(joinedResult, aes(y=missingness_rate, x=ethnicity)) +
  geom_boxplot() +
  stat_summary(fun.data=get_boxplot_fun_data, geom="text", position=position_dodge(width=0.9), vjust=-0.5) +
  scale_y_continuous(limits=c(0, NA), labels=percent_format()) +
  ylab("Missingness Rate") +
  xlab("Sequencing Platform") +
  ggtitle("Genome-Specific Missingness") +
  theme(axis.text.x=element_text(angle=50, hjust=1))
```

<img src="./platinum_genomes/figure/sampleMissingnessSummary-1.png" title="plot of chunk sampleMissingnessSummary" alt="plot of chunk sampleMissingnessSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(joinedResult) +
  geom_point(aes(x=call_set_name, y=missingness_rate, color=sex)) +
  scale_x_discrete(expand=c(0.05, 1)) +
  scale_y_continuous(limits=c(0, NA), labels=percent_format()) +
  xlab("Sample") +
  ylab("Missingness Rate") +
  ggtitle("Scatter Plot: Genome-Specific Missingness")
if(nrow(result) <= 20) {
  p + theme(axis.text.x=element_text(angle=50, hjust=1))
} else {
  p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())
}
```

<img src="./platinum_genomes/figure/sampleMissingness-1.png" title="plot of chunk sampleMissingness" alt="plot of chunk sampleMissingness" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
sampleResults <- full_join(sampleResults, result)
```

## Singleton Rate

For each genome, count the number of variants shared by no other member of the cohort.  Too many private calls for a particular individual may indicate a problem.


```r
result <- DisplayAndDispatchQuery("./sql/private-variants.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Compute private variants counts for each sample.
WITH filtered_called_alleles AS (
  SELECT
    reference_name,
    start,
    reference_bases AS ref,
    alt,
    call.call_set_name,
    (SELECT COUNT(CAST(gt = alt_offset+1 AS INT64)) FROM call.genotype gt) AS allele_cnt
  FROM
    `genomics-public-data.platinum_genomes.variants` v, v.call call, v.alternate_bases alt WITH OFFSET alt_offset
  WHERE
    # Skip homozygous reference calls, no-calls, and non-passing variants.
    EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
),

grouped_alleles AS (
  SELECT
    reference_name,
    start,
    ref,
    alt,
    STRING_AGG(call_set_name) AS call_set_name,
    COUNT(call_set_name) AS num_samples_with_variant
  FROM filtered_called_alleles
  GROUP BY
    reference_name,
    start,
    ref,
    alt
)

SELECT
  call_set_name,
  COUNT(call_set_name) AS private_variant_count
FROM grouped_alleles
WHERE
  num_samples_with_variant = 1
GROUP BY
  call_set_name
ORDER BY
  private_variant_count DESC
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:08:05 2016 -->
<table border=1>
<tr> <th> call_set_name </th> <th> private_variant_count </th>  </tr>
  <tr> <td> NA12892 </td> <td align="right"> 303216 </td> </tr>
  <tr> <td> NA12890 </td> <td align="right"> 292892 </td> </tr>
  <tr> <td> NA12889 </td> <td align="right"> 283839 </td> </tr>
  <tr> <td> NA12891 </td> <td align="right"> 268484 </td> </tr>
  <tr> <td> NA12877 </td> <td align="right"> 24783 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 23264 </td> </tr>
   </table>

Let's join this with the sample information and visualize the results:

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
ggplot(joinedResult, aes(y=private_variant_count, x=ethnicity, fill=ethnicity)) +
  geom_boxplot() +
  stat_summary(fun.data=get_boxplot_fun_data, geom="text", position=position_dodge(width=0.9), vjust=-0.5) +
  scale_y_continuous(labels=comma, expand = c(0.3, 0)) +
  ylab("Number of Singletons") +
  xlab("Ethnicity") +
  ggtitle("Box plot: Count of singletons per genome by ethnicity") +
  theme(axis.text.x=element_text(angle=50, hjust=1))
```

<img src="./platinum_genomes/figure/singletonsSummary-1.png" title="plot of chunk singletonsSummary" alt="plot of chunk singletonsSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(joinedResult) +
  geom_point(aes(x=call_set_name, y=private_variant_count, color=sex)) +
  scale_x_discrete(expand=c(0.05, 1)) +
  scale_y_continuous(labels=comma) +
  xlab("Sample") +
  ylab("Number of Singletons") +
  ggtitle("Scatter Plot: Count of Singletons Per Genome")
if(nrow(result) <= 20) {
  p + theme(axis.text.x=element_text(angle=50, hjust=1))
} else {
  p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())
}
```

<img src="./platinum_genomes/figure/singletons-1.png" title="plot of chunk singletons" alt="plot of chunk singletons" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
sampleResults <- full_join(sampleResults, result)
```

## Heterozygosity Rate

For each genome, determine the number of heterozygous variants.


```r
result <- DisplayAndDispatchQuery("./sql/heterozygous-calls-by-sample.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Count the number of heterozygous variants per sample.
WITH filtered_snp_calls AS (
  SELECT
    call.call_set_name,
    call.genotype[ORDINAL(1)] AS first_allele,
    call.genotype[ORDINAL(2)] AS second_allele
  FROM
    `genomics-public-data.platinum_genomes.variants` v, v.call call
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    # Skip homozygous reference calls, no-calls, non-passing variants, and non-diploid calls.
    AND EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
    AND ARRAY_LENGTH(call.genotype) = 2
)

SELECT
  call_set_name,
  SUM(CAST((first_allele != second_allele) AS INT64)) AS heterozygous_variant_count
FROM filtered_snp_calls
GROUP BY
  call_set_name
ORDER BY
  call_set_name
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:08:11 2016 -->
<table border=1>
<tr> <th> call_set_name </th> <th> heterozygous_variant_count </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 1916100 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 1961370 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 1959626 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 1982368 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 1984345 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 1919017 </td> </tr>
   </table>

Let's join this with the sample information and visualize the results:

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
ggplot(joinedResult, aes(y=heterozygous_variant_count, x=ethnicity, fill=ethnicity)) +
  geom_boxplot() +
  stat_summary(fun.data=get_boxplot_fun_data, geom="text", position=position_dodge(width=0.9), vjust=-0.5) +
  scale_y_continuous(labels=comma, expand = c(0.3, 0)) +
  ylab("Number of Heterozyous Variants") +
  xlab("Ethnicity") +
  ggtitle("Box plot: Count of heterozygous variants per genome by ethnicity") +
  theme(axis.text.x=element_text(angle=50, hjust=1))
```

<img src="./platinum_genomes/figure/heterozygousSummary-1.png" title="plot of chunk heterozygousSummary" alt="plot of chunk heterozygousSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(joinedResult) +
  geom_point(aes(x=call_set_name, y=heterozygous_variant_count, color=sex)) +
  scale_x_discrete(expand=c(0.05, 1)) +
  scale_y_continuous(labels=comma) +
  xlab("Sample") +
  ylab("Number of Heterozygous Variants") +
  ggtitle("Scatter Plot: Count of Heterozygous Variants Per Genome")
if(nrow(result) <= 20) {
  p + theme(axis.text.x=element_text(angle=50, hjust=1))
} else {
  p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())
}
```

<img src="./platinum_genomes/figure/heterozygous-1.png" title="plot of chunk heterozygous" alt="plot of chunk heterozygous" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
sampleResults <- full_join(sampleResults, result)
```

## Homozygosity Rate

For each genome, calculate the fraction of homozygous positions per chromosome.  This is useful to identify uniparental disomy (UPD) or large stretches of homozygosity.


```r
result <- DisplayAndDispatchQuery("./sql/homozygous-variant-rate-by-sample-and-reference.sql",
                                  project=project,
                                  replacements=queryReplacements,
                                  max_pages=Inf)
```

```
# Compute the ratio of homozygous vs. heterozygous variant calls for each individual.
WITH filtered_snp_calls AS (
  SELECT
    reference_name,
    reference_bases AS ref,
    ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat,
    call.call_set_name,
    call.genotype[ORDINAL(1)] AS first_allele,
    call.genotype[ORDINAL(2)] AS second_allele
  FROM
    `genomics-public-data.platinum_genomes.variants` v, v.call call
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    # Skip homozygous reference calls, no-calls, non-passing variants, and non-diploid calls.
    AND EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
    AND ARRAY_LENGTH(call.genotype) = 2
),

variant_counts AS (
  SELECT
    call_set_name,
    reference_name,
    SUM(CAST(first_allele = 1 AND second_allele = 1 AS INT64)) AS HOM_ALT,
    SUM(CAST(first_allele = 1 OR second_allele = 1 AS INT64))  AS HAS_ALT,
    COUNT(call_set_name) AS N_SITES
  FROM filtered_snp_calls
  GROUP BY
    call_set_name,
    reference_name
)

SELECT
  call_set_name,
  reference_name,
  HOM_ALT,
  HAS_ALT,
  N_SITES,
  ROUND((HOM_ALT) / (HAS_ALT), 5) AS F
FROM variant_counts
ORDER BY
  call_set_name,
  reference_name
```
Number of rows returned by this query: **408**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:08:16 2016 -->
<table border=1>
<tr> <th> call_set_name </th> <th> reference_name </th> <th> HOM_ALT </th> <th> HAS_ALT </th> <th> N_SITES </th> <th> F </th>  </tr>
  <tr> <td> NA12877 </td> <td> chr1 </td> <td align="right"> 103349 </td> <td align="right"> 247958 </td> <td align="right"> 247958 </td> <td align="right"> 0.42 </td> </tr>
  <tr> <td> NA12877 </td> <td> chr10 </td> <td align="right"> 66163 </td> <td align="right"> 161993 </td> <td align="right"> 161993 </td> <td align="right"> 0.41 </td> </tr>
  <tr> <td> NA12877 </td> <td> chr11 </td> <td align="right"> 67922 </td> <td align="right"> 164243 </td> <td align="right"> 164243 </td> <td align="right"> 0.41 </td> </tr>
  <tr> <td> NA12877 </td> <td> chr12 </td> <td align="right"> 60750 </td> <td align="right"> 154954 </td> <td align="right"> 154954 </td> <td align="right"> 0.39 </td> </tr>
  <tr> <td> NA12877 </td> <td> chr13 </td> <td align="right"> 55000 </td> <td align="right"> 126927 </td> <td align="right"> 126927 </td> <td align="right"> 0.43 </td> </tr>
  <tr> <td> NA12877 </td> <td> chr14 </td> <td align="right"> 40460 </td> <td align="right"> 101858 </td> <td align="right"> 101858 </td> <td align="right"> 0.40 </td> </tr>
   </table>
Let's join this with the sample information and visualize the results:

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
  ggplot(joinedResult, aes(y=F, x=reference_name, color=sex)) +
  geom_boxplot() +
  facet_grid(sex ~ .) +
  scale_y_continuous(labels=comma) +
  ylab("Fraction of Homozygous Variants") +
  xlab("Reference Name") +
  ggtitle("Fraction of Homozygous Variants Per Genome") +
  theme(axis.text.x=element_text(angle=50, hjust=1))
```

<img src="./platinum_genomes/figure/homozygositySummary-1.png" title="plot of chunk homozygositySummary" alt="plot of chunk homozygositySummary" style="display: block; margin: auto;" />


```r
sampleReferenceResults <- result
```

## Inbreeding Coefficient

For each genome, compare the expected and observed rates of homozygosity.


```r
result <- DisplayAndDispatchQuery("./sql/homozygosity-coefficient.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Compute the expected and observed homozygosity rate for each individual.
WITH variant_calls AS (
  SELECT
    call.call_set_name,
    call.genotype[ORDINAL(1)] = call.genotype[ORDINAL(2)] AS O_HOM,
    (SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt WHERE gt >= 0)) FROM v.call) AS AN,
    (SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt WHERE gt = 1)) FROM v.call) AS AC
  FROM
    `google.com:biggene.platinum_genomes.expanded_variants` v, v.call call
  WHERE
    # Only include biallelic snps within autosomes. (Concise but inexact regexp used for brevity.)
    REGEXP_CONTAINS(reference_name, r'^(chr)?([1-2])?[0-9]$')
    AND reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(v.alternate_bases) = 1
    AND v.alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    # Skip no-calls, non-passing variants, and non-diploid calls.
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
    AND ARRAY_LENGTH(call.genotype) = 2
),

grouped_values AS (
  SELECT
    call_set_name,
    SUM(CAST(O_HOM AS INT64)) AS O_HOM,
    SUM(1.0 - (2.0 * (AC/AN) * (1.0 - (AC/AN)) * (AN / (AN - 1.0)))) AS E_HOM,
    COUNT(call_set_name) AS N_SITES
  FROM variant_calls
  WHERE
    # Ensure allelic frequency is in the range (0, 1).
    AN > 0 AND AC < AN
  GROUP BY
    call_set_name
)

SELECT
  call_set_name,
  O_HOM,
  ROUND(E_HOM, 2) AS E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM grouped_values
ORDER BY
  call_set_name
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:08:26 2016 -->
<table border=1>
<tr> <th> call_set_name </th> <th> O_HOM </th> <th> E_HOM </th> <th> N_SITES </th> <th> F </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 5107357 </td> <td align="right"> 5131169.71 </td> <td align="right"> 7022727 </td> <td align="right"> -0.01 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 5079568 </td> <td align="right"> 5098324.62 </td> <td align="right"> 6984556 </td> <td align="right"> -0.01 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 4937093 </td> <td align="right"> 4987386.15 </td> <td align="right"> 6841573 </td> <td align="right"> -0.03 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 4903339 </td> <td align="right"> 4979031.62 </td> <td align="right"> 6832036 </td> <td align="right"> -0.04 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 4877364 </td> <td align="right"> 4958460.55 </td> <td align="right"> 6808732 </td> <td align="right"> -0.04 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 5152560 </td> <td align="right"> 5176429.00 </td> <td align="right"> 7070614 </td> <td align="right"> -0.01 </td> </tr>
   </table>

Let's join this with the sample information and visualize the results:

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
limits <- c(min(result$O_HOM, result$E_HOM),
            max(result$O_HOM, result$E_HOM))
ggplot(result) +
  geom_point(aes(x=O_HOM, y=E_HOM, label=call_set_name), alpha=1/1.5) +
  geom_abline(color="darkslateblue") +
  scale_x_continuous(limits=limits, labels=comma) +
  scale_y_continuous(limits=limits, labels=comma) +
  xlab("Observed Homozygous Variants") +
  ylab("Expected Homozygous Variants") +
  ggtitle("Homozygosity")
```

<img src="./platinum_genomes/figure/homozygosityCoeff-1.png" title="plot of chunk homozygosityCoeff" alt="plot of chunk homozygosityCoeff" style="display: block; margin: auto;" />

And with labels:

```r
ggplot(result) +
  geom_text(aes(x=O_HOM, y=E_HOM, label=call_set_name), alpha=1/1.5) +
  geom_abline(color="darkslateblue") +
  scale_x_continuous(limits=limits, labels=comma, expand=c(0.05, 5)) +
  scale_y_continuous(limits=limits, labels=comma) +
  xlab("Observed Homozygous Variants") +
  ylab("Expected Homozygous Variants") +
  ggtitle("Homozygosity")
```

<img src="./platinum_genomes/figure/homozygosityCoeffLabelled-1.png" title="plot of chunk homozygosityCoeffLabelled" alt="plot of chunk homozygosityCoeffLabelled" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
sampleResults <- full_join(sampleResults, result)
```

## Ti/Tv Ratio per Chromosome

For each genome, determine the Ti/Tv ratio per chromosome.

```r
result <- DisplayAndDispatchQuery("./sql/ti-tv-by-sample-and-reference.sql",
                                  project=project,
                                  replacements=queryReplacements,
                                  max_pages=Inf)
```

```
# Compute the transition/transversion ratio per sample and reference name.
WITH filtered_snp_calls AS (
  SELECT
    reference_name,
    call.call_set_name,
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
    reference_name,
    call_set_name,
    SUM(CAST(mutation IN ('A->G', 'G->A', 'C->T', 'T->C') AS INT64)) AS transitions,
    SUM(CAST(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                          'A->T', 'T->A', 'C->G', 'G->C') AS INT64)) AS transversions
  FROM filtered_snp_calls
  GROUP BY
    reference_name,
    call_set_name
)

SELECT
  reference_name,
  call_set_name,
  transitions,
  transversions,
  transitions/transversions AS titv
FROM mutation_type_counts
WHERE
  transversions > 0
ORDER BY
  titv DESC,
  call_set_name
```
Number of rows returned by this query: **408**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:08:31 2016 -->
<table border=1>
<tr> <th> reference_name </th> <th> call_set_name </th> <th> transitions </th> <th> transversions </th> <th> titv </th>  </tr>
  <tr> <td> chr22 </td> <td> NA12887 </td> <td align="right"> 29028 </td> <td align="right"> 11960 </td> <td align="right"> 2.43 </td> </tr>
  <tr> <td> chr22 </td> <td> NA12893 </td> <td align="right"> 28724 </td> <td align="right"> 11864 </td> <td align="right"> 2.42 </td> </tr>
  <tr> <td> chr22 </td> <td> NA12886 </td> <td align="right"> 28476 </td> <td align="right"> 11767 </td> <td align="right"> 2.42 </td> </tr>
  <tr> <td> chr22 </td> <td> NA12892 </td> <td align="right"> 30420 </td> <td align="right"> 12616 </td> <td align="right"> 2.41 </td> </tr>
  <tr> <td> chr22 </td> <td> NA12879 </td> <td align="right"> 28914 </td> <td align="right"> 11995 </td> <td align="right"> 2.41 </td> </tr>
  <tr> <td> chr22 </td> <td> NA12889 </td> <td align="right"> 29342 </td> <td align="right"> 12183 </td> <td align="right"> 2.41 </td> </tr>
   </table>

Let's join this with the sample information and visualize the results:

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
ggplot(joinedResult, aes(y=titv, x=reference_name, color=sex)) +
  geom_boxplot() +
  facet_grid(~ ethnicity)+
  scale_y_continuous(labels=comma) +
  ylab("Ti/Tv ratio") +
  xlab("Chromosome") +
  ggtitle("Ti/Tv ratio per genome") +
  theme(axis.text.x=element_text(angle=50, hjust=1))
```

<img src="./platinum_genomes/figure/titvSummary-1.png" title="plot of chunk titvSummary" alt="plot of chunk titvSummary" style="display: block; margin: auto;" />


```r
sampleReferenceResults <- full_join(sampleReferenceResults, result)
```

## Sex Inference

For each genome, compare the sex from the sample information to the heterozygosity rate on the chromosome X calls.

In the query that follows we specifically examine the percentage of SNP variants that are heterozygous but note that the Inbreeding Coefficient query above can also be used as a sex check when run upon only chromosome X omitting the pseudoautosomal regions.  For more detail, see the [comparison](./comparison/QC-Comparison.md) against results from other tools.


```r
result <- DisplayAndDispatchQuery("./sql/check-sex.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Compute the homozygous and heterozygous variant counts for each individual
# within chromosome X to help determine whether the sex phenotype value is
# correct for each individual.
WITH filtered_snp_calls AS (
  SELECT
    call.call_set_name,
    CAST((SELECT LOGICAL_AND(gt > 0) FROM UNNEST(call.genotype) gt) AS INT64) AS hom_AA,
    CAST(EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
      AND EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt = 0) AS INT64) AS het_RA
  FROM
    `genomics-public-data.platinum_genomes.variants` v, v.call call
  WHERE
    reference_name IN ('chrX', 'X')
    # Locations of PAR1 and PAR2 on GRCh37.
    AND start NOT BETWEEN 59999 AND 2699519
    AND start NOT BETWEEN 154931042 AND 155260559
    # Only include biallelic snps.
    AND reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    # Skip non-passing calls.
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
)

SELECT
  call_set_name,
  ROUND(SUM(het_RA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_het_alt_in_snvs,
  ROUND(SUM(hom_AA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_hom_alt_in_snvs,
  SUM(hom_AA) AS hom_AA_count,
  SUM(het_RA) AS het_RA_count
FROM filtered_snp_calls
GROUP BY
  call_set_name
ORDER BY
  call_set_name
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.3 by xtable 1.8-2 package -->
<!-- Wed Nov 23 13:08:37 2016 -->
<table border=1>
<tr> <th> call_set_name </th> <th> perct_het_alt_in_snvs </th> <th> perct_hom_alt_in_snvs </th> <th> hom_AA_count </th> <th> het_RA_count </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 0.01 </td> <td align="right"> 0.99 </td> <td align="right"> 65922 </td> <td align="right"> 655 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 0.61 </td> <td align="right"> 0.39 </td> <td align="right"> 35769 </td> <td align="right"> 56382 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 0.59 </td> <td align="right"> 0.41 </td> <td align="right"> 37747 </td> <td align="right"> 55144 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 0.58 </td> <td align="right"> 0.42 </td> <td align="right"> 39094 </td> <td align="right"> 53670 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 0.57 </td> <td align="right"> 0.43 </td> <td align="right"> 39969 </td> <td align="right"> 52976 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 0.01 </td> <td align="right"> 0.99 </td> <td align="right"> 67242 </td> <td align="right"> 763 </td> </tr>
   </table>

Let's join this with the sample information and visualize the results:

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
ggplot(joinedResult) +
  geom_boxplot(aes(x=sex, y=perct_het_alt_in_snvs, fill=sex)) +
  scale_y_continuous(labels = percent_format()) +
  xlab("Sex") +
  ylab("Heterozygosity Rate ") +
  ggtitle("Box Plot: Heterozygosity Rate on the X Chromosome")
```

<img src="./platinum_genomes/figure/sexCheckSummary-1.png" title="plot of chunk sexCheckSummary" alt="plot of chunk sexCheckSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(joinedResult) +
  geom_point(aes(x=call_set_name, y=perct_het_alt_in_snvs, color=sex)) +
  scale_x_discrete(expand=c(0.05, 1)) +
  scale_y_continuous(labels = percent_format()) +
  xlab("Sample") +
  ylab("Heterozygosity Rate ") +
  ggtitle("Scatter Plot: Heterozygosity Rate on the X Chromosome")
if(nrow(result) <= 20) {
  p + theme(axis.text.x=element_text(angle=50, hjust=1))
} else {
  p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())
}
```

<img src="./platinum_genomes/figure/sexCheck-1.png" title="plot of chunk sexCheck" alt="plot of chunk sexCheck" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
sampleResults <- full_join(sampleResults, result)
```

## Ethnicity Inference

For each genome, compare the ethncity from the sample information to the clustering in this analysis.

For this check, we:

* use the intersection of common variants found in both 1,000 Genomes phase 1 variants and Platinum Genomes
* compute PCA on those variants in common between the two data
* examine whether the individuals in Platinum Genomes cluster with other samples of the same ethnicity

See the Google Genomics [2-way PCA cookbook entry](http://googlegenomics.readthedocs.org/en/latest/use_cases/compute_principal_coordinate_analysis/2-way-pca.html) for the details as to how to run this pipeline.

Note that this `n^2` analysis is a cluster compute job instead of a BigQuery query.

### Results


```r
# Read in the demographic information for 1,000 Genomes.
sampleData1kg <- read.csv("http://storage.googleapis.com/genomics-public-data/1000-genomes/other/sample_info/sample_info.csv")
sampleInfo1kg <- dplyr::select(sampleData1kg, call_set_name=Sample, sex=Gender, ethnicity=Super_Population)

# Update our sample information for Platinum Genomes as "Unknown" since this is what we are trying to check.
sampleInfoToCheck <- mutate(sampleInfo, ethnicity="Unknown")

# Note that 5 samples are in both datasets, so those will be plotted twice with different symbols.
pcaPlatinumX1kg <- inner_join(pca, rbind(sampleInfoToCheck, sampleInfo1kg), by=c("call_call_set_name" = "call_set_name"))
pcaPlatinumX1kg <- mutate(pcaPlatinumX1kg, unknown=(ethnicity == "Unknown"))
```


```r
ggplot(pcaPlatinumX1kg) +
  geom_point(aes(x=PC1, y=PC2,
                 color=ethnicity,
                 shape=ethnicity,
                 size=unknown)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  scale_shape_manual(values=c(3, 3, 3, 3, 19)) +
  scale_size_manual(values=c(2,4)) +
  ggtitle("2-way Principal Coordinate Analysis upon Platinum Genomes and 1,000 Genomes")
```

<img src="./platinum_genomes/figure/pca-with-ethnicity-1.png" title="plot of chunk pca-with-ethnicity" alt="plot of chunk pca-with-ethnicity" style="display: block; margin: auto;" />

## Genome Similarity

Perform a simplistic similarity check on each pair of genomes to identify any mislabled or cross-contaminated samples.  See the Google Genomics [Identity-By-State cookbook entry](http://googlegenomics.readthedocs.org/en/latest/use_cases/compute_identity_by_state/index.html) for the details as to how to run this pipeline.

Note that this `n^2` analysis is a cluster compute job instead of a BigQuery query.

### Results


```r
ggplot(ibs) +
  geom_tile(aes(x=sample1, y=sample2, fill=ibsScore), colour="white") +
  scale_fill_gradient(low="white", high="steelblue",
                      na.value="black", trans="log",
                      guide=guide_colourbar(title= "IBS Score")) +
  theme(axis.text.x=element_text(angle=50, hjust=1)) +
  xlab("Sample 1") +
  ylab("Sample 2") +
  ggtitle("Identity By State (IBS) Heat Map")
```

<img src="./platinum_genomes/figure/ibs-1.png" title="plot of chunk ibs" alt="plot of chunk ibs" style="display: block; margin: auto;" />

# Removing Genomes from the Cohort

To remove a genome from BigQuery only:

* Re-export the table to BigQuery using the `--call-set-ids` flag on the `gcloud alpha genomics variantsets export` command.

To exclude a genome from data returned by the Genomics API:

* See the `callSetIds` property on the [variants search](https://cloud.google.com/genomics/reference/rest/v1/variants/search) method.

To entirely remove a genome from a variant set in the Genomics API:

* See the [callsets delete](https://cloud.google.com/genomics/reference/rest/v1/callsets/delete) method.
* To delete a callset using a command line tool, see the `gcloud alpha genomics callsets delete` command.
* *Note:* deletion cannot be undone.

# Summary

Accumulated results for per sample analyses:

```r
dim(sampleResults)
```

```
[1] 17 14
```

```r
summary(sampleResults)
```

```
 call_set_name      number_of_calls      no_calls        
 Length:17          Min.   :3825240   Min.   : 84981715  
 Class :character   1st Qu.:3870044   1st Qu.: 86275658  
 Mode  :character   Median :3877100   Median : 88847153  
                    Mean   :3875038   Mean   : 93487248  
                    3rd Qu.:3888561   3rd Qu.:101300696  
                    Max.   :3903770   Max.   :102372220  
 missingness_rate  private_variant_count heterozygous_variant_count
 Min.   :0.02968   Min.   : 12846        Min.   :1885606           
 1st Qu.:0.03013   1st Qu.: 16293        1st Qu.:1919017           
 Median :0.03103   Median : 18669        Median :1939962           
 Mean   :0.03265   Mean   : 81112        Mean   :1949622           
 3rd Qu.:0.03538   3rd Qu.: 24783        3rd Qu.:1980057           
 Max.   :0.03576   Max.   :303216        Max.   :2016053           
     O_HOM             E_HOM            N_SITES              F           
 Min.   :4877364   Min.   :4958461   Min.   :6808732   Min.   :-0.04611  
 1st Qu.:4943497   1st Qu.:4998364   1st Qu.:6862189   1st Qu.:-0.03635  
 Median :4980269   Median :5034754   Median :6902259   Median :-0.02712  
 Mean   :4994753   Mean   :5041605   Mean   :6917916   Mean   :-0.02509  
 3rd Qu.:5039706   3rd Qu.:5072849   3rd Qu.:6980214   3rd Qu.:-0.01259  
 Max.   :5152560   Max.   :5176429   Max.   :7070614   Max.   : 0.00028  
 perct_het_alt_in_snvs perct_hom_alt_in_snvs  hom_AA_count  
 Min.   :0.0080        Min.   :0.3820        Min.   :34943  
 1st Qu.:0.0120        1st Qu.:0.4060        1st Qu.:37747  
 Median :0.0160        Median :0.9840        Median :64397  
 Mean   :0.2866        Mean   :0.7134        Mean   :52428  
 3rd Qu.:0.5940        3rd Qu.:0.9880        3rd Qu.:65922  
 Max.   :0.6180        Max.   :0.9920        Max.   :67242  
  het_RA_count  
 Min.   :  556  
 1st Qu.:  769  
 Median : 1073  
 Mean   :26351  
 3rd Qu.:55144  
 Max.   :56543  
```

```r
write.csv(sampleResults, file=file.path(kResultsDir, "sampleResults.csv"))
```

Accumulated results for per sample, reference analyses:

```r
dim(sampleReferenceResults)
```

```
[1] 408   9
```

```r
summary(sampleReferenceResults)
```

```
 call_set_name      reference_name        HOM_ALT          HAS_ALT      
 Length:408         Length:408         Min.   :   123   Min.   :   135  
 Class :character   Class :character   1st Qu.: 35368   1st Qu.: 84462  
 Mode  :character   Mode  :character   Median : 55347   Median :131487  
                                       Mean   : 54947   Mean   :136181  
                                       3rd Qu.: 72902   3rd Qu.:187800  
                                       Max.   :109168   Max.   :269962  
    N_SITES             F           transitions     transversions  
 Min.   :   135   Min.   :0.3122   Min.   :    90   Min.   :   45  
 1st Qu.: 84462   1st Qu.:0.3802   1st Qu.: 59466   1st Qu.:25027  
 Median :131487   Median :0.3986   Median : 88706   Median :42931  
 Mean   :136181   Mean   :0.4236   Mean   : 92279   Mean   :43904  
 3rd Qu.:187800   3rd Qu.:0.4167   3rd Qu.:127154   3rd Qu.:61680  
 Max.   :269962   Max.   :1.0000   Max.   :182576   Max.   :87391  
      titv      
 Min.   :1.369  
 1st Qu.:2.053  
 Median :2.111  
 Mean   :2.109  
 3rd Qu.:2.164  
 Max.   :2.427  
```

```r
write.csv(sampleReferenceResults, file=file.path(kResultsDir, "sampleReferenceResults.csv"))
```


```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.6 (El Capitan)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitr_1.13           ggplot2_2.1.0        bigrquery_0.3.0.9000
## [4] scales_0.4.0         dplyr_0.5.0          RCurl_1.95-4.8      
## [7] bitops_1.0-6         xtable_1.8-2        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.7      magrittr_1.5     munsell_0.4.3    lattice_0.20-33 
##  [5] colorspace_1.2-6 R6_2.1.2         stringr_1.0.0    httr_1.2.1      
##  [9] plyr_1.8.3       tools_3.2.3      grid_3.2.3       nlme_3.1-128    
## [13] gtable_0.2.0     mgcv_1.8-12      DBI_0.5-1        htmltools_0.3.5 
## [17] openssl_0.9.5    lazyeval_0.2.0   assertthat_0.1   digest_0.6.9    
## [21] tibble_1.2       Matrix_1.2-6     formatR_1.4      reshape2_1.4.1  
## [25] curl_2.2         evaluate_0.9     rmarkdown_0.9.6  labeling_0.3    
## [29] stringi_1.0-1    jsonlite_1.1
```
--------------------------------------------------------
_Next_: [Part 4: Variant-Level QC](./Variant-Level-QC.md)
