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

* [Genome Call Rate](#genome-call-rate)
* [Missingness Rate](#missingness-rate)
* [Singleton Rate](#singleton-rate)
* [Heterozygosity Rate](#heterozygosity-rate)
* [Homozygosity Rate](#homozygosity-rate)
* [Inbreeding Coefficient](#inbreeding-coefficient)
* [Sex Inference](#sex-inference)
* [Ethnicity Inference](#ethnicity-inference)
* [Genome Similarity](#genome-similarity)

By default this codelab runs upon the Illumina Platinum Genomes Variants. Update the table and change the source of sample information here if you wish to run the queries against a different dataset.

```r
queryReplacements <- list("_GENOME_CALL_TABLE_"="genomics-public-data:platinum_genomes.variants",
                          "_MULTISAMPLE_VARIANT_TABLE_"="google.com:biggene:platinum_genomes.expanded_variants")

sampleData <- read.csv("http://storage.googleapis.com/genomics-public-data/platinum-genomes/other/platinum_genomes_sample_info.csv")
sampleInfo <- dplyr::select(sampleData, call_call_set_name=Catalog.ID, sex=Gender)

ibs <- read.table("./data/platinum-genomes-ibs.tsv",
                  col.names=c("sample1", "sample2", "ibsScore", "similar", "observed"))

# Read in the results of the 2-way PCA over BRCA1.
pca <- read.table("./data/platinum-genomes-X-1kg-brca1-pca.tsv",
                  col.names=c("call_call_set_name", "PC1", "PC2", "count"))

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpDataset.R")
```

## Genome Call Rate

For each genome, count the number of calls (both variants and non-variant calls, if applicable).  Any genomes whose count is far away from the mean may indicate a problem such as sample quality or identical data loaded multiple times.


```r
result <- DisplayAndDispatchQuery("./sql/genome-calls.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Count the number of calls per genome.
SELECT
  call.call_set_name,
  COUNT(call.call_set_name) AS number_of_calls,
FROM
  [genomics-public-data:platinum_genomes.variants]
GROUP BY
  call.call_set_name
ORDER BY
  call.call_set_name
  
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.2 by xtable 1.7-4 package -->
<!-- Mon Feb  8 16:00:53 2016 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> number_of_calls </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 51612762 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 50722005 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 50994084 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 51707412 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 53298483 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 51187390 </td> </tr>
   </table>

And visualizing the results:

```r
ggplot(result, aes(x=number_of_calls)) +
  geom_histogram(color="black", fill="#FF6666") +
  scale_x_continuous(labels=comma) +
  xlab("Number of Calls") +
  ylab("Sample Count") +
  ggtitle("Histogram: Count of Calls Per Genome")
```

<img src="figure/genomeCallsSummary-1.png" title="plot of chunk genomeCallsSummary" alt="plot of chunk genomeCallsSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(result) +
  geom_point(aes(x=call_call_set_name, y=number_of_calls)) +
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

<img src="figure/genomeCalls-1.png" title="plot of chunk genomeCalls" alt="plot of chunk genomeCalls" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
allResults <- result
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
SELECT
  call.call_set_name,
  no_calls,
  all_calls,
  (no_calls/all_calls) AS missingness_rate
FROM (
  SELECT
    call.call_set_name,
    SUM(IF(has_no_calls, delta, 0)) AS no_calls,
    SUM(delta) AS all_calls
  FROM (
    SELECT
      END - start AS delta,
      call.call_set_name,
      SOME(call.genotype == -1) WITHIN call AS has_no_calls,
    FROM
      [genomics-public-data:platinum_genomes.variants]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
    )
  GROUP BY
    call.call_set_name)
ORDER BY
  call.call_set_name
```

```
Warning: NAs introduced by coercion to integer range
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.2 by xtable 1.7-4 package -->
<!-- Mon Feb  8 16:00:55 2016 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> no_calls </th> <th> all_calls </th> <th> missingness_rate </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 41927032 </td> <td align="right">  </td> <td align="right"> 0.01 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 58122228 </td> <td align="right">  </td> <td align="right"> 0.02 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 59224162 </td> <td align="right">  </td> <td align="right"> 0.02 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 58539440 </td> <td align="right">  </td> <td align="right"> 0.02 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 58261455 </td> <td align="right">  </td> <td align="right"> 0.02 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 42249595 </td> <td align="right">  </td> <td align="right"> 0.01 </td> </tr>
   </table>

For some datasets, we see numeric overflow when bigrquery converts 64-bit integer results from BigQuery to 32-bit R integers in the dataframe. For this query, the particular column with the overflow is not used in our downstream analysis in R, so we can omit it.

```r
.Machine$integer.max
```

```
## [1] 2147483647
```

```r
result <- dplyr::select(result, -all_calls)
```

And visualizing the results:

```r
ggplot(result, aes(x=missingness_rate)) +
  geom_histogram(color="black", fill="#FF6666") +
  scale_x_continuous(limits=c(0, NA), labels=percent_format()) +
  xlab("Missingness Rate") +
  ylab("Sample Count") +
  ggtitle("Histogram: Genome-Specific Missingness")
```

<img src="figure/sampleMissingnessSummary-1.png" title="plot of chunk sampleMissingnessSummary" alt="plot of chunk sampleMissingnessSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(result) +
  geom_point(aes(x=call_call_set_name, y=missingness_rate)) +
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

<img src="figure/sampleMissingness-1.png" title="plot of chunk sampleMissingness" alt="plot of chunk sampleMissingness" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
allResults <- full_join(allResults, result)
```

```
## Joining by: "call_call_set_name"
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
SELECT
  call.call_set_name,
  COUNT(call.call_set_name) AS private_variant_count,
FROM (
  SELECT
    reference_name,
    start,
    GROUP_CONCAT(CASE WHEN cnt = 1 THEN 'S'
      WHEN cnt = 2 THEN 'D'
      ELSE STRING(cnt) END) AS SINGLETON_DOUBLETON,
    reference_bases,
    alternate_bases,
    GROUP_CONCAT(call.call_set_name) AS call.call_set_name,
    GROUP_CONCAT(genotype) AS genotype,
    COUNT(call.call_set_name) AS num_samples_with_variant
  FROM (
    SELECT
      reference_name,
      start,
      reference_bases,
      alternate_bases,
      alt_num,
      call.call_set_name,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
      SUM(call.genotype == alt_num) WITHIN call AS cnt,
    FROM (
        FLATTEN((
          SELECT
            reference_name,
            start,
            reference_bases,
            alternate_bases,
            POSITION(alternate_bases) AS alt_num,
            call.call_set_name,
            call.genotype,
          FROM
            [genomics-public-data:platinum_genomes.variants]
          # Optionally add a clause here to limit the query to a particular
          # region of the genome.
          #_WHERE_
          OMIT call IF EVERY(call.genotype <= 0)
        ), alternate_bases)
        )
    OMIT RECORD IF EVERY(alternate_bases IS NULL) OR EVERY(alternate_bases = "<NON_REF>")
    HAVING
      cnt > 0
      )
    GROUP EACH BY
      reference_name,
      start,
      reference_bases,
      alternate_bases
  HAVING
    num_samples_with_variant = 1
    )
GROUP BY
  call.call_set_name
ORDER BY
  private_variant_count DESC
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.2 by xtable 1.7-4 package -->
<!-- Mon Feb  8 16:00:58 2016 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> private_variant_count </th>  </tr>
  <tr> <td> NA12890 </td> <td align="right"> 418760 </td> </tr>
  <tr> <td> NA12892 </td> <td align="right"> 415630 </td> </tr>
  <tr> <td> NA12889 </td> <td align="right"> 415513 </td> </tr>
  <tr> <td> NA12891 </td> <td align="right"> 394767 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 86565 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 82767 </td> </tr>
   </table>

And visualizing the results:

```r
ggplot(result, aes(x=private_variant_count)) +
  geom_histogram(color="black", fill="#FF6666") +
  scale_x_continuous(labels=comma) +
  xlab("Number of Singletons") +
  ylab("Sample Count") +
  ggtitle("Histogram: Count of Singletons Per Genome")
```

<img src="figure/singletonsSummary-1.png" title="plot of chunk singletonsSummary" alt="plot of chunk singletonsSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(result) +
  geom_point(aes(x=call_call_set_name, y=private_variant_count)) +
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

<img src="figure/singletons-1.png" title="plot of chunk singletons" alt="plot of chunk singletons" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
allResults <- full_join(allResults, result)
```

```
## Joining by: "call_call_set_name"
```

## Heterozygosity Rate

For each genome, determine the number of heterozygous variants.


```r
result <- DisplayAndDispatchQuery("./sql/heterozygous-calls-count.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Count the number of heterozygous variants per sample.
SELECT
  call.call_set_name,
  SUM(first_allele != second_allele) AS heterozygous_variant_count
FROM (
  SELECT
    reference_name,
    start,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    call.call_set_name,
    NTH(1, call.genotype) WITHIN call AS first_allele,
    NTH(2, call.genotype) WITHIN call AS second_allele,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts
  FROM
    [genomics-public-data:platinum_genomes.variants]
  # Skip no-calls and single-allele genotypes
  OMIT call IF SOME(call.genotype < 0) OR (2 > COUNT(call.genotype))
  HAVING
    num_alts = 1
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
  )
GROUP BY
  call.call_set_name
ORDER BY
  call.call_set_name
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.2 by xtable 1.7-4 package -->
<!-- Mon Feb  8 16:01:01 2016 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> heterozygous_variant_count </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 3410507 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 3471185 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 3526898 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 3560275 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 3586659 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 3396755 </td> </tr>
   </table>

And visualizing the results:

```r
ggplot(result, aes(x=heterozygous_variant_count)) +
  geom_histogram(color="black", fill="#FF6666") +
  scale_x_continuous(labels=comma) +
  xlab("Number of Heterozygous Variants") +
  ylab("Sample Count") +
  ggtitle("Histogram: Count of Heterozygous Variants Per Genome")
```

<img src="figure/heterozygousSummary-1.png" title="plot of chunk heterozygousSummary" alt="plot of chunk heterozygousSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(result) +
  geom_point(aes(x=call_call_set_name, y=heterozygous_variant_count)) +
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

<img src="figure/heterozygous-1.png" title="plot of chunk heterozygous" alt="plot of chunk heterozygous" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
allResults <- full_join(allResults, result)
```

```
## Joining by: "call_call_set_name"
```

## Homozygosity Rate

For each genome, calculate the fraction of homozygous positions per chromosome.  This is useful to identify uniparental disomy (UPD) or large stretches of homozygosity.


```r
result <- DisplayAndDispatchQuery("./sql/homozygous-variants-alt.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Compute the fraction of homozygous calls for each individual.

SELECT
call.call_set_name,
reference_name,
HOM,
HET_HOM,
N_SITES,
ROUND((HOM) / (HET_HOM), 5) AS F
FROM (
  SELECT
  call.call_set_name,
  reference_name,
  SUM(first_allele = 1 and second_allele = 1) AS HOM,
  SUM(first_allele + second_allele > 0)  AS HET_HOM,
  COUNT(call.call_set_name) AS N_SITES,
  FROM (
    SELECT
    reference_name,
    start,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    call.call_set_name,
    NTH(1, call.genotype) WITHIN call AS first_allele,
    NTH(2, call.genotype) WITHIN call AS second_allele,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    SUM(call.genotype > 0) WITHIN RECORD AS called_alt_allele_count,
    FROM
    [google.com:biggene:platinum_genomes.expanded_variants]
    # Optionally add a clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
    # Skip no calls and haploid sites
    OMIT call IF SOME(call.genotype < 0) OR EVERY(call.genotype = 0) OR (2 != COUNT(call.genotype))
    HAVING
    # Skip 1/2 genotypes.
    num_alts = 1
    # Only use SNPs since non-variant segments are only included for SNPs.
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
  )
  GROUP BY
  call.call_set_name,
  reference_name
)
ORDER BY
call.call_set_name,
reference_name
```
Number of rows returned by this query: **425**.

Let's join this with the sample information:


```r
joinedResult <- inner_join(result, sampleInfo)
```

And visualizing the results:

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

<img src="figure/homozygosityAltSummary-1.png" title="plot of chunk homozygosityAltSummary" alt="plot of chunk homozygosityAltSummary" style="display: block; margin: auto;" />

## Inbreeding Coefficient

For each genome, compare the expected and observed rates of homozygosity.


```r
result <- DisplayAndDispatchQuery("./sql/homozygous-variants.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Compute the expected and observed homozygosity rate for each individual.
SELECT
  call.call_set_name,
  O_HOM,
  ROUND(E_HOM, 2) as E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM (
  SELECT
    call.call_set_name,
    SUM(first_allele = second_allele) AS O_HOM,
    SUM(1.0 - (2.0 * freq * (1.0 - freq) * (called_allele_count / (called_allele_count - 1.0)))) AS E_HOM,
    COUNT(call.call_set_name) AS N_SITES,
  FROM (
    SELECT
      reference_name,
      start,
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      call.call_set_name,
      NTH(1, call.genotype) WITHIN call AS first_allele,
      NTH(2, call.genotype) WITHIN call AS second_allele,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      SUM(call.genotype >= 0) WITHIN RECORD AS called_allele_count,
      IF((SUM(1 = call.genotype) > 0),
        SUM(call.genotype = 1)/SUM(call.genotype >= 0),
        -1)  WITHIN RECORD AS freq
    FROM
      [google.com:biggene:platinum_genomes.expanded_variants]
    # Optionally add a clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
    # Skip no calls and haploid sites
    OMIT call IF SOME(call.genotype < 0) OR (2 != COUNT(call.genotype))
    HAVING
      # Skip 1/2 genotypes.
      num_alts = 1
      # Only use SNPs since non-variant segments are only included for SNPs.
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      # Skip records where all samples have the same allele.
      AND freq > 0 AND freq < 1 
      )
  GROUP BY
    call.call_set_name
    )
ORDER BY
  call.call_set_name
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.2 by xtable 1.7-4 package -->
<!-- Mon Feb  8 16:01:08 2016 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> O_HOM </th> <th> E_HOM </th> <th> N_SITES </th> <th> F </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 6135459 </td> <td align="right"> 6770355.20 </td> <td align="right"> 9546033 </td> <td align="right"> -0.23 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 6044715 </td> <td align="right"> 6748644.58 </td> <td align="right"> 9515963 </td> <td align="right"> -0.25 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 5962198 </td> <td align="right"> 6729953.98 </td> <td align="right"> 9489139 </td> <td align="right"> -0.28 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 5932991 </td> <td align="right"> 6732824.08 </td> <td align="right"> 9493300 </td> <td align="right"> -0.29 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 5909067 </td> <td align="right"> 6734507.42 </td> <td align="right"> 9495760 </td> <td align="right"> -0.30 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 6147961 </td> <td align="right"> 6769854.25 </td> <td align="right"> 9544771 </td> <td align="right"> -0.22 </td> </tr>
   </table>

And visualizing the results:

```r
limits <- c(min(result$O_HOM, result$E_HOM),
            max(result$O_HOM, result$E_HOM))
ggplot(result) +
  geom_point(aes(x=O_HOM, y=E_HOM, label=call_call_set_name), alpha=1/1.5) +
  geom_abline(color="darkslateblue") +
  scale_x_continuous(limits=limits, labels=comma) +
  scale_y_continuous(limits=limits, labels=comma) +
  xlab("Observed Homozygous Variants") +
  ylab("Expected Homozygous Variants") +
  ggtitle("Homozygosity")
```

<img src="figure/homozygosity-1.png" title="plot of chunk homozygosity" alt="plot of chunk homozygosity" style="display: block; margin: auto;" />

And with labels:

```r
ggplot(result) +
  geom_text(aes(x=O_HOM, y=E_HOM, label=call_call_set_name), alpha=1/1.5) +
  geom_abline(color="darkslateblue") +
  scale_x_continuous(limits=limits, labels=comma, expand=c(0.05, 5)) +
  scale_y_continuous(limits=limits, labels=comma) +
  xlab("Observed Homozygous Variants") +
  ylab("Expected Homozygous Variants") +
  ggtitle("Homozygosity")
```

<img src="figure/homozygosityLabelled-1.png" title="plot of chunk homozygosityLabelled" alt="plot of chunk homozygosityLabelled" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
allResults <- full_join(allResults, result)
```

```
## Joining by: "call_call_set_name"
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
# Compute the the homozygous and heterozygous variant counts for each individual
# within chromosome X to help determine whether the sex phenotype value is
# correct for each individual.
SELECT
  call.call_set_name,
  ROUND(SUM(het_RA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_het_alt_in_snvs,
  ROUND(SUM(hom_AA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_hom_alt_in_snvs,
  SUM(hom_AA) AS hom_AA_count,
  SUM(het_RA) AS het_RA_count,
  SUM(hom_RR) AS hom_RR_count,
FROM (
  SELECT
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    call.call_set_name,
    SOME(call.genotype = 0) AND NOT SOME(call.genotype > 0) WITHIN call AS hom_RR,
    SOME(call.genotype > 0) AND NOT SOME(call.genotype = 0) WITHIN call AS hom_AA,
    SOME(call.genotype > 0) AND SOME(call.genotype = 0) WITHIN call AS het_RA
  FROM
    [google.com:biggene:platinum_genomes.expanded_variants]
  WHERE
    (reference_name = 'chrX' OR reference_name = 'X')
    AND start NOT BETWEEN 59999 AND 2699519
    AND start NOT BETWEEN 154931042 AND 155260559
  HAVING
    # Skip 1/2 genotypes _and non-SNP variants
    num_alts = 1
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
  )
GROUP BY
  call.call_set_name
ORDER BY
  call.call_set_name
```
Number of rows returned by this query: **17**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.2.2 by xtable 1.7-4 package -->
<!-- Mon Feb  8 16:01:11 2016 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> perct_het_alt_in_snvs </th> <th> perct_hom_alt_in_snvs </th> <th> hom_AA_count </th> <th> het_RA_count </th> <th> hom_RR_count </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 0.32 </td> <td align="right"> 0.68 </td> <td align="right"> 79739 </td> <td align="right"> 37299 </td> <td align="right"> 212773 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 0.71 </td> <td align="right"> 0.29 </td> <td align="right"> 43666 </td> <td align="right"> 106358 </td> <td align="right"> 183525 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 0.70 </td> <td align="right"> 0.30 </td> <td align="right"> 45655 </td> <td align="right"> 105692 </td> <td align="right"> 180162 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 0.69 </td> <td align="right"> 0.31 </td> <td align="right"> 47261 </td> <td align="right"> 105744 </td> <td align="right"> 178206 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 0.69 </td> <td align="right"> 0.31 </td> <td align="right"> 47446 </td> <td align="right"> 105364 </td> <td align="right"> 178591 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 0.31 </td> <td align="right"> 0.69 </td> <td align="right"> 78815 </td> <td align="right"> 34893 </td> <td align="right"> 214852 </td> </tr>
   </table>

Let's join this with the sample information:

```r
joinedResult <- inner_join(result, sampleInfo)
```

And visualize the results:

```r
ggplot(joinedResult) +
  geom_boxplot(aes(x=sex, y=perct_het_alt_in_snvs, fill=sex)) +
  scale_y_continuous(labels = percent_format()) +
  xlab("Sex") +
  ylab("Heterozygosity Rate ") +
  ggtitle("Box Plot: Heterozygosity Rate on the X Chromosome")
```

<img src="figure/sexCheckSummary-1.png" title="plot of chunk sexCheckSummary" alt="plot of chunk sexCheckSummary" style="display: block; margin: auto;" />


```r
p <- ggplot(joinedResult) +
  geom_point(aes(x=call_call_set_name, y=perct_het_alt_in_snvs, color=sex)) +
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

<img src="figure/sexCheck-1.png" title="plot of chunk sexCheck" alt="plot of chunk sexCheck" style="display: block; margin: auto;" />

Let's accumulate our sample-specific results for later use.

```r
allResults <- full_join(allResults, result)
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
sampleInfo1kg <- dplyr::select(sampleData1kg, call_call_set_name=Sample, sex=Gender, ethnicity=Super_Population)

# Update our sample information for Platinum Genomes as "Unknown" since this is what we are trying to check.
sampleInfoToCheck <- mutate(sampleInfo, ethnicity="Unknown")

# Note that 5 samples are in both datasets, so those will be plotted twice with different symbols.
pcaPlatinumX1kg <- inner_join(pca, rbind(sampleInfoToCheck, sampleInfo1kg))
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

<img src="figure/pca-with-ethnicity-1.png" title="plot of chunk pca-with-ethnicity" alt="plot of chunk pca-with-ethnicity" style="display: block; margin: auto;" />

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

<img src="figure/ibs-1.png" title="plot of chunk ibs" alt="plot of chunk ibs" style="display: block; margin: auto;" />

# Removing Genomes from the Cohort

To only remove a genome from BigQuery only:

* Re-export the table to BigQuery using the `--call-set-ids` flag on the `gcloud alpha genomics variantsets export` command.

To exclude a genome from data returned by the Genomics API:

* See the `callSetIds` property on the [variants search](https://cloud.google.com/genomics/reference/rest/v1/variants/search) method.

To entirely remove a genome from a variant set in the Genomics API:

* See the [callsets delete](https://cloud.google.com/genomics/reference/rest/v1/callsets/delete) method.
* To delete a callset using a command line tool, see the `gcloud alpha genomics callsets delete` command.
* *Note:* deletion cannot be undone.

# Summary

Let's wrap up with a quick comparison using the variables we've collected for each sample.

```r
plot(dplyr::select(allResults, number_of_calls, missingness_rate,
                   private_variant_count, heterozygous_variant_count, O_HOM))
```

<img src="figure/summary-1.png" title="plot of chunk summary" alt="plot of chunk summary" style="display: block; margin: auto;" />

If we see any relationships that we do not expect, it may be worth a closer look.  We could also compare additional columns (only a few were included to keep this plot scrutable).

--------------------------------------------------------
_Next_: [Part 4: Variant-Level QC](./Variant-Level-QC.md)
