<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->

<!-- Copyright 2014 Google Inc. All rights reserved. -->

<!-- Licensed under the Apache License, Version 2.0 (the "License"); -->
<!-- you may not use this file except in compliance with the License. -->
<!-- You may obtain a copy of the License at -->

<!--     http://www.apache.org/licenses/LICENSE-2.0 -->

<!-- Unless required by applicable law or agreed to in writing, software -->
<!-- distributed under the License is distributed on an "AS IS" BASIS, -->
<!-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. -->
<!-- See the License for the specific language governing permissions and -->
<!-- limitations under the License. -->

Performing QC on gVCF Data using Google Genomics
================================================

The following example makes use of [Illumina Platinum Genomes](http://www.illumina.com/platinumgenomes/).  For more detail about how this data was loaded into the Google Genomics API, please see [Google Genomics Public Data](https://cloud.google.com/genomics/data/platinum-genomes).

Setting Up and Describing the Data
----------------------------------





```r
# Setup for BigQuery access
require(bigrquery)
require(xtable)
require(RCurl)
require(dplyr)

project <- "gbsc-gcp-project-mvp"                   # put your projectID here

DisplayAndDispatchQuery <- function(queryUri, replacements=list()) {
  if(grepl("^https.*", queryUri)) {
    querySql <- getURL(queryUri, ssl.verifypeer=FALSE)
  } else {
    querySql <- readChar(queryUri, nchars=1e6)
  }
  for(replacement in names(replacements)) {
    querySql <- gsub(replacement, replacements[[replacement]], querySql, fixed=TRUE)
  }
  cat(querySql)
  query_exec(querySql, project)
}

table_replacement <- list("_THE_TABLE_"="genomics-public-data:platinum_genomes.variants",
                          "_THE_EXPANDED_TABLE_"="google.com:biggene:platinum_genomes.expanded_variants")
```

Let's take a look at a few of the [variants within BRCA1 via BigQuery](https://github.com/googlegenomics/getting-started-bigquery/blob/master/RMarkdown/literate-programming-demo.md#data-visualization)

```r
result <- DisplayAndDispatchQuery("https://raw.githubusercontent.com/googlegenomics/getting-started-bigquery/master/sql/variant-level-data-for-brca1.sql",
                                  replacements=table_replacement)
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
HAVING
  alternate_bases IS NOT NULL
ORDER BY
  start,
  alternate_bases
```
Number of rows returned by this query: 335.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:31 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> quality </th> <th> filter </th> <th> names </th> <th> num_samples </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td align="right"> 733.47 </td> <td> PASS </td> <td>  </td> <td align="right">   7 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td align="right"> 41196822 </td> <td> CT </td> <td> C </td> <td align="right"> 63.74 </td> <td> LowQD </td> <td>  </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td align="right"> 41196823 </td> <td> CTT </td> <td> C,CT </td> <td align="right"> 314.59 </td> <td> PASS </td> <td>  </td> <td align="right">   3 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196840 </td> <td align="right"> 41196841 </td> <td> G </td> <td> T </td> <td align="right"> 85.68 </td> <td> TruthSensitivityTranche99.90to100.00,LowQD </td> <td>  </td> <td align="right">   2 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td align="right"> 41197274 </td> <td> C </td> <td> A </td> <td align="right"> 1011.08 </td> <td> PASS </td> <td>  </td> <td align="right">   7 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197938 </td> <td align="right"> 41197939 </td> <td> A </td> <td> AT </td> <td align="right"> 86.95 </td> <td> LowQD </td> <td>  </td> <td align="right">   3 </td> </tr>
   </table>

And then let's take a look at the domain and range of values for alternate_bases:

```r
result <- DisplayAndDispatchQuery("./sql/characterize-alts.sql",
                                  replacements=table_replacement)
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
OMIT
  RECORD IF EVERY(alternate_bases IS NULL)
GROUP BY
  alt_contains_no_special_characters
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:33 2015 -->
<table border=1>
<tr> <th> number_of_variant_records </th> <th> alt_contains_no_special_characters </th> <th> max_ref_len </th> <th> max_alt_len </th>  </tr>
  <tr> <td align="right"> 12634588 </td> <td> TRUE </td> <td align="right">  56 </td> <td align="right">  47 </td> </tr>
   </table>
We see from the query results that there are no special charaters in alternate_bases and the maximum length is ~50 base pairs.

And finally let's take a look at the domain and range of values for genotype:

```r
result <- DisplayAndDispatchQuery("./sql/genotypes-brca1.sql",
                                  replacements=table_replacement)
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

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:34 2015 -->
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
We see from the query results the variety of genotypes within BRCA1.

To summarize attributes of this particular dataset that we need to consider when performing QC:
* It is in gVCF format which adds complexity above and beyond [similar examples for the 1,000 Genomes dataset](https://github.com/googlegenomics/bigquery-examples/blob/master/1000genomes/sql/README.md).
* It is comprised only of SNPs and INDELs (contains no structural variants).
* The values for `alternate_bases` are just comprised of the letters A,C,G,T (e.g., contains no `<DEL>` values).
* It contains some single-allele and 1/2 genotypes.

Working with gVCF Data
======================

Data in gVCF format can be challenging to query.  For more detail see this [comparison](https://github.com/googlegenomics/bigquery-examples/tree/master/pgp/data-stories/schema-comparisons#motivation).

In the analyses below, sometimes we work with the original gVCF data and sometimes we work with data converted from gVCF to VCF to make querying easier.  

For more detail about the conversion process see the [PlatinumGenomes gVCF-to-VCF codelab](https://github.com/deflaux/codelabs/tree/qc-codelab/Python/PlatinumGenomes-gVCF-to-VCF).

Sample-Level QC
===============

Check Singletons
----------------

```r
result <- DisplayAndDispatchQuery("./sql/private-variants-brca1.sql",
                                  replacements=table_replacement)
```

```
# Private variants within BRCA1.
SELECT
  reference_name AS CHROM,
  start AS POS,
  GROUP_CONCAT(CASE WHEN cnt = 1 THEN 'S'
    WHEN cnt = 2 THEN 'D'
    ELSE STRING(cnt) END) AS SINGLETON_DOUBLETON,
  reference_bases AS REF,
  alternate_bases AS ALT,
  GROUP_CONCAT(call.call_set_name) AS INDV,
  GROUP_CONCAT(genotype) AS genotype,
  SUM(num_samples_with_variant) AS num_samples_with_variant
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
    COUNT(call.call_set_name) WITHIN RECORD AS num_samples_with_variant
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
        WHERE
          reference_name = 'chr17'
          AND start BETWEEN 41196311
          AND 41277499
        OMIT
          call IF EVERY(call.genotype = -1)
          ),
        alternate_bases)
      )
  OMIT
    RECORD IF alternate_bases IS NULL
  HAVING
    cnt > 0
    )
GROUP BY
  chrom,
  pos,
  ref,
  alt
HAVING
  num_samples_with_variant = 1
ORDER BY
  POS,
  INDV
```
Number of rows returned by this query: 63.

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:38 2015 -->
<table border=1>
<tr> <th> CHROM </th> <th> POS </th> <th> SINGLETON_DOUBLETON </th> <th> REF </th> <th> ALT </th> <th> INDV </th> <th> genotype </th> <th> num_samples_with_variant </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td> S </td> <td> CT </td> <td> C </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197958 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12884 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198195 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12893 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198220 </td> <td> S </td> <td> T </td> <td> G </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198226 </td> <td> S </td> <td> A </td> <td> C </td> <td> NA12890 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204831 </td> <td> D </td> <td> ACACACACACT </td> <td> A </td> <td> NA12890 </td> <td> "1,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204839 </td> <td> D </td> <td> ACT </td> <td> A </td> <td> NA12892 </td> <td> "1,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204862 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12882 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41208206 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12880 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210033 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12888 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211352 </td> <td> S </td> <td> T </td> <td> G </td> <td> NA12892 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211485 </td> <td> S </td> <td> CACA </td> <td> CACAACA </td> <td> NA12878 </td> <td> "1,2" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41217551 </td> <td> S </td> <td> A </td> <td> C </td> <td> NA12889 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218851 </td> <td> S </td> <td> G </td> <td> T </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218884 </td> <td> S </td> <td> A </td> <td> G </td> <td> NA12893 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219261 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12889 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219852 </td> <td> S </td> <td> AT </td> <td> A </td> <td> NA12890 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219872 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12887 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219873 </td> <td> S </td> <td> G </td> <td> T </td> <td> NA12887 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219905 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12880 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41223537 </td> <td> S </td> <td> G </td> <td> GAATGTTCACTGTAACAATGCTTGT </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224871 </td> <td> S </td> <td> A </td> <td> C </td> <td> NA12893 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41227082 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12887 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41227082 </td> <td> S </td> <td> C </td> <td> CGGAA </td> <td> NA12892 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229351 </td> <td> S </td> <td> C </td> <td> T </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229761 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12878 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232485 </td> <td> S </td> <td> T </td> <td> G </td> <td> NA12881 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41236026 </td> <td> S </td> <td> A </td> <td> C </td> <td> NA12893 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41238125 </td> <td> S </td> <td> T </td> <td> TACACAC </td> <td> NA12880 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239705 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12891 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242077 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12880 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250046 </td> <td> S </td> <td> C </td> <td> T </td> <td> NA12888 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250220 </td> <td> S </td> <td> AT </td> <td> A </td> <td> NA12892 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252576 </td> <td> S </td> <td> A </td> <td> G </td> <td> NA12891 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252645 </td> <td> S </td> <td> ATGT </td> <td> A </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252648 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12878 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252675 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12890 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252683 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12878 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252692 </td> <td> S </td> <td> ATAAT </td> <td> A </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252693 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12877 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12877 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12878 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254373 </td> <td> S </td> <td> C </td> <td> T </td> <td> NA12892 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254393 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12884 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254979 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12887 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256088 </td> <td> S </td> <td> AAAAAAAAAGAAAAG </td> <td> A </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256091 </td> <td> S </td> <td> A </td> <td> G </td> <td> NA12878 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256094 </td> <td> S </td> <td> A </td> <td> G </td> <td> NA12878 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258182 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12887 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259078 </td> <td> S </td> <td> A </td> <td> ATT </td> <td> NA12892 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260721 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12891 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260984 </td> <td> S </td> <td> G </td> <td> C </td> <td> NA12892 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263116 </td> <td> S </td> <td> C </td> <td> CA </td> <td> NA12878 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263429 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12891 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264110 </td> <td> S </td> <td> C </td> <td> CT </td> <td> NA12891 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270777 </td> <td> S </td> <td> C </td> <td> CT </td> <td> NA12878 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41271292 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12887 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td> S </td> <td> G </td> <td> C </td> <td> NA12889 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273699 </td> <td> S </td> <td> C </td> <td> CAAAA </td> <td> NA12888 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273774 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12892 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274876 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12883 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275080 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12887 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276517 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12880 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
   </table>

Compare to [brca1.singletons](./data/singletons/brca1.singletons) which has 85 some of which are for 0/0 genotypes from reference matching blocks (see the [vcftools command line](./data/singletons/brca1.log) used to create this file).


```r
expectedResult <- read.table("./data/singletons/brca1.singletons", header=TRUE)
# Convert to zero-based coordinates
expectedResult <- mutate(expectedResult, POS = POS - 1)
# Clean colnames to match
colnames(expectedResult) <- gsub('\\.+', '_', colnames(expectedResult))
```

How many singletons do the two results have in common?

```r
nrow(inner_join(result, expectedResult))
```

```
## Joining by: c("CHROM", "POS", "SINGLETON_DOUBLETON", "INDV")
```

```
## [1] 75
```

Which singletons were only identified by BigQuery?

```r
onlyBQ <- anti_join(result, expectedResult)
```

```
## Joining by: c("CHROM", "POS", "SINGLETON_DOUBLETON", "INDV")
```

```r
print(xtable(onlyBQ), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:38 2015 -->
<table border=1>
<tr> <th> CHROM </th> <th> POS </th> <th> SINGLETON_DOUBLETON </th> <th> REF </th> <th> ALT </th> <th> INDV </th> <th> genotype </th> <th> num_samples_with_variant </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211485 </td> <td> S </td> <td> CACA </td> <td> CACAACA </td> <td> NA12878 </td> <td> "1,2" </td> <td align="right">   1 </td> </tr>
   </table>

Which singletons were only identified by vcftools?

```r
onlyVcftools <- anti_join(expectedResult, result)
```

```
## Joining by: c("CHROM", "POS", "SINGLETON_DOUBLETON", "INDV")
```

```r
print(xtable(onlyVcftools), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:38 2015 -->
<table border=1>
<tr> <th> CHROM </th> <th> POS </th> <th> SINGLETON_DOUBLETON </th> <th> ALLELE </th> <th> INDV </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694.00 </td> <td> S </td> <td> AAT </td> <td> NA12886 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841.00 </td> <td> S </td> <td> T </td> <td> NA12888 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196320.00 </td> <td> D </td> <td> T </td> <td> NA12886 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196319.00 </td> <td> D </td> <td> T </td> <td> NA12886 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196318.00 </td> <td> D </td> <td> G </td> <td> NA12886 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196317.00 </td> <td> D </td> <td> T </td> <td> NA12886 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196316.00 </td> <td> D </td> <td> G </td> <td> NA12886 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196315.00 </td> <td> D </td> <td> A </td> <td> NA12886 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196314.00 </td> <td> D </td> <td> A </td> <td> NA12886 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196313.00 </td> <td> D </td> <td> G </td> <td> NA12886 </td> </tr>
   </table>

Retrieving the gVCF data for the singletons identified only by vcftools:

```r
having <- paste("start = ", onlyVcftools$POS,
                sep="", collapse=" OR ")
result <- DisplayAndDispatchQuery("./sql/examine-data.sql",
                                  replacements=c(table_replacement,
                                                 "_HAVING_"=having))
```

```
# Examine the data for particular calls.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS gt,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
  GROUP_CONCAT(STRING(call.pl)) WITHIN call AS likelihood,
FROM
  [genomics-public-data:platinum_genomes.variants]
WHERE
  reference_name = 'chr17'
HAVING
  start = 41252694 OR start = 41204841 OR start = 41196320 OR start = 41196319 OR start = 41196318 OR start = 41196317 OR start = 41196316 OR start = 41196315 OR start = 41196314 OR start = 41196313
ORDER BY
  start,
  end,
  call.call_set_name
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:40 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> call_call_set_name </th> <th> gt </th> <th> quality </th> <th> filter </th> <th> likelihood </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196313 </td> <td align="right"> 41196746 </td> <td> G </td> <td>  </td> <td> NA12886 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12877 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 793,44,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12878 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 1025,63,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12879 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 1175,48,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12880 </td> <td> -1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12881 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 725,63,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12882 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 917,54,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12883 </td> <td> -1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12884 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 943,81,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12885 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 809,69,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12886 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 785,38,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12887 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 754,39,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12888 </td> <td> 0,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 359,0,25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12889 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 815,54,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12890 </td> <td> -1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12891 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 877,78,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12892 </td> <td> -1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12893 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 859,72,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td> T </td> <td> NA12877 </td> <td> 0,1 </td> <td align="right"> 143.44 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 173,0,874 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12878 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12880 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12881 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12883 </td> <td> 0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12884 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12889 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12890 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12892 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252697 </td> <td> AAT </td> <td> A </td> <td> NA12884 </td> <td> 1,1 </td> <td align="right"> 648.06 </td> <td> LowGQX </td> <td> 690,17,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252697 </td> <td> AAT </td> <td> A </td> <td> NA12886 </td> <td> 0,1 </td> <td align="right"> 648.06 </td> <td> LowGQX </td> <td> 823,0,38 </td> </tr>
   </table>

It appears that they correspond either to:
* a reference-matching block, so not actually a singleton and just perhaps violating an assumption in the vcftools code
* or a non-singleon variant, perhaps due to a problem in converting the gVCF data to all-positions VCF via gvcftools?

Check Individual Heterozygosity
-----------------------------------


```r
result <- DisplayAndDispatchQuery("./sql/homozygous-variants-brca1.sql",
                                  replacements=table_replacement)
```

```
# Individual Homozygosity
SELECT
  INDV,
  O_HOM,
  ROUND(E_HOM, 2) as E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM (
  SELECT
    call.call_set_name AS INDV,
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
      NTH(1,
        call.genotype) WITHIN call AS first_allele,
      NTH(2,
        call.genotype) WITHIN call AS second_allele,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      SUM(call.genotype >= 0) WITHIN RECORD AS called_allele_count,
      IF((SUM(1 = call.genotype) > 0),
        SUM(call.genotype = 1)/SUM(call.genotype >= 0),
        -1)  WITHIN RECORD AS freq
    FROM
      [google.com:biggene:platinum_genomes.expanded_variants]
    WHERE
      reference_name = 'chr17'
      AND start BETWEEN 41196311
      AND 41277499
    OMIT
      call IF SOME(call.genotype < 0)
      OR (2 > COUNT(call.genotype))
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      )
  GROUP BY
    INDV
    )
ORDER BY
  INDV
```
Number of rows returned by this query: 17.

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:42 2015 -->
<table border=1>
<tr> <th> INDV </th> <th> O_HOM </th> <th> E_HOM </th> <th> N_SITES </th> <th> F </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right"> 252 </td> <td align="right"> 233.49 </td> <td align="right"> 274 </td> <td align="right"> 0.46 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right"> 100 </td> <td align="right"> 210.84 </td> <td align="right"> 268 </td> <td align="right"> -1.94 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right"> 246 </td> <td align="right"> 233.49 </td> <td align="right"> 274 </td> <td align="right"> 0.31 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right"> 101 </td> <td align="right"> 195.83 </td> <td align="right"> 266 </td> <td align="right"> -1.35 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right"> 237 </td> <td align="right"> 233.49 </td> <td align="right"> 274 </td> <td align="right"> 0.09 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right"> 252 </td> <td align="right"> 233.49 </td> <td align="right"> 274 </td> <td align="right"> 0.46 </td> </tr>
  <tr> <td> NA12883 </td> <td align="right">  94 </td> <td align="right"> 186.30 </td> <td align="right"> 254 </td> <td align="right"> -1.36 </td> </tr>
  <tr> <td> NA12884 </td> <td align="right"> 246 </td> <td align="right"> 232.52 </td> <td align="right"> 272 </td> <td align="right"> 0.34 </td> </tr>
  <tr> <td> NA12885 </td> <td align="right"> 251 </td> <td align="right"> 233.49 </td> <td align="right"> 274 </td> <td align="right"> 0.43 </td> </tr>
  <tr> <td> NA12886 </td> <td align="right"> 250 </td> <td align="right"> 232.52 </td> <td align="right"> 272 </td> <td align="right"> 0.44 </td> </tr>
  <tr> <td> NA12887 </td> <td align="right">  85 </td> <td align="right"> 203.06 </td> <td align="right"> 264 </td> <td align="right"> -1.94 </td> </tr>
  <tr> <td> NA12888 </td> <td align="right">  95 </td> <td align="right"> 205.64 </td> <td align="right"> 267 </td> <td align="right"> -1.80 </td> </tr>
  <tr> <td> NA12889 </td> <td align="right">  98 </td> <td align="right"> 205.87 </td> <td align="right"> 267 </td> <td align="right"> -1.76 </td> </tr>
  <tr> <td> NA12890 </td> <td align="right"> 242 </td> <td align="right"> 222.27 </td> <td align="right"> 271 </td> <td align="right"> 0.40 </td> </tr>
  <tr> <td> NA12891 </td> <td align="right"> 241 </td> <td align="right"> 229.03 </td> <td align="right"> 269 </td> <td align="right"> 0.30 </td> </tr>
  <tr> <td> NA12892 </td> <td align="right">  92 </td> <td align="right"> 204.72 </td> <td align="right"> 266 </td> <td align="right"> -1.84 </td> </tr>
  <tr> <td> NA12893 </td> <td align="right"> 238 </td> <td align="right"> 230.91 </td> <td align="right"> 271 </td> <td align="right"> 0.18 </td> </tr>
   </table>

Compare to [brca1.het](./data/heterozygous/brca1.het) (see the [vcftools command line](./data/heterozygous/brca1.log) used to create this file).



```r
expectedResult <- read.table("./data/heterozygous/brca1.het", header=TRUE)
# Clean colnames to match
colnames(expectedResult) <- gsub('\\.+$', '', colnames(expectedResult))
colnames(expectedResult) <- gsub('\\.+', '_', colnames(expectedResult))
```


```r
joinedResult <- inner_join(expectedResult, result, by=c("INDV"))
print(xtable(joinedResult[,order(colnames(joinedResult))]), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:42 2015 -->
<table border=1>
<tr> <th> E_HOM.x </th> <th> E_HOM.y </th> <th> F.x </th> <th> F.y </th> <th> INDV </th> <th> N_SITES.x </th> <th> N_SITES.y </th> <th> O_HOM.x </th> <th> O_HOM.y </th>  </tr>
  <tr> <td align="right"> 185.60 </td> <td align="right"> 233.49 </td> <td align="right"> 0.71 </td> <td align="right"> 0.46 </td> <td> NA12877 </td> <td align="right"> 254 </td> <td align="right"> 274 </td> <td align="right"> 234 </td> <td align="right"> 252 </td> </tr>
  <tr> <td align="right"> 195.70 </td> <td align="right"> 210.84 </td> <td align="right"> -1.28 </td> <td align="right"> -1.94 </td> <td> NA12878 </td> <td align="right"> 278 </td> <td align="right"> 268 </td> <td align="right">  90 </td> <td align="right"> 100 </td> </tr>
  <tr> <td align="right"> 186.70 </td> <td align="right"> 233.49 </td> <td align="right"> 0.57 </td> <td align="right"> 0.31 </td> <td> NA12879 </td> <td align="right"> 256 </td> <td align="right"> 274 </td> <td align="right"> 226 </td> <td align="right"> 246 </td> </tr>
  <tr> <td align="right"> 195.50 </td> <td align="right"> 195.83 </td> <td align="right"> -1.31 </td> <td align="right"> -1.35 </td> <td> NA12880 </td> <td align="right"> 277 </td> <td align="right"> 266 </td> <td align="right">  89 </td> <td align="right"> 101 </td> </tr>
  <tr> <td align="right"> 185.10 </td> <td align="right"> 233.49 </td> <td align="right"> 0.47 </td> <td align="right"> 0.09 </td> <td> NA12881 </td> <td align="right"> 253 </td> <td align="right"> 274 </td> <td align="right"> 217 </td> <td align="right"> 237 </td> </tr>
  <tr> <td align="right"> 186.10 </td> <td align="right"> 233.49 </td> <td align="right"> 0.68 </td> <td align="right"> 0.46 </td> <td> NA12882 </td> <td align="right"> 255 </td> <td align="right"> 274 </td> <td align="right"> 233 </td> <td align="right"> 252 </td> </tr>
  <tr> <td align="right"> 197.40 </td> <td align="right"> 186.30 </td> <td align="right"> -1.19 </td> <td align="right"> -1.36 </td> <td> NA12883 </td> <td align="right"> 285 </td> <td align="right"> 254 </td> <td align="right">  93 </td> <td align="right">  94 </td> </tr>
  <tr> <td align="right"> 187.20 </td> <td align="right"> 232.52 </td> <td align="right"> 0.60 </td> <td align="right"> 0.34 </td> <td> NA12884 </td> <td align="right"> 257 </td> <td align="right"> 272 </td> <td align="right"> 229 </td> <td align="right"> 246 </td> </tr>
  <tr> <td align="right"> 186.50 </td> <td align="right"> 233.49 </td> <td align="right"> 0.65 </td> <td align="right"> 0.43 </td> <td> NA12885 </td> <td align="right"> 256 </td> <td align="right"> 274 </td> <td align="right"> 232 </td> <td align="right"> 251 </td> </tr>
  <tr> <td align="right"> 186.10 </td> <td align="right"> 232.52 </td> <td align="right"> 0.65 </td> <td align="right"> 0.44 </td> <td> NA12886 </td> <td align="right"> 255 </td> <td align="right"> 272 </td> <td align="right"> 231 </td> <td align="right"> 250 </td> </tr>
  <tr> <td align="right"> 195.10 </td> <td align="right"> 203.06 </td> <td align="right"> -1.44 </td> <td align="right"> -1.94 </td> <td> NA12887 </td> <td align="right"> 277 </td> <td align="right"> 264 </td> <td align="right">  77 </td> <td align="right">  85 </td> </tr>
  <tr> <td align="right"> 196.90 </td> <td align="right"> 205.64 </td> <td align="right"> -1.34 </td> <td align="right"> -1.80 </td> <td> NA12888 </td> <td align="right"> 280 </td> <td align="right"> 267 </td> <td align="right">  86 </td> <td align="right">  95 </td> </tr>
  <tr> <td align="right"> 195.10 </td> <td align="right"> 205.87 </td> <td align="right"> -1.35 </td> <td align="right"> -1.76 </td> <td> NA12889 </td> <td align="right"> 275 </td> <td align="right"> 267 </td> <td align="right">  87 </td> <td align="right">  98 </td> </tr>
  <tr> <td align="right"> 184.60 </td> <td align="right"> 222.27 </td> <td align="right"> 0.59 </td> <td align="right"> 0.40 </td> <td> NA12890 </td> <td align="right"> 253 </td> <td align="right"> 271 </td> <td align="right"> 225 </td> <td align="right"> 242 </td> </tr>
  <tr> <td align="right"> 181.60 </td> <td align="right"> 229.03 </td> <td align="right"> 0.55 </td> <td align="right"> 0.30 </td> <td> NA12891 </td> <td align="right"> 250 </td> <td align="right"> 269 </td> <td align="right"> 219 </td> <td align="right"> 241 </td> </tr>
  <tr> <td align="right"> 196.90 </td> <td align="right"> 204.72 </td> <td align="right"> -1.32 </td> <td align="right"> -1.84 </td> <td> NA12892 </td> <td align="right"> 282 </td> <td align="right"> 266 </td> <td align="right">  85 </td> <td align="right">  92 </td> </tr>
  <tr> <td align="right"> 183.80 </td> <td align="right"> 230.91 </td> <td align="right"> 0.49 </td> <td align="right"> 0.18 </td> <td> NA12893 </td> <td align="right"> 252 </td> <td align="right"> 271 </td> <td align="right"> 217 </td> <td align="right"> 238 </td> </tr>
   </table>

The logic in the query looks similar to vcftools [output_het method](http://sourceforge.net/p/vcftools/code/HEAD/tree/trunk/cpp/variant_file_output.cpp#l165) but there is clearly a difference.  TODO: investigate the difference further.

Cohort Level QC
===============

Check Hardy-Weinberg Equilibrium
-----------------------------------

```r
result <- DisplayAndDispatchQuery("./sql/hardy-weinberg-brca1-expanded.sql",
                                  replacements=table_replacement)
```

```
# The following query computes the Hardy-Weinberg equilibrium for BRCA1 variants.
SELECT
  CHR,
  POS,
  ref,
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
  IF((POW(OBS_HOM1 - E_HOM1, 2) / E_HOM1)
  + (POW(OBS_HET - E_HET, 2) / E_HET)
  + (POW(OBS_HOM2 - E_HOM2, 2) / E_HOM2)
  > 5.991, "TRUE", "FALSE") AS PVALUE_SIG

FROM (
    SELECT
      CHR,
      POS,
      ref,
      alt,
      OBS_HOM1,
      OBS_HET,
      OBS_HOM2,

      # Expected AA
      # p^2
      # ((COUNT(AA) + (COUNT(Aa)/2) /
      #  SAMPLE_COUNT) ^ 2) * SAMPLE_COUNT
      ROUND(POW((OBS_HOM1 + (OBS_HET/2)) /
        SAMPLE_COUNT, 2) * SAMPLE_COUNT, 2)
        AS E_HOM1,

      # Expected Aa
      # 2pq
      # 2 * (COUNT(AA) + (COUNT(Aa)/2) / SAMPLE_COUNT) *
      # (COUNT(aa) + (COUNT(Aa)/2) / SAMPLE_COUNT)
      # * SAMPLE_COUNT
      ROUND(2 * ((OBS_HOM1 + (OBS_HET/2)) / SAMPLE_COUNT) *
        ((OBS_HOM2 + (OBS_HET/2)) / SAMPLE_COUNT)
        * SAMPLE_COUNT, 2)
        AS E_HET,

      # Expected aa
      # q^2
      # (COUNT(aa) + (COUNT(Aa)/2) /
      #  SAMPLE_COUNT) ^ 2 * SAMPLE_COUNT
      ROUND(POW((OBS_HOM2 + (OBS_HET/2)) /
        SAMPLE_COUNT, 2) * SAMPLE_COUNT, 2)
        AS E_HOM2,

  FROM (
    SELECT
      reference_name AS CHR,
      start AS POS,
      reference_bases AS ref,
      alternate_bases AS alt,
      HOM_REF AS OBS_HOM1,
      HET AS OBS_HET,
      HOM_ALT AS OBS_HOM2,
      HOM_REF + HET + HOM_ALT AS SAMPLE_COUNT,
    FROM (
      SELECT
        reference_name,
        start,
        END,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
        SUM(EVERY(1 = call.genotype)) WITHIN call AS HOM_ALT,
        SUM(SOME(0 = call.genotype)
          AND SOME(1 = call.genotype)) WITHIN call AS HET,
      FROM
        [google.com:biggene:platinum_genomes.expanded_variants]
      WHERE
        reference_name = 'chr17'
        AND start BETWEEN 41196311
        AND 41277499
      HAVING
        # Skip 1/2 genotypes
        num_alts = 1
        )))
ORDER BY
  CHR,
  POS,
  ref
```
Number of rows returned by this query: 333.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:46 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ref </th> <th> alt </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th> <th> ChiSq </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   7 </td> <td align="right">   0 </td> <td align="right"> 10.72 </td> <td align="right"> 5.56 </td> <td align="right"> 0.72 </td> <td align="right"> 1.14 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td> CT </td> <td> C </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> <td align="right"> 1.00 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196840 </td> <td> G </td> <td> T </td> <td align="right">  15 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 15.06 </td> <td align="right"> 1.88 </td> <td align="right"> 0.06 </td> <td align="right"> 0.07 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   7 </td> <td align="right">   0 </td> <td align="right"> 10.72 </td> <td align="right"> 5.56 </td> <td align="right"> 0.72 </td> <td align="right"> 1.14 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197938 </td> <td> A </td> <td> AT </td> <td align="right">   0 </td> <td align="right">   3 </td> <td align="right">   0 </td> <td align="right"> 0.75 </td> <td align="right"> 1.50 </td> <td align="right"> 0.75 </td> <td align="right"> 3.00 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197957 </td> <td> G </td> <td> T </td> <td align="right">   5 </td> <td align="right">  12 </td> <td align="right">   0 </td> <td align="right"> 7.12 </td> <td align="right"> 7.76 </td> <td align="right"> 2.12 </td> <td align="right"> 5.07 </td> <td> FALSE </td> </tr>
   </table>

Compare to [brca1.hwe](./data/hwe/brca1.hwe) (see the [vcftools command line](./data/hwe/brca1.log) used to create this file).


```r
require(dplyr)
df <- read.table("./data/hwe/brca1.hwe", header=TRUE)
obsSplitCol <- "OBS.HOM1.HET.HOM2."
obsTemp <- read.table(text=as.character(df[, obsSplitCol]), sep = "/")
names(obsTemp) <- c("OBS_HOM1", "OBS_HET", "OBS_HOM2")
eSplitCol <- "E.HOM1.HET.HOM2."
eTemp <- read.table(text=as.character(df[, eSplitCol]), sep = "/")
names(eTemp) <- c("E_HOM1", "E_HET", "E_HOM2")
expectedResult <- cbind(cbind(df[setdiff(names(df), c(obsSplitCol,eSplitCol))], obsTemp), eTemp)
# Convert to zero-based coordinates
expectedResult <- mutate(expectedResult, POS = POS - 1)
```

How many results do the two results have in common?

```r
nrow(inner_join(result, expectedResult, by=c("CHR", "POS", "OBS_HOM1", "OBS_HET", "OBS_HOM2")))
```

```
## [1] 305
```

Which results were only identified by BigQuery?

```r
onlyBQ <- anti_join(result, expectedResult, , by=c("CHR", "POS", "OBS_HOM1", "OBS_HET", "OBS_HOM2"))
print(xtable(arrange(onlyBQ, CHR, POS)), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:46 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ref </th> <th> alt </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th> <th> ChiSq </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   7 </td> <td align="right">   0 </td> <td align="right"> 10.72 </td> <td align="right"> 5.56 </td> <td align="right"> 0.72 </td> <td align="right"> 1.14 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td> CT </td> <td> C </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> <td align="right"> 1.00 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204837 </td> <td> A </td> <td> T </td> <td align="right">  14 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 14.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204839 </td> <td> A </td> <td> T </td> <td align="right">  14 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 14.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td> T </td> <td> A </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">  12 </td> <td align="right"> 0.02 </td> <td align="right"> 0.96 </td> <td align="right"> 12.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211485 </td> <td> CACA </td> <td> C </td> <td align="right">   0 </td> <td align="right">   7 </td> <td align="right">   0 </td> <td align="right"> 1.75 </td> <td align="right"> 3.50 </td> <td align="right"> 1.75 </td> <td align="right"> 7.00 </td> <td> TRUE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211485 </td> <td> C </td> <td> CACA </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right"> 0.13 </td> <td align="right"> 0.75 </td> <td align="right"> 1.13 </td> <td align="right"> 0.23 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41214209 </td> <td> A </td> <td> T </td> <td align="right">  16 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 16.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41214210 </td> <td> A </td> <td> C </td> <td align="right">  16 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 16.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td> T </td> <td> TAA </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">   3 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right"> 3.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td> T </td> <td> TA </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">   3 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right"> 3.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226740 </td> <td> T </td> <td> G </td> <td align="right">  16 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 16.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239915 </td> <td> T </td> <td> A </td> <td align="right">  16 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 16.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242077 </td> <td> G </td> <td> A </td> <td align="right">  13 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 13.02 </td> <td align="right"> 0.96 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252648 </td> <td> T </td> <td> A </td> <td align="right">  13 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 13.02 </td> <td align="right"> 0.96 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252693 </td> <td> T </td> <td> A </td> <td align="right">  16 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 16.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td> A </td> <td> T </td> <td align="right">  16 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 16.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252695 </td> <td> A </td> <td> T </td> <td align="right">   2 </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right"> 3.50 </td> <td align="right"> 7.00 </td> <td align="right"> 3.50 </td> <td align="right"> 2.57 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696 </td> <td> T </td> <td> C </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 1.13 </td> <td align="right"> 0.75 </td> <td align="right"> 0.13 </td> <td align="right"> 0.23 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696 </td> <td> T </td> <td> A </td> <td align="right">   1 </td> <td align="right">   9 </td> <td align="right">   3 </td> <td align="right"> 2.33 </td> <td align="right"> 6.35 </td> <td align="right"> 4.33 </td> <td align="right"> 2.27 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256091 </td> <td> A </td> <td> G </td> <td align="right">  16 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 16.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256094 </td> <td> A </td> <td> G </td> <td align="right">  16 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 16.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256097 </td> <td> G </td> <td> A </td> <td align="right">  13 </td> <td align="right">   3 </td> <td align="right">   0 </td> <td align="right"> 13.14 </td> <td align="right"> 2.72 </td> <td align="right"> 0.14 </td> <td align="right"> 0.17 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256100 </td> <td> A </td> <td> G </td> <td align="right">  14 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 14.06 </td> <td align="right"> 1.88 </td> <td align="right"> 0.06 </td> <td align="right"> 0.07 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256102 </td> <td> G </td> <td> A </td> <td align="right">  12 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 12.25 </td> <td align="right"> 3.50 </td> <td align="right"> 0.25 </td> <td align="right"> 0.33 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41271293 </td> <td> A </td> <td> G </td> <td align="right">  16 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 16.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   6 </td> <td align="right">   0 </td> <td align="right"> 10.56 </td> <td align="right"> 4.88 </td> <td align="right"> 0.56 </td> <td align="right"> 0.85 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 10.02 </td> <td align="right"> 0.95 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
   </table>

Note vcftools appears to skip variants with single allele genotypes:
```
zgrep 41242078 platinum_genomes_brca1_expanded_merged.vcf.gz 
chr17  41242078  .  G	A	143	LowGQX;TruthSensitivityTranche99.90to100.00;LowQD;SiteConflict	BLOCKAVG_min30p3a;MQ=57;MQ0=0;BaseQRankSum=0.781;Dels=0.3;FS=1.561;HRun=11;HaplotypeScore=77.7361;MQRankSum=0.093;QD=2.01;ReadPosRankSum=-2.871;SB=-45.67;VQSLOD=-1.8762;culprit=QD;set=FilteredInAll;DP=425;AF=0.5;AN=25;AC=1	GT:DP:GQX:MQ:AD:GQ:PL:VF	0/0:57:99:59:.:.:.:.	0:27:25:57:26:25.38:.:.	0/0:51:99:57:.:.:.:.	0/1:50:99:59:42,8:99:173,0,1238:0.16	0/0:46:99:59:.:.:.:.	0/0:44:99:60:.:.:.:.	.:46:.:59:40,6:.:.:.	0/0:40:85:59:.:.:.:.	0/0:40:85:59:.:.:.:.	0/0:63:99:58:.:.:.:.	0:42:2:58:37:1.58:.:.	0:33:0:57:29:0.03:.:.	.:44:.:58:31,12:.:.:.	0/0:44:90:58:.:.:.:.	0/0:40:87:58:.:.:.:.	.:44:.:57:39,5:.:.:.	0/0:55:99:59:.:.:.:.
```

Which results were only identified by vcftools?

```r
onlyVcftools <- anti_join(expectedResult, result, , by=c("CHR", "POS", "OBS_HOM1", "OBS_HET", "OBS_HOM2"))
print(xtable(arrange(onlyVcftools, CHR, POS)), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:46 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ChiSq </th> <th> P </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407.00 </td> <td align="right"> 1.39 </td> <td align="right"> 0.53 </td> <td align="right">   8 </td> <td align="right">   7 </td> <td align="right">   0 </td> <td align="right"> 8.82 </td> <td align="right"> 5.37 </td> <td align="right"> 0.82 </td> </tr>
   </table>

Retrieving the gVCF data for the results identified only by vcftools:

```r
having <- paste("start <= ", onlyVcftools$POS,
                "AND",
                "end >= ", onlyVcftools$POS+1)
result <- DisplayAndDispatchQuery("./sql/examine-data.sql",
                                  replacements=c(table_replacement,
                                                 "_HAVING_"=having))
```

```
# Examine the data for particular calls.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS gt,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
  GROUP_CONCAT(STRING(call.pl)) WITHIN call AS likelihood,
FROM
  [genomics-public-data:platinum_genomes.variants]
WHERE
  reference_name = 'chr17'
HAVING
  start <=  41196407 AND end >=  41196408
ORDER BY
  start,
  end,
  call.call_set_name
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:48 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> call_call_set_name </th> <th> gt </th> <th> quality </th> <th> filter </th> <th> likelihood </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196196 </td> <td align="right"> 41196429 </td> <td> A </td> <td>  </td> <td> NA12891 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196228 </td> <td align="right"> 41196606 </td> <td> T </td> <td>  </td> <td> NA12882 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196313 </td> <td align="right"> 41196746 </td> <td> G </td> <td>  </td> <td> NA12886 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196337 </td> <td align="right"> 41196620 </td> <td> T </td> <td>  </td> <td> NA12881 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196339 </td> <td align="right"> 41196489 </td> <td> C </td> <td>  </td> <td> NA12893 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196349 </td> <td align="right"> 41196417 </td> <td> A </td> <td>  </td> <td> NA12877 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196355 </td> <td align="right"> 41196477 </td> <td> A </td> <td>  </td> <td> NA12879 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196376 </td> <td align="right"> 41196621 </td> <td> T </td> <td>  </td> <td> NA12890 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196390 </td> <td align="right"> 41196469 </td> <td> C </td> <td>  </td> <td> NA12884 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196396 </td> <td align="right"> 41196814 </td> <td> C </td> <td>  </td> <td> NA12885 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td> NA12878 </td> <td> 0,1 </td> <td align="right"> 733.47 </td> <td> PASS </td> <td> 763,0,946 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td> NA12880 </td> <td> 0,1 </td> <td align="right"> 733.47 </td> <td> PASS </td> <td> 1089,0,1199 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td> NA12883 </td> <td> 0,1 </td> <td align="right"> 733.47 </td> <td> PASS </td> <td> 1372,0,733 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td> NA12887 </td> <td> 0,1 </td> <td align="right"> 733.47 </td> <td> PASS </td> <td> 1140,0,707 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td> NA12888 </td> <td> 0,1 </td> <td align="right"> 733.47 </td> <td> PASS </td> <td> 1169,0,1265 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td> NA12889 </td> <td> 0,1 </td> <td align="right"> 733.47 </td> <td> PASS </td> <td> 1042,0,911 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td> NA12892 </td> <td> 0,1 </td> <td align="right"> 733.47 </td> <td> PASS </td> <td> 964,0,1187 </td> </tr>
   </table>

It appears that with BigQuery we are computing HWE for all the same variants as vcftools and the expected and Chi-Squared values are only slightly different.

See also: the [gVCF version of this query](./sql/hardy-weinberg-brca1.sql), which is close but only works for SNPs and needs a RIGHT OUTER JOIN to compute values for variants for which all the samples have the variant.

===============

Check Transition-Transversion Ratio
-----------------------------------

```r
result <- DisplayAndDispatchQuery("./sql/ti-tv-ratio.sql",
                                  replacements=table_replacement)
```

```
# Compute the Ti/Tv ratio of the 1,000 Genomes dataset.
SELECT
  transitions,
  transversions,
  transitions/transversions AS titv,
  COUNT
FROM (
  SELECT
    SUM(IF(mutation IN ('A->G',
          'G->A',
          'C->T',
          'T->C'),
        INTEGER(num_snps),
        INTEGER(0))) AS transitions,
    SUM(IF(mutation IN ('A->C',
          'C->A',
          'G->T',
          'T->G',
          'A->T',
          'T->A',
          'C->G',
          'G->C'),
        INTEGER(num_snps),
        INTEGER(0))) AS transversions,
        COUNT(mutation) AS COUNT
  FROM (
    SELECT
      CONCAT(reference_bases,
        CONCAT(STRING('->'),
          alternate_bases)) AS mutation,
      COUNT(alternate_bases) AS num_snps,
    FROM
      [genomics-public-data:platinum_genomes.variants]
    WHERE
      reference_name = 'chr17'
      AND start BETWEEN 41196311
      AND 41277499     
      AND LENGTH(alternate_bases) == 1
      AND LENGTH(reference_bases) == 1
    GROUP BY
      mutation,
    ORDER BY
      mutation))
```
The result:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:51 2015 -->
<table border=1>
<tr> <th> transitions </th> <th> transversions </th> <th> titv </th> <th> COUNT </th>  </tr>
  <tr> <td align="right"> 143 </td> <td align="right"> 132 </td> <td align="right"> 1.08 </td> <td align="right">  12 </td> </tr>
   </table>


Let's compare this to what we get from vcftools.  For information about the vcftools command see the [log](./titv/platinum_genomes_brca1_expanded_merged.log).  

```r
expectedResult <- read.table("./data/titv/platinum_genomes_brca1_expanded_merged.TsTv.summary", header=TRUE)
```
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:51 2015 -->
<table border=1>
<tr> <th> MODEL </th> <th> COUNT </th>  </tr>
  <tr> <td> AC </td> <td align="right">  48 </td> </tr>
  <tr> <td> AG </td> <td align="right">  75 </td> </tr>
  <tr> <td> AT </td> <td align="right">  37 </td> </tr>
  <tr> <td> CG </td> <td align="right">  11 </td> </tr>
  <tr> <td> CT </td> <td align="right">  66 </td> </tr>
  <tr> <td> GT </td> <td align="right">  34 </td> </tr>
  <tr> <td> Ts </td> <td align="right"> 141 </td> </tr>
  <tr> <td> Tv </td> <td align="right"> 130 </td> </tr>
   </table>
We can see that with BigQuery we get 143 transition mutations, and 132 transversion mutations.  Using vcftools we get two less of each category, 141 transitions and 130 transversions.  


First, we need to get the specific mutations from BigQuery.

```r
result <- DisplayAndDispatchQuery("./sql/ti-tv-variants.sql",
                                  replacements=table_replacement)
```

```
#Show the location of each mutation
SELECT
      reference_name,
      start,
      reference_bases,
      alternate_bases
    FROM
      [genomics-public-data:platinum_genomes.variants]
    WHERE
      reference_name = 'chr17'
      AND start BETWEEN 41196311
      AND 41277499     
      AND LENGTH(alternate_bases) == 1
      AND LENGTH(reference_bases) == 1
    GROUP BY
      reference_name,
      start,
      end,
      reference_bases,
      alternate_bases
    ORDER BY
      start
```
Here's the first few variants reported by BigQuery:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:52 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> reference_bases </th> <th> alternate_bases </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td> G </td> <td> A </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196840 </td> <td> G </td> <td> T </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td> C </td> <td> A </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197957 </td> <td> G </td> <td> T </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197958 </td> <td> A </td> <td> T </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198182 </td> <td> A </td> <td> C </td> </tr>
   </table>
Next, import the [transition and transversion mutations](./data/titv/brca.titv.csv) from the raw vcf file.  These were pulled out of the original vcf using a [custom perl script](./titv/pull_titv).  vcftools does not output a new vcf file with only the transitions and transversions so we need to use a proxy method.

```r
expectedResult <- read.csv("./data/titv/brca.titv.csv", header=FALSE)
# Set column names
names(expectedResult) <- c("reference_name","start","reference_bases","alternate_bases")
# Convert to zero-based coordinates
expectedResult <- mutate(expectedResult, start = start - 1)
```
Which variants were only identified by BigQuery?

```r
onlyBQ <- anti_join(result, expectedResult)
```

```
## Joining by: c("reference_name", "start", "reference_bases", "alternate_bases")
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon Jan 12 18:08:52 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> reference_bases </th> <th> alternate_bases </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td> G </td> <td> C </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td> G </td> <td> A </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696 </td> <td> T </td> <td> C </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696 </td> <td> T </td> <td> A </td> </tr>
   </table>
The four variants identified only by BigQuery are from two positions, each having multiple alternate alleles.  The perl script used to identify positions with transitions and transversion did not account for this, it is likely vcftools does not either.  Because vcftools does not output we cannot say for sure whether these are the 4 variants that vcftools missed, but it is a safe assumption given that we have two additional transitions and two transversions from these positions (which matches the discrepencey we originally had) as well as a logical reason for a bug in vcftools.

Let's double check that no variants were identified only by vcftools.

```r
nrow(anti_join(expectedResult, result))
```

```
## Joining by: c("reference_name", "start", "reference_bases", "alternate_bases")
```

```
## [1] 0
```
