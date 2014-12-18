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
project <- "stanford.edu:gbsc-stanford-google"                   # put your projectID here
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
table_replacement <- list("_THE_TABLE_"="genomics-public-data:platinum_genomes.variants")
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
<<<<<<< HEAD
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Wed Dec 17 16:05:56 2014 -->
=======
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec 15 18:28:01 2014 -->
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
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

<<<<<<< HEAD
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Wed Dec 17 16:05:58 2014 -->
=======
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec 15 18:28:03 2014 -->
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
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

<<<<<<< HEAD
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Wed Dec 17 16:06:00 2014 -->
=======
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec 15 18:28:04 2014 -->
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
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

<<<<<<< HEAD
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Wed Dec 17 16:06:03 2014 -->
=======
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec 15 18:28:08 2014 -->
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
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
require(dplyr)
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
onlyBQ
```

```
##   CHROM      POS SINGLETON_DOUBLETON  REF     ALT    INDV genotype
## 1 chr17 41211485                   S CACA CACAACA NA12878    "1,2"
##   num_samples_with_variant
## 1                        1
```

Which singletons were only identified by vcftools?

```r
onlyVcftools <- anti_join(expectedResult, result)
```

```
## Joining by: c("CHROM", "POS", "SINGLETON_DOUBLETON", "INDV")
```

```r
onlyVcftools
```

```
##    CHROM      POS SINGLETON_DOUBLETON ALLELE    INDV
## 1  chr17 41252694                   S    AAT NA12886
## 2  chr17 41204841                   S      T NA12888
## 3  chr17 41196320                   D      T NA12886
## 4  chr17 41196319                   D      T NA12886
## 5  chr17 41196318                   D      G NA12886
## 6  chr17 41196317                   D      T NA12886
## 7  chr17 41196315                   D      A NA12886
## 8  chr17 41196316                   D      G NA12886
## 9  chr17 41196314                   D      A NA12886
## 10 chr17 41196313                   D      G NA12886
```

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
  start = 41252694 OR start = 41204841 OR start = 41196320 OR start = 41196319 OR start = 41196318 OR start = 41196317 OR start = 41196315 OR start = 41196316 OR start = 41196314 OR start = 41196313
ORDER BY
  start,
  end,
  call.call_set_name
```

<<<<<<< HEAD
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Wed Dec 17 16:06:04 2014 -->
=======
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec 15 18:28:10 2014 -->
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> call_call_set_name </th> <th> gt </th> <th> quality </th> <th> filter </th> <th> likelihood </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196313 </td> <td align="right"> 41196746 </td> <td> G </td> <td>  </td> <td> NA12886 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
<<<<<<< HEAD
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12878 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 1025,63,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12881 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 725,63,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12892 </td> <td> -1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12882 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 917,54,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12893 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 859,72,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12890 </td> <td> -1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12884 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 943,81,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12886 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 785,38,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12891 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 877,78,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12877 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 793,44,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12889 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 815,54,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12885 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 809,69,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12887 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 754,39,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12883 </td> <td> -1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12879 </td> <td> 1,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 1175,48,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12880 </td> <td> -1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204841 </td> <td align="right"> 41204842 </td> <td> T </td> <td> A </td> <td> NA12888 </td> <td> 0,1 </td> <td align="right"> 759.64 </td> <td> TruthSensitivityTranche99.90to100.00 </td> <td> 359,0,25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12889 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12883 </td> <td> 0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12880 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12884 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12890 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12878 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12892 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td>  </td> <td> NA12881 </td> <td> 0,0 </td> <td align="right"> 48.05 </td> <td> LowGQX </td> <td>  </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td> T </td> <td> NA12877 </td> <td> 0,1 </td> <td align="right"> 143.44 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 173,0,874 </td> </tr>
=======
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
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252697 </td> <td> AAT </td> <td> A </td> <td> NA12886 </td> <td> 0,1 </td> <td align="right"> 648.06 </td> <td> LowGQX </td> <td> 823,0,38 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252697 </td> <td> AAT </td> <td> A </td> <td> NA12884 </td> <td> 1,1 </td> <td align="right"> 648.06 </td> <td> LowGQX </td> <td> 690,17,0 </td> </tr>
   </table>

It appears that they correspond either to:
* a reference-matching block, so not actually a singleton and just perhaps violating an assumption in the vcftools code
* or a non-singleon variant, perhaps due to a problem in converting the gVCF data to all-positions VCF via gvcftools?

Check Hardy-Weinberg Equilibrium
-----------------------------------

```r
result <- DisplayAndDispatchQuery("./sql/hardy-weinberg-brca1.sql",
                                  replacements=table_replacement)
```

```
SELECT
  vars.reference_name AS CHR,
  vars.start AS POS,
  reference_bases AS ref,
  alternate_bases AS alt,
  SUM(refs.HOM_REF) + vars.HOM_REF AS OBS_HOM1,
  vars.HET AS OBS_HET,
  vars.HOM_ALT AS OBS_HOM2,
  
  # variables not supported :(
  #@sample_count := (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT), 
  
  # Expected AA
  # p^2
  # ((COUNT(AA) + (COUNT(Aa)/2) / 
  #  TOTAL_SAMPLE_COUNT) ^ 2) * 
  #  TOTAL_SAMPLE_COUNT
  
  POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)
    AS EXP_HOM1,
  
  # Expected Aa
  # 2pq
  # 2 *
  # ((COUNT(AA) + (COUNT(Aa)/2) / 
  #  TOTAL_SAMPLE_COUNT) ^ 2) * 
  # (COUNT(aa) + (COUNT(Aa)/2) / 
  #  TOTAL_SAMPLE_COUNT) ^ 2 *
  # TOTAL_SAMPLE_COUNT
  
  2 *
  ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
  ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)
    AS EXP_HET,
    
  # Expected aa
  # q^2
  # (COUNT(aa) + (COUNT(Aa)/2) / 
  #  TOTAL_SAMPLE_COUNT) ^ 2 *
  #  TOTAL_SAMPLE_COUNT
  
  POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT) 
    AS EXP_HOM2,
  
  ########################
  # Chi square
  # SUM(((Observed - Expected)^2) / Expected )
    
    # AA chisq value
    # Numerator 
    POW(
    # Observed
    (SUM(refs.HOM_REF) + vars.HOM_REF) 
      - 
    # Expected
    (POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2) 
    /
    # Denominator 
    # Expected  
      (POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
      (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
      (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    
    +
    
    # Aa chisq value
    # Numerator
    POW(
    # Observed
    vars.HET 
      -
    # Expected
    (2 *
    ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2)
    / 
    # Denominator
    # Expected
    (2 *
    ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
      
    +
      
    #aa chisq value
    # Numerator
    POW(
    # Observed
    vars.HOM_ALT - 
    # Expected
    (POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2)
    / 
    # Denominator
    # Expected
    (POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    
    AS CHISQ,
    #########################
    
    #########################
    # p-value boolean
    # Need to recalculate chisq value and determine if it is greater than 
    # the cutoff for significance with two degrees of freedom, 5.991.
    IF ( 
    POW(
    # Observed
    (SUM(refs.HOM_REF) + vars.HOM_REF) 
      - 
    # Expected
    (POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2) 
    /
    # Denominator 
    # Expected  
      (POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
      (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
      (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    
    +
    
    # Aa chisq value
    # Numerator
    POW(
    # Observed
    vars.HET 
      -
    # Expected
    (2 *
    ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2)
    / 
    # Denominator
    # Expected
    (2 *
    ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
      
    +
      
    #aa chisq value
    # Numerator
    POW(
    # Observed
    vars.HOM_ALT - 
    # Expected
    (POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2)
    / 
    # Denominator
    # Expected
    (POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))

    > 5.991, "TRUE", "FALSE") AS PVALUE,
    #########################
  
FROM (
<<<<<<< HEAD
	    # Constrain the left hand side of the _join to reference-matching blocks.
	  SELECT
	    reference_name,
	    start,
	    END,
	    SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
	  FROM
	    [genomics-public-data:platinum_genomes.variants]
	  WHERE
	    reference_name = 'chr17'
	  OMIT
	    RECORD IF EVERY(alternate_bases IS NOT NULL)
	    ) AS refs
	JOIN (
		  SELECT
		    reference_name,
		    start,
		    END,
		    reference_bases,
		    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
		    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
		    SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
		    SUM(EVERY(1 = call.genotype)) WITHIN call AS HOM_ALT,
		    SUM(SOME(0 = call.genotype) AND SOME(1 = call.genotype)) WITHIN call AS HET,
        
		  FROM
		    [genomics-public-data:platinum_genomes.variants]
		  WHERE
		    reference_name = 'chr17'
		    AND start BETWEEN 41196311
		    AND 41277499
		#  OMIT call IF 2 != COUNT(call.genotype)
		  HAVING
		    # Skip ref-matching blocks, 1/2 genotypes, and non-SNP variants
		    num_alts = 1
		    AND reference_bases IN ('A','C','G','T')
		    AND alternate_bases IN ('A','C','G','T')
		    ) AS vars
		  # The _join criteria _is complicated since we are trying to see if a variant
		  # overlaps a reference-matching interval.
		ON
		  vars.reference_name = refs.reference_name
		WHERE
		  refs.start <= vars.start
		  AND refs.END >= vars.start+1
		GROUP BY
		  CHR,
		  POS,
		  ref,
		  alt,
		  vars.HOM_REF,
		  OBS_HET,
		  OBS_HOM2,
      vars.HET,
      vars.HOM_ALT
		ORDER BY
		  CHR,
		  POS,
		  ref,
		  alt
=======
    # Constrain the left hand side of the _join to reference-matching blocks.
  SELECT
    reference_name,
    start,
    END,
    SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
  FROM
    [genomics-public-data:platinum_genomes.variants]
  WHERE
    reference_name = 'chr17'
  OMIT
    RECORD IF EVERY(alternate_bases IS NOT NULL)
    ) AS refs
JOIN (
  SELECT
    reference_name,
    start,
    END,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
    SUM(EVERY(1 = call.genotype)) WITHIN call AS HOM_ALT,
    SUM(SOME(0 = call.genotype) AND SOME(1 = call.genotype)) WITHIN call AS HET,
  FROM
    [genomics-public-data:platinum_genomes.variants]
  WHERE
    reference_name = 'chr17'
    AND start BETWEEN 41196311
    AND 41277499
#  OMIT call IF 2 != COUNT(call.genotype)
  HAVING
    # Skip ref-matching blocks, 1/2 genotypes, and non-SNP variants
    num_alts = 1
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
    ) AS vars
  # The _join criteria _is complicated since we are trying to see if a variant
  # overlaps a reference-matching interval.
ON
  vars.reference_name = refs.reference_name
WHERE
  refs.start <= vars.start
  AND refs.END >= vars.start+1
GROUP BY
  CHR,
  POS,
  ref,
  alt,
  vars.HOM_REF,
  OBS_HET,
  OBS_HOM2
ORDER BY
  CHR,
  POS,
  ref,
  alt
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
```
Number of rows returned by this query: 271.

Displaying the first few results:
<<<<<<< HEAD
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Wed Dec 17 16:06:07 2014 -->
=======
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec 15 18:28:14 2014 -->
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ref </th> <th> alt </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> EXP_HOM1 </th> <th> EXP_HET </th> <th> EXP_HOM2 </th> <th> CHISQ </th> <th> PVALUE </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   7 </td> <td align="right">   0 </td> <td align="right"> 10.72 </td> <td align="right"> 5.56 </td> <td align="right"> 0.72 </td> <td align="right"> 1.14 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196840 </td> <td> G </td> <td> T </td> <td align="right">  15 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 15.06 </td> <td align="right"> 1.88 </td> <td align="right"> 0.06 </td> <td align="right"> 0.07 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   7 </td> <td align="right">   0 </td> <td align="right"> 10.72 </td> <td align="right"> 5.56 </td> <td align="right"> 0.72 </td> <td align="right"> 1.14 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197957 </td> <td> G </td> <td> T </td> <td align="right">   5 </td> <td align="right">  12 </td> <td align="right">   0 </td> <td align="right"> 7.12 </td> <td align="right"> 7.76 </td> <td align="right"> 2.12 </td> <td align="right"> 5.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197958 </td> <td> A </td> <td> T </td> <td align="right">  16 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 16.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198182 </td> <td> A </td> <td> C </td> <td align="right">  11 </td> <td align="right">   6 </td> <td align="right">   0 </td> <td align="right"> 11.53 </td> <td align="right"> 4.94 </td> <td align="right"> 0.53 </td> <td align="right"> 0.78 </td> <td> FALSE </td> </tr>
   </table>

Compare to [brca1.hwe](./data/singletons/brca1.hwe) (see the [vcftools command line](./data/hwe/brca1.log) used to create this file).


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
nrow(inner_join(result, expectedResult))
```

```
## Joining by: c("CHR", "POS", "OBS_HOM1", "OBS_HET", "OBS_HOM2")
```

```
## [1] 249
```

Which results were only identified by BigQuery?

```r
onlyBQ <- anti_join(result, expectedResult, by=c("CHR", "POS"))
onlyBQ
```

```
<<<<<<< HEAD
##      CHR      POS ref alt OBS_HOM1 OBS_HET OBS_HOM2  EXP_HOM1   EXP_HET
## 1  chr17 41256102   G   A       12       4        0 12.250000 3.5000000
## 2  chr17 41256100   A   G       14       2        0 14.062500 1.8750000
## 3  chr17 41256097   G   A       13       3        0 13.140625 2.7187500
## 4  chr17 41256094   A   G       16       1        0 16.014706 0.9705882
## 5  chr17 41256091   A   G       16       1        0 16.014706 0.9705882
## 6  chr17 41252696   T   A        1       9        3  2.326923 6.3461538
## 7  chr17 41252696   T   C        1       1        0  1.125000 0.7500000
## 8  chr17 41252695   A   T        2      10        2  3.500000 7.0000000
## 9  chr17 41252693   T   A       16       1        0 16.014706 0.9705882
## 10 chr17 41252648   T   A       13       1        0 13.017857 0.9642857
## 11 chr17 41271293   A   G       16       0        0 16.000000 0.0000000
## 12 chr17 41242077   G   A       13       1        0 13.017857 0.9642857
## 13 chr17 41239915   T   A       16       0        0 16.000000 0.0000000
## 14 chr17 41226740   T   G       16       0        0 16.000000 0.0000000
## 15 chr17 41273094   G   A       10       6        0 10.562500 4.8750000
## 16 chr17 41273094   G   C       10       1        0 10.022727 0.9545455
## 17 chr17 41214210   A   C       16       0        0 16.000000 0.0000000
## 18 chr17 41214209   A   T       16       0        0 16.000000 0.0000000
## 19 chr17 41204837   A   T       14       0        0 14.000000 0.0000000
##      EXP_HOM2      CHISQ PVALUE
## 1  0.25000000 0.32653061  FALSE
## 2  0.06250000 0.07111111  FALSE
## 3  0.14062500 0.17122473  FALSE
## 4  0.01470588 0.01561065  FALSE
## 5  0.01470588 0.01561065  FALSE
## 6  4.32692308 2.27338843  FALSE
## 7  0.12500000 0.22222222  FALSE
## 8  3.50000000 2.57142857  FALSE
## 9  0.01470588 0.01561065  FALSE
## 10 0.01785714 0.01920439  FALSE
## 11 0.00000000         NA  FALSE
## 12 0.01785714 0.01920439  FALSE
## 13 0.00000000         NA  FALSE
## 14 0.00000000         NA  FALSE
## 15 0.56250000 0.85207101  FALSE
## 16 0.02272727 0.02494331  FALSE
## 17 0.00000000         NA  FALSE
## 18 0.00000000         NA  FALSE
## 19 0.00000000         NA  FALSE
=======
##      CHR      POS ref alt OBS_HOM1 OBS_HET OBS_HOM2
## 1  chr17 41273094   G   A       10       6        0
## 2  chr17 41273094   G   C       10       1        0
## 3  chr17 41271293   A   G       16       0        0
## 4  chr17 41256102   G   A       12       4        0
## 5  chr17 41256100   A   G       14       2        0
## 6  chr17 41256097   G   A       13       3        0
## 7  chr17 41256094   A   G       16       1        0
## 8  chr17 41256091   A   G       16       1        0
## 9  chr17 41252695   A   T        2      10        2
## 10 chr17 41252693   T   A       16       1        0
## 11 chr17 41252648   T   A       13       1        0
## 12 chr17 41242077   G   A       13       1        0
## 13 chr17 41252696   T   A        1       9        3
## 14 chr17 41252696   T   C        1       1        0
## 15 chr17 41226740   T   G       16       0        0
## 16 chr17 41214210   A   C       16       0        0
## 17 chr17 41214209   A   T       16       0        0
## 18 chr17 41239915   T   A       16       0        0
## 19 chr17 41204837   A   T       14       0        0
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
```

Note vcftools appears to skip variants with single allele genotypes:
```
zgrep 41242078 platinum_genomes_brca1_expanded_merged.vcf.gz 
chr17  41242078  .  G	A	143	LowGQX;TruthSensitivityTranche99.90to100.00;LowQD;SiteConflict	BLOCKAVG_min30p3a;MQ=57;MQ0=0;BaseQRankSum=0.781;Dels=0.3;FS=1.561;HRun=11;HaplotypeScore=77.7361;MQRankSum=0.093;QD=2.01;ReadPosRankSum=-2.871;SB=-45.67;VQSLOD=-1.8762;culprit=QD;set=FilteredInAll;DP=425;AF=0.5;AN=25;AC=1	GT:DP:GQX:MQ:AD:GQ:PL:VF	0/0:57:99:59:.:.:.:.	0:27:25:57:26:25.38:.:.	0/0:51:99:57:.:.:.:.	0/1:50:99:59:42,8:99:173,0,1238:0.16	0/0:46:99:59:.:.:.:.	0/0:44:99:60:.:.:.:.	.:46:.:59:40,6:.:.:.	0/0:40:85:59:.:.:.:.	0/0:40:85:59:.:.:.:.	0/0:63:99:58:.:.:.:.	0:42:2:58:37:1.58:.:.	0:33:0:57:29:0.03:.:.	.:44:.:58:31,12:.:.:.	0/0:44:90:58:.:.:.:.	0/0:40:87:58:.:.:.:.	.:44:.:57:39,5:.:.:.	0/0:55:99:59:.:.:.:.
```

Which results were only identified by vcftools?

```r
onlyVcftools <- anti_join(expectedResult, result, by=c("CHR", "POS"))
onlyVcftools
```

```
##      CHR      POS     ChiSq        P OBS_HOM1 OBS_HET OBS_HOM2 E_HOM1
## 1  chr17 41270777  1.000000 1.000000        0       1        0   0.25
## 2  chr17 41268207  7.000000 0.037296        0       7        0   1.75
## 3  chr17 41267517  5.000000 0.126984        0       5        0   1.25
<<<<<<< HEAD
## 4  chr17 41264754  7.000000 0.037296        0       7        0   1.75
## 5  chr17 41264750  7.000000 0.037296        0       7        0   1.75
=======
## 4  chr17 41266406  4.000000 0.314286        0       4        0   1.00
## 5  chr17 41264754  7.000000 0.037296        0       7        0   1.75
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
## 6  chr17 41264742  7.000000 0.037296        0       7        0   1.75
## 7  chr17 41264739  7.000000 0.037296        0       7        0   1.75
## 8  chr17 41264110  1.000000 1.000000        0       1        0   0.25
## 9  chr17 41259078  1.000000 1.000000        0       1        0   0.25
## 10 chr17 41256088  1.000000 1.000000        0       1        0   0.25
## 11 chr17 41256073  5.000000 0.126984        0       5        0   1.25
<<<<<<< HEAD
## 12 chr17 41258134  3.000000 0.400000        0       3        0   0.75
## 13 chr17 41254964  3.644628 0.176471        0       7        2   1.36
## 14 chr17 41252590  3.000000 0.400000        0       3        0   0.75
## 15 chr17 41250677  7.000000 0.037296        0       7        0   1.75
## 16 chr17 41250220  1.000000 1.000000        0       1        0   0.25
## 17 chr17 41266406  4.000000 0.314286        0       4        0   1.00
## 18 chr17 41249362  7.000000 0.037296        0       7        0   1.75
## 19 chr17 41242074  6.000000 0.090909        0       6        0   1.50
## 20 chr17 41241567  4.000000 0.314286        0       4        0   1.00
## 21 chr17 41239914  7.000000 0.037296        0       7        0   1.75
## 22 chr17 41232343       NaN 1.000000        0       0       17   0.00
=======
## 12 chr17 41254964  3.644628 0.176471        0       7        2   1.36
## 13 chr17 41252692  1.000000 1.000000        0       1        0   0.25
## 14 chr17 41252590  3.000000 0.400000        0       3        0   0.75
## 15 chr17 41250677  7.000000 0.037296        0       7        0   1.75
## 16 chr17 41249362  7.000000 0.037296        0       7        0   1.75
## 17 chr17 41242074  6.000000 0.090909        0       6        0   1.50
## 18 chr17 41241567  4.000000 0.314286        0       4        0   1.00
## 19 chr17 41250220  1.000000 1.000000        0       1        0   0.25
## 20 chr17 41239914  7.000000 0.037296        0       7        0   1.75
## 21 chr17 41232343       NaN 1.000000        0       0       17   0.00
## 22 chr17 41230104  5.000000 0.126984        0       5        0   1.25
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
## 23 chr17 41229776       NaN 1.000000        0       0       17   0.00
## 24 chr17 41229759  7.000000 0.037296        0       7        0   1.75
## 25 chr17 41226735  7.000000 0.037296        0       7        0   1.75
## 26 chr17 41225653 10.000000 0.006906        0      10        0   2.50
## 27 chr17 41223537  1.000000 1.000000        0       1        0   0.25
## 28 chr17 41219906 17.000000 0.000056        0      17        0   4.25
<<<<<<< HEAD
## 29 chr17 41219852  1.000000 1.000000        0       1        0   0.25
## 30 chr17 41218817 17.000000 0.000056        0      17        0   4.25
## 31 chr17 41252645  1.000000 1.000000        0       1        0   0.25
## 32 chr17 41214208  7.000000 0.037296        0       7        0   1.75
## 33 chr17 41252692  1.000000 1.000000        0       1        0   0.25
## 34 chr17 41247121  7.000000 0.037296        0       7        0   1.75
## 35 chr17 41213601  5.000000 0.126984        0       5        0   1.25
## 36 chr17 41230104  5.000000 0.126984        0       5        0   1.25
## 37 chr17 41208190  2.000000 1.000000        0       2        0   0.50
## 38 chr17 41206760  6.000000 0.090909        0       6        0   1.50
## 39 chr17 41204835       NaN 1.000000        0       0        2   0.00
=======
## 29 chr17 41218817 17.000000 0.000056        0      17        0   4.25
## 30 chr17 41219852  1.000000 1.000000        0       1        0   0.25
## 31 chr17 41214208  7.000000 0.037296        0       7        0   1.75
## 32 chr17 41213601  5.000000 0.126984        0       5        0   1.25
## 33 chr17 41264750  7.000000 0.037296        0       7        0   1.75
## 34 chr17 41252645  1.000000 1.000000        0       1        0   0.25
## 35 chr17 41208190  2.000000 1.000000        0       2        0   0.50
## 36 chr17 41206760  6.000000 0.090909        0       6        0   1.50
## 37 chr17 41247121  7.000000 0.037296        0       7        0   1.75
## 38 chr17 41204835       NaN 1.000000        0       0        2   0.00
## 39 chr17 41256088  1.000000 1.000000        0       1        0   0.25
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
## 40 chr17 41204831       NaN 1.000000        0       0        1   0.00
## 41 chr17 41200703  6.000000 0.090909        0       6        0   1.50
## 42 chr17 41197938  3.000000 0.400000        0       3        0   0.75
##    E_HET E_HOM2
## 1   0.50   0.25
## 2   3.50   1.75
## 3   2.50   1.25
## 4   3.50   1.75
## 5   3.50   1.75
## 6   3.50   1.75
## 7   3.50   1.75
## 8   0.50   0.25
## 9   0.50   0.25
## 10  0.50   0.25
## 11  2.50   1.25
<<<<<<< HEAD
## 12  1.50   0.75
## 13  4.28   3.36
## 14  1.50   0.75
## 15  3.50   1.75
## 16  0.50   0.25
## 17  2.00   1.00
## 18  3.50   1.75
## 19  3.00   1.50
## 20  2.00   1.00
## 21  3.50   1.75
## 22  0.00  17.00
=======
## 12  4.28   3.36
## 13  0.50   0.25
## 14  1.50   0.75
## 15  3.50   1.75
## 16  3.50   1.75
## 17  3.00   1.50
## 18  2.00   1.00
## 19  0.50   0.25
## 20  3.50   1.75
## 21  0.00  17.00
## 22  2.50   1.25
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
## 23  0.00  17.00
## 24  3.50   1.75
## 25  3.50   1.75
## 26  5.00   2.50
## 27  0.50   0.25
## 28  8.50   4.25
<<<<<<< HEAD
## 29  0.50   0.25
## 30  8.50   4.25
## 31  0.50   0.25
## 32  3.50   1.75
## 33  0.50   0.25
## 34  3.50   1.75
## 35  2.50   1.25
## 36  2.50   1.25
## 37  1.00   0.50
## 38  3.00   1.50
## 39  0.00   2.00
=======
## 29  8.50   4.25
## 30  0.50   0.25
## 31  3.50   1.75
## 32  2.50   1.25
## 33  3.50   1.75
## 34  0.50   0.25
## 35  1.00   0.50
## 36  3.00   1.50
## 37  3.50   1.75
## 38  0.00   2.00
## 39  0.50   0.25
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
## 40  0.00   1.00
## 41  3.00   1.50
## 42  1.50   0.75
```

Retrieving the gVCF data for the results identified only by vcftools:

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
<<<<<<< HEAD
  start = 41270777 OR start = 41268207 OR start = 41267517 OR start = 41264754 OR start = 41264750 OR start = 41264742 OR start = 41264739 OR start = 41264110 OR start = 41259078 OR start = 41256088 OR start = 41256073 OR start = 41258134 OR start = 41254964 OR start = 41252590 OR start = 41250677 OR start = 41250220 OR start = 41266406 OR start = 41249362 OR start = 41242074 OR start = 41241567 OR start = 41239914 OR start = 41232343 OR start = 41229776 OR start = 41229759 OR start = 41226735 OR start = 41225653 OR start = 41223537 OR start = 41219906 OR start = 41219852 OR start = 41218817 OR start = 41252645 OR start = 41214208 OR start = 41252692 OR start = 41247121 OR start = 41213601 OR start = 41230104 OR start = 41208190 OR start = 41206760 OR start = 41204835 OR start = 41204831 OR start = 41200703 OR start = 41197938
=======
  start = 41270777 OR start = 41268207 OR start = 41267517 OR start = 41266406 OR start = 41264754 OR start = 41264742 OR start = 41264739 OR start = 41264110 OR start = 41259078 OR start = 41258134 OR start = 41256073 OR start = 41254964 OR start = 41252692 OR start = 41252590 OR start = 41250677 OR start = 41249362 OR start = 41242074 OR start = 41241567 OR start = 41250220 OR start = 41239914 OR start = 41232343 OR start = 41230104 OR start = 41229776 OR start = 41229759 OR start = 41226735 OR start = 41225653 OR start = 41223537 OR start = 41219906 OR start = 41218817 OR start = 41219852 OR start = 41214208 OR start = 41213601 OR start = 41264750 OR start = 41252645 OR start = 41208190 OR start = 41206760 OR start = 41247121 OR start = 41204835 OR start = 41256088 OR start = 41204831 OR start = 41200703 OR start = 41197938
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
ORDER BY
  start,
  end,
  call.call_set_name
```

Let's filter out indels and reference-matching blocks from this result:

```r
result <- filter(result, reference_bases %in% c('A','C','G','T') & alternate_bases %in% c('A','C','G','T'))
```

<<<<<<< HEAD
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Wed Dec 17 16:06:10 2014 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> END </th> <th> reference_bases </th> <th> alternate_bases </th> <th> call_call_set_name </th> <th> gt </th> <th> quality </th> <th> filter </th> <th> likelihood </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12885 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 250,0,708 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12882 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 213,0,1195 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12877 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 70,0,881 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12881 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 290,0,988 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12892 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 138,0,863 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12878 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 234,0,1058 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12893 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 354,0,1202 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12890 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 225,0,842 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12884 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 87,0,427 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12886 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 397,0,757 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12888 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 195,0,867 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12891 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 312,0,806 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12880 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 35,0,777 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12879 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 176,0,721 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12883 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 322,0,982 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12887 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 228,0,813 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12889 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 194,0,1006 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12893 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 293,0,667 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12878 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 261,0,326 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12891 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 235,0,204 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12882 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 208,0,366 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12877 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 211,0,369 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12889 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 165,0,353 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12885 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 210,0,370 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12887 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 176,0,226 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12883 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 222,0,298 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12879 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 134,0,489 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12880 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 266,0,171 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12888 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 295,0,336 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12886 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 305,0,557 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12884 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 224,0,196 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12890 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 161,0,404 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12881 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 227,0,492 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12892 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 302,0,260 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12883 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1733,129,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12879 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2753,211,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12880 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2543,199,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12887 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2200,169,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12885 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2283,172,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12889 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2528,193,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12877 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2002,153,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12882 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1880,144,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12891 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2059,156,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12888 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1834,141,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12886 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2491,193,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12884 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1712,132,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12890 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1787,138,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12893 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2434,184,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12878 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2157,163,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12892 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2246,172,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12881 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1811,144,0 </td> </tr>
=======
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec 15 18:28:18 2014 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> call_call_set_name </th> <th> gt </th> <th> quality </th> <th> filter </th> <th> likelihood </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12877 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 70,0,881 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12878 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 234,0,1058 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12879 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 176,0,721 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12880 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 35,0,777 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12881 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 290,0,988 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12882 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 213,0,1195 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12883 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 322,0,982 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12884 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 87,0,427 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12885 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 250,0,708 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12886 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 397,0,757 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12887 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 228,0,813 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12888 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 195,0,867 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12889 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 194,0,1006 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12890 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 225,0,842 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12891 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 312,0,806 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12892 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 138,0,863 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218817 </td> <td align="right"> 41218818 </td> <td> A </td> <td> C </td> <td> NA12893 </td> <td> 0,1 </td> <td align="right"> 40.36 </td> <td> TruthSensitivityTranche99.00to99.90,LowQD </td> <td> 354,0,1202 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12877 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 211,0,369 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12878 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 261,0,326 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12879 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 134,0,489 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12880 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 266,0,171 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12881 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 227,0,492 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12882 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 208,0,366 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12883 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 222,0,298 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12884 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 224,0,196 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12885 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 210,0,370 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12886 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 305,0,557 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12887 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 176,0,226 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12888 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 295,0,336 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12889 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 165,0,353 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12890 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 161,0,404 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12891 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 235,0,204 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12892 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 302,0,260 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219906 </td> <td align="right"> 41219907 </td> <td> T </td> <td> A </td> <td> NA12893 </td> <td> 0,1 </td> <td align="right"> 180.52 </td> <td> TruthSensitivityTranche99.00to99.90 </td> <td> 293,0,667 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12877 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2002,153,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12878 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2157,163,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12879 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2753,211,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12880 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2543,199,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12881 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1811,144,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12882 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1880,144,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12883 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1733,129,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12884 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1712,132,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12885 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2283,172,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12886 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2491,193,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12887 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2200,169,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12888 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1834,141,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12889 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2528,193,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12890 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 1787,138,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12891 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2059,156,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12892 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2246,172,0 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343 </td> <td align="right"> 41232344 </td> <td> G </td> <td> C </td> <td> NA12893 </td> <td> 1,1 </td> <td align="right"> 1968.62 </td> <td> PASS </td> <td> 2434,184,0 </td> </tr>
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
   </table>

It appears that the difference in results only returned by vcftools correspond either to:
* indels, which we are not examining here
* or variants for which all samples have the alternate for one or both alleles -> a RIGHT OUTER JOIN is needed in this query

TODO(deflaux):
* find a a way to work around the lack of RIGHT OUTER JOIN
* add Chi-Squared test to query from [this sample](https://github.com/googlegenomics/bigquery-examples/tree/master/1000genomes/data-stories/reproducing-hardy-weinberg-equilibrium)

Check Individual Heterozygosity
-----------------------------------

More calculations needed here . . . these are just counts.


```r
result <- DisplayAndDispatchQuery("./sql/homozygous-variants-brca1.sql",
                                  replacements=table_replacement)
```

```
# Individual Homozygosity
SELECT
  call.call_set_name AS INDV,
  SUM(first_allele = second_allele) AS O_HOM,
  COUNT(call.call_set_name) AS N_SITES,
FROM (
  SELECT
    call.call_set_name,
    NTH(1,
      call.genotype) WITHIN call AS first_allele,
    NTH(2,
      call.genotype) WITHIN call AS second_allele,
    COUNT(call.genotype) WITHIN call AS ploidy
  FROM
    [genomics-public-data:platinum_genomes.variants]
  WHERE
    reference_name = 'chr17'
    AND start BETWEEN 41196311
    AND 41277499
  OMIT RECORD IF EVERY(alternate_bases IS NULL)
  HAVING
    ploidy = 2
    )
GROUP BY
  INDV
ORDER BY
  INDV
```
Number of rows returned by this query: 17.

<<<<<<< HEAD
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Wed Dec 17 16:06:12 2014 -->
=======
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Mon Dec 15 18:28:20 2014 -->
>>>>>>> e06aff5a2555ca5c10d68d8d62a8aae9747dc65d
<table border=1>
<tr> <th> INDV </th> <th> O_HOM </th> <th> N_SITES </th>  </tr>
  <tr> <td> NA12877 </td> <td align="right">   3 </td> <td align="right">  27 </td> </tr>
  <tr> <td> NA12878 </td> <td align="right">   3 </td> <td align="right"> 198 </td> </tr>
  <tr> <td> NA12879 </td> <td align="right">   5 </td> <td align="right">  37 </td> </tr>
  <tr> <td> NA12880 </td> <td align="right">   3 </td> <td align="right"> 193 </td> </tr>
  <tr> <td> NA12881 </td> <td align="right">   3 </td> <td align="right">  42 </td> </tr>
  <tr> <td> NA12882 </td> <td align="right">   6 </td> <td align="right">  31 </td> </tr>
  <tr> <td> NA12883 </td> <td align="right">   3 </td> <td align="right"> 197 </td> </tr>
  <tr> <td> NA12884 </td> <td align="right">   6 </td> <td align="right">  35 </td> </tr>
  <tr> <td> NA12885 </td> <td align="right">   4 </td> <td align="right">  31 </td> </tr>
  <tr> <td> NA12886 </td> <td align="right">   4 </td> <td align="right">  29 </td> </tr>
  <tr> <td> NA12887 </td> <td align="right">   4 </td> <td align="right"> 211 </td> </tr>
  <tr> <td> NA12888 </td> <td align="right">   3 </td> <td align="right"> 205 </td> </tr>
  <tr> <td> NA12889 </td> <td align="right">   6 </td> <td align="right"> 198 </td> </tr>
  <tr> <td> NA12890 </td> <td align="right">   3 </td> <td align="right">  33 </td> </tr>
  <tr> <td> NA12891 </td> <td align="right">   5 </td> <td align="right">  37 </td> </tr>
  <tr> <td> NA12892 </td> <td align="right">   6 </td> <td align="right"> 209 </td> </tr>
  <tr> <td> NA12893 </td> <td align="right">   6 </td> <td align="right">  41 </td> </tr>
   </table>


