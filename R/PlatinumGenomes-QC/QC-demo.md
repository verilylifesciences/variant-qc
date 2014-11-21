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
project <- "genomics-public-data"                   # put your projectID here
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
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Thu Nov 20 18:54:14 2014 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> quality </th> <th> filter </th> <th> names </th> <th> num_samples </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td align="right"> 733.47 </td> <td> PASS </td> <td>  </td> <td align="right">   7 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td align="right"> 41196822 </td> <td> CT </td> <td> C </td> <td align="right"> 63.74 </td> <td> LowQD </td> <td>  </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td align="right"> 41196823 </td> <td> CTT </td> <td> C,CT </td> <td align="right"> 314.59 </td> <td> PASS </td> <td>  </td> <td align="right">   3 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196840 </td> <td align="right"> 41196841 </td> <td> G </td> <td> T </td> <td align="right"> 85.68 </td> <td> TruthSensitivityTranche99.90to100.00,LowQD </td> <td>  </td> <td align="right">   2 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td align="right"> 41197274 </td> <td> C </td> <td> A </td> <td align="right"> 1011.08 </td> <td> PASS </td> <td>  </td> <td align="right">   7 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197938 </td> <td align="right"> 41197939 </td> <td> A </td> <td> AT </td> <td align="right"> 86.95 </td> <td> LowQD </td> <td>  </td> <td align="right">   3 </td> </tr>
   </table>

Notes about this particular dataset:
* It is in gVCF format which adds complexity above and beyond [similar examples for the 1,000 Genomes dataset](https://github.com/googlegenomics/bigquery-examples/blob/master/1000genomes/sql/README.md).
* It is comprised only of SNPs and INDELs (contains no structural variants).
* The values for `alternate_bases` are just comprised of the letters A,C,G,T (e.g., contains no `<DEL>` values).
* It contains some contradictory calls.


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

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Thu Nov 20 18:54:16 2014 -->
<table border=1>
<tr> <th> number_of_variant_records </th> <th> alt_contains_no_special_characters </th> <th> max_ref_len </th> <th> max_alt_len </th>  </tr>
  <tr> <td align="right"> 12634588 </td> <td> TRUE </td> <td align="right">  56 </td> <td align="right">  47 </td> </tr>
   </table>
We see from the query results that there are no special charaters in alternate_bases and the maximum length is ~50 base pairs.


Here is one specific example of a contradictory call.

```r
result <- DisplayAndDispatchQuery("./sql/contradictory-data.sql",
                                  replacements=table_replacement)
```

```
# An example of contradictory data.
SELECT
  reference_name,
  start,
  END,
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
  AND start <= 41198773
  AND END >= 41198774
HAVING
  call.call_set_name = 'NA12883'
ORDER BY
  start,
  END
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Thu Nov 20 18:54:18 2014 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> END </th> <th> reference_bases </th> <th> alternate_bases </th> <th> call_call_set_name </th> <th> gt </th> <th> quality </th> <th> filter </th> <th> likelihood </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198773 </td> <td align="right"> 41198774 </td> <td> C </td> <td> CA </td> <td> NA12883 </td> <td> 0,1 </td> <td align="right"> 317.76 </td> <td> PASS </td> <td> 165,0,334 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198773 </td> <td align="right"> 41198774 </td> <td> C </td> <td>  </td> <td> NA12883 </td> <td> 0,0 </td> <td align="right"> 0.00 </td> <td> PASS </td> <td>  </td> </tr>
   </table>
We see that NA12883 matches the reference for both alleles at position 41198773 and also has an insertion is one allele.  Which is the correct call?

Check Singletons
----------------

```r
result <- DisplayAndDispatchQuery("./sql/private-variants-brca1.sql",
                                  replacements=table_replacement)
```

```
# Private variants within BRCA1
SELECT
  reference_name AS CHROM,
  start AS POS,
  CASE WHEN cnt = 1 THEN 'S'
  WHEN cnt = 2 THEN 'D'
  ELSE STRING(cnt) END AS SINGLETON_DOUBLETON,
  reference_bases AS REF,
  alternate_bases AS ALLELE,
  call.call_set_name AS INDV,
  alt_num,
  genotype,
  cnt
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
          ),
        alternate_bases)
      )
  OMIT
    RECORD IF alternate_bases IS NULL
  HAVING
    num_samples_with_variant = 1
    AND cnt > 0
    )
ORDER BY
  POS, INDV
```
Number of rows returned by this query: 63.

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Thu Nov 20 18:54:23 2014 -->
<table border=1>
<tr> <th> CHROM </th> <th> POS </th> <th> SINGLETON_DOUBLETON </th> <th> REF </th> <th> ALLELE </th> <th> INDV </th> <th> alt_num </th> <th> genotype </th> <th> cnt </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196820 </td> <td> S </td> <td> CT </td> <td> C </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197958 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12884 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198195 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12893 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198220 </td> <td> S </td> <td> T </td> <td> G </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198226 </td> <td> S </td> <td> A </td> <td> C </td> <td> NA12890 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204831 </td> <td> D </td> <td> ACACACACACT </td> <td> A </td> <td> NA12890 </td> <td align="right">   1 </td> <td> 1,1 </td> <td align="right">   2 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204839 </td> <td> D </td> <td> ACT </td> <td> A </td> <td> NA12892 </td> <td align="right">   1 </td> <td> 1,1 </td> <td align="right">   2 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204862 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12882 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41208206 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12880 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210033 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12888 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211352 </td> <td> S </td> <td> T </td> <td> G </td> <td> NA12892 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211485 </td> <td> S </td> <td> CACA </td> <td> CACAACA </td> <td> NA12878 </td> <td align="right">   2 </td> <td> 1,2 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211485 </td> <td> S </td> <td> CACA </td> <td> C </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 1,2 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41217551 </td> <td> S </td> <td> A </td> <td> C </td> <td> NA12889 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218851 </td> <td> S </td> <td> G </td> <td> T </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218884 </td> <td> S </td> <td> A </td> <td> G </td> <td> NA12893 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219261 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12889 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219852 </td> <td> S </td> <td> AT </td> <td> A </td> <td> NA12890 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219872 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12887 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219873 </td> <td> S </td> <td> G </td> <td> T </td> <td> NA12887 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219905 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12880 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41223537 </td> <td> S </td> <td> G </td> <td> GAATGTTCACTGTAACAATGCTTGT </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224871 </td> <td> S </td> <td> A </td> <td> C </td> <td> NA12893 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41227082 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12887 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41227082 </td> <td> S </td> <td> C </td> <td> CGGAA </td> <td> NA12892 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229351 </td> <td> S </td> <td> C </td> <td> T </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229761 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232485 </td> <td> S </td> <td> T </td> <td> G </td> <td> NA12881 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41236026 </td> <td> S </td> <td> A </td> <td> C </td> <td> NA12893 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41238125 </td> <td> S </td> <td> T </td> <td> TACACAC </td> <td> NA12880 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239705 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12891 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250046 </td> <td> S </td> <td> C </td> <td> T </td> <td> NA12888 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250220 </td> <td> S </td> <td> AT </td> <td> A </td> <td> NA12892 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252576 </td> <td> S </td> <td> A </td> <td> G </td> <td> NA12891 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252645 </td> <td> S </td> <td> ATGT </td> <td> A </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252648 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252675 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12890 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252683 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252692 </td> <td> S </td> <td> ATAAT </td> <td> A </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252693 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12877 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12877 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254373 </td> <td> S </td> <td> C </td> <td> T </td> <td> NA12892 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254393 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12884 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254979 </td> <td> S </td> <td> A </td> <td> T </td> <td> NA12887 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256088 </td> <td> S </td> <td> AAAAAAAAAGAAAAG </td> <td> A </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256091 </td> <td> S </td> <td> A </td> <td> G </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256094 </td> <td> S </td> <td> A </td> <td> G </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258182 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12887 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259078 </td> <td> S </td> <td> A </td> <td> ATT </td> <td> NA12892 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260721 </td> <td> S </td> <td> T </td> <td> A </td> <td> NA12891 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260984 </td> <td> S </td> <td> G </td> <td> C </td> <td> NA12892 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263116 </td> <td> S </td> <td> C </td> <td> CA </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263429 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12891 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264110 </td> <td> S </td> <td> C </td> <td> CT </td> <td> NA12891 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270777 </td> <td> S </td> <td> C </td> <td> CT </td> <td> NA12878 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41271292 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12887 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td> S </td> <td> G </td> <td> C </td> <td> NA12889 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273699 </td> <td> S </td> <td> C </td> <td> CAAAA </td> <td> NA12888 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273774 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12892 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274876 </td> <td> S </td> <td> C </td> <td> A </td> <td> NA12883 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275080 </td> <td> S </td> <td> G </td> <td> A </td> <td> NA12887 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276517 </td> <td> S </td> <td> T </td> <td> C </td> <td> NA12880 </td> <td align="right">   1 </td> <td> 0,1 </td> <td align="right">   1 </td> </tr>
   </table>

Check Hardy-Weinberg Equilibrium
-----------------------------------

More calculations needed here . . . this is just allelic frequency.


```r
result <- DisplayAndDispatchQuery("./sql/allelic-frequency-brca1.sql",
                                  replacements=table_replacement)
```

```
# The following query computes the allelic frequency for BRCA1 SNPs.
SELECT
  vars.reference_name,
  vars.start,
  vars.end,
  reference_bases,
  alternate_bases,
  (ref_count + SUM(called_allele_count))/(ref_count + SUM(called_allele_count) + alt_count) AS ref_freq,
  alt_count/(ref_count + SUM(called_allele_count) + alt_count) AS alt_freq,
  ref_count + alt_count + SUM(called_allele_count) AS num_sample_alleles,
  ref_count,
  alt_count,
  SUM(called_allele_count) AS called_allele_count,
FROM (
    # Constrain the left hand side of the _join to reference-matching blocks.
  SELECT
    reference_name,
    start,
    END,
    SUM(0 = call.genotype) WITHIN RECORD AS called_allele_count,
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
    SUM(0 = call.genotype) WITHIN RECORD AS ref_count,
    SUM(1 = call.genotype) WITHIN RECORD AS alt_count,
  FROM
    [genomics-public-data:platinum_genomes.variants]
  WHERE
    reference_name = 'chr17'
    AND start BETWEEN 41196311
    AND 41277499
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
  vars.reference_name,
  vars.start,
  vars.end,
  reference_bases,
  alternate_bases,
  ref_count,
  alt_count
ORDER BY
  vars.reference_name,
  vars.start,
  reference_bases
```
Number of rows returned by this query: 271.

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Thu Nov 20 18:54:27 2014 -->
<table border=1>
<tr> <th> vars_reference_name </th> <th> vars_start </th> <th> vars_end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> ref_freq </th> <th> alt_freq </th> <th> num_sample_alleles </th> <th> ref_count </th> <th> alt_count </th> <th> called_allele_count </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td align="right"> 41196408 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196840 </td> <td align="right"> 41196841 </td> <td> G </td> <td> T </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td align="right"> 41197274 </td> <td> C </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197957 </td> <td align="right"> 41197958 </td> <td> G </td> <td> T </td> <td align="right"> 0.65 </td> <td align="right"> 0.35 </td> <td align="right">  34 </td> <td align="right">  12 </td> <td align="right">  12 </td> <td align="right">  10 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197958 </td> <td align="right"> 41197959 </td> <td> A </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198182 </td> <td align="right"> 41198183 </td> <td> A </td> <td> C </td> <td align="right"> 0.82 </td> <td align="right"> 0.18 </td> <td align="right">  34 </td> <td align="right">   6 </td> <td align="right">   6 </td> <td align="right">  22 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198186 </td> <td align="right"> 41198187 </td> <td> A </td> <td> C </td> <td align="right"> 0.74 </td> <td align="right"> 0.26 </td> <td align="right">  34 </td> <td align="right">   9 </td> <td align="right">   9 </td> <td align="right">  16 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198195 </td> <td align="right"> 41198196 </td> <td> T </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198220 </td> <td align="right"> 41198221 </td> <td> T </td> <td> G </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198226 </td> <td align="right"> 41198227 </td> <td> A </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198231 </td> <td align="right"> 41198232 </td> <td> A </td> <td> C </td> <td align="right"> 0.74 </td> <td align="right"> 0.26 </td> <td align="right">  34 </td> <td align="right">   9 </td> <td align="right">   9 </td> <td align="right">  16 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198264 </td> <td align="right"> 41198265 </td> <td> G </td> <td> T </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198620 </td> <td align="right"> 41198621 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198769 </td> <td align="right"> 41198770 </td> <td> A </td> <td> C </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198773 </td> <td align="right"> 41198774 </td> <td> C </td> <td> A </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41199912 </td> <td align="right"> 41199913 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200108 </td> <td align="right"> 41200109 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200536 </td> <td align="right"> 41200537 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201701 </td> <td align="right"> 41201702 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202634 </td> <td align="right"> 41202635 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202687 </td> <td align="right"> 41202688 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203324 </td> <td align="right"> 41203325 </td> <td> T </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203590 </td> <td align="right"> 41203591 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203671 </td> <td align="right"> 41203672 </td> <td> G </td> <td> T </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204376 </td> <td align="right"> 41204377 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204389 </td> <td align="right"> 41204390 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204796 </td> <td align="right"> 41204797 </td> <td> A </td> <td> C </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204837 </td> <td align="right"> 41204838 </td> <td> A </td> <td> T </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right">  28 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204839 </td> <td align="right"> 41204840 </td> <td> A </td> <td> T </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right">  28 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204855 </td> <td align="right"> 41204856 </td> <td> G </td> <td> A </td> <td align="right"> 0.74 </td> <td align="right"> 0.26 </td> <td align="right">  34 </td> <td align="right">   9 </td> <td align="right">   9 </td> <td align="right">  16 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204862 </td> <td align="right"> 41204863 </td> <td> T </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41205771 </td> <td align="right"> 41205772 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41205940 </td> <td align="right"> 41205941 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41206055 </td> <td align="right"> 41206056 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41208206 </td> <td align="right"> 41208207 </td> <td> T </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41208207 </td> <td align="right"> 41208208 </td> <td> C </td> <td> T </td> <td align="right"> 0.68 </td> <td align="right"> 0.32 </td> <td align="right">  34 </td> <td align="right">  11 </td> <td align="right">  11 </td> <td align="right">  12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41209577 </td> <td align="right"> 41209578 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210033 </td> <td align="right"> 41210034 </td> <td> C </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210034 </td> <td align="right"> 41210035 </td> <td> T </td> <td> A </td> <td align="right"> 0.85 </td> <td align="right"> 0.15 </td> <td align="right">  34 </td> <td align="right">   5 </td> <td align="right">   5 </td> <td align="right">  24 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210035 </td> <td align="right"> 41210036 </td> <td> C </td> <td> A </td> <td align="right"> 0.68 </td> <td align="right"> 0.32 </td> <td align="right">  34 </td> <td align="right">  11 </td> <td align="right">  11 </td> <td align="right">  12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210062 </td> <td align="right"> 41210063 </td> <td> T </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210395 </td> <td align="right"> 41210396 </td> <td> A </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211352 </td> <td align="right"> 41211353 </td> <td> T </td> <td> G </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211652 </td> <td align="right"> 41211653 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212168 </td> <td align="right"> 41212169 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212180 </td> <td align="right"> 41212181 </td> <td> G </td> <td> A </td> <td align="right"> 0.74 </td> <td align="right"> 0.26 </td> <td align="right">  34 </td> <td align="right">   9 </td> <td align="right">   9 </td> <td align="right">  16 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212337 </td> <td align="right"> 41212338 </td> <td> A </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212546 </td> <td align="right"> 41212547 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212804 </td> <td align="right"> 41212805 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213625 </td> <td align="right"> 41213626 </td> <td> G </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213659 </td> <td align="right"> 41213660 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213747 </td> <td align="right"> 41213748 </td> <td> T </td> <td> C </td> <td align="right"> 0.74 </td> <td align="right"> 0.26 </td> <td align="right">  34 </td> <td align="right">   9 </td> <td align="right">   9 </td> <td align="right">  16 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213892 </td> <td align="right"> 41213893 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213995 </td> <td align="right"> 41213996 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41214209 </td> <td align="right"> 41214210 </td> <td> A </td> <td> T </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right">  26 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41214210 </td> <td align="right"> 41214211 </td> <td> A </td> <td> C </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right">  26 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41215824 </td> <td align="right"> 41215825 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216204 </td> <td align="right"> 41216205 </td> <td> G </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216865 </td> <td align="right"> 41216866 </td> <td> G </td> <td> T </td> <td align="right"> 0.74 </td> <td align="right"> 0.26 </td> <td align="right">  34 </td> <td align="right">   9 </td> <td align="right">   9 </td> <td align="right">  16 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216932 </td> <td align="right"> 41216933 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41217551 </td> <td align="right"> 41217552 </td> <td> A </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41217873 </td> <td align="right"> 41217874 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218332 </td> <td align="right"> 41218333 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218571 </td> <td align="right"> 41218572 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218814 </td> <td align="right"> 41218815 </td> <td> A </td> <td> C </td> <td align="right"> 0.82 </td> <td align="right"> 0.18 </td> <td align="right">  34 </td> <td align="right">   6 </td> <td align="right">   6 </td> <td align="right">  22 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218821 </td> <td align="right"> 41218822 </td> <td> T </td> <td> C </td> <td align="right"> 0.68 </td> <td align="right"> 0.32 </td> <td align="right">  34 </td> <td align="right">  11 </td> <td align="right">  11 </td> <td align="right">  12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218851 </td> <td align="right"> 41218852 </td> <td> G </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218855 </td> <td align="right"> 41218856 </td> <td> C </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218858 </td> <td align="right"> 41218859 </td> <td> G </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218884 </td> <td align="right"> 41218885 </td> <td> A </td> <td> G </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219261 </td> <td align="right"> 41219262 </td> <td> A </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219290 </td> <td align="right"> 41219291 </td> <td> G </td> <td> T </td> <td align="right"> 0.76 </td> <td align="right"> 0.24 </td> <td align="right">  34 </td> <td align="right">   8 </td> <td align="right">   8 </td> <td align="right">  18 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219340 </td> <td align="right"> 41219341 </td> <td> G </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219559 </td> <td align="right"> 41219560 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219779 </td> <td align="right"> 41219780 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219803 </td> <td align="right"> 41219804 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219872 </td> <td align="right"> 41219873 </td> <td> A </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219873 </td> <td align="right"> 41219874 </td> <td> G </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219892 </td> <td align="right"> 41219893 </td> <td> A </td> <td> C </td> <td align="right"> 0.65 </td> <td align="right"> 0.35 </td> <td align="right">  34 </td> <td align="right">  12 </td> <td align="right">  12 </td> <td align="right">  10 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219900 </td> <td align="right"> 41219901 </td> <td> T </td> <td> G </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219903 </td> <td align="right"> 41219904 </td> <td> A </td> <td> G </td> <td align="right"> 0.56 </td> <td align="right"> 0.44 </td> <td align="right">  34 </td> <td align="right">  15 </td> <td align="right">  15 </td> <td align="right">   4 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219905 </td> <td align="right"> 41219906 </td> <td> G </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219929 </td> <td align="right"> 41219930 </td> <td> T </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220222 </td> <td align="right"> 41220223 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222461 </td> <td align="right"> 41222462 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222722 </td> <td align="right"> 41222723 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41223093 </td> <td align="right"> 41223094 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224832 </td> <td align="right"> 41224833 </td> <td> G </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224867 </td> <td align="right"> 41224868 </td> <td> A </td> <td> C </td> <td align="right"> 0.76 </td> <td align="right"> 0.24 </td> <td align="right">  34 </td> <td align="right">   8 </td> <td align="right">   8 </td> <td align="right">  18 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224871 </td> <td align="right"> 41224872 </td> <td> A </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224891 </td> <td align="right"> 41224892 </td> <td> G </td> <td> T </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224907 </td> <td align="right"> 41224908 </td> <td> A </td> <td> T </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225578 </td> <td align="right"> 41225579 </td> <td> C </td> <td> A </td> <td align="right"> 0.68 </td> <td align="right"> 0.32 </td> <td align="right">  34 </td> <td align="right">  11 </td> <td align="right">  11 </td> <td align="right">  12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225780 </td> <td align="right"> 41225781 </td> <td> A </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225782 </td> <td align="right"> 41225783 </td> <td> A </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225838 </td> <td align="right"> 41225839 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226600 </td> <td align="right"> 41226601 </td> <td> G </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226674 </td> <td align="right"> 41226675 </td> <td> A </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226740 </td> <td align="right"> 41226741 </td> <td> T </td> <td> G </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right">  26 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41227082 </td> <td align="right"> 41227083 </td> <td> C </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41227103 </td> <td align="right"> 41227104 </td> <td> A </td> <td> G </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41227105 </td> <td align="right"> 41227106 </td> <td> G </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229351 </td> <td align="right"> 41229352 </td> <td> C </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229365 </td> <td align="right"> 41229366 </td> <td> G </td> <td> T </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229385 </td> <td align="right"> 41229386 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229761 </td> <td align="right"> 41229762 </td> <td> T </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229762 </td> <td align="right"> 41229763 </td> <td> T </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229772 </td> <td align="right"> 41229773 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229811 </td> <td align="right"> 41229812 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229856 </td> <td align="right"> 41229857 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229907 </td> <td align="right"> 41229908 </td> <td> T </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230227 </td> <td align="right"> 41230228 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230335 </td> <td align="right"> 41230336 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230375 </td> <td align="right"> 41230376 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230523 </td> <td align="right"> 41230524 </td> <td> T </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230536 </td> <td align="right"> 41230537 </td> <td> A </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230989 </td> <td align="right"> 41230990 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231220 </td> <td align="right"> 41231221 </td> <td> A </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231515 </td> <td align="right"> 41231516 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231901 </td> <td align="right"> 41231902 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232485 </td> <td align="right"> 41232486 </td> <td> T </td> <td> G </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232566 </td> <td align="right"> 41232567 </td> <td> A </td> <td> C </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232697 </td> <td align="right"> 41232698 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41234469 </td> <td align="right"> 41234470 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41235798 </td> <td align="right"> 41235799 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41236026 </td> <td align="right"> 41236027 </td> <td> A </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41237952 </td> <td align="right"> 41237953 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41238125 </td> <td align="right"> 41238126 </td> <td> T </td> <td> C </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239471 </td> <td align="right"> 41239472 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239490 </td> <td align="right"> 41239491 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239627 </td> <td align="right"> 41239628 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239680 </td> <td align="right"> 41239681 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239698 </td> <td align="right"> 41239699 </td> <td> G </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239705 </td> <td align="right"> 41239706 </td> <td> G </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239915 </td> <td align="right"> 41239916 </td> <td> T </td> <td> A </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right">  26 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41240276 </td> <td align="right"> 41240277 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241389 </td> <td align="right"> 41241390 </td> <td> C </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241502 </td> <td align="right"> 41241503 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241568 </td> <td align="right"> 41241569 </td> <td> T </td> <td> C </td> <td align="right"> 0.85 </td> <td align="right"> 0.15 </td> <td align="right">  34 </td> <td align="right">   5 </td> <td align="right">   5 </td> <td align="right">  24 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241581 </td> <td align="right"> 41241582 </td> <td> G </td> <td> T </td> <td align="right"> 0.85 </td> <td align="right"> 0.15 </td> <td align="right">  34 </td> <td align="right">   5 </td> <td align="right">   5 </td> <td align="right">  24 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242077 </td> <td align="right"> 41242078 </td> <td> G </td> <td> A </td> <td align="right"> 0.96 </td> <td align="right"> 0.04 </td> <td align="right">  25 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  23 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242284 </td> <td align="right"> 41242285 </td> <td> T </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41243189 </td> <td align="right"> 41243190 </td> <td> T </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41243999 </td> <td align="right"> 41244000 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41244434 </td> <td align="right"> 41244435 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41244935 </td> <td align="right"> 41244936 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245236 </td> <td align="right"> 41245237 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245465 </td> <td align="right"> 41245466 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247603 </td> <td align="right"> 41247604 </td> <td> A </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248163 </td> <td align="right"> 41248164 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248392 </td> <td align="right"> 41248393 </td> <td> C </td> <td> A </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248411 </td> <td align="right"> 41248412 </td> <td> G </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248413 </td> <td align="right"> 41248414 </td> <td> G </td> <td> A </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248415 </td> <td align="right"> 41248416 </td> <td> G </td> <td> A </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248419 </td> <td align="right"> 41248420 </td> <td> G </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248483 </td> <td align="right"> 41248484 </td> <td> G </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41249093 </td> <td align="right"> 41249094 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250046 </td> <td align="right"> 41250047 </td> <td> C </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250063 </td> <td align="right"> 41250064 </td> <td> G </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250065 </td> <td align="right"> 41250066 </td> <td> A </td> <td> G </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250259 </td> <td align="right"> 41250260 </td> <td> C </td> <td> G </td> <td align="right"> 0.74 </td> <td align="right"> 0.26 </td> <td align="right">  34 </td> <td align="right">   9 </td> <td align="right">   9 </td> <td align="right">  16 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250851 </td> <td align="right"> 41250852 </td> <td> G </td> <td> T </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250922 </td> <td align="right"> 41250923 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250986 </td> <td align="right"> 41250987 </td> <td> A </td> <td> C </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250989 </td> <td align="right"> 41250990 </td> <td> A </td> <td> C </td> <td align="right"> 0.71 </td> <td align="right"> 0.29 </td> <td align="right">  34 </td> <td align="right">  10 </td> <td align="right">  10 </td> <td align="right">  14 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251494 </td> <td align="right"> 41251495 </td> <td> C </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251645 </td> <td align="right"> 41251646 </td> <td> T </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251930 </td> <td align="right"> 41251931 </td> <td> G </td> <td> A </td> <td align="right"> 0.71 </td> <td align="right"> 0.29 </td> <td align="right">  34 </td> <td align="right">  10 </td> <td align="right">  10 </td> <td align="right">  14 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252574 </td> <td align="right"> 41252575 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252576 </td> <td align="right"> 41252577 </td> <td> A </td> <td> G </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252611 </td> <td align="right"> 41252612 </td> <td> T </td> <td> A </td> <td align="right"> 0.82 </td> <td align="right"> 0.18 </td> <td align="right">  34 </td> <td align="right">   6 </td> <td align="right">   6 </td> <td align="right">  22 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252612 </td> <td align="right"> 41252613 </td> <td> T </td> <td> C </td> <td align="right"> 0.82 </td> <td align="right"> 0.18 </td> <td align="right">  34 </td> <td align="right">   6 </td> <td align="right">   6 </td> <td align="right">  22 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252648 </td> <td align="right"> 41252649 </td> <td> T </td> <td> A </td> <td align="right"> 0.96 </td> <td align="right"> 0.04 </td> <td align="right">  27 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252649 </td> <td align="right"> 41252650 </td> <td> T </td> <td> C </td> <td align="right"> 0.93 </td> <td align="right"> 0.07 </td> <td align="right">  28 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  24 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252657 </td> <td align="right"> 41252658 </td> <td> T </td> <td> C </td> <td align="right"> 0.89 </td> <td align="right"> 0.11 </td> <td align="right">  28 </td> <td align="right">   1 </td> <td align="right">   3 </td> <td align="right">  24 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252675 </td> <td align="right"> 41252676 </td> <td> T </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  32 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252683 </td> <td align="right"> 41252684 </td> <td> T </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  32 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252693 </td> <td align="right"> 41252694 </td> <td> T </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  33 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  31 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694 </td> <td align="right"> 41252695 </td> <td> A </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  33 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  31 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252695 </td> <td align="right"> 41252696 </td> <td> A </td> <td> T </td> <td align="right"> 0.50 </td> <td align="right"> 0.50 </td> <td align="right">  28 </td> <td align="right">  10 </td> <td align="right">  14 </td> <td align="right">   4 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696 </td> <td align="right"> 41252697 </td> <td> T </td> <td> A </td> <td align="right"> 0.42 </td> <td align="right"> 0.58 </td> <td align="right">  26 </td> <td align="right">   9 </td> <td align="right">  15 </td> <td align="right">   2 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696 </td> <td align="right"> 41252697 </td> <td> T </td> <td> C </td> <td align="right"> 0.75 </td> <td align="right"> 0.25 </td> <td align="right">   4 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">   2 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254173 </td> <td align="right"> 41254174 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254373 </td> <td align="right"> 41254374 </td> <td> C </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254392 </td> <td align="right"> 41254393 </td> <td> G </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254393 </td> <td align="right"> 41254394 </td> <td> A </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254404 </td> <td align="right"> 41254405 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254485 </td> <td align="right"> 41254486 </td> <td> T </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254979 </td> <td align="right"> 41254980 </td> <td> A </td> <td> T </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255101 </td> <td align="right"> 41255102 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255110 </td> <td align="right"> 41255111 </td> <td> A </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255482 </td> <td align="right"> 41255483 </td> <td> C </td> <td> A </td> <td align="right"> 0.68 </td> <td align="right"> 0.32 </td> <td align="right">  34 </td> <td align="right">  11 </td> <td align="right">  11 </td> <td align="right">  12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256086 </td> <td align="right"> 41256087 </td> <td> G </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256091 </td> <td align="right"> 41256092 </td> <td> A </td> <td> G </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  33 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  31 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256094 </td> <td align="right"> 41256095 </td> <td> A </td> <td> G </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  33 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  31 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256097 </td> <td align="right"> 41256098 </td> <td> G </td> <td> A </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  32 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256100 </td> <td align="right"> 41256101 </td> <td> A </td> <td> G </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  32 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256102 </td> <td align="right"> 41256103 </td> <td> G </td> <td> A </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  32 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  24 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41257133 </td> <td align="right"> 41257134 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41257457 </td> <td align="right"> 41257458 </td> <td> A </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258042 </td> <td align="right"> 41258043 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258152 </td> <td align="right"> 41258153 </td> <td> G </td> <td> A </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258182 </td> <td align="right"> 41258183 </td> <td> C </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259048 </td> <td align="right"> 41259049 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259064 </td> <td align="right"> 41259065 </td> <td> T </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259099 </td> <td align="right"> 41259100 </td> <td> G </td> <td> T </td> <td align="right"> 0.68 </td> <td align="right"> 0.32 </td> <td align="right">  34 </td> <td align="right">  11 </td> <td align="right">  11 </td> <td align="right">  12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259100 </td> <td align="right"> 41259101 </td> <td> A </td> <td> T </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259112 </td> <td align="right"> 41259113 </td> <td> G </td> <td> A </td> <td align="right"> 0.82 </td> <td align="right"> 0.18 </td> <td align="right">  34 </td> <td align="right">   6 </td> <td align="right">   6 </td> <td align="right">  22 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259417 </td> <td align="right"> 41259418 </td> <td> G </td> <td> T </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260351 </td> <td align="right"> 41260352 </td> <td> C </td> <td> A </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260693 </td> <td align="right"> 41260694 </td> <td> A </td> <td> C </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260721 </td> <td align="right"> 41260722 </td> <td> T </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260722 </td> <td align="right"> 41260723 </td> <td> C </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260761 </td> <td align="right"> 41260762 </td> <td> T </td> <td> G </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260807 </td> <td align="right"> 41260808 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260984 </td> <td align="right"> 41260985 </td> <td> G </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261053 </td> <td align="right"> 41261054 </td> <td> T </td> <td> C </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261054 </td> <td align="right"> 41261055 </td> <td> A </td> <td> C </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261056 </td> <td align="right"> 41261057 </td> <td> C </td> <td> A </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261057 </td> <td align="right"> 41261058 </td> <td> T </td> <td> C </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261084 </td> <td align="right"> 41261085 </td> <td> G </td> <td> A </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261085 </td> <td align="right"> 41261086 </td> <td> T </td> <td> C </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261232 </td> <td align="right"> 41261233 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41262105 </td> <td align="right"> 41262106 </td> <td> T </td> <td> G </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41262357 </td> <td align="right"> 41262358 </td> <td> C </td> <td> T </td> <td align="right"> 0.82 </td> <td align="right"> 0.18 </td> <td align="right">  34 </td> <td align="right">   6 </td> <td align="right">   6 </td> <td align="right">  22 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263043 </td> <td align="right"> 41263044 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263116 </td> <td align="right"> 41263117 </td> <td> C </td> <td> A </td> <td align="right"> 0.82 </td> <td align="right"> 0.18 </td> <td align="right">  34 </td> <td align="right">   6 </td> <td align="right">   6 </td> <td align="right">  22 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263429 </td> <td align="right"> 41263430 </td> <td> G </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263565 </td> <td align="right"> 41263566 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264145 </td> <td align="right"> 41264146 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264363 </td> <td align="right"> 41264364 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264738 </td> <td align="right"> 41264739 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264748 </td> <td align="right"> 41264749 </td> <td> C </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264749 </td> <td align="right"> 41264750 </td> <td> A </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264752 </td> <td align="right"> 41264753 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41265775 </td> <td align="right"> 41265776 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41266423 </td> <td align="right"> 41266424 </td> <td> G </td> <td> T </td> <td align="right"> 0.85 </td> <td align="right"> 0.15 </td> <td align="right">  34 </td> <td align="right">   5 </td> <td align="right">   5 </td> <td align="right">  24 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41267049 </td> <td align="right"> 41267050 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41268205 </td> <td align="right"> 41268206 </td> <td> A </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270228 </td> <td align="right"> 41270229 </td> <td> T </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270276 </td> <td align="right"> 41270277 </td> <td> C </td> <td> T </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270462 </td> <td align="right"> 41270463 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270665 </td> <td align="right"> 41270666 </td> <td> C </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270793 </td> <td align="right"> 41270794 </td> <td> G </td> <td> T </td> <td align="right"> 0.88 </td> <td align="right"> 0.12 </td> <td align="right">  34 </td> <td align="right">   4 </td> <td align="right">   4 </td> <td align="right">  26 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41271292 </td> <td align="right"> 41271293 </td> <td> G </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41271293 </td> <td align="right"> 41271294 </td> <td> A </td> <td> G </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right">  27 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  27 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td align="right"> 41273095 </td> <td> G </td> <td> C </td> <td align="right"> 0.95 </td> <td align="right"> 0.05 </td> <td align="right">  22 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td align="right"> 41273095 </td> <td> G </td> <td> A </td> <td align="right"> 0.81 </td> <td align="right"> 0.19 </td> <td align="right">  32 </td> <td align="right">   6 </td> <td align="right">   6 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273347 </td> <td align="right"> 41273348 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273378 </td> <td align="right"> 41273379 </td> <td> G </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273536 </td> <td align="right"> 41273537 </td> <td> A </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273697 </td> <td align="right"> 41273698 </td> <td> C </td> <td> A </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273698 </td> <td align="right"> 41273699 </td> <td> T </td> <td> A </td> <td align="right"> 0.91 </td> <td align="right"> 0.09 </td> <td align="right">  34 </td> <td align="right">   3 </td> <td align="right">   3 </td> <td align="right">  28 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273699 </td> <td align="right"> 41273700 </td> <td> C </td> <td> A </td> <td align="right"> 0.62 </td> <td align="right"> 0.38 </td> <td align="right">  34 </td> <td align="right">  13 </td> <td align="right">  13 </td> <td align="right">   8 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273725 </td> <td align="right"> 41273726 </td> <td> G </td> <td> A </td> <td align="right"> 0.62 </td> <td align="right"> 0.38 </td> <td align="right">  34 </td> <td align="right">  13 </td> <td align="right">  13 </td> <td align="right">   8 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273739 </td> <td align="right"> 41273740 </td> <td> C </td> <td> A </td> <td align="right"> 0.85 </td> <td align="right"> 0.15 </td> <td align="right">  34 </td> <td align="right">   5 </td> <td align="right">   5 </td> <td align="right">  24 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273746 </td> <td align="right"> 41273747 </td> <td> C </td> <td> A </td> <td align="right"> 0.68 </td> <td align="right"> 0.32 </td> <td align="right">  34 </td> <td align="right">  11 </td> <td align="right">  11 </td> <td align="right">  12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273748 </td> <td align="right"> 41273749 </td> <td> G </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273774 </td> <td align="right"> 41273775 </td> <td> C </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273793 </td> <td align="right"> 41273794 </td> <td> T </td> <td> A </td> <td align="right"> 0.94 </td> <td align="right"> 0.06 </td> <td align="right">  34 </td> <td align="right">   2 </td> <td align="right">   2 </td> <td align="right">  30 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274777 </td> <td align="right"> 41274778 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274876 </td> <td align="right"> 41274877 </td> <td> C </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274905 </td> <td align="right"> 41274906 </td> <td> G </td> <td> A </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275080 </td> <td align="right"> 41275081 </td> <td> G </td> <td> A </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275150 </td> <td align="right"> 41275151 </td> <td> G </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275391 </td> <td align="right"> 41275392 </td> <td> G </td> <td> T </td> <td align="right"> 0.62 </td> <td align="right"> 0.38 </td> <td align="right">  34 </td> <td align="right">  13 </td> <td align="right">  13 </td> <td align="right">   8 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275644 </td> <td align="right"> 41275645 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276246 </td> <td align="right"> 41276247 </td> <td> A </td> <td> G </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276347 </td> <td align="right"> 41276348 </td> <td> T </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276517 </td> <td align="right"> 41276518 </td> <td> T </td> <td> C </td> <td align="right"> 0.97 </td> <td align="right"> 0.03 </td> <td align="right">  34 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">  32 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41277186 </td> <td align="right"> 41277187 </td> <td> G </td> <td> C </td> <td align="right"> 0.79 </td> <td align="right"> 0.21 </td> <td align="right">  34 </td> <td align="right">   7 </td> <td align="right">   7 </td> <td align="right">  20 </td> </tr>
   </table>

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

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Thu Nov 20 18:54:29 2014 -->
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


