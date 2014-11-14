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

Performing QC using Google Genomics
===================================

The following example makes use of [Illumina Platinum Genomes](http://www.illumina.com/platinumgenomes/).  For more detail about how this data was loaded into the Google Genomics API, please see [Google Genomics Public Data](https://cloud.google.com/genomics/data/platinum-genomes).

DISCLAIMER: this is a work-in-progress, don't expect this to make sense just yet :-)

For provenance details of the results being reproduced in this codelab, see <coming soon>.



Setting Up
-----------

```r
# Setup for BigQuery access
require(bigrquery)
require(RCurl)
project <- "genomics-public-data"                   # put your projectID here
DisplayAndDispatchQuery <- function(queryUri, replacements=list()) {
  if(grepl("^https.*", queryUri)) {
    querySql <- getURL(queryUri, ssl.verifypeer=FALSE)
  } else {
    querySql <- readChar(queryUri, nchars=1e6)
  }
  for(replacement in names(replacements)) {
    querySql <- sub(replacement, replacements[[replacement]], querySql, fixed=TRUE)
  }
  cat(querySql)
  query_exec(querySql, project)
}
table_replacement <- list("_THE_TABLE_"="genomics-public-data:platinum_genomes.variants")
```

Comparing Source Data
--------------------- 
Number of rows in our multi-sample VCF [platinum_genomes_brca1_merged.vcf.gz](./data/platinum_genomes_brca1_merged.vcf.gz)
```
> zcat platinum_genomes_brca1_merged.vcf.gz | grep -v '#' | cut -f 5 | grep -v '\.' | wc -l
332
```

Compare to [BRCA1 variants in BigQuery](https://github.com/googlegenomics/getting-started-bigquery/blob/master/RMarkdown/literate-programming-demo.md#data-visualization)

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

The differences are due to:
* how the alternates for the same genomic position are nested
* any gVCF reference-matching blocks overlapping the start of BRCA1 are included by the region filter for bcftools merge
* Google Genomics uses a 0-based coordinate system but VCF files are 1-based

Singletons
----------

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
<!-- Thu Nov 13 16:22:43 2014 -->
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

Compare to [brca1.singletons](./data/singletons/brca1.singletons) which has 12,384 entries.  The majority of those are for 0/0 genotypes from reference matching blocks.  Should those actually be considered singletons?
