# Part 2: Data Conversion

Data from source files in [genome VCF](https://sites.google.com/site/gvcftools/home/about-gvcf/gvcf-conventions) (gVCF) format or in Complete Genomics format can be challenging to query due to the presence of non-variant segment records.  For more detail see this [comparison](https://github.com/googlegenomics/bigquery-examples/tree/master/pgp/data-stories/schema-comparisons#motivation).

In this portion of the codelab we will use a cluster computing job to transform data with non-variant segments to variant-only data with calls from non-variant-segments merged into the variants with which they overlap. This is currently done only for SNP variants. Indels and structural variants are left as-is.

## Running the Cluster Compute Job

We will use either data from variant set [3049512673186936334](https://cloud.google.com/genomics/data/platinum-genomes) or BigQuery table genomics-public-data:platinum_genomes.variants as the data source for this job.

More coming soon!  Its a combination of:
* **Java** https://github.com/googlegenomics/dataflow-java/blob/master/src/main/java/com/google/cloud/genomics/dataflow/pipelines/ConvertNonVariantSegmentsToVariants.java which will be moved into codelabs/Java
* **Python** https://github.com/deflaux/codelabs/tree/qc-codelab/Python/PlatinumGenomes-variant-transformation with an update for PySpark in addition to the current support for Hadoop Streaming

## Results

We have now created table [google.com:biggene:platinum_genomes.expanded_variants](https://bigquery.cloud.google.com/table/google.com:biggene:platinum_genomes.expanded_variants?pli=1)

In the analyses that follow in this codelab, sometimes we work with the original data and sometimes we work with the expanded variant-only data.

## Optional: modify the Cluster Compute Job

More coming soon!  

Some ideas:
* Note that for a large cohort with many, many more rare variants we may wish to instead modify the logic here to summarize the number of calls that match the reference for each variant instead of adding those individual calls to the record.
* merge calls in records with 1/2 genotypes with the same variant found in records with 0/1 genotypes

Reference Name | Start     | End       | Reference Bases | Alternate Bases
---------------|-----------|-----------|-----------------|-----------------
chr6           | 120458771 | 120458773 |TA               |TAA			
chr6           | 120458771 | 120458773 |TA               |TAA,T
 
--------------------------------------------------------
_Next_: [Part 3: Sample-Level QC](./Sample-Level-QC.md)

