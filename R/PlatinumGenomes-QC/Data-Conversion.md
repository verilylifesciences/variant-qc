# Part 2: Data Conversion

Data from source files in [genome VCF](https://sites.google.com/site/gvcftools/home/about-gvcf/gvcf-conventions) (gVCF) format or in Complete Genomics format can be challenging to query due to the presense of non-variant segment records.  For more detail see this [comparison](https://github.com/googlegenomics/bigquery-examples/tree/master/pgp/data-stories/schema-comparisons#motivation).

In the analyses within this codelab, sometimes we work with the original data and sometimes we work with variant-only data which was converted to include reference-matching calls from overlapping non-variant segment records to make querying easier.  

More coming soon!  Its a combination of:
* https://github.com/googlegenomics/dataflow-java/blob/master/src/main/java/com/google/cloud/genomics/dataflow/pipelines/ConvertNonVariantSegmentsToVariants.java which will be moved into codelabs/Java
* https://github.com/deflaux/codelabs/tree/qc-codelab/Python/PlatinumGenomes-gVCF-to-VCF with an update for PySpark in addition to the current support for Hadoop Streaming

--------------------------------------------------------
_Next_: [Part 3: Sample-Level QC](./Sample-Level-QC.md)
