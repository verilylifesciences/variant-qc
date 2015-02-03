# Data Conversion

Data in [genome VCF](https://sites.google.com/site/gvcftools/home/about-gvcf/gvcf-conventions) (gVCF) format can be challenging to query.  For more detail see this [comparison](https://github.com/googlegenomics/bigquery-examples/tree/master/pgp/data-stories/schema-comparisons#motivation).

In the analyses within this codelab, sometimes we work with the original gVCF data and sometimes we work with data converted from gVCF to VCF to make querying easier.  

More coming soon!  Its a combination of:
* https://github.com/googlegenomics/dataflow-java/blob/master/src/main/java/com/google/cloud/genomics/dataflow/pipelines/ConvertGvcfToVcf.java which will be moved into codelabs/Java
* https://github.com/deflaux/codelabs/tree/qc-codelab/Python/PlatinumGenomes-gVCF-to-VCF with an update for PySpark in addition to the current support for Hadoop Streaming

