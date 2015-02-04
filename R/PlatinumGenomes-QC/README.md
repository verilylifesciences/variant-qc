# Quality Control using Google Genomics

In this codelab you will use [Google Genomics](https://cloud.google.com/genomics/), [Google BigQuery](https://cloud.google.com/bigquery/what-is-bigquery), [Google Cloud Dataflow](https://cloud.google.com/dataflow/), [Apache Spark](http://spark.apache.org/), and [R](http://www.r-project.org/) to perform quality control checks on the [Illumina Platinum Genomes dataset](https://cloud.google.com/genomics/data/platinum-genomes).  Specifically, you will:

1. examine the format, domain and range of the Illumina Platinum Genomes data
2. convert data from Genome VCF format to VCF format via your choice of
 + Java: a Google Cloud Dataflow job
 + Python: an Apache PySpark job or a Hadoop Streaming job
3. perform sample-level Quality Control checks such as
 + homozygosity rate via BigQuery
 + ethnicity using Principal Coordinate Analysis (either from scratch via Apache Spark or using pre-computed results)
 + relatedness among individuals using Identity-By-State results (either from scratch via Google Cloud Dataflow or using pre-computed results)
4. perform variant-level Quality Control checks via BigQuery such as
 + Hardy-Weinberg Equilibrium via BigQuery

## Instructions
1. Make sure you have completed the necessary [prerequisites](../README.md).
2. Read the rendered version of the RMarkdown for a particular analysis.
 + [Data Overview](./Data-Overview.md)
 + [Data Conversion](./Data-Conversion.md)
 + [Sample-Level QC](./Sample-Level-QC.md)
 + [Variant-Level QC](./Variant-Level-QC.md)
3. Execute the code *chunk by chunk* in RStudio or *line by line* in R.
4. Optionally, [read more](./comparison/QC-Comparison.md) about how these results compare to those from tools such as [vcftools](http://vcftools.sourceforge.net/).

_This is very much a work-in-progress._  Please submit pull requests and file issues for feedback and corrections!
