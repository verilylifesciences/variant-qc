# Quality Control using Google Genomics

In this codelab you will use [Google Genomics](https://cloud.google.com/genomics/), [Google BigQuery](https://cloud.google.com/bigquery/what-is-bigquery), [Apache Beam](https://beam.apache.org/), [Apache Spark](http://spark.apache.org/), and [R](http://www.r-project.org/) to perform quality control checks on the [Illumina Platinum Genomes dataset](https://cloud.google.com/genomics/data/platinum-genomes).  Specifically, you will:

1. Examine the format, domain and range of the Illumina Platinum Genomes data.
2. Convert data with non-variant segment records to variant-only data which includes reference-matching calls from overlapping non-variant segment records by running an Apache Beam pipeline on [Google Cloud Dataflow](https://cloud.google.com/dataflow/).
3. Perform sample-level Quality Control checks such as:
    + homozygosity rate via BigQuery
    + ethnicity using Principal Coordinate Analysis (either by running this Apache Spark job from scratch or reading in pre-computed results)
    + relatedness among individuals using Identity-By-State results (either by running the Apache Beam pipeline from scratch or reading in pre-computed results)
4. Perform variant-level Quality Control checks via BigQuery such as:
    + Hardy-Weinberg Equilibrium via BigQuery

## Instructions
You can read or execute the analysis:

* **Read** the rendered version of the RMarkdown for a particular analysis.
    + [Part 1: Data Overview](./Data-Overview.md)
    + [Part 2: Data Transformation](./Data-Transformation.md)
    + [Part 3: Sample-Level QC](./Sample-Level-QC.md)
    + [Part 4: Variant-Level QC](./Variant-Level-QC.md)
* **Execute** the code *chunk by chunk* in RStudio or *line by line* in R.
    1. Clone this github repository.
    2. Make sure you have completed the necessary [prerequisites](../README.md#prerequisites).
    3. Run the RMarkdown file for a particular analysis.
      + [Part 1: Data Overview](./Data-Overview.Rmd)
      + [Part 2: Data Transformation](./Data-Transformation.md)
      + [Part 3: Sample-Level QC](./Sample-Level-QC.Rmd)
      + [Part 4: Variant-Level QC](./Variant-Level-QC.Rmd)
    4. If you prefer, `Data-Overview.R` can be created from [Data-Overview.Rmd](./Data-Overview.Rmd) via
```
library(knitr)
purl("./Data-Overview.Rmd", documentation=2)
```

* **Optionally**, [read more](./comparison/QC-Comparison.md) about how these results compare to those from tools such as [vcftools](http://vcftools.sourceforge.net/).

For more detail, see [Cloud-based interactive analytics for terabytes of genomic variants data](https://academic.oup.com/bioinformatics/article-abstract/doi/10.1093/bioinformatics/btx468/4036385/Cloud-based-Interactive-Analytics-for-Terabytes-of?redirectedFrom=fulltext).
