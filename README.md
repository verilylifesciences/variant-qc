### Disclaimer

This is not an official Verily product.

variant-qc
==========

This repository contains code to perform cohort-level quality control checks on human genomic variants. Cloud technology is used to perform queries in parallel. For prior work, see [Cloud-based interactive analytics for terabytes of genomic variants data](https://academic.oup.com/bioinformatics/article/33/23/3709/4036385).

View output from these queries run on public data
--------------------------------------------------

Before running the queries yourself, you can see the results on a few public datasets:

* [QC overview reports](https://console.cloud.google.com/storage/genomics-public-data/platinum-genomes/reports) on [DeepVariant Platinum Genomes](https://cloud.google.com/genomics/docs/public-datasets/illumina-platinum-genomes) data
* [QC overview reports](https://console.cloud.google.com/storage/genomics-public-data/simons-genome-diversity-project/reports) on [Simons Genome Diversity Project](https://cloud.google.com/genomics/docs/public-datasets/simons) data
* [QC overview reports](https://console.cloud.google.com/storage/genomics-public-data/1000-genomes-phase-3/reports) on [1000 Genomes](https://cloud.google.com/genomics/docs/public-datasets/1000-genomes) data
* [example ad hoc explorations of QC results](./notebooks)

Run these queries on your own data
----------------------------------

#### Load data to BigQuery

The queries in this repository assume that the VCFs were loaded to BigQuery using [Variant Transforms](https://cloud.google.com/genomics/docs/how-tos/load-variants) with the [MOVE_TO_CALLS merge strategy](https://github.com/googlegenomics/gcp-variant-transforms/blob/master/docs/variant_merging.md#move_to_calls-strategy) included.

Using the MOVE_TO_CALLS merge strategy will produce a core set of columns common to all tables created from VCFs and calls for the exact same (`reference_name`, `start_position`, `end_position`, `reference_bases`, and all `alternate_bases`) grouped together in a single row.

We recommend loading single-sample VCFs into a "genome call table" and also the multisample VCF into a "multisample-variants table".

If you do not have a multisample VCF, you could:

* use https://github.com/gatk-workflows/gatk4-germline-snps-indels#joint-discovery-gatk- to create one
* use https://github.com/verilylifesciences/joint-genotype to create one
* or skip the queries that require knowing how many samples match the reference such as Hardy-Weinberg Equilibrium

#### Predict ancestry

If your sample information does not already include ancestry, you can predict the ancestry for each genome using [Genomic ancestry inference with deep learning](https://cloud.google.com/blog/big-data/2017/09/genomic-ancestry-inference-with-deep-learning).

#### Run the QC overview reports

Run the [RMarkdown parameterized reports](https://bookdown.org/yihui/rmarkdown/parameterized-reports.html) to get an overview of your data.

#### Drill down on results

Drill down further on results by creating additional plots and/or performing additional queries. For example, these queries can be used from the context of Jupyter [notebooks](./notebooks), and then additional queries or other queries can be used to further explain the results for a particular dataset.

Technologies used
-----------------

The methods make use of:

  * Standard SQL via [BigQuery](https://cloud.google.com/bigquery/docs/)
  * [Apache Beam](https://beam.apache.org/) via [Dataflow](https://cloud.google.com/dataflow/docs/)
  * [TensorFlow](https://www.tensorflow.org/) via [Cloud Machine Learning Engine](https://cloud.google.com/ml-engine/docs/)

Each technology has introductory material that may help you when working with the code in this repository.
