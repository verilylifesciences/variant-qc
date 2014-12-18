Using Hadoop Streaming to convert gVCF to VCF
==============================================

The following example makes use of [Illumina Platinum Genomes](http://www.illumina.com/platinumgenomes/).  For more detail about how this data was loaded into the Google Genomics API, please see [Google Genomics Public Data](https://cloud.google.com/genomics/data/platinum-genomes).

This job creates a new BigQuery table containing data that is reshaped, but the schema is exactly the same.  It adds the reference-matching calls to the relevant SNP variant records --> essentially adding redundant data in an effort to enable easier querying. See [gvcf_expander.py](./gvcf_expander.py) for more detail including unit tests.

[Hadoop on Google Cloud Platform](https://cloud.google.com/solutions/hadoop/click-to-deploy) supports using BigQuery as a source and/or a sink for jobs.  For more detail, see the [BigQuery Connector](https://cloud.google.com/hadoop/bigquery-connector) and specific details for [Hadoop Streaming](https://groups.google.com/forum/#!topic/gcp-hadoop-announce/bzji9yjj304)

First copy over the code and schema to the Hadoop master:
```
gcutil hadoop-m gvcf* platinum_genomes.variants.schema .
```

Then ssh to the master and kick off the job:
```
export OUTPUT_PROJECT=YOUR-PROJECT
export OUTPUT_DATASET=YOUR-DATASET
export OUTPUT_TABLE=YOUR-TABLE

export SCHEMA=`cat platinum_genomes.variants.schema`
hadoop jar /home/hadoop/hadoop-install/contrib/streaming/hadoop-streaming-1.2.1.jar \
      -D mapred.bq.input.project.id=genomics-public-data \
      -D mapred.bq.input.dataset.id=platinum_genomes \
      -D mapred.bq.input.table.id=variants \
      -D mapred.bq.output.project.id=$OUTPUT_PROJECT \
      -D mapred.bq.output.dataset.id=$OUTPUT_DATASET \
      -D mapred.bq.output.table.id=$OUTPUT_TABLE \
      -D mapred.bq.output.table.schema=${SCHEMA} \
      -D mapred.output.committer.class=com.google.cloud.hadoop.io.bigquery.mapred.BigQueryMapredOutputCommitter \
      -inputformat com.google.cloud.hadoop.io.bigquery.mapred.BigQueryMapredInputFormat \
      -input requiredButUnused \
      -outputformat com.google.cloud.hadoop.io.bigquery.mapred.BigQueryMapredOutputFormat \
      -output requiredButUnused \
      -file gvcf_expander.py \
      -mapper gvcf-expand-mapper.py -file gvcf-expand-mapper.py  \
      -reducer gvcf-expand-reducer.py -file gvcf-expand-reducer.py
```

Lastly, note that the BigQuery connector when used with Hadoop Streaming does not clean up its temporary files automatically.
```
export CONFIGBUCKET=YOUR-HADOOP-BUCKET
gsutil rm -r gs://$CONFIGBUCKET/hadoop/tmp/bigquery/job_*
```

Caveat
------
Note: using BigQuery as the sink appears to be problematic.  For now, send the output to Google Cloud Storage instead.

Edit [gvcf-expand-reducer.py](./gvcf-expand-reducer.py) to set `BIG_QUERY_SINK=False` and then:
```
hadoop jar /home/hadoop/hadoop-install/contrib/streaming/hadoop-streaming-1.2.1.jar \
-D mapred.bq.input.project.id=genomics-public-data \
-D mapred.bq.input.dataset.id=platinum_genomes \
-D mapred.bq.input.table.id=variants \
-inputformat com.google.cloud.hadoop.io.bigquery.mapred.BigQueryMapredInputFormat \
-input requiredButUnused  \
-file gvcf_expander.py \
-mapper gvcf-expand-mapper.py -file gvcf-expand-mapper.py \
-reducer gvcf-expand-reducer.py -file gvcf-expand-reducer.py \
-output gs://YOUR-BUCKET/YOUR_OUTPUT_PREFIX
```

Finally, load the json from the Google Cloud Storage destination to BigQuery using schema [platinum_genomes.variants.schema](./platinum_genomes.variants.schema) 

Appendix
--------

Here's how we grabbed the schema from the source table.  Do something similar when running this job against a different table.
```
bq --project genomics-public-data show --format json platinum_genomes.variants | python -c "import json,sys ; print \"'%s'\" % (json.dumps(json.loads(sys.stdin.readline())['schema']['fields']).replace(\"'\", \"_\"))" > platinum_genomes.variants.schema
```
