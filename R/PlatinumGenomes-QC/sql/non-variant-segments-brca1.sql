# Retrieve non-variant segments for BRCA1.
SELECT
  call.call_set_name,
  (SELECT STRING_AGG(CAST(gt AS STRING)) from UNNEST(call.genotype) gt) AS genotype,
  reference_name,
  start,
  `end`,
  reference_bases AS ref,
  ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat
FROM
  `@GENOME_CALL_TABLE` v, v.call call
WHERE
  reference_name IN ('17', 'chr17')
  AND start BETWEEN 41196311 AND 41277499 # per GRCh37
  # Skip all variant sites.
  AND NOT EXISTS (SELECT alt FROM UNNEST(v.alternate_bases) alt WHERE alt NOT IN ("<NON_REF>", "<*>"))
ORDER BY
  start,
  call.call_set_name
LIMIT
  10000
