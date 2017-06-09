#standardSQL
--
-- Retrieve variant-level information for BRCA1 variants.
--
SELECT
  reference_name,
  start,
  `end`,
  reference_bases AS ref,
  ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat,
  quality,
  ARRAY_TO_STRING(v.filter, ',') AS filters,
  ARRAY_TO_STRING(v.names, ',') AS names,
  ARRAY_LENGTH(v.call) AS num_samples
FROM
  `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v
WHERE
  reference_name IN ('17', 'chr17')
  AND start BETWEEN 41196311 AND 41277499 # per GRCh37
  # Skip non-variant segments.
  AND EXISTS (SELECT alt FROM UNNEST(v.alternate_bases) alt WHERE alt NOT IN ("<NON_REF>", "<*>"))
ORDER BY
  start,
  alt_concat
