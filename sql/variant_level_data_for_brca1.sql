#standardSQL
--
-- Retrieve variant-level information for BRCA1 variants.
--
SELECT
  reference_name,
  start_position,
  end_position,
  reference_bases AS ref,
  ARRAY_TO_STRING(ARRAY(SELECT a.alt FROM UNNEST(alternate_bases) AS a), ',') AS alt_concat,
  quality,
  ARRAY_TO_STRING(filter, ',') AS filters,
  ARRAY_TO_STRING(names, ',') AS names,
  ARRAY_LENGTH(call) AS num_samples
FROM
  `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}`
WHERE
  reference_name IN ('17', 'chr17')
  AND start_position BETWEEN {{ BRCA1_START }} AND {{ BRCA1_END }}
  # Skip non-variant segments.
  AND EXISTS (SELECT alt FROM UNNEST(alternate_bases) WHERE alt NOT IN ("<NON_REF>", "<*>"))
ORDER BY
  start_position,
  alt_concat
