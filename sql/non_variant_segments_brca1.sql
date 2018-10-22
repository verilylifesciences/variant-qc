#standardSQL
--
-- Retrieve non-variant segments for BRCA1.
--
SELECT
  name,
  (SELECT STRING_AGG(CAST(g AS STRING)) from UNNEST(c.genotype) AS g) AS genotype,
  reference_name,
  start_position,
  end_position,
  reference_bases AS ref,
  ARRAY_TO_STRING(ARRAY(SELECT a.alt FROM UNNEST(alternate_bases) AS a), ',') AS alt_concat
FROM
  `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.call) AS c
WHERE
  reference_name IN ('17', 'chr17')
  AND start_position BETWEEN {{ BRCA1_START }} AND {{ BRCA1_END }}
  # Skip all variant sites.
  AND NOT EXISTS (SELECT a.alt FROM UNNEST(alternate_bases) AS a WHERE alt NOT IN ("<NON_REF>", "<*>"))
ORDER BY
  start_position,
  name
LIMIT
  10000
