#standardSQL
--
-- Count the number of variant calls per genome.
--
SELECT
  call.call_set_name,
  COUNT(call.call_set_name) AS number_of_calls
FROM
  `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v, v.call call
WHERE
  # Skip homozygous reference calls, no-calls, and non-passing variants.
  EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
  AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
  AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
GROUP BY
  call_set_name
ORDER BY
  call_set_name
