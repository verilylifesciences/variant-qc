#standardSQL
--
-- Count the number of variant calls per genome.
--
SELECT
  c.name,
  reference_name,
  COUNT(c.name) AS number_of_calls
FROM
  `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.call) AS c
WHERE
  # Skip homozygous reference calls, no-calls, and non-passing variants.
  EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g > 0)
  AND NOT EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g < 0)
  # Include only high quality calls.
  AND {{ HIGH_QUALITY_CALLS_FILTER }}
GROUP BY
  name,
  reference_name
ORDER BY
  name,
  reference_name

