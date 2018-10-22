#standardSQL
--
-- Retrieve heterozygous haplotype calls on chromosomes X and Y for male samples.
--
SELECT
  reference_name,
  start_position,
  end_position,
  reference_bases,
  ARRAY_TO_STRING(ARRAY(SELECT a.alt FROM UNNEST(alternate_bases) AS a), ',') AS alt_concat,
  c.name,
  (SELECT STRING_AGG(CAST(g AS STRING)) from UNNEST(c.genotype) AS g) AS genotype
FROM
  `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.call) AS c
WHERE
  reference_name IN ('chrX', 'X', 'chrY', 'Y')
  AND name IN ({{ MALE_SAMPLES_QUERY }})
  AND c.genotype[SAFE_ORDINAL(1)] != c.genotype[SAFE_ORDINAL(2)]
  # Skip no-calls.
  AND NOT EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g < 0)
ORDER BY reference_name, start_position, alt_concat, name
LIMIT 1000
