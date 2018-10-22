#standardSQL
--
-- Compute the ratio of homozygous vs. heterozygous variant calls for each individual.
--
WITH filtered_snp_calls AS (
  SELECT
    reference_name,
    reference_bases AS ref,
    ARRAY_TO_STRING(ARRAY(SELECT a.alt FROM UNNEST(alternate_bases) AS a), ',') AS alt_concat,
    c.name,
    c.genotype[ORDINAL(1)] AS first_allele,
    c.genotype[ORDINAL(2)] AS second_allele
  FROM
    `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.call) AS c
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND alternate_bases[ORDINAL(1)].alt IN ('A','C','G','T')
    AND (ARRAY_LENGTH(alternate_bases) = 1
      OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)].alt = '<*>'))
    # Skip homozygous reference calls, no-calls, and non-diploid calls.
    AND EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g > 0)
    AND NOT EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g < 0)
    AND ARRAY_LENGTH(c.genotype) = 2
    # Include only high quality calls.
    AND {{ HIGH_QUALITY_CALLS_FILTER }}
),

variant_counts AS (
  SELECT
    name,
    reference_name,
    SUM(CAST(first_allele = 1 AND second_allele = 1 AS INT64)) AS HOM_ALT,
    SUM(CAST(first_allele = 1 OR second_allele = 1 AS INT64))  AS HAS_ALT,
    COUNT(name) AS N_SITES
  FROM filtered_snp_calls
  GROUP BY
    name,
    reference_name
)

SELECT
  name,
  reference_name,
  HOM_ALT,
  HAS_ALT,
  N_SITES,
  ROUND((HOM_ALT) / (HAS_ALT), 5) AS F
FROM variant_counts
ORDER BY
  name,
  reference_name
