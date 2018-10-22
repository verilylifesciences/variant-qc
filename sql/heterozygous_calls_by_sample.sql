#standardSQL
--
-- Count the number of heterozygous variants per sample.
--
WITH filtered_snp_calls AS (
  SELECT
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
)

SELECT
  name,
  SUM(CAST((first_allele != second_allele) AS INT64)) AS heterozygous_variant_count
FROM filtered_snp_calls
GROUP BY
  name
ORDER BY
  name
