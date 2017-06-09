#standardSQL
--
-- Compute the ratio of homozygous vs. heterozygous variant calls for each individual.
--
WITH filtered_snp_calls AS (
  SELECT
    reference_name,
    reference_bases AS ref,
    ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat,
    call.call_set_name,
    call.genotype[ORDINAL(1)] AS first_allele,
    call.genotype[ORDINAL(2)] AS second_allele
  FROM
    `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v, v.call call
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    # Skip homozygous reference calls, no-calls, non-passing variants, and non-diploid calls.
    AND EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
    AND ARRAY_LENGTH(call.genotype) = 2
),

variant_counts AS (
  SELECT
    call_set_name,
    reference_name,
    SUM(CAST(first_allele = 1 AND second_allele = 1 AS INT64)) AS HOM_ALT,
    SUM(CAST(first_allele = 1 OR second_allele = 1 AS INT64))  AS HAS_ALT,
    COUNT(call_set_name) AS N_SITES
  FROM filtered_snp_calls
  GROUP BY
    call_set_name,
    reference_name
)

SELECT
  call_set_name,
  reference_name,
  HOM_ALT,
  HAS_ALT,
  N_SITES,
  ROUND((HOM_ALT) / (HAS_ALT), 5) AS F
FROM variant_counts
ORDER BY
  call_set_name,
  reference_name
