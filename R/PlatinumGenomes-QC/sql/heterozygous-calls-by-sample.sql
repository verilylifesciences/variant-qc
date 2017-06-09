#standardSQL
--
-- Count the number of heterozygous variants per sample.
--
WITH filtered_snp_calls AS (
  SELECT
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
)

SELECT
  call_set_name,
  SUM(CAST((first_allele != second_allele) AS INT64)) AS heterozygous_variant_count
FROM filtered_snp_calls
GROUP BY
  call_set_name
ORDER BY
  call_set_name
