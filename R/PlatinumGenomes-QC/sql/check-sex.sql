#standardSQL
--
-- Compute the homozygous and heterozygous variant counts for each individual
-- within chromosome X to help determine whether the sex phenotype value is
-- correct for each individual.
--
WITH filtered_snp_calls AS (
  SELECT
    call.call_set_name,
    CAST((SELECT LOGICAL_AND(gt > 0) FROM UNNEST(call.genotype) gt) AS INT64) AS hom_AA,
    CAST(EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
      AND EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt = 0) AS INT64) AS het_RA
  FROM
    `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v, v.call call
  WHERE
    reference_name IN ('chrX', 'X')
    # Locations of PAR1 and PAR2 on GRCh37.
    AND start NOT BETWEEN 59999 AND 2699519
    AND start NOT BETWEEN 154931042 AND 155260559
    # Only include biallelic snps.
    AND reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    # Skip non-passing calls.
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
)

SELECT
  call_set_name,
  ROUND(SUM(het_RA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_het_alt_in_snvs,
  ROUND(SUM(hom_AA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_hom_alt_in_snvs,
  SUM(hom_AA) AS hom_AA_count,
  SUM(het_RA) AS het_RA_count
FROM filtered_snp_calls
GROUP BY
  call_set_name
ORDER BY
  call_set_name
