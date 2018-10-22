#standardSQL
--
-- Compute the homozygous and heterozygous variant counts for each individual
-- within chromosome X to help determine whether the sex phenotype value is
-- correct for each individual.
--
WITH filtered_snp_calls AS (
  SELECT
    c.name,
    CAST((SELECT LOGICAL_AND(g > 0) FROM UNNEST(c.genotype) AS g) AS INT64) AS hom_AA,
    CAST(EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g > 0)
      AND EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g = 0) AS INT64) AS het_RA
  FROM
    `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.call) AS c
  WHERE
    reference_name IN ('chrX', 'X')
    AND start_position NOT BETWEEN {{ PAR1_START }} AND {{ PAR1_END }}
    AND start_position NOT BETWEEN {{ PAR2_START }} AND {{ PAR2_END }}
    # Only include biallelic snps.
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases[ORDINAL(1)].alt IN ('A','C','G','T')
    AND (ARRAY_LENGTH(alternate_bases) = 1
      OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)].alt = '<*>'))
    # Include only high quality calls.
    AND {{ HIGH_QUALITY_CALLS_FILTER }}
)

SELECT
  name,
  ROUND(SAFE_DIVIDE(SUM(het_RA), SUM(hom_AA) + SUM(het_RA)), 3) AS perct_het_alt_in_snvs,
  ROUND(SAFE_DIVIDE(SUM(hom_AA), SUM(hom_AA) + SUM(het_RA)), 3) AS perct_hom_alt_in_snvs,
  SUM(hom_AA) AS hom_AA_count,
  SUM(het_RA) AS het_RA_count
FROM filtered_snp_calls
GROUP BY
  name
ORDER BY
  name
