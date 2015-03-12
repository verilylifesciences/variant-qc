# Compute the the homozygous and heterozygous variant counts for each individual
# within chromosome X to help determine whether the gender phenotype value is
# correct for each individual.
SELECT
  call.call_set_name,
  ROUND(SUM(het_RA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_het_alt_in_snvs,
  ROUND(SUM(hom_AA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_hom_alt_in_snvs,
  SUM(hom_AA) AS hom_AA_count,
  SUM(het_RA) AS het_RA_count,
  SUM(hom_RR) AS hom_RR_count,
FROM (
  SELECT
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    call.call_set_name,
    SOME(call.genotype = 0) AND NOT SOME(call.genotype > 0) WITHIN call AS hom_RR,
    SOME(call.genotype > 0) AND NOT SOME(call.genotype = 0) WITHIN call AS hom_AA,
    SOME(call.genotype > 0) AND SOME(call.genotype = 0) WITHIN call AS het_RA
  FROM
    [_THE_EXPANDED_TABLE_]
  WHERE
    reference_name = 'chrX'
    AND start NOT BETWEEN 59999 AND 2699519
    AND start NOT BETWEEN 154931042 AND 155260559
  HAVING
    # Skip 1/2 genotypes _and non-SNP variants
    num_alts = 1
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
  )
GROUP BY
  call.call_set_name
ORDER BY
  call.call_set_name
