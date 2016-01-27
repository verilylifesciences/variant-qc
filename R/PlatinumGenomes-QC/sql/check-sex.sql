# Compute homozygosity (inbreeding) rate on chromosome X for each individual,
# to help determine the correctness of the sex phenotype value.
SELECT
  call.call_set_name,
  O_HOM,
  ROUND(E_HOM, 2) as E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM (
  SELECT
    call.call_set_name,
    SUM(first_allele = second_allele) AS O_HOM,
    SUM(1.0 - (2.0 * freq * (1.0 - freq) * (called_allele_count / (called_allele_count - 1.0)))) AS E_HOM,
    COUNT(call.call_set_name) AS N_SITES,
  FROM (
    SELECT
      reference_name,
      start,
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      call.call_set_name,
      NTH(1, call.genotype) WITHIN call AS first_allele,
      NTH(2, call.genotype) WITHIN call AS second_allele,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      SUM(call.genotype >= 0) WITHIN RECORD AS called_allele_count,
      IF((SUM(1 = call.genotype) > 0),
        SUM(call.genotype = 1)/SUM(call.genotype >= 0),
        -1)  WITHIN RECORD AS freq
    FROM
      [_MULTISAMPLE_VARIANT_TABLE_]
    WHERE
      (reference_name = 'chrX' OR reference_name = 'X')
      # Omit pseudoautosomal regions.
      AND start NOT BETWEEN 59999 AND 2699519
      AND start NOT BETWEEN 154931042 AND 155260559
    # Skip no calls and haploid sites
    OMIT call IF SOME(call.genotype < 0) OR (2 != COUNT(call.genotype))
    HAVING
      # Skip 1/2 genotypes.
      num_alts = 1
      # Only use SNPs since non-variant segments are only included for SNPs.
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      # Skip records where all samples have the same allele.
      AND freq > 0 AND freq < 1 
      )
  GROUP BY
    call.call_set_name
    )
ORDER BY
  call.call_set_name
