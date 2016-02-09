# Compute the fraction of homozygous calls for each individual.

SELECT
call.call_set_name,
reference_name,
HOM,
HET_HOM,
N_SITES,
ROUND((HOM) / (HET_HOM), 5) AS F
FROM (
  SELECT
  call.call_set_name,
  reference_name,
  SUM(first_allele = 1 and second_allele = 1) AS HOM,
  SUM(first_allele + second_allele > 0)  AS HET_HOM,
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
    SUM(call.genotype > 0) WITHIN RECORD AS called_alt_allele_count,
    FROM
    [_GENOME_CALL_TABLE_]
    # Optionally add a clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
    # Skip no calls and haploid sites
    OMIT call IF SOME(call.genotype < 0) OR EVERY(call.genotype = 0) OR (2 != COUNT(call.genotype))
    HAVING
    # Skip 1/2 genotypes.
    num_alts = 1
    # Only use SNPs since non-variant segments are only included for SNPs.
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
  )
  GROUP BY
  call.call_set_name,
  reference_name
)
ORDER BY
call.call_set_name,
reference_name

