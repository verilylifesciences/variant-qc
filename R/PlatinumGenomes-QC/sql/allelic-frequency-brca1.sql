# The following query computes the allelic frequency for BRCA1 SNPs.
SELECT
  vars.reference_name,
  vars.start,
  reference_bases AS ref,
  alternate_bases AS alt,
  (ref_count + SUM(called_allele_count))/(ref_count + SUM(called_allele_count) + alt_count) AS ref_freq,
  alt_count/(ref_count + SUM(called_allele_count) + alt_count) AS alt_freq,
  ref_count,
  alt_count,
  SUM(called_allele_count) AS called_allele_count,
FROM (
    # Constrain the left hand side of the _join to reference-matching blocks.
  SELECT
    reference_name,
    start,
    END,
    SUM(0 = call.genotype) WITHIN RECORD AS called_allele_count,
  FROM
    [_THE_TABLE_]
  WHERE
    reference_name = 'chr17'
  OMIT RECORD IF EVERY(alternate_bases IS NOT NULL)
    ) AS refs
JOIN (
  SELECT
    reference_name,
    start,
    END,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    SUM(0 = call.genotype) WITHIN RECORD AS ref_count,
    SUM(1 = call.genotype) WITHIN RECORD AS alt_count,
  FROM
    [_THE_TABLE_]
  WHERE
    reference_name = 'chr17'
    AND start BETWEEN 41196311
    AND 41277499
  HAVING
    # Skip ref-matching blocks, 1/2 genotypes, and non-SNP variants
    num_alts = 1
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
    ) AS vars
  # The _join criteria _is complicated since we are trying to see if a variant
  # overlaps a reference-matching interval.
ON
  vars.reference_name = refs.reference_name
WHERE
  refs.start <= vars.start
  AND refs.END >= vars.start+1
GROUP BY
  vars.reference_name,
  vars.start,
  vars.end,
  ref,
  alt,
  ref_count,
  alt_count
ORDER BY
  vars.reference_name,
  vars.start,
  ref
