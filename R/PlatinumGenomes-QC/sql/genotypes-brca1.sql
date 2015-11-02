# Query to show the variety of genotypes within BRCA1, such as
# single allele genotypes.
SELECT
  genotype,
  COUNT(genotype) AS genotype_count
FROM (
  SELECT
    GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  FROM
  [_GENOME_CALL_TABLE_]
  WHERE
    reference_name CONTAINS '17' # To match both 'chr17' and '17'
    AND start BETWEEN 41196311
    AND 41277499
    )
GROUP BY
  genotype
ORDER BY
  genotype_count DESC
