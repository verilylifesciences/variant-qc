# Individual Homozygosity
SELECT
  call.call_set_name AS INDV,
  SUM(first_allele = second_allele) AS O_HOM,
  COUNT(call.call_set_name) AS N_SITES,
FROM (
  SELECT
    call.call_set_name,
    NTH(1,
      call.genotype) WITHIN call AS first_allele,
    NTH(2,
      call.genotype) WITHIN call AS second_allele,
    COUNT(call.genotype) WITHIN call AS ploidy
  FROM
    [_THE_TABLE_]
  WHERE
    reference_name = 'chr17'
    AND start BETWEEN 41196311
    AND 41277499
  OMIT RECORD IF EVERY(alternate_bases IS NULL)
  HAVING
    ploidy = 2
    )
GROUP BY
  INDV
ORDER BY
  INDV
