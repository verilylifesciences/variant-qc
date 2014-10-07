# Count alternate alleles for each sample.
SELECT
  Sample,
  SUM(IF((first_allele > 0
        AND second_allele = 0)
      OR (first_allele = 0
        AND second_allele > 0),
      1,
      0))
  AS single,
  SUM(IF(first_allele > 0
      AND second_allele > 0,
      1,
      0))
  AS double,
FROM
  (
  SELECT
    reference_name,
    call.call_set_name AS Sample,
    NTH(1,
      call.genotype) WITHIN call AS first_allele,
    NTH(2,
      call.genotype) WITHIN call AS second_allele
  FROM
    [genomics-public-data:1000_genomes.variants])
OMIT
  RECORD IF
  reference_name IN ("X",
    "Y",
    "MT")
GROUP BY
  Sample
ORDER BY
  Sample
