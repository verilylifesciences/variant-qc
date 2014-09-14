# Sample alternate allele counts within this bi-allelic data.
SELECT
  genotype.sample_id AS Sample,
  SUM(IF((genotype.first_allele > 0
        AND genotype.second_allele = 0)
      OR (genotype.first_allele = 0
        AND genotype.second_allele > 0),
      1,
      0))
  AS single,
  SUM(IF(genotype.first_allele > 0
      AND genotype.second_allele > 0,
      1,
      0)) AS double,
FROM
  [google.com:biggene:1000genomes.variants1kG]
OMIT
  RECORD IF contig IN ("X",
    "Y")
GROUP BY
  Sample
ORDER BY
  Sample