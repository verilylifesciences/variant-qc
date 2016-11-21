# Query to show the variety of genotypes.
SELECT
  genotype,
  COUNT(genotype) AS genotype_count
FROM (
  SELECT
  (SELECT STRING_AGG(CAST(gt AS STRING)) from UNNEST(call.genotype) gt) AS genotype
  FROM
  `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v, v.call call)
GROUP BY
  genotype
ORDER BY
  genotype_count DESC,
  genotype
