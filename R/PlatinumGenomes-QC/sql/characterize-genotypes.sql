#standardSQL
--
-- Query to show the variety of genotypes.
--
SELECT
  genotype,
  --- Cast to a float since the values can sometimes be larger than int32.
  --- Environments such as R do not have native support for int64.
  CAST(COUNT(genotype) AS FLOAT64) AS genotype_count
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
