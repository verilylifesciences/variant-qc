#standardSQL
--
-- Retrieve heterozygous haplotype calls on chromosomes X and Y.
--
SELECT
  reference_name,
  start,
  `end`,
  reference_bases,
  ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat,
  call.call_set_name,
  (SELECT STRING_AGG(CAST(gt AS STRING)) from UNNEST(call.genotype) gt) AS genotype
FROM
  `@@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v, v.call call
WHERE
  reference_name IN ('chrX', 'X', 'chrY', 'Y')
  AND call_set_name IN (@MALE_SAMPLE_IDS)
  AND (SELECT LOGICAL_OR(gt = 0) AND LOGICAL_OR(gt = 1) FROM UNNEST(call.genotype) gt)
-- Optionally add a clause here to constrain the results.
