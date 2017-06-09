#standardSQL
--
-- Compute the ratio of no-calls for each variant.
--
WITH variant_missingness AS (
  SELECT
    reference_name,
    start,
    `end`,
    reference_bases,
    ARRAY_TO_STRING(v.alternate_bases, ',') AS alt_concat,
    ARRAY_LENGTH(v.refMatchCallsets) * 2
      + (SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt)) FROM v.call) AS all_calls,
    (SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt WHERE gt < 0)) FROM v.call) AS no_calls
  FROM
    `@MULTISAMPLE_VARIANT_TABLE` v
)

SELECT
  reference_name,
  start,
  `end`,
  reference_bases,
  alt_concat,
  no_calls,
  all_calls,
  no_calls/all_calls AS missingness_rate
FROM variant_missingness
WHERE
  all_calls > 0
-- Optionally add a clause here to constrain the results.
