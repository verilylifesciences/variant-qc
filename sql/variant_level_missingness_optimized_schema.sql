#standardSQL
--
-- Compute the ratio of no-calls for each variant.
--
WITH variant_missingness AS (
  SELECT
    reference_name,
    start_position,
    end_position,
    reference_bases,
    ARRAY_TO_STRING(ARRAY(SELECT a.alt FROM UNNEST(alternate_bases) AS a), ',') AS alt_concat,
    ARRAY_LENGTH(refMatchCallsets) * 2
      + (SELECT SUM((SELECT COUNT(g) FROM UNNEST(c.genotype) AS g)) FROM UNNEST(call) AS c) AS all_calls,
    (SELECT SUM((SELECT COUNT(g) FROM UNNEST(c.genotype) AS g WHERE g < 0)) FROM UNNEST(call) AS c) AS no_calls
  FROM
    `{{ MULTISAMPLE_VARIANT_TABLE }}`
)

SELECT
  reference_name,
  start_position,
  end_position,
  reference_bases,
  alt_concat,
  no_calls,
  all_calls,
  no_calls/all_calls AS missingness_rate
FROM variant_missingness
WHERE
  all_calls > 0
ORDER BY missingness_rate DESC, reference_name, start_position, reference_bases, alt_concat
LIMIT 1000
