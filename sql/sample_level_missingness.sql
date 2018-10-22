#standardSQL
--
-- Compute the ratio of positions corresponding to no-calls versus all positions
-- called (reference, variant, and no-calls).
--
WITH deltas AS (
  SELECT
    end_position - start_position AS delta,
    c.name,
    EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g < 0)
      # Include low quality calls.
      OR {{ LOW_QUALITY_CALLS_FILTER }} AS has_no_calls
  FROM
    `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.call) AS c
),

positions_called AS (
  SELECT
    name,
    SUM(IF(has_no_calls, delta, 0)) AS no_calls,
    SUM(delta) AS all_calls
  FROM deltas
  GROUP BY
    name
)

SELECT
  name,
  no_calls,
  --- Cast to a float since the values can sometimes be larger than int32.
  --- Environments such as R do not have native support for int64.
  CAST(all_calls AS FLOAT64) AS all_calls,
  (no_calls/all_calls) AS missingness_rate
FROM positions_called
ORDER BY
  name
