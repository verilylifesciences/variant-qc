# Compute the ratio of positions corresponding to no-calls versus all positions
# called (reference, variant, and no-calls).
WITH deltas AS (
  SELECT
    `end` - start AS delta,
    call.call_set_name,
    EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
      OR
      EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.')) AS has_no_calls
  FROM
    `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v, v.call call
),

positions_called AS (
  SELECT
    call_set_name,
    SUM(IF(has_no_calls, delta, 0)) AS no_calls,
    SUM(delta) AS all_calls
  FROM deltas
  GROUP BY
    call_set_name
)

SELECT
  call_set_name,
  no_calls,
  all_calls,
  (no_calls/all_calls) AS missingness_rate
FROM positions_called
ORDER BY
  call_set_name
