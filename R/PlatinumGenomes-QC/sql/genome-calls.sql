# Count the number of calls per genome.
SELECT
  call.call_set_name,
  COUNT(call.call_set_name) AS number_of_calls,
FROM
  [_GENOME_CALL_TABLE_]
GROUP BY
  call.call_set_name
