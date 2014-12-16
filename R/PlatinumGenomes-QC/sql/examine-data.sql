# Examine the data for particular calls.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS gt,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
  GROUP_CONCAT(STRING(call.pl)) WITHIN call AS likelihood,
FROM
  [_THE_TABLE_]
WHERE
  reference_name = 'chr17'
HAVING
  _HAVING_
ORDER BY
  start,
  end,
  call.call_set_name
