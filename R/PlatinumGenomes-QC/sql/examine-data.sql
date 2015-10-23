# Examine the data for particular calls.
SELECT
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
FROM
  [_GENOME_CALL_TABLE_]
WHERE
  reference_name = 'chr17'
HAVING
  _HAVING_
ORDER BY
  start,
  end,
  call.call_set_name
