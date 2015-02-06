# Retrieve non-variant segments for BRCA1, flattening by sample.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
FROM
  [_THE_TABLE_]
WHERE
  reference_name = 'chr17'
  AND start BETWEEN 41196311
  AND 41277499
OMIT RECORD IF SOME(alternate_bases IS NOT NULL)
ORDER BY
  start,
  call.call_set_name
