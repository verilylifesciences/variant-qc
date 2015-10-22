# Retrieve non-variant segments for BRCA1, implicitly flattening by sample.
SELECT
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
FROM
  [_GENOME_CALL_TABLE_]
WHERE
  reference_name CONTAINS '17' # To match both 'chr17' and '17'
  AND start BETWEEN 41196311
  AND 41277499
# In some datasets, alternate_bases will be empty (therefore NULL) for non-variant segments.
# In other datasets, alternate_bases will have the value "<NON_REF>" for non-variant segments.
OMIT RECORD IF (SOME(alternate_bases IS NOT NULL) AND SOME(alternate_bases != "<NON_REF>"))
ORDER BY
  start,
  call.call_set_name
LIMIT
  10000
