# Retrieve variant-level information for BRCA1 variants.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
  GROUP_CONCAT(names) WITHIN RECORD AS names,
  COUNT(call.call_set_name) WITHIN RECORD AS num_samples,
FROM
  [_GENOME_CALL_TABLE_]
WHERE
  reference_name CONTAINS '17' # To match both 'chr17' and '17'
  AND start BETWEEN 41196311
  AND 41277499
# In some datasets, alternate_bases will be empty (therefore NULL) for non-variant segments.
# In other datasets, alternate_bases will have the value "<NON_REF>" for non-variant segments.
OMIT RECORD IF EVERY(alternate_bases IS NULL) OR EVERY(alternate_bases = "<NON_REF>")
ORDER BY
  start,
  alternate_bases
