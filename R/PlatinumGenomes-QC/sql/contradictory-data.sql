# An example of contradictory data.
SELECT
  reference_name,
  start,
  END,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS gt,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
  GROUP_CONCAT(STRING(call.pl)) WITHIN call AS likelihood,
FROM
  [genomics-public-data:platinum_genomes.variants]
WHERE
  reference_name = 'chr17'
  AND start <= 41198773
  AND END >= 41198774
HAVING
  call.call_set_name = 'NA12883'
ORDER BY
  start,
  END
