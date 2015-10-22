# Retrieve heterozygous haplotype calls on chromosomes X and Y.
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
  reference_name CONTAINS 'X' OR reference_name CONTAINS 'Y'
OMIT
  call if (2 > COUNT(call.genotype))
  OR EVERY(call.genotype <= 0)
  OR EVERY(call.genotype = 1)
HAVING call.call_set_name IN (_MALE_SAMPLE_IDS_)
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_
