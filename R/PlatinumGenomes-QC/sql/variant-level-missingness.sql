# Compute the ratio no-calls for each variant.
SELECT
  reference_name,
  start,
  END,
  reference_bases,
  alternate_bases,
  no_calls,
  all_calls,
  (no_calls/all_calls) AS missingness_rate
FROM (
  SELECT
    reference_name,
    start,
    END,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    SUM(call.genotype == -1) WITHIN RECORD AS no_calls,
    COUNT(call.genotype) WITHIN RECORD AS all_calls,
  FROM
      [_MULTISAMPLE_VARIANT_TABLE_]
  # Optionally add clause here to limit the query to a particular
  # region of the genome.
  #_WHERE_
  HAVING
    # Only use SNPs since non-variant segments are only included for SNPs.
    reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
  )
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_
