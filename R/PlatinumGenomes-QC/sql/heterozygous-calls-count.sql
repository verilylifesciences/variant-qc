  SELECT
    call.call_set_name,
    SUM(first_allele != second_allele) AS O_HET
  FROM (
    SELECT
      reference_name,
      start,
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      call.call_set_name,
      NTH(1, call.genotype) WITHIN call AS first_allele,
      NTH(2, call.genotype) WITHIN call AS second_allele,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts
    FROM
      [_THE_EXPANDED_TABLE_]
    OMIT call IF SOME(call.genotype < 0) OR (2 > COUNT(call.genotype)) 
    HAVING
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
    )
  GROUP BY
    call.call_set_name
