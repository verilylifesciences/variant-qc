# Transition/Transversion Ratio by Depth of Coverage.
SELECT
  call.call_set_name,
  (transitions/transversions) AS titv_ratio,
  call.DP AS average_depth,
FROM (
  SELECT
    call.call_set_name,
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
    call.DP
  FROM (
    SELECT
      reference_bases,
      alternate_bases,
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      call.DP
    FROM (
      SELECT
        call.call_set_name,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        call.genotype,
        call.DP,
      FROM
        [_MULTISAMPLE_VARIANT_TABLE_]
      # Optionally add clause here to limit the query to a particular
      # region of the genome.
      #_WHERE_
    )
    OMIT call if call.DP IS NULL
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T'))
  GROUP BY
    call.call_set_name,
    call.DP)
WHERE
  transversions > 0
GROUP BY
  call.call_set_name,
  titv_ratio,
  average_depth
ORDER BY
  call.call_set_name,
  average_depth
