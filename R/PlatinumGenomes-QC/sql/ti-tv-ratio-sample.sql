SELECT
  reference_name,
  call_call_set_name,
  transitions,
  transversions,
  transitions/transversions AS titv,
FROM (
  SELECT
    reference_name,
    call_call_set_name,
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
  FROM (
    SELECT
      reference_name,
      call_call_set_name,
      reference_bases,
      alternate_bases,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      num_alts
    FROM
      (SELECT
        reference_name,
        call.call_set_name as call_call_set_name,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        
      FROM
        [_GENOME_CALL_TABLE_]
        OMIT call IF EVERY(call.genotype <= 0) OR SOME(call.FILTER NOT IN ('PASS'))
      # Optionally add clause here to limit the query to a particular
      # region of the genome.
      #_WHERE_
      HAVING
        # Skip 1/2 genotypes
        num_alts = 1
        AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')))
     
  GROUP BY
    call_call_set_name, reference_name)
ORDER BY
  titv desc
