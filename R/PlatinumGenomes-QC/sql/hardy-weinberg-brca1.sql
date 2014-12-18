# The following query computes the Hardy-Weinberg equilibrium for BRCA1 SNPs.
SELECT
  vars.reference_name AS CHR,
  vars.start AS POS,
  reference_bases AS ref,
  alternate_bases AS alt,
  SUM(refs.HOM_REF) + vars.HOM_REF AS OBS_HOM1,
  vars.HET AS OBS_HET,
  vars.HOM_ALT AS OBS_HOM2,
  # Expected AA
  # p^2
  # ((COUNT(AA) + (COUNT(Aa)/2) / 
  #  TOTAL_SAMPLE_COUNT) ^ 2) * 
  #  TOTAL_SAMPLE_COUNT
  
  POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)
    AS EXP_HOM1,
  
  # Expected Aa
  # 2pq
  # 2 *
  # ((COUNT(AA) + (COUNT(Aa)/2) / 
  #  TOTAL_SAMPLE_COUNT) ^ 2) * 
  # (COUNT(aa) + (COUNT(Aa)/2) / 
  #  TOTAL_SAMPLE_COUNT) ^ 2 *
  # TOTAL_SAMPLE_COUNT
  
  2 *
  ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
  ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)
    AS EXP_HET,
    
  # Expected aa
  # q^2
  # (COUNT(aa) + (COUNT(Aa)/2) / 
  #  TOTAL_SAMPLE_COUNT) ^ 2 *
  #  TOTAL_SAMPLE_COUNT
  
  POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT) 
    AS EXP_HOM2,
  
  ########################
  # Chi square
  # SUM(((Observed - Expected)^2) / Expected )
    
    # AA chisq value
    # Numerator 
    POW(
    # Observed
    (SUM(refs.HOM_REF) + vars.HOM_REF) 
      - 
    # Expected
    (POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2) 
    /
    # Denominator 
    # Expected  
      (POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
      (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
      (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    
    +
    
    # Aa chisq value
    # Numerator
    POW(
    # Observed
    vars.HET 
      -
    # Expected
    (2 *
    ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2)
    / 
    # Denominator
    # Expected
    (2 *
    ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
      
    +
      
    #aa chisq value
    # Numerator
    POW(
    # Observed
    vars.HOM_ALT - 
    # Expected
    (POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2)
    / 
    # Denominator
    # Expected
    (POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    
    AS CHISQ,
    #########################
    
    #########################
    # p-value boolean
    # Need to recalculate chisq value and determine if it is greater than 
    # the cutoff for significance with two degrees of freedom, 5.991.
    IF ( 
    POW(
    # Observed
    (SUM(refs.HOM_REF) + vars.HOM_REF) 
      - 
    # Expected
    (POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2) 
    /
    # Denominator 
    # Expected  
      (POW(((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
      (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2) *
      (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    
    +
    
    # Aa chisq value
    # Numerator
    POW(
    # Observed
    vars.HET 
      -
    # Expected
    (2 *
    ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2)
    / 
    # Denominator
    # Expected
    (2 *
    ((SUM(refs.HOM_REF) + vars.HOM_REF + (vars.HET / 2))  /  
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    ((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)) *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
      
    +
      
    #aa chisq value
    # Numerator
    POW(
    # Observed
    vars.HOM_ALT - 
    # Expected
    (POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))
    , 2)
    / 
    # Denominator
    # Expected
    (POW(((vars.HOM_ALT + (vars.HET / 2)) / 
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT)), 2)  *
    (SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT))

    > 5.991, "TRUE", "FALSE") AS PVALUE,
    #########################
    
FROM (
    # Constrain the left hand side of the _join to reference-matching blocks.
  SELECT
    reference_name,
    start,
    END,
    SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
  FROM
    [_THE_TABLE_]
  WHERE
    reference_name = 'chr17'
  OMIT
    RECORD IF EVERY(alternate_bases IS NOT NULL)
    ) AS refs
JOIN (
  SELECT
    reference_name,
    start,
    END,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
    SUM(EVERY(1 = call.genotype)) WITHIN call AS HOM_ALT,
    SUM(SOME(0 = call.genotype) AND SOME(1 = call.genotype)) WITHIN call AS HET,
  FROM
    [_THE_TABLE_]
  WHERE
    reference_name = 'chr17'
    AND start BETWEEN 41196311
    AND 41277499
#  OMIT call IF 2 != COUNT(call.genotype)
  HAVING
    # Skip ref-matching blocks, 1/2 genotypes, and non-SNP variants
    num_alts = 1
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases IN ('A','C','G','T')
    ) AS vars
  # The _join criteria _is complicated since we are trying to see if a variant
  # overlaps a reference-matching interval.
ON
  vars.reference_name = refs.reference_name
WHERE
  refs.start <= vars.start
  AND refs.END >= vars.start+1
GROUP BY
  CHR,
  POS,
  ref,
  alt,
  vars.HOM_REF,
  OBS_HET,
  OBS_HOM2
ORDER BY
  CHR,
  POS,
  ref,
  alt
