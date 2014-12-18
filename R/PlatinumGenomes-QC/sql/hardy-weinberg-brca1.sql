# The following query computes the Hardy-Weinberg equilibrium for BRCA1 SNPs.
SELECT 
  calcs.CHR AS CHR,
  calcs.POS AS POS,
  calcs.ref AS ref,
  calcs.alt AS alt,
  calcs.OBS_HOM1 AS OBS_HOM1,
  calcs.OBS_HET AS OBS_HET,
  calcs.OBS_HOM2 AS OBS_HOM2,
  calcs.EXP_HOM1 AS EXP_HOM1,
  calcs.EXP_HET AS EXP_HET,
  calcs.EXP_HOM2 AS EXP_HOM2,
  
  # Chi Squared Calculation
  # SUM(((Observed - Expected)^2) / Expected )
  ROUND((POW(calcs.OBS_HOM1 - calcs.EXP_HOM1, 2) / calcs.EXP_HOM1)
  + (POW(calcs.OBS_HET - calcs.EXP_HET, 2) / calcs.EXP_HET)
  + (POW(calcs.OBS_HOM2 - calcs.EXP_HOM2, 2) / calcs.EXP_HOM2), 3)
  AS CHI_SQ,
  
  # Determine if Chi Sq value is significant
  IF((POW(calcs.OBS_HOM1 - calcs.EXP_HOM1, 2) / calcs.EXP_HOM1)
  + (POW(calcs.OBS_HET - calcs.EXP_HET, 2) / calcs.EXP_HET)
  + (POW(calcs.OBS_HOM2 - calcs.EXP_HOM2, 2) / calcs.EXP_HOM2) 
  > 5.991, "TRUE", "FALSE") AS PVALUE_SIG
  
FROM (
    SELECT
      vals.CHR AS CHR,
      vals.POS AS POS,
      vals.ref AS ref,
      vals.alt AS alt,
      vals.OBS_HOM1 AS OBS_HOM1,
      vals.OBS_HET AS OBS_HET,
      vals.OBS_HOM2 AS OBS_HOM2,
    
      # Expected AA
      # p^2
      # ((COUNT(AA) + (COUNT(Aa)/2) / 
      #  SAMPLE_COUNT) ^ 2) * SAMPLE_COUNT
      ROUND(POW((vals.OBS_HOM1 + (vals.OBS_HET/2)) /
        vals.SAMPLE_COUNT, 2) * vals.SAMPLE_COUNT, 2)
        AS EXP_HOM1,
    
      # Expected Aa
      # 2pq
      # 2 * (COUNT(AA) + (COUNT(Aa)/2) / SAMPLE_COUNT) * 
      # (COUNT(aa) + (COUNT(Aa)/2) / SAMPLE_COUNT) 
      # * SAMPLE_COUNT
      ROUND(2 * ((vals.OBS_HOM1 + (vals.OBS_HET/2)) / vals.SAMPLE_COUNT) *
        ((vals.OBS_HOM2 + (vals.OBS_HET/2)) / vals.SAMPLE_COUNT) 
        * vals.SAMPLE_COUNT, 2)
        AS EXP_HET,
    
      # Expected aa
      # q^2
      # (COUNT(aa) + (COUNT(Aa)/2) / 
      #  SAMPLE_COUNT) ^ 2 * SAMPLE_COUNT    
      ROUND(POW((vals.OBS_HOM2 + (vals.OBS_HET/2)) /
        vals.SAMPLE_COUNT, 2) * vals.SAMPLE_COUNT, 2)
        AS EXP_HOM2,
      
    FROM (
        SELECT
          vars.reference_name AS CHR,
          vars.start AS POS,
          reference_bases AS ref,
          alternate_bases AS alt,
          SUM(refs.HOM_REF) + vars.HOM_REF AS OBS_HOM1,
          vars.HET AS OBS_HET,
          vars.HOM_ALT AS OBS_HOM2, 
          SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT AS SAMPLE_COUNT,
        
        FROM (
              # Constrain the left hand side of the _join to reference-matching blocks.
            SELECT
              reference_name,
              start,
              END,
              SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
            FROM
              [genomics-public-data:platinum_genomes.variants]
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
                [genomics-public-data:platinum_genomes.variants]
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
              OBS_HOM2,
              vars.HET,
              vars.HOM_ALT
            ORDER BY
              CHR,
              POS,
              ref,
              alt ) AS vals ) AS calcs
  
