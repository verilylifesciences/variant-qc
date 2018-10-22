#standardSQL
--
-- The following query computes the Hardy-Weinberg equilibrium for variants.
--
WITH variants AS (
  SELECT
    reference_name,
    start_position,
    end_position,
    reference_bases,
    alternate_bases[ORDINAL(1)].alt AS alt,
    -- Within each variant count the number of HOM_REF/HOM_ALT/HET samples.
    (SELECT SUM(CAST((SELECT LOGICAL_AND(g = 0)
      FROM UNNEST(genotype) AS g) AS INT64)) FROM UNNEST(call)) AS HOM_REF,
    (SELECT SUM(CAST((SELECT LOGICAL_AND(g = 1)
      FROM UNNEST(genotype) AS g) AS INT64)) FROM UNNEST(call)) AS HOM_ALT,
    (SELECT SUM(CAST((SELECT LOGICAL_OR(g = 0) AND LOGICAL_OR(g = 1)
      FROM UNNEST(genotype) AS g) AS INT64)) FROM UNNEST(call)) AS HET
  FROM
    `{{ MULTISAMPLE_VARIANT_TABLE }}`
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND alternate_bases[ORDINAL(1)].alt IN ('A','C','G','T')
    AND (ARRAY_LENGTH(alternate_bases) = 1
      OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)].alt = '<*>'))
),

observations AS (
  SELECT
    reference_name,
    start_position,
    reference_bases,
    alt,
    HOM_REF AS OBS_HOM1,
    HET AS OBS_HET,
    HOM_ALT AS OBS_HOM2,
    HOM_REF + HET + HOM_ALT AS SAMPLE_COUNT
  FROM variants
),

expectations AS (
  SELECT
    reference_name,
    start_position,
    reference_bases,
    alt,
    OBS_HOM1,
    OBS_HET,
    OBS_HOM2,

    # Expected AA
    # p^2
    # ((COUNT(AA) + (COUNT(Aa)/2) /
    #  SAMPLE_COUNT) ^ 2) * SAMPLE_COUNT
    POW((OBS_HOM1 + (OBS_HET/2)) /
      SAMPLE_COUNT, 2) * SAMPLE_COUNT
      AS E_HOM1,

    # Expected Aa
    # 2pq
    # 2 * (COUNT(AA) + (COUNT(Aa)/2) / SAMPLE_COUNT) *
    # (COUNT(aa) + (COUNT(Aa)/2) / SAMPLE_COUNT)
    # * SAMPLE_COUNT
    2 * ((OBS_HOM1 + (OBS_HET/2)) / SAMPLE_COUNT) *
      ((OBS_HOM2 + (OBS_HET/2)) / SAMPLE_COUNT)
      * SAMPLE_COUNT
      AS E_HET,

    # Expected aa
    # q^2
    # (COUNT(aa) + (COUNT(Aa)/2) /
    #  SAMPLE_COUNT) ^ 2 * SAMPLE_COUNT
    POW((OBS_HOM2 + (OBS_HET/2)) /
      SAMPLE_COUNT, 2) * SAMPLE_COUNT
      AS E_HOM2

  FROM observations
  WHERE SAMPLE_COUNT > 0
)

SELECT
  reference_name,
  start_position,
  reference_bases,
  alt,
  OBS_HOM1,
  OBS_HET,
  OBS_HOM2,
  E_HOM1,
  E_HET,
  E_HOM2,

  # Chi Squared Calculation
  # SUM(((Observed - Expected)^2) / Expected )
  ROUND((POW(OBS_HOM1 - E_HOM1, 2) / E_HOM1)
  + (POW(OBS_HET - E_HET, 2) / E_HET)
  + (POW(OBS_HOM2 - E_HOM2, 2) / E_HOM2), 6)
  AS ChiSq,

  # Determine if Chi Sq value is significant
  # The chi-squared score corresponding to a nominal P-value of 0.05
  # for a table with 2 degrees of freedom is 5.991.
  (POW(OBS_HOM1 - E_HOM1, 2) / E_HOM1)
  + (POW(OBS_HET - E_HET, 2) / E_HET)
  + (POW(OBS_HOM2 - E_HOM2, 2) / E_HOM2) > 5.991 AS PVALUE_SIG

FROM expectations
WHERE
  E_HOM1 > 0 AND E_HET > 0 AND E_HOM2 > 0
ORDER BY ChiSq DESC, reference_name, start_position, alt
LIMIT 1000
