#standardSQL
--
-- The following query computes the Hardy-Weinberg equilibrium for variants.
--
WITH variants AS (
  SELECT
    reference_name,
    start,
    `end`,
    reference_bases,
    v.alternate_bases[ORDINAL(1)] AS alt,
    -- Within each variant count the number of HOM_REF/HOM_ALT/HET samples.
    (SELECT SUM(CAST((SELECT LOGICAL_AND(gt = 0)
      FROM UNNEST(call.genotype) gt) AS INT64)) FROM v.call) AS HOM_REF,
    (SELECT SUM(CAST((SELECT LOGICAL_AND(gt = 1)
      FROM UNNEST(call.genotype) gt) AS INT64)) FROM v.call) AS HOM_ALT,
    (SELECT SUM(CAST((SELECT LOGICAL_OR(gt = 0) AND LOGICAL_OR(gt = 1)
      FROM UNNEST(call.genotype) gt) AS INT64)) FROM v.call) AS HET
  FROM
    `@MULTISAMPLE_VARIANT_TABLE` v
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
),

observations AS (
  SELECT
    reference_name,
    start,
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
    start,
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
  start,
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
  IF((POW(OBS_HOM1 - E_HOM1, 2) / E_HOM1)
  + (POW(OBS_HET - E_HET, 2) / E_HET)
  + (POW(OBS_HOM2 - E_HOM2, 2) / E_HOM2)
  > 5.991, "TRUE", "FALSE") AS PVALUE_SIG

FROM expectations
WHERE
  E_HOM1 > 0 AND E_HET > 0 AND E_HOM2 > 0
-- Optionally add a clause here to constrain the results.
