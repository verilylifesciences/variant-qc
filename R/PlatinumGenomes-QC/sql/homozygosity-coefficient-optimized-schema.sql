#standardSQL
--
-- Compute the expected and observed homozygosity rate for each individual.
-- This query is for the large cohort optimized schema.
--
WITH variant_calls AS (
  SELECT
    call.call_set_name,
    call.genotype[ORDINAL(1)] = 1 AND call.genotype[ORDINAL(2)] = 1 AS O_HOM,
    1.0 - (2.0 * alt[ORDINAL(1)].AF * (1.0 - alt[ORDINAL(1)].AF) * (AN / (AN - 1.0))) AS E_HOM
  FROM
    `@MULTISAMPLE_VARIANT_TABLE` v, v.call call
  WHERE
    # Only include biallelic snps within autosomes.
    REGEXP_CONTAINS(reference_name, r'^(chr)?([1-2])?[0-9]$')
    AND reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alt) = 1
    AND alt[ORDINAL(1)].alternate_bases IN ('A','C','G','T')
    AND alt[ORDINAL(1)].AF > 0
    AND alt[ORDINAL(1)].AF < 1
    # Skip no-calls, non-passing variants, and non-diploid calls.
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
    AND ARRAY_LENGTH(call.genotype) = 2
),

reference_calls AS (
  SELECT
    call_set_name,
    TRUE AS O_HOM,
    1.0 - (2.0 * alt[ORDINAL(1)].AF * (1.0 - alt[ORDINAL(1)].AF) * (AN / (AN - 1.0))) AS E_HOM
  FROM
    `@MULTISAMPLE_VARIANT_TABLE` v, v.refMatchCallsets call_set_name
  WHERE
    # Only include biallelic snps within autosomes. (Concise but inexact regexp used for brevity.)
    REGEXP_CONTAINS(reference_name, r'^(chr)?([1-2])?[0-9]$')
    AND reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alt) = 1
    AND alt[ORDINAL(1)].alternate_bases IN ('A','C','G','T')
    # Ensure allelic frequency is in the range (0, 1).
    AND alt[ORDINAL(1)].AF > 0
    AND alt[ORDINAL(1)].AF < 1
),

grouped_values AS (
  SELECT
    call_set_name,
    SUM(CAST(O_HOM AS INT64)) AS O_HOM,
    SUM(E_HOM) AS E_HOM,
    COUNT(call_set_name) AS N_SITES
  FROM
    (SELECT * FROM variant_calls
    UNION ALL
    SELECT * FROM reference_calls)
  GROUP BY
    call_set_name
)

SELECT
  call_set_name,
  O_HOM,
  ROUND(E_HOM, 2) AS E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM grouped_values
ORDER BY
  call_set_name
