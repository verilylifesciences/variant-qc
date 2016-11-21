# Compute the expected and observed homozygosity rate for each individual.
WITH variant_calls AS (
  SELECT
    call.call_set_name,
    call.genotype[ORDINAL(1)] = call.genotype[ORDINAL(2)] AS O_HOM,
    (SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt WHERE gt >= 0)) FROM v.call) AS AN,
    (SELECT SUM((SELECT COUNT(gt) FROM UNNEST(call.genotype) gt WHERE gt = 1)) FROM v.call) AS AC
  FROM
    `@MULTISAMPLE_VARIANT_TABLE` v, v.call call
  WHERE
    # Only include biallelic snps within autosomes. (Concise but inexact regexp used for brevity.)
    REGEXP_CONTAINS(reference_name, r'^(chr)?([1-2])?[0-9]$')
    AND reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(v.alternate_bases) = 1
    AND v.alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
    # Skip no-calls, non-passing variants, and non-diploid calls.
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
    AND ARRAY_LENGTH(call.genotype) = 2
),

grouped_values AS (
  SELECT
    call_set_name,
    SUM(CAST(O_HOM AS INT64)) AS O_HOM,
    SUM(1.0 - (2.0 * (AC/AN) * (1.0 - (AC/AN)) * (AN / (AN - 1.0)))) AS E_HOM,
    COUNT(call_set_name) AS N_SITES
  FROM variant_calls
  WHERE
    # Ensure allelic frequency is in the range (0, 1).
    AN > 0 AND AC < AN
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
