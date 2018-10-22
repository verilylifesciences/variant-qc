#standardSQL
--
-- Compute the expected and observed homozygosity rate for each individual.
-- This query is for the large cohort optimized schema.
--
WITH variant_calls AS (
  SELECT
    c.name,
    c.genotype[ORDINAL(1)] = 1 AND c.genotype[ORDINAL(2)] = 1 AS O_HOM,
    1.0 - (2.0 * alternate_bases[ORDINAL(1)].AF * (1.0 - alternate_bases[ORDINAL(1)].AF) * (AN / (AN - 1.0))) AS E_HOM
  FROM
    `{{ MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.call) AS c
  WHERE
    # Only include biallelic snps within autosomes.
    REGEXP_CONTAINS(reference_name, r'^(chr)?([1-2])?[0-9]$')
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases[ORDINAL(1)].alt IN ('A','C','G','T')
    AND (ARRAY_LENGTH(alternate_bases) = 1
      OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)].alt = '<*>'))
    AND alternate_bases[ORDINAL(1)].AF > 0
    AND alternate_bases[ORDINAL(1)].AF < 1
    # Skip no-calls and non-diploid calls.
    AND NOT EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g < 0)
    AND ARRAY_LENGTH(c.genotype) = 2
    # Include only high quality calls.
    AND {{ HIGH_QUALITY_CALLS_FILTER }}
),

reference_calls AS (
  SELECT
    name,
    TRUE AS O_HOM,
    1.0 - (2.0 * alternate_bases[ORDINAL(1)].AF * (1.0 - alternate_bases[ORDINAL(1)].AF) * (AN / (AN - 1.0))) AS E_HOM
  FROM
    `{{ MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.refMatchCallsets) AS name
  WHERE
    # Only include biallelic snps within autosomes. (Concise but inexact regexp used for brevity.)
    REGEXP_CONTAINS(reference_name, r'^(chr)?([1-2])?[0-9]$')
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases[ORDINAL(1)].alt IN ('A','C','G','T')
    AND (ARRAY_LENGTH(alternate_bases) = 1
      OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)].alt = '<*>'))
    # Ensure allelic frequency is in the range (0, 1).
    AND alternate_bases[ORDINAL(1)].AF > 0
    AND alternate_bases[ORDINAL(1)].AF < 1
),

grouped_values AS (
  SELECT
    name,
    SUM(CAST(O_HOM AS INT64)) AS O_HOM,
    SUM(E_HOM) AS E_HOM,
    COUNT(name) AS N_SITES
  FROM
    (SELECT * FROM variant_calls
    UNION ALL
    SELECT * FROM reference_calls)
  GROUP BY
    name
)

SELECT
  name,
  O_HOM,
  ROUND(E_HOM, 2) AS E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM grouped_values
ORDER BY
  name
