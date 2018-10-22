#standardSQL
--
-- Compute the expected and observed homozygosity rate for each individual.
--
WITH filtered_variant_calls AS (
  SELECT
    ARRAY(SELECT STRUCT(c.name, c.genotype) FROM UNNEST(call) AS c
          WHERE
            # Skip no-calls and non-diploid calls.
            NOT EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g < 0)
            AND ARRAY_LENGTH(c.genotype) = 2
            # Include only high quality calls.
            AND {{ HIGH_QUALITY_CALLS_FILTER }}
    ) AS call
  FROM
    `{{ MULTISAMPLE_VARIANT_TABLE }}`
  WHERE
    # Only include biallelic snps within autosomes. (Concise but inexact regexp used for brevity.)
    REGEXP_CONTAINS(reference_name, r'^(chr)?([1-2])?[0-9]$')
    AND reference_bases IN ('A','C','G','T')
    AND alternate_bases[ORDINAL(1)].alt IN ('A','C','G','T')
    AND (ARRAY_LENGTH(alternate_bases) = 1
      OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)].alt = '<*>'))
),
allele_counts_per_sample_variant AS (
  SELECT
    c.name,
    c.genotype[ORDINAL(1)] = c.genotype[ORDINAL(2)] AS O_HOM,
    -- Compute AN and AC for all samples in the variant record.
    (SELECT SUM((SELECT COUNT(g) FROM UNNEST(vc.genotype) AS g WHERE g >= 0)) FROM v.call AS vc) AS AN,
    (SELECT SUM((SELECT COUNT(g) FROM UNNEST(vc.genotype) AS g WHERE g = 1)) FROM v.call AS vc) AS AC
  FROM
    filtered_variant_calls AS v, UNNEST(v.call) AS c
),
grouped_values AS (
  SELECT
    name,
    SUM(CAST(O_HOM AS INT64)) AS O_HOM,
    SUM(1.0 - (2.0 * (AC/AN) * (1.0 - (AC/AN)) * (AN / (AN - 1.0)))) AS E_HOM,
    COUNT(name) AS N_SITES
  FROM allele_counts_per_sample_variant
  WHERE
    # Ensure allelic frequency is in the range (0, 1).
    AN > 0 AND AC < AN
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
