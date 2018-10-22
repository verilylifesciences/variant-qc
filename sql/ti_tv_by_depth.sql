#standardSQL
--
-- Transition/Transversion Ratio by Depth of Coverage.
--
WITH filtered_snp_calls AS (
  SELECT
    c.name,
    c.DP AS depth,
    CONCAT(reference_bases, '->', alternate_bases[ORDINAL(1)].alt) AS mutation
  FROM
    `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.call) AS c
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND alternate_bases[ORDINAL(1)].alt IN ('A','C','G','T')
    AND (ARRAY_LENGTH(alternate_bases) = 1
      OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)].alt = '<*>'))
    # Skip homozygous reference calls and no-calls.
    AND EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g > 0)
    AND NOT EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g < 0)
    # Include only high quality calls.
    AND {{ HIGH_QUALITY_CALLS_FILTER }}
),

mutation_type_counts AS (
  SELECT
    name,
    depth,
    SUM(CAST(mutation IN ('A->G', 'G->A', 'C->T', 'T->C') AS INT64)) AS transitions,
    SUM(CAST(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                          'A->T', 'T->A', 'C->G', 'G->C') AS INT64)) AS transversions,
    COUNT(mutation) AS num_variants_in_group
  FROM filtered_snp_calls
  WHERE
    depth IS NOT NULL
    AND depth > 0
  GROUP BY
    name,
    depth
)

SELECT
  name,
  depth,
  transitions,
  transversions,
  transitions/transversions AS titv,
  num_variants_in_group
FROM mutation_type_counts
WHERE
  transversions > 0 AND transitions > 0
ORDER BY
  name,
  depth
