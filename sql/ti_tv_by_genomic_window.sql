#standardSQL
--
-- Compute the Ti/Tv ratio for variants within genomic windows.
--
WITH filtered_snp_calls AS (
  SELECT
    reference_name,
    CAST(FLOOR(start_position / {{ WINDOW_SIZE }}) AS INT64) AS genomic_window,
    CONCAT(reference_bases, '->', alternate_bases[ORDINAL(1)].alt) AS mutation
  FROM
    `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}`
  WHERE
    # Only include biallelic snps with at least one passing variant call.
    reference_bases IN ('A','C','G','T')
    AND alternate_bases[ORDINAL(1)].alt IN ('A','C','G','T')
    AND (ARRAY_LENGTH(alternate_bases) = 1
      OR (ARRAY_LENGTH(alternate_bases) = 2 AND alternate_bases[ORDINAL(2)].alt = '<*>'))
    AND EXISTS (SELECT g FROM UNNEST(call) AS c, UNNEST(c.genotype) AS g WHERE g > 0)
    # Include variants with at least one high quality call.
    AND {{ HIGH_QUALITY_VARIANTS_FILTER }}
),

mutation_type_counts AS (
  SELECT
    reference_name,
    genomic_window,
    SUM(CAST(mutation IN ('A->G', 'G->A', 'C->T', 'T->C') AS INT64)) AS transitions,
    SUM(CAST(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                          'A->T', 'T->A', 'C->G', 'G->C') AS INT64)) AS transversions,
    COUNT(mutation) AS num_variants_in_group
  FROM filtered_snp_calls
  GROUP BY
    reference_name,
    genomic_window
)

SELECT
  reference_name,
  genomic_window * {{ WINDOW_SIZE }} AS window_start,
  transitions,
  transversions,
  transitions/transversions AS titv,
  num_variants_in_group
FROM mutation_type_counts
WHERE
  transversions > 0
ORDER BY
  reference_name,
  window_start
