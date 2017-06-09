#standardSQL
--
-- Compute the Ti/Tv ratio for variants binned by alternate allele count.
--
WITH filtered_snp_calls AS (
  SELECT
    ( -- Compute the allele count for this site of variation.
      SELECT COUNTIF(gt = 1)
      FROM UNNEST(v.call) AS call, UNNEST(call.genotype) AS gt
      WHERE
        # Skip homozygous reference calls, no-calls, and non-passing variants.
        EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
        AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
        AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
    ) AS AC,
    CONCAT(reference_bases, '->', alternate_bases[ORDINAL(1)]) AS mutation
  FROM
    `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v
  WHERE
    # Only include biallelic snps.
    reference_bases IN ('A','C','G','T')
    AND ARRAY_LENGTH(alternate_bases) = 1
    AND alternate_bases[ORDINAL(1)] IN ('A','C','G','T')
),

mutation_type_counts AS (
  SELECT
    AC,
    SUM(CAST(mutation IN ('A->G', 'G->A', 'C->T', 'T->C') AS INT64)) AS transitions,
    SUM(CAST(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                          'A->T', 'T->A', 'C->G', 'G->C') AS INT64)) AS transversions,
    COUNT(mutation) AS num_variants_in_group
  FROM filtered_snp_calls
  GROUP BY
    AC
)

SELECT
  AC,
  transitions,
  transversions,
  transitions/transversions AS titv,
  num_variants_in_group
FROM mutation_type_counts
WHERE
  transversions > 0 AND AC > 0
ORDER BY
  AC DESC
