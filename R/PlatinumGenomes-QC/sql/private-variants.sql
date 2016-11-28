# Compute private variants counts for each sample.
WITH filtered_called_alleles AS (
  SELECT
    reference_name,
    start,
    reference_bases AS ref,
    alt,
    call.call_set_name,
    (SELECT COUNT(CAST(gt = alt_offset+1 AS INT64)) FROM call.genotype gt) AS allele_cnt
  FROM
    `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v, v.call call, v.alternate_bases alt WITH OFFSET alt_offset
  WHERE
    # Skip homozygous reference calls, no-calls, and non-passing variants.
    EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt > 0)
    AND NOT EXISTS (SELECT gt FROM UNNEST(call.genotype) gt WHERE gt < 0)
    AND NOT EXISTS (SELECT ft FROM UNNEST(call.FILTER) ft WHERE ft NOT IN ('PASS', '.'))
),

grouped_alleles AS (
  SELECT
    reference_name,
    start,
    ref,
    alt,
    STRING_AGG(call_set_name) AS call_set_name,
    COUNT(call_set_name) AS num_samples_with_variant
  FROM filtered_called_alleles
  GROUP BY
    reference_name,
    start,
    ref,
    alt
)

SELECT
  call_set_name,
  COUNT(call_set_name) AS private_variant_count
FROM grouped_alleles
WHERE
  num_samples_with_variant = 1
GROUP BY
  call_set_name
ORDER BY
  private_variant_count DESC
