#standardSQL
--
-- Compute private variants counts for each sample.
--
WITH filtered_called_alleles AS (
  SELECT
    reference_name,
    start_position,
    reference_bases AS ref,
    a.alt,
    c.name
  FROM
    `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v,
    UNNEST(v.call) AS c,
    UNNEST(v.alternate_bases) AS a WITH OFFSET alt_offset
  WHERE
    EXISTS (SELECT g FROM UNNEST(c.genotype) AS g WHERE g = alt_offset+1)
    # Include only high quality calls.
    AND {{ HIGH_QUALITY_CALLS_FILTER }}
),

grouped_alleles AS (
  SELECT
    reference_name,
    start_position,
    ref,
    alt,
    STRING_AGG(name) AS name,
    COUNT(name) AS num_samples_with_variant
  FROM filtered_called_alleles
  GROUP BY
    reference_name,
    start_position,
    ref,
    alt
)

SELECT
  name,
  COUNT(name) AS private_variant_count
FROM grouped_alleles
WHERE
  num_samples_with_variant = 1
GROUP BY
  name
ORDER BY
  private_variant_count DESC
