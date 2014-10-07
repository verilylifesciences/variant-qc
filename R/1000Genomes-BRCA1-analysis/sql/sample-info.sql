# Get ethnicity, gender and family relationship information for the samples.
SELECT
  Sample,
  Gender,
  Family_ID,
  Relationship,
  Population,
  Population_Description,
  Super_Population,
  Super_Population_Description
FROM
  [genomics-public-data:1000_genomes.sample_info]
WHERE
  In_Phase1_Integrated_Variant_Set = TRUE
ORDER BY
  Sample
