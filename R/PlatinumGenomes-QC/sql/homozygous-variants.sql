# Compute the expected and observed homozygosity rate for each individual.
# FROM DOCUMENTATION IN VCFTOOLS
# P(Homo) = F + (1-F)P(Homo by chance)
# P(Homo by chance) = p^2+q^2 for a biallelic locus.
# For an individual with N genotyped loci, we
#   1. count the total observed number of loci which are homozygous (O),
#   2. calculate the total expected number of loci homozygous by chance (E)
# Then, using the method of moments, we have
#    O = NF + (1-F)E
# Which rearranges to give
#    F = (O-E)/(N-E)

SELECT
  call_call_set_name,
  O_HOM,
  ROUND(E_HOM, 2) AS E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM (
  SELECT
    call_call_set_name,
    SUM(O_HOM) AS O_HOM,
    SUM(E_HOM) AS E_HOM,
    COUNT(call_call_set_name) AS N_SITES,
  FROM js(
  // Input Table
  (SELECT
      reference_bases,
      alt.alternate_bases,
      AN,
      alt.AF,
      refMatchCallsets,
      call.call_set_name,
      call.genotype,
    FROM
      [_MULTISAMPLE_VARIANT_TABLE_]
    # Optionally add a clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
   ),
  // Input Columns
  reference_bases, alt.alternate_bases, AN, alt.AF, refMatchCallsets, call.call_set_name, call.genotype,
  // Output Schema
  "[{name: 'call_call_set_name', type: 'string'},
    {name: 'O_HOM', type: 'boolean'},
    {name: 'E_HOM', type: 'float'}]",
  // Function
  "function(row, emit) {
    // Only operate on bi-allelic variants.
    if (1 != row.alt.length) return;

    // Skip records where all samples have the same allele.
    if (row.alt[0].AF == 0 || row.alt[0].AF == 1) return;

    // Only operate on SNPs
    if (!['A', 'C', 'G', 'T'].includes(row.reference_bases)) return;
    if (!['A', 'C', 'G', 'T'].includes(row.alt[0].alternate_bases)) return;

    // Compute the expected homozygosity at this site.
    var expectedHom =
          1.0 - (2.0 * row.alt[0].AF * (1.0 - row.alt[0].AF) * (row.AN / (row.AN - 1.0)));

    // Emit the observed homozygosity for each variant call.
    for(var i = 0; i < row.call.length; i++) {
      // Skip genotypes with anything other than two values.
      if (2 != row.call[i].genotype.length) continue;
      // Skip genotypes with no-calls.
      if (row.call[i].genotype.includes(-1)) continue;

      emit({call_call_set_name: row.call[i].call_set_name,
            O_HOM: (row.call[i].genotype[0] == 1 && row.call[i].genotype[1] == 1),
            E_HOM: expectedHom});
    }

    // Emit the observed homozygosity for each reference-match call.
    for(var i = 0; i < row.refMatchCallsets.length; i++) {
      emit({call_call_set_name: row.refMatchCallsets[i],
            O_HOM: true,
            E_HOM: expectedHom});
    }
  }")
  GROUP BY
    call_call_set_name
  )
ORDER BY
  call_call_set_name
