This directory contains the inputs and outputs of PLINK check-sex analysis.

#### Inputs:
  * [chrX-all-samples-merged.bim](chrX-all-samples-merged.bim) and [chrX-all-samples-merged.bed](chrX-all-samples-merged.bed):
    * These files were transformed from a VCF file that was created using
      the query included below. The primary reason for using a customized input
      file was to ensure PLINK used the same input as BigQuery. For instance,
      the 'expanded table' used in the query only contains non-variant segments
      for SNPs (i.e. indels matching the reference are not included). Thus,
      indels needed to be excluded in the provided inputs.
      Also, to obtain accurate results in PLINK, the input files must
      contain merged info about all samples for each SNP. Thus, using
      separate VCF files (each containing info for a specific individual)
      was not feasible.

  * [chrX-all-samples-merged.fam](chrX-all-samples-merged.fam):
    * This file contains pedigree information for all 17 individuals.

#### Outputs:
  * [sex-check-with-fam](sex-check-with-fam.sexcheck):
    * PLINK output when using actual pedigree info.

  * [sex-check-no-fam](sex-check-no-fam.sexcheck):
    * PLINK output when using a dummy .fam file where all individuals have
      unknown parents. This was used to obtain comparable results to
      BiqQuery since it also ignores pedigree info.

#### Query used to generate the input VCF file:
    SELECT reference_name, start, reference_bases, alternate_bases,
    GROUP_CONCAT(calls, ',') AS all_calls
    FROM (
    SELECT
      reference_name,
      start,
      reference_bases,
      alternate_bases,
      call.call_set_name + ':' + string(IFNULL(first_allele, -1)) + '/' + string(IFNULL(second_allele, -1)) as calls
    FROM (
        SELECT
          reference_name,
          start,
          reference_bases,
          GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
          call.call_set_name,
          NTH(1, call.genotype) WITHIN call AS first_allele,
          NTH(2, call.genotype) WITHIN call AS second_allele,
        FROM
          [google.com:biggene:platinum_genomes.expanded_variants]
        WHERE
          reference_name = 'chrX' AND
          (start NOT BETWEEN 59999 AND 2699519) AND
          (start NOT BETWEEN 154931042 AND 155260559)
          AND reference_bases IN ('A', 'C', 'T', 'G')
        HAVING
           alternate_bases IN ('A', 'C', 'T', 'G')
    ))
    GROUP BY reference_name, start, reference_bases, alternate_bases
