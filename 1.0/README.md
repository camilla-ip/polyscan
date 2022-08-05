# polyscan 1.0

## Usage

```
Usage: polyscan.py [9 arguments]

For each sample DNA, scan for the specified polymorphisms and output a table of polymorphisms in each sample.

Tabulate the number of samples, grouped by HCC or non-HCC diagnosis, that contain each
polymorphism, and the p-value for the Fisher's Exact Test for the HCC vs non-HCC counts.

Arguments
---------
            refname : The reference genome name used in input file IDs (e.g., A_X02763)
  refprot_fastafile : Amino-acid sequences with id=REFNAME_PROTNAME [FASTA]
  refnucl_fastafile : Nucleotide genome with id=REFNAME [FASTA]
  polytable_tsvfile : Table of polymorphisms with 1 header row and fields
                      poly_ prefix then [_group, _name, _desc, _level, _gene, _type, _refval _refpos, _altval]
                      Multiple rows for the same (poly_group, poly_name) means the criteria specified in all rows
                      must be present for the (poly_group, poly_name) variant to be considered present[TSV]
sampletable_tsvfile : Table of samples with 1 header row and fields
                      [dataid, samplename, studynumber, diagnosis, HCC] [TSV]
          aaaln_dir : Directory containing DATAID_PROTNAME_AA_align_REFNAME.fasta files
                      Each file contains an alignment of the reference protein (>REFNAME_PROTNAME) and
                      the sample protein (>DATAID_bam_ref_PROTNAME).
                      If a file for a particular PROTNAME does not exist, the protein could not be recovered
                      from the de novo assembly for that sample.
          ntaln_dir : Directory containing DATAID_nt_aligned.fasta files
                      Each file contains an alignment of the reference genome (>REFNAME) and
                      the de novo assembly for the sample (>DATAID)
        outdata_dir : Output directory
     outdata_prefix : Prefix for output files
```

## Installation

```
cd WORKING_DIR
git clone https://github.com/camilla-ip/polyscan.git
```

## Examples

The example contains a sub-directories containing the 'input' files and the expected output files in 'output-correct'.

Executing run_examples.sh will create files in 'output-test' and a .summary.txt file recording any difference between the 'correct' and 'test' output.

```
cd WORKING_DIR/polyscan/1.0/examples
run_examples.sh
```

<h2>INPUT</h2>

1. HBV reference strain (X02763) proteins as amino acids [FASTA]
2. HBV reference strain (X02763) proteins as nucleotides [FASTA]
3. polymorphism-list.tsv : Table of polymorphisms to scan for [TSV]
4. sample_list.tsv : Table of samples and their diagnosis (i.e., HCC vs not-HCC) [TSV]
5. rawdata/protein_alignments/DATAID_GENENAME_AA_align_A_X02763.fasta (264 files) [FASTA]

<h2>OUTPUT</h2>

1. sample_polymorphisms.tsv : Table of sample of polymorphisms [TSV]
2. summary_polymorphisms_with.tsv : Table 2 with 6 columns [TSV]
3. summary_polymorphisms_with_pvalue.tsv : As above, with additional column for pvalue [TSV]

<h2>DATA FORMATS</h2>

In the following explanations, character sequence ' | ' (space pipe space) is used to show where the 'tab' character is placed in tab-separated-variable (TSV) file examples.

<h3>polymorphism-list.tsv</h3>

A TSV file 1 header row, with columns:
```
poly_group | poly_name | poly_desc | poly_level | poly_gene | poly_type | poly_refval | poly_refpos | poly_altval
HBsAg | W172* | aa | Surface_S | SNP | W | 172 | *
```

where
- Name should not contain any spaces or commas
- Level must be { aa, nt }
- Type must be { SNP, DELat, DELany }
  - For DELat polymorphisms, RefVal and AltVal are '-', and RefPos is a range of "N-M" where N are 1-based positions, where the deletion is for positions N to M (inclusive)
  - For DELany polymorphisms, RefVal, RefPos, and AltVal are all '-' because there just has to be some deletions in the sample sequence somewhere in the alignment with the reference protein sequence.
- RefVal must be a valid capitalized 1-letter amino acid code, or the '-' character to denote 'not applicable'
- RefPos must be a valid positive integer (1-based coordinates)
- AltVal must be a valid capitalized 1-letter amino acid code,  '-' to denote 'not applicable', or '\*' for a stop codon
- multiple rows for a Name indicates that all polymormphisms for this site must be present for this sample (i.e., we AND all criteria together for each sample being scanned)

<h3>sample_list.tsv</h3>

A TSV file with 1 header row, with columns:
```
dataid | samplename_clin_final | Study_number | Diagnosis | HCc
```

Example row:
```
wtchgD00007502 | B1-1_11611878_SV | 11611878 (SV) | Hcc | 1
```

where
- HCC must be { 1, 0 }
- AltVal with comma-separated values means that AltVal can have any of those values

<h3>Protein alignments</h3>

Files are named DATAID_GENENAME_AA_align_A_X02763.fasta where GENENAME is one of:
1. Core
2. Pol_RH
3. Pol_RT
4. Pol_TP
5. Pol_spacer
6. PreCore
7. Surface_S
8. Surface_preS1
9. Surface_preS2
10. X

where
- the reference sequence name is A_X02763_GENENAME (e.g., A_X02763_X)
- the DATAID_bam_ref_GENENAME (e.g., wtchgD10001024_bam_ref_X)
- the FASTA file may contain blank lines, which should be ignored

<h3>sample_polymorphisms.tsv</h3>

TSV file with 1 header row, and 4 columns plus one column per polymorphism name:
- dataid
- samplename_clin_final
- Study_number
- HCC [1, 0]
- [poly_name]+

<h3>summary_polymorphisms_with_pvalue.tsv</h3>

TSV file with 1 header row, and 7 columns:
- poly_group
- poly_name
- freq_in_HCC
- pct_in_HCC
- freq_in_nonHCC
- pct_in_nonHCC
- p-value

## Run-time warnings

polyscan.py prints warnings when the polymorphism in the reference sequence is not the same as the value in the 'poly_refval' listed in the polymorphism table.

For example, in example 1, the reference sequence has amino acid 'T' as position 13 of the Pre-Core protein-coding gene. The polymorphism-list-example1.tsv file states that the 'poly_refval' is amino-acid 'S' and the 'poly_altval' is amino-acid 'T'.

The program prints the following warning:

```
Warning: Polytable SNP refval disagrees with refprot FASTA (Pre-Core, S13T): table=S fasta=T
```

However, this is nothing to worry about. The warning is printed to alert the user to possible data-entry issues. But there is no reason why the reference genome does not contain the 'mutation' rather than the 'wild-type' value at any given position.

## Licence

polyscan is released under the [MIT Licence](https://github.com/camilla-ip/software/blob/main/LICENSE).
