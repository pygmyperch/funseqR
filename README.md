# funseqR

<br/> Functional annotation of sequence data for non-model species in R <br/>

# How to Install

The preferred way to install this package is using devtools:

First install devtools

``` r
install.packages("devtools")
```

Then install funseqR

``` r
devtools::install_github("pygmyperch/funseqR")
```

A quick overview of some of the key functions:

-   `read_vcf`: reads a VCF (Variant Call Format) file into R using the vcfR package

-   `vcf2bed`: generates a data frame in BED format and (optional) .bed file from a vcfR object

-   `extract_flanking_sequences`: extracts flanking sequences for SNPs from a reference genome

-   `dwnld_ncbi_db`: generates a script for downloading the specified database using one of several methods

-   `perform_blast`: performs either a blastn or blastx search using BLAST+

-   `process_blast_results`: Process BLAST results and retrieve UniProt annotations including GO terms and KEGG references

-   `summarize_annotations`:

-   `xxx`:

-   `xxx`:

<br/>

Feel free to use or modify these scripts for your own needs. The code is intended to be functional and straightforward, though improvements and optimizations are always welcome so please feel free to submit pull requests.
