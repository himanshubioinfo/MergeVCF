# MergeVCF

This script takes two vcf files at a time and merges them into one output VCF file with:
- common and unique metadata as well as variants are present only once in the output file
- label each variant according to its source file, so that its origin can be traced

Usage (terminal) "python MergeCVF.py file1.vcf file2.vcf"
