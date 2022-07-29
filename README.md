[![PyPI version](https://badge.fury.io/py/ScanCE.svg)](https://badge.fury.io/py/ScanCE) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-360/)

# ScanCE

A computational workflow for cryptic exon identification. Cryptic exon (CE) is one of the non-canonical splicing events, which is defined as a novel exon within an annotated intron.


Prerequisites
----------------
`ScanCE` runs under Python 3.7+ and is available through python3-pip.

# Installation

The recommended way to install ScanCE is using pip:

```console
$ pip install ScanCE
```

# Usage

ScanCE has two modes, `Scan_ce_strict` and `Scan_ce_loose`. 
For `Scan_ce_strict` model, we want to find novel junctions that satisfy the conditions as shown in the following image. When extracting novel spliced junctions, we only extract one end at the intron boundary and one end inside the intron.
<img  width="650" src="https://github.com/ylab-hi/ScanCE/blob/master/ScanCE_strict.png">

For `Scan_ce_loose` model, in addition to the cryptic exons found in strict mode, the cryptic exons shown below will also be found. When extracting novel spliced junctions, we only guarantee that one end of them falls in the intron, and the other end may be inside the exon.
<img  width="650" src="https://github.com/ylab-hi/ScanCE/blob/master/ScanCE_loose.png">

There are two reqiured inputs:
(1) a BAM alignment file of short-reads.
(2) sorted and bgzip'd gene annotation files. To ensure as much as possible that the cryptic exons found are novel unannotated exons, we recommend using multiple annotation files. For example, for human RNA sequencing, we used annotation files from genecode and NCBI. Examples are as follows:
```console
$ ScanCE Scan_ce_strict -i LNCaP_CCLE_chr21_test.bam -r1 gencode.v38.annotation_sorted_chr21.gtf.gz -r2 GRCh38_latest_genomic.sorted_chr21.gff.gz
```







