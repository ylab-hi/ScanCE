# ScanCE

**Unified *de novo* detection and quantification of cryptic exons from short-read, long-read, and single-cell RNA-seq**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

Cryptic exons (CEs) are unannotated exonic sequences within annotated introns that become aberrantly included in mature mRNA when repressive RNA-binding proteins (e.g., TDP-43) are depleted. CEs are implicated in neurodegeneration (ALS/FTD) and cancer.

**ScanCE** is the first purpose-built computational tool for *de novo* CE detection. It operates on a single BAM/CRAM file — no case-control or multi-sample design is required.

## Features

- **Three sequencing modes**: short-read (Illumina), long-read (ONT/PacBio), and single-cell (scISO-seq/MAS-seq)
- **Dual-annotation cross-referencing**: validates CE candidates against both GENCODE and NCBI RefSeq to reduce false positives
- **Multi-exon CE detection**: identifies complex CEs comprising multiple novel exonic segments within a single intron
- **Streaming architecture**: constant ~46 MB peak memory regardless of sequencing depth
- **PSI quantification**: percent-spliced-in calculation for all modes; phased single-molecule evidence for long-read data

## Installation

### From PyPI (recommended)

```bash
pip install ScanCE
```

### From source

```bash
git clone https://github.com/ylab-hi/ScanCE.git
cd ScanCE
pip install .
```

### Dependencies

- Python >= 3.8
- [pysam](https://github.com/pysam-developers/pysam) >= 0.19.0
- [gffutils](https://github.com/daler/gffutils) >= 0.12

## Reference file preparation

ScanCE requires three reference files. Prepare them once and point to them via `config.ini`.

### 1. Reference genome FASTA

Download from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) or [Ensembl](https://www.ensembl.org/):

```bash
# Example: hg38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

### 2. GENCODE GTF annotation (sorted + indexed)

```bash
# Download GENCODE v38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip gencode.v38.annotation.gtf.gz

# Sort, compress, and index
(grep "^#" gencode.v38.annotation.gtf; \
 grep -v "^#" gencode.v38.annotation.gtf | sort -k1,1 -k4,4n) \
 | bgzip > gencode.v38.annotation_sorted.gtf.gz
tabix -p gff gencode.v38.annotation_sorted.gtf.gz

# Build gffutils database
python -c "
import gffutils
gffutils.create_db(
    'gencode.v38.annotation_sorted.gtf.gz',
    'gencode.v38.annotation_sorted.gtf.gz.db',
    force=True, merge_strategy='merge')
"
```

### 3. NCBI RefSeq GFF3 annotation (sorted + indexed)

```bash
# Download NCBI RefSeq
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.gff.gz

# Sort, compress, and index
(grep "^#" GCF_000001405.40_GRCh38.p14_genomic.gff; \
 grep -v "^#" GCF_000001405.40_GRCh38.p14_genomic.gff | sort -k1,1 -k4,4n) \
 | bgzip > GRCh38_latest_genomic.sorted.gff.gz
tabix -p gff GRCh38_latest_genomic.sorted.gff.gz
```

### 4. Configure paths

```bash
cp config.ini.example config.ini
# Edit config.ini with your actual paths
```

## Quick start

```bash
# Short-read mode (Illumina)
scanCE -i sample.bam --mode sr --ce_type single -o sample.sr.single.ce

# Long-read mode (ONT/PacBio)
scanCE -i sample.bam --mode lr --ce_type single -o sample.lr.single.ce

# Multi-exon CE detection
scanCE -i sample.bam --mode lr --ce_type multi -o sample.lr.multi.ce

# Single-cell mode (PacBio CCS/MAS-seq)
scanCE -i cell.bam --mode sc --ce_type single --primary_only -o cell.sc.single.ce
```

## Usage

```
usage: scanCE [-h] -i INPUT --mode {sr,lr,sc} [--ce_type {single,multi}]
              [-o OUTPUT] [-m MAPQ] [-a AO] [-s {no,fr-firststrand,fr-secondstrand}]
              [-p PSI] [--min_junction_reads MIN_JUNCTION_READS]
              [--stringency {loose,strict}] [--cell_id CELL_ID]
              [--primary_only] [-v]

Required arguments:
  -i, --input           Input BAM/CRAM file (must be indexed)
  --mode {sr,lr,sc}     sr=Illumina bulk, lr=ONT/PacBio bulk, sc=single-cell long-read

Optional arguments:
  --ce_type {single,multi}
                        single=single-exon CE, multi=multi-exon CE (default: multi)
  -o, --output          Output filename (default: <sample>.<mode>.<ce_type>.ce)
  -m, --mapq            Min MAPQ (default: sr=50, lr/sc=0)
  -a, --ao              Min reads for each outer junction (default: 1)
  -s, --stranded        SR strand mode: no, fr-firststrand, fr-secondstrand (default: no)
  -p, --psi             Min PSI threshold (default: 0.0)
  --min_junction_reads  Min internal junction support for multi mode (default: 2)
  --stringency {loose,strict}
                        loose=endpoint anywhere within exon; strict=exact boundary (default: loose)
  --primary_only        Use primary alignments only (recommended for PacBio CCS)
  --cell_id             Cell ID for sc mode (auto-inferred from BAM filename if not given)
```

## Output format

ScanCE outputs a tab-delimited `.ce` file. Column definitions vary by mode and CE type:

### Short-read, single-exon (`sr` + `single`)

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| D | Donor site (annotated intron start) |
| A | Acceptor site (annotated intron end) |
| ce_start | Cryptic exon start coordinate |
| ce_end | Cryptic exon end coordinate |
| ao1 | Junction reads at 5' CE boundary |
| ao2 | Junction reads at 3' CE boundary |
| ao3 | Canonical skip junction reads |
| a_count | Internal spanning reads |
| PSI | Percent spliced-in |
| strand | Strand (+/-) |
| gene_id | Ensembl gene ID |
| gene_name | Gene symbol |

### Long-read, single-exon (`lr` + `single`)

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| D | Donor site (annotated intron start) |
| A | Acceptor site (annotated intron end) |
| ce_start | Cryptic exon start coordinate |
| ce_end | Cryptic exon end coordinate |
| ao | Full-spanning reads (crossing both CE junctions) |
| ao1 | Junction reads at 5' CE boundary |
| ao2 | Junction reads at 3' CE boundary |
| ao_canon | Canonical skip junction reads |
| PSI | Percent spliced-in: ao / (ao + ao_canon) |
| strand | Strand (+/-) |
| gene_id | Ensembl gene ID |
| gene_name | Gene symbol |

### Multi-exon mode (`multi`)

Multi-exon output includes additional columns `ce_start_1`, `ce_end_1`, `ce_start_2`, `ce_end_2` for the two sub-exon coordinates, replacing `ce_start`/`ce_end`.

### Single-cell mode (`sc`)

Single-cell output prepends a `cell_id` column to the corresponding `lr` format.

## Algorithm

ScanCE executes a five-step pipeline:

1. **Junction extraction**: Parses CIGAR strings for splice junctions. In LR/SC mode, deletions >= 30 bp are reclassified as introns with adjacent-deletion boundary correction.
2. **Strand assignment**: XS tag (SR) or ts tag (LR/SC) for strand inference. Dual-strand strategy for unstranded libraries.
3. **Dual-annotation cross-referencing**: Matches junctions against both GENCODE and NCBI RefSeq. Only junctions absent from both databases are retained as CE candidates.
4. **CE assembly and typology**: Pairs donor and acceptor junctions anchored in the same annotated intron. Classifies into single-exon or multi-exon CEs.
5. **Quantification**: PSI computation with mode-specific evidence (junction counts for SR; phased single-molecule reads for LR/SC).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
