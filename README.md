[![PyPI version](https://badge.fury.io/py/ScanCE.svg)](https://badge.fury.io/py/ScanCE) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-360/)

# ScanCE

A computational workflow for cryptic exon (CE) identification. Cryptic exon is defined as a novel exon within an annotated intron, identified from RNA-seq data.

## Prerequisites

`ScanCE` runs under Python 3.7+ and requires:

```
pysam
gffutils
```

## Installation

```console
$ pip install ScanCE
```

## Usage

ScanCE v2.0 introduces a unified command `scan_ce` that supports short-read (SR), long-read (LR), and single-cell long-read (SC) data in a single tool.

```
ScanCE <command> [options]

Commands:
  scan_ce              Unified CE detection for SR/LR/SC (v2, recommended)
  Scan_ce_loose        Short-read CE detection (loose, v1)
  Scan_ce_strict       Short-read CE detection (strict, v1)
  Scan_ce_lr_loose     Long-read CE detection (loose, v1)
  Scan_ce_lr_strict    Long-read CE detection (strict, v1)
```

---

## `scan_ce` — Unified Mode (v2, Recommended)

### Key arguments

| Argument | Description |
|----------|-------------|
| `-i` | Input BAM file (must be indexed) |
| `--mode` | `sr` (Illumina bulk), `lr` (ONT/PacBio bulk), `sc` (single-cell long-read) |
| `--ce_type` | `single` (single-exon CE), `multi` (multi-exon CE) |
| `-m` | Min MAPQ (default: 50 for SR, 0 for LR/SC) |
| `-o` | Output filename |
| `--primary_only` | Use primary alignments only (recommended for PacBio CCS) |
| `--min_junction_reads` | Min reads for internal junction in multi mode (default: 2) |
| `--cell_id` | Cell ID for SC mode (auto-inferred from BAM filename if not set) |

### Short-read bulk RNA-seq (Illumina)

```console
$ ScanCE scan_ce --mode sr --ce_type single -i sample.bam -m 50 -o sample.sr.single.ce
$ ScanCE scan_ce --mode sr --ce_type multi  -i sample.bam -m 50 -o sample.sr.multi.ce
```

### Long-read bulk RNA-seq (ONT / PacBio)

```console
$ ScanCE scan_ce --mode lr --ce_type single -i sample.bam -m 0 -o sample.lr.single.ce
$ ScanCE scan_ce --mode lr --ce_type multi  -i sample.bam -m 0 -o sample.lr.multi.ce
```

### Single-cell long-read (PacBio CCS / MAS-seq)

```console
$ ScanCE scan_ce --mode sc --ce_type single -i cell.bam -m 0 --primary_only -o cell.sc.single.ce
```

### Configuration

`scan_ce` reads annotation paths from `config.ini` in the same directory:

```ini
[sorted GENCODE annotation]
annotation = /path/to/gencode.annotation_sorted.gtf.gz
```

### Output columns

**SR / single-exon:**
`chrom  D  A  ce_start  ce_end  ao1  ao2  ao3  a_count  PSI  strand  gene_id  gene_name`

**LR / single-exon:**
`chrom  D  A  ce_start  ce_end  ao  ao1  ao2  strand  gene_id  gene_name`

**LR / multi-exon:**
`chrom  D  A  ce_start_1  ce_end_1  ce_start_2  ce_end_2  ao  ao1  ao2  ao3  strand  gene_id  gene_name`

**SC / single-exon:**
`cell_id  chrom  D  A  ce_start  ce_end  ao  ao1  ao2  strand  gene_id  gene_name`

> `D`/`A`: donor/acceptor positions of the outer canonical intron
> `ao`: spanning reads (same read covers both outer junctions, LR/SC only)
> `ao1`/`ao2`: reads supporting each outer cryptic junction
> `ao3`: reads supporting the canonical (skipping) junction
> `PSI`: Percent Spliced In (SR only)

---

## v1 Legacy Modes

For `Scan_ce_strict`: novel junctions where one end is at the intron boundary and the other is inside the intron.

<img width="650" src="https://github.com/ylab-hi/ScanCE/blob/master/ScanCE_strict.png">

For `Scan_ce_loose`: additionally includes CEs where one end falls within an exon.

<img width="650" src="https://github.com/ylab-hi/ScanCE/blob/master/ScanCE_loose.png">

```console
$ ScanCE Scan_ce_strict -i sample.bam -r1 gencode.v38.annotation_sorted.gtf.gz -r2 GRCh38_latest_genomic.sorted.gff.gz
$ ScanCE Scan_ce_loose  -i sample.bam -r1 gencode.v38.annotation_sorted.gtf.gz -r2 GRCh38_latest_genomic.sorted.gff.gz
```
