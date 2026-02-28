[![PyPI version](https://badge.fury.io/py/ScanCE.svg)](https://badge.fury.io/py/ScanCE) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-360/)

# ScanCE

A computational workflow for cryptic exon (CE) identification. A cryptic exon is a novel exon located within an annotated intron, detected from RNA-seq data.

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

## Overview

```
ScanCE <command> [options]

Commands:
  scan_ce              Unified CE detection for SR / LR / SC  (v2, recommended)
  Scan_ce_loose        Short-read, loose stringency  (v1 legacy)
  Scan_ce_strict       Short-read, strict stringency (v1 legacy)
  Scan_ce_lr_loose     Long-read,  loose stringency  (v1 legacy)
  Scan_ce_lr_strict    Long-read,  strict stringency (v1 legacy)
```

---

## `scan_ce` — Unified Command (v2, Recommended)

`scan_ce` combines all sequencing types and stringency modes into a single tool.

### Arguments

| Argument | Choices / Default | Description |
|----------|-------------------|-------------|
| `-i` | — | Input BAM file (must be indexed) |
| `--mode` | `sr` / `lr` / `sc` | Sequencing type: Illumina bulk / ONT·PacBio bulk / single-cell long-read |
| `--ce_type` | `single` / `multi` (default: `multi`) | Single-exon CE (no internal junction) or multi-exon CE (has internal junction) |
| `--stringency` | `loose` / `strict` (default: `loose`) | Junction matching stringency (see below) |
| `-m` | int (default: 50 for SR, 0 for LR/SC) | Min MAPQ filter |
| `-a` | int (default: 1) | Min reads for each outer junction (ao1, ao2) |
| `-o` | — | Output filename |
| `-s` | `no` / `fr-firststrand` / `fr-secondstrand` | Strand library type for SR (default: `no`, uses XS tag) |
| `-p` | float (default: 0.0) | Min PSI threshold (SR mode only) |
| `--min_junction_reads` | int (default: 2) | Min reads for internal junction in `multi` mode |
| `--primary_only` | flag | Use primary alignments only (recommended for PacBio CCS) |
| `--cell_id` | string | Cell ID for SC mode (auto-inferred from BAM filename if not set) |

### Stringency modes

The `--stringency` parameter controls how cryptic junction endpoints are matched to annotated exons:

| Mode | Donor match | Acceptor match |
|------|-------------|----------------|
| `loose` | Junction start falls **anywhere within** an annotated exon | Junction end falls **anywhere within** an annotated exon |
| `strict` | Junction start falls **exactly at the exon end** boundary | Junction end falls **exactly at the exon start** boundary |

**loose** (default) — finds more CEs, including those where the cryptic splice site is internal to a known exon:

<img width="650" src="https://github.com/ylab-hi/ScanCE/blob/master/ScanCE_loose.png">

**strict** — only finds CEs where the cryptic junction shares an exact boundary with a known exon:

<img width="650" src="https://github.com/ylab-hi/ScanCE/blob/master/ScanCE_strict.png">

### CE type

| Type | Description |
|------|-------------|
| `single` | No junctions inside the CE region — the cryptic exon is a single novel exon |
| `multi` | At least one junction inside the CE region — the CE spans multiple novel exons |

### Usage examples

**Short-read bulk RNA-seq (Illumina):**
```console
$ ScanCE scan_ce --mode sr --ce_type single --stringency loose -i sample.bam -m 50 -o sample.sr.single.loose.ce
$ ScanCE scan_ce --mode sr --ce_type single --stringency strict -i sample.bam -m 50 -o sample.sr.single.strict.ce
$ ScanCE scan_ce --mode sr --ce_type multi  --stringency loose -i sample.bam -m 50 -o sample.sr.multi.loose.ce
```

**Long-read bulk RNA-seq (ONT / PacBio):**
```console
$ ScanCE scan_ce --mode lr --ce_type single --stringency loose -i sample.bam -m 0 -o sample.lr.single.loose.ce
$ ScanCE scan_ce --mode lr --ce_type multi  --stringency loose -i sample.bam -m 0 -o sample.lr.multi.loose.ce
```

**Single-cell long-read (PacBio CCS / MAS-seq):**
```console
$ ScanCE scan_ce --mode sc --ce_type single --stringency loose -i cell.bam -m 0 --primary_only -o cell.sc.single.ce
```

### Configuration

`scan_ce` reads the GENCODE annotation path from `config.ini` placed in the same directory as `ScanCE_v3.py`:

```ini
[sorted GENCODE annotation]
annotation = /path/to/gencode.vXX.annotation_sorted.gtf.gz
```

The annotation file must be bgzip-compressed and tabix-indexed (`.tbi`), and a `gffutils` database (`.db`) must exist alongside it.

### Output columns

| Mode / CE type | Columns |
|----------------|---------|
| SR / single | `chrom  D  A  ce_start  ce_end  ao1  ao2  ao3  a_count  PSI  strand  gene_id  gene_name` |
| SR / multi | `chrom  D  A  ce_start_1  ce_end_1  ce_start_2  ce_end_2  ao1  ao2  ao3  strand  gene_id  gene_name` |
| LR / single | `chrom  D  A  ce_start  ce_end  ao  ao1  ao2  strand  gene_id  gene_name` |
| LR / multi | `chrom  D  A  ce_start_1  ce_end_1  ce_start_2  ce_end_2  ao  ao1  ao2  ao3  strand  gene_id  gene_name` |
| SC / single | `cell_id  chrom  D  A  ce_start  ce_end  ao  ao1  ao2  strand  gene_id  gene_name` |
| SC / multi | `cell_id  chrom  D  A  ce_start_1  ce_end_1  ce_start_2  ce_end_2  ao  ao1  ao2  ao3  strand  gene_id  gene_name` |

**Column definitions:**

| Column | Description |
|--------|-------------|
| `D` / `A` | Donor / acceptor positions of the outer canonical (skipping) intron |
| `ce_start` / `ce_end` | Start and end of the cryptic exon |
| `ao` | Spanning reads: same read covers both outer junctions (LR/SC only) |
| `ao1` / `ao2` | Reads supporting each outer cryptic junction |
| `ao3` | Reads supporting the canonical (exon-skipping) junction |
| `a_count` | Reads fully contained within the CE region (SR only) |
| `PSI` | Percent Spliced In = (ao1 + ao2 + a_count) / (ao1 + ao2 + ao3 + a_count) (SR only) |

---

## v1 Legacy Commands

The original v1 scripts are still available for backward compatibility. They require annotation files passed directly via `-r1` and `-r2`.

```console
$ ScanCE Scan_ce_strict -i sample.bam -r1 gencode.v38.annotation_sorted.gtf.gz -r2 GRCh38_latest_genomic.sorted.gff.gz
$ ScanCE Scan_ce_loose  -i sample.bam -r1 gencode.v38.annotation_sorted.gtf.gz -r2 GRCh38_latest_genomic.sorted.gff.gz
```
