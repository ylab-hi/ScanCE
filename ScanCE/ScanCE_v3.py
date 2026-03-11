#!/usr/bin/env python3
"""
ScanCE v3.0 — Unified Cryptic Exon Detector

Extended from Scan_ce_hpc_longread_loose_multi_v2.py.

New features:
  --mode   sr  : Short-read (Illumina bulk), uses XS tag for strand inference, computes PSI
           lr  : Long-read (ONT/PacBio bulk), computes ao_canon/PSI
           sc  : Single-cell long-read (PacBio CCS/MAS-seq), LR logic + cell_id + primary_only

  --ce_type single : Single-exon CE (no internal junction in CE region)
            multi  : Multi-exon CE (has internal junction, original v2 behavior)
"""

import sys
import os
import argparse
import pysam
from collections import Counter, defaultdict
import gffutils
from configparser import ConfigParser

__version__ = 'v3.0'


# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="ScanCE v3.0: Unified Cryptic Exon Detector",
        epilog="Based on Scan_ce_hpc_longread_loose_multi_v2.py. "
               "Modes: sr (short-read), lr (long-read), sc (single-cell long-read).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input BAM/CRAM file (must be indexed).')
    parser.add_argument('--mode', required=True, choices=['sr', 'lr', 'sc'],
                        help='sr=Illumina bulk, lr=ONT/PacBio bulk, sc=single-cell long-read.')
    parser.add_argument('--ce_type', choices=['single', 'multi'], default='multi',
                        help='single=single-exon CE (no internal junction), '
                             'multi=multi-exon CE (has internal junction, original v2 behavior).')
    parser.add_argument('-o', '--output', default=None,
                        help='Output filename. Default: <sample>.<mode>.<ce_type>.ce')
    parser.add_argument('-m', '--mapq', type=int, default=None,
                        help='Min MAPQ. Auto-default: sr=50, lr/sc=0.')
    parser.add_argument('-a', '--ao', type=int, default=1,
                        help='Min reads for each outer junction (ao1, ao2).')
    parser.add_argument('-s', '--stranded',
                        choices=['no', 'fr-firststrand', 'fr-secondstrand'], default='no',
                        help='SR strand mode (uses XS tag when "no"). LR/SC always use ts tag.')
    parser.add_argument('-p', '--psi', type=float, default=0.0,
                        help='Min PSI threshold (all modes).')
    parser.add_argument('--min_junction_reads', type=int, default=2,
                        help='Threshold for internal junction support (multi mode): '
                             'keeps junctions with count > this value. '
                             'Default 2 means >=3 reads. '
                             'Use 1 for >=2 reads, or 0 (sc multi recommended) for >=1 read.')
    parser.add_argument('--stringency', choices=['loose', 'strict'], default='loose',
                        help='loose=cryptic junction endpoint anywhere within an annotated exon; '
                             'strict=endpoint must be exactly at the exon boundary.')
    # SC-specific
    parser.add_argument('--cell_id', default=None,
                        help='Cell ID (sc mode). Auto-inferred from BAM filename if not given.')
    parser.add_argument('--primary_only', action='store_true',
                        help='Only use primary alignments. Recommended for PacBio CCS (lr/sc).')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    args = parser.parse_args()
    if args.mapq is None:
        args.mapq = 50 if args.mode == 'sr' else 0
    return args


def config_getter(config_file='config.ini'):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    config = ConfigParser(os.environ)
    config.read(os.path.join(this_dir, config_file))
    annotation_ref = config.get('sorted GENCODE annotation', 'annotation')
    refseq_ref = config.get('NCBI RefSeq annotation', 'refseq')
    return annotation_ref, refseq_ref


# ─────────────────────────────────────────────────────────────────────────────
# Junction extraction (v2 original + SR extension)
# ─────────────────────────────────────────────────────────────────────────────

def find_introns(read_iterator, mode='lr', stranded='no'):
    """
    Extract splice junctions from aligned reads.

    lr/sc : Same as v2 -- ts tag for strand inference, large deletions (>=30bp)
            treated as introns, adjacent deletion boundary correction.
    sr    : XS tag (or library strandedness) for strand inference, N-CIGAR only.

    Returns
    -------
    introns : Counter  {(junc_start, junc_end, strand): count}
    reads   : dict     {(junc_start, junc_end, strand): [read_name, ...]}
    """
    BAM_CREF_SKIP = 3
    match = {0, 7, 8}   # M, =, X

    introns = Counter()
    reads   = defaultdict(list)

    for r in read_iterator:
        if r.cigartuples is None:
            continue
        base_position = r.pos

        # ── Short-read: XS tag ───────────────────────────────────────────
        if mode == 'sr':
            strands_to_use = []
            if stranded == 'no':
                if r.has_tag('XS'):
                    strands_to_use = [r.get_tag('XS')]
                else:
                    # No XS tag: count both strands; downstream GTF annotation
                    # determines the correct strand direction
                    strands_to_use = ['+', '-']
            elif stranded == 'fr-firststrand':
                strand = ('+' if r.is_reverse else '-') if r.is_read1 \
                    else ('-' if r.is_reverse else '+')
                strands_to_use = [strand]
            elif stranded == 'fr-secondstrand':
                strand = ('-' if r.is_reverse else '+') if r.is_read1 \
                    else ('+' if r.is_reverse else '-')
                strands_to_use = [strand]
            if not strands_to_use:
                continue

            for tag, nt in r.cigartuples:
                if tag in match:
                    base_position += nt
                elif tag == BAM_CREF_SKIP:
                    junc_start = base_position
                    base_position += nt
                    junc_end = base_position
                    for s in strands_to_use:
                        introns[(junc_start, junc_end, s)] += 1
                        reads[(junc_start, junc_end, s)].append(r.query_name)
                elif tag == 2:
                    base_position += nt

        # ── Long-read / single-cell: ts tag (v2 original logic) ──────────
        else:
            cigar = r.cigartuples
            n_ops = len(cigar)
            for i, (tag, nt) in enumerate(cigar):
                if tag in match or (tag == 2 and nt < 30):
                    base_position += nt
                elif tag == BAM_CREF_SKIP or (tag == 2 and nt >= 30):
                    junc_start = base_position
                    base_position += nt
                    junc_end = base_position
                    try:
                        ts = r.get_tag('ts')
                        if ts == '+':
                            strand = '-' if r.is_reverse else '+'
                        elif ts == '-':
                            strand = '+' if r.is_reverse else '-'
                        else:
                            continue
                    except KeyError:
                        continue
                    # v2 Bug-fix: adjacent deletion boundary correction (with bounds check)
                    if i > 0 and cigar[i - 1][0] == 2:
                        junc_start -= cigar[i - 1][1]
                    if i + 1 < n_ops and cigar[i + 1][0] == 2:
                        junc_end += cigar[i + 1][1]
                    introns[(junc_start, junc_end, strand)] += 1
                    reads[(junc_start, junc_end, strand)].append(r.query_name)

    return introns, reads


# ─────────────────────────────────────────────────────────────────────────────
# CE detection (v2 original logic + mode/ce_type branches)
# ─────────────────────────────────────────────────────────────────────────────

def _iter_reads(bamfile, chrm, mapq, primary_only):
    """Stream BAM reads, skipping low-quality and non-primary alignments."""
    for r in bamfile.fetch(chrm):
        if r.mapping_quality < mapq:
            continue
        if primary_only and (r.is_secondary or r.is_supplementary):
            continue
        yield r


def ce_caller(bamfile, referencename, refseq_ref, chrm,
              mode='lr', stranded='no', mapq=0,
              primary_only=False, min_junc=2, ce_type='multi',
              stringency='loose'):
    """
    Detect cryptic exons on a single chromosome. Core algorithm is the same
    as v2, with the following additions:
      - SR mode (PSI, a_count)
      - LR/SC mode (ao_canon + PSI quantification)
      - Single CE type (no internal junction)
      - primary_only filtering (SC/LR)
      - Streaming read iteration to avoid loading all reads into memory
    """
    known_splices_D = defaultdict(list)
    known_splices_A = defaultdict(list)

    gtf  = pysam.TabixFile(referencename, parser=pysam.asGTF())
    gtf2 = pysam.TabixFile(refseq_ref, parser=pysam.asGTF())
    db = gffutils.FeatureDB(referencename + '.db')

    # ── Stream reads without accumulating into a list (memory optimization) ──
    introns, reads = find_introns(
        _iter_reads(bamfile, chrm, mapq, primary_only),
        mode=mode, stranded=stranded,
    )

    # ── Match junctions against annotation (v2 original logic) ───────────
    for intron in introns:
        intron_start    = intron[0] - 1
        intron_end      = intron[1] + 2
        intron_witnesses = introns[intron]
        try:
            intersection = list(gtf.fetch(chrm, intron_start, intron_end + 5))
        except ValueError:
            continue

        for feature in intersection:
            if feature.feature != 'exon':
                continue
            region_start  = feature.start
            region_end    = feature.end
            transcript_id = feature.transcript_id
            gene_name     = feature.gene_name
            gene_id       = feature.gene_id
            strand        = feature.strand

            # Donor match
            # loose: junction start anywhere within exon; strict: exactly at exon end
            donor_range = range(region_end, region_end + 1) if stringency == 'strict' \
                else range(region_start, region_end + 1)
            if intron_start + 1 in donor_range:
                try:
                    a = list(db.children(db[transcript_id],
                                         featuretype='exon', order_by='start'))
                    for idx, exon in enumerate(a):
                        if exon.start > intron_start + 2:
                            ref_intron_end_D   = exon.start
                            ref_intron_start_D = a[idx - 1].end
                            if ref_intron_end_D > intron_end - 2 + 1:
                                known_splices_D[(chrm, ref_intron_start_D,
                                                ref_intron_end_D, strand)].append(
                                    (intron_start + 1, intron_end - 2 + 1,
                                     intron_witnesses, transcript_id,
                                     strand, gene_id, gene_name))
                                break
                            if ref_intron_end_D <= intron_end - 2 + 1:
                                break
                except Exception:
                    continue

            # Acceptor match
            # loose: junction end anywhere within exon (or up to 5 bp before exon start,
            #        to tolerate HISAT2 ±1-2 bp de-novo junction coordinate shift);
            # strict: exactly at exon start
            acceptor_range = range(region_start, region_start + 1) if stringency == 'strict' \
                else range(region_start - 5, region_end + 1)
            if intron_end - 2 in acceptor_range:
                try:
                    b = list(db.children(db[transcript_id],
                                         featuretype='exon', order_by='end'))[::-1]
                    for idx, exon in enumerate(b):
                        if exon.end < intron_end - 2 + 1:
                            ref_intron_start_A = exon.end
                            ref_intron_end_A   = b[idx - 1].start
                            if ref_intron_start_A < intron_start + 1:
                                known_splices_A[(chrm, ref_intron_start_A,
                                                ref_intron_end_A, strand)].append(
                                    (intron_start + 1, intron_end - 2 + 1,
                                     intron_witnesses, transcript_id,
                                     strand, gene_id, gene_name))
                                break
                            if ref_intron_start_A >= intron_start + 1:
                                break
                except Exception:
                    continue

    # ── CE candidate pairing (v2 original logic + branches) ──────────────
    ce = []
    chrm_gencode = ['chr1',  'chr2',  'chr3',  'chr4',  'chr5',
                    'chr6',  'chr7',  'chr8',  'chr9',  'chr10',
                    'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                    'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                    'chr21', 'chr22', 'chrX',  'chrY',  'chrM']
    chrm_ncbi    = ['NC_000001.11', 'NC_000002.12', 'NC_000003.12', 'NC_000004.12',
                    'NC_000005.10', 'NC_000006.12', 'NC_000007.14', 'NC_000008.11',
                    'NC_000009.12', 'NC_000010.11', 'NC_000011.10', 'NC_000012.12',
                    'NC_000013.11', 'NC_000014.9',  'NC_000015.10', 'NC_000016.10',
                    'NC_000017.11', 'NC_000018.10', 'NC_000019.10', 'NC_000020.11',
                    'NC_000021.9',  'NC_000022.11', 'NC_000023.11', 'NC_000024.10',
                    'NC_012920.1']
    dic_gencodetoncbi = dict(zip(chrm_gencode, chrm_ncbi))

    for (chrms, ref_intron_start, ref_intron_end, strand) in known_splices_D:
        for i in known_splices_D[(chrms, ref_intron_start, ref_intron_end, strand)]:
            for x in known_splices_A[(chrms, ref_intron_start, ref_intron_end, strand)]:
                if x[0] - i[1] <= 20 or i[3] != x[3]:
                    continue

                # Check that CE candidate region has no known exon overlap
                overlap_exon = []
                try:
                    for f in gtf.fetch(chrms, i[1], x[0]):
                        if f.feature == 'exon':
                            overlap_exon.append(f.gene_id)
                except ValueError:
                    pass
                try:
                    for f in gtf2.fetch(dic_gencodetoncbi[chrms], i[1], x[0]):
                        if f.feature == 'exon':
                            overlap_exon.append(f.contig)
                except (ValueError, KeyError):
                    pass
                if overlap_exon:
                    continue

                # ── Collect internal junctions (v2 bug-fix: use list to avoid key conflicts)
                junction_inside = []
                for intron in introns:
                    if intron[0] > i[1] and intron[1] < x[0]:
                        junction_inside.append(
                            (intron[0], intron[1], introns[intron]))

                # ce_type branching
                if ce_type == 'single':
                    # Single-exon: require no significant internal junction
                    # (filter out noise junctions with count <= min_junc, same as multi mode)
                    junction_inside = [(s, e, c) for s, e, c in junction_inside
                                       if c > min_junc]
                    if len(junction_inside) != 0:
                        continue
                else:
                    # Multi-exon (v2 original behavior): require internal junction
                    # meeting minimum support threshold
                    junction_inside = [(s, e, c) for s, e, c in junction_inside
                                       if c > min_junc]
                    if len(junction_inside) == 0:
                        continue

                # ── AO computation: SR vs LR/SC ─────────────────────────
                if mode == 'sr':
                    ao1 = i[2]
                    ao2 = x[2]
                    ao3 = introns.get((ref_intron_start, ref_intron_end, strand), 0)
                    a_count = sum(
                        1 for r in bamfile.fetch(chrms, i[1], x[0])
                        if r.mapping_quality >= mapq
                        and r.reference_start is not None
                        and r.reference_end  is not None
                        and r.reference_start >= i[1]
                        and r.reference_end   <= x[0]
                    )
                    total = ao1 + ao2 + ao3 + a_count
                    psi   = round((ao1 + ao2 + a_count) / total, 4) if total > 0 else 0.0

                    if ce_type == 'single':
                        ce.append({
                            'chrom': chrms, 'D': ref_intron_start, 'A': ref_intron_end,
                            'ce_start': i[1], 'ce_end': x[0],
                            'ao1': ao1, 'ao2': ao2, 'ao3': ao3,
                            'a_count': a_count, 'psi': psi,
                            'strand': strand, 'gene_id': i[5], 'gene_name': i[6],
                        })
                    else:
                        for (js, je, jc) in junction_inside:
                            ce.append({
                                'chrom': chrms, 'D': ref_intron_start, 'A': ref_intron_end,
                                'ce_start_1': i[1],  'ce_end_1': js,
                                'ce_start_2': je,    'ce_end_2': x[0],
                                'ao1': ao1, 'ao2': jc, 'ao3': ao2,
                                'strand': strand, 'gene_id': i[5], 'gene_name': i[6],
                            })

                elif mode == 'lr':
                    # LR: require same read spanning both outer junctions
                    list_sameread = list(
                        set(reads[(i[0], i[1] - 1, i[4])]).intersection(
                            reads[(x[0], x[1] - 1, x[4])]))
                    ao = len(list_sameread)
                    if ao == 0:
                        continue
                    # PSI = ao (CE-spanning reads) / [ao + canonical junction reads]
                    ao_canon = introns.get(
                        (ref_intron_start, ref_intron_end, strand), 0)
                    psi = round(ao / (ao + ao_canon), 4) \
                        if (ao + ao_canon) > 0 else 0.0

                    if ce_type == 'single':
                        ce.append({
                            'chrom': chrms, 'D': ref_intron_start, 'A': ref_intron_end,
                            'ce_start': i[1], 'ce_end': x[0],
                            'ao': ao, 'ao1': i[2], 'ao2': x[2],
                            'ao_canon': ao_canon, 'psi': psi,
                            'strand': strand, 'gene_id': i[5], 'gene_name': i[6],
                        })
                    else:
                        # v2 Bug-fix: ce_start_2/ce_end_2 order corrected
                        for (js, je, jc) in junction_inside:
                            ce.append({
                                'chrom': chrms, 'D': ref_intron_start, 'A': ref_intron_end,
                                'ce_start_1': i[1],  'ce_end_1': js,
                                'ce_start_2': je,    'ce_end_2': x[0],
                                'ao': ao, 'ao1': i[2], 'ao2': jc, 'ao3': x[2],
                                'ao_canon': ao_canon, 'psi': psi,
                                'strand': strand, 'gene_id': i[5], 'gene_name': i[6],
                            })

                else:
                    # SC: same long read spans both outer junctions + PSI quantification
                    list_sameread = list(
                        set(reads[(i[0], i[1] - 1, i[4])]).intersection(
                            reads[(x[0], x[1] - 1, x[4])]))
                    ao = len(list_sameread)
                    if ao == 0:
                        continue
                    # PSI = ao (CE-spanning reads) / [ao + canonical junction reads]
                    ao_canon = introns.get(
                        (ref_intron_start, ref_intron_end, strand), 0)
                    psi = round(ao / (ao + ao_canon), 4) \
                        if (ao + ao_canon) > 0 else 0.0

                    if ce_type == 'single':
                        ce.append({
                            'chrom': chrms, 'D': ref_intron_start, 'A': ref_intron_end,
                            'ce_start': i[1], 'ce_end': x[0],
                            'ao': ao, 'ao1': i[2], 'ao2': x[2],
                            'ao_canon': ao_canon, 'psi': psi,
                            'strand': strand, 'gene_id': i[5], 'gene_name': i[6],
                        })
                    else:
                        for (js, je, jc) in junction_inside:
                            ce.append({
                                'chrom': chrms, 'D': ref_intron_start, 'A': ref_intron_end,
                                'ce_start_1': i[1],  'ce_end_1': js,
                                'ce_start_2': je,    'ce_end_2': x[0],
                                'ao': ao, 'ao1': i[2], 'ao2': jc, 'ao3': x[2],
                                'ao_canon': ao_canon, 'psi': psi,
                                'strand': strand, 'gene_id': i[5], 'gene_name': i[6],
                            })

    return ce


# ─────────────────────────────────────────────────────────────────────────────
# Filtering and deduplication
# ─────────────────────────────────────────────────────────────────────────────

def filter_ce(ces, ao_min=1, psi_min=0.0, mode='lr', ce_type='multi'):
    ce_filtered = set()
    for ce in ces:
        if ce['ao1'] < ao_min or ce['ao2'] < ao_min:
            continue
        if mode in ('sr', 'sc', 'lr') and ce.get('psi', 0.0) < psi_min:
            continue

        if ce_type == 'single':
            if mode == 'sr':
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start'], ce['ce_end'],
                    ce['ao1'], ce['ao2'], ce['ao3'], ce['a_count'], ce['psi'],
                    ce['strand'], ce['gene_id'], ce['gene_name'],
                ))
            elif mode == 'sc':
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start'], ce['ce_end'],
                    ce['ao'], ce['ao1'], ce['ao2'], ce['ao_canon'], ce['psi'],
                    ce['strand'], ce['gene_id'], ce['gene_name'],
                ))
            else:  # lr
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start'], ce['ce_end'],
                    ce['ao'], ce['ao1'], ce['ao2'],
                    ce['ao_canon'], ce['psi'],
                    ce['strand'], ce['gene_id'], ce['gene_name'],
                ))
        else:  # multi
            if mode == 'sr':
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start_1'], ce['ce_end_1'],
                    ce['ce_start_2'], ce['ce_end_2'],
                    ce['ao1'], ce['ao2'], ce['ao3'],
                    ce['strand'], ce['gene_id'], ce['gene_name'],
                ))
            elif mode == 'sc':
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start_1'], ce['ce_end_1'],
                    ce['ce_start_2'], ce['ce_end_2'],
                    ce['ao'], ce['ao1'], ce['ao2'], ce['ao3'],
                    ce['ao_canon'], ce['psi'],
                    ce['strand'], ce['gene_id'], ce['gene_name'],
                ))
            else:  # lr
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start_1'], ce['ce_end_1'],
                    ce['ce_start_2'], ce['ce_end_2'],
                    ce['ao'], ce['ao1'], ce['ao2'], ce['ao3'],
                    ce['ao_canon'], ce['psi'],
                    ce['strand'], ce['gene_id'], ce['gene_name'],
                ))

    return ce_filtered


# ─────────────────────────────────────────────────────────────────────────────
# Output headers
# ─────────────────────────────────────────────────────────────────────────────

HEADERS = {
    ('sr', 'single'): 'chrom\tD\tA\tce_start\tce_end\tao1\tao2\tao3\ta_count\tPSI\tstrand\tgene_id\tgene_name',
    ('sr', 'multi'):  'chrom\tD\tA\tce_start_1\tce_end_1\tce_start_2\tce_end_2\tao1\tao2\tao3\tstrand\tgene_id\tgene_name',
    ('lr', 'single'): 'chrom\tD\tA\tce_start\tce_end\tao\tao1\tao2\tao_canon\tPSI\tstrand\tgene_id\tgene_name',
    ('lr', 'multi'):  'chrom\tD\tA\tce_start_1\tce_end_1\tce_start_2\tce_end_2\tao\tao1\tao2\tao3\tao_canon\tPSI\tstrand\tgene_id\tgene_name',
    # SC: ao = spanning reads across both outer junctions; ao1/ao2 = per-junction reads;
    #     ao_canon = canonical junction reads (skipping CE); PSI = ao / (ao + ao_canon)
    ('sc', 'single'): 'cell_id\tchrom\tD\tA\tce_start\tce_end\tao\tao1\tao2\tao_canon\tPSI\tstrand\tgene_id\tgene_name',
    ('sc', 'multi'):  'cell_id\tchrom\tD\tA\tce_start_1\tce_end_1\tce_start_2\tce_end_2\tao\tao1\tao2\tao3\tao_canon\tPSI\tstrand\tgene_id\tgene_name',
}


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    annotation_ref, refseq_ref = config_getter('config.ini')

    try:
        bamfile = pysam.AlignmentFile(args.input, 'rb', require_index=True)
    except FileNotFoundError:
        print('Error: BAM file not found: {}'.format(args.input))
        sys.exit(1)
    except ValueError:
        print('Error: BAM index (.bai) is required.')
        sys.exit(1)
    except Exception as e:
        print('Error opening BAM file: {}'.format(e))
        sys.exit(1)

    # SC mode: infer cell_id from BAM filename if not provided
    cell_id = None
    if args.mode == 'sc':
        cell_id = args.cell_id or \
            os.path.basename(bamfile.filename.decode('UTF-8')).split('.')[0]
        print('Cell ID         : {}'.format(cell_id))

    print('Mode            : {}'.format(args.mode))
    print('CE type         : {}'.format(args.ce_type))
    print('MAPQ >=         : {}'.format(args.mapq))
    print('ao_min          : {}'.format(args.ao))
    if args.mode in ('sr', 'sc', 'lr'):
        print('PSI >=          : {}'.format(args.psi))
    if args.mode == 'sr':
        print('Stranded        : {}'.format(args.stranded))
    if args.mode in ('lr', 'sc'):
        print('Primary only    : {}'.format(args.primary_only))
    if args.ce_type == 'multi':
        print('Min junc reads  : {}'.format(args.min_junction_reads))
    print('Stringency      : {}'.format(args.stringency))
    sys.stdout.flush()

    chrms = ['chr1',  'chr2',  'chr3',  'chr4',  'chr5',
             'chr6',  'chr7',  'chr8',  'chr9',  'chr10',
             'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
             'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22', 'chrX',  'chrY',  'chrM']

    ce_total = []
    for chrm in chrms:
        print('Finding cryptic exon in {}'.format(chrm))
        sys.stdout.flush()
        ce = ce_caller(
            bamfile, annotation_ref, refseq_ref, chrm,
            mode=args.mode,
            stranded=args.stranded,
            mapq=args.mapq,
            primary_only=args.primary_only,
            min_junc=args.min_junction_reads,
            ce_type=args.ce_type,
            stringency=args.stringency,
        )
        ce_total.extend(ce)

    ce_filtered = filter_ce(ce_total,
                            ao_min=args.ao,
                            psi_min=args.psi,
                            mode=args.mode,
                            ce_type=args.ce_type)

    # Output filename
    if args.output:
        out_file_name = args.output
    else:
        prefix = cell_id if args.mode == 'sc' \
            else os.path.basename(bamfile.filename.decode('UTF-8')).split('.')[0]
        out_file_name = '{}.{}.{}.ce'.format(prefix, args.mode, args.ce_type)

    print('Finished. {} events -> {}'.format(len(ce_filtered), out_file_name))

    header = HEADERS[(args.mode, args.ce_type)]

    with open(out_file_name, 'w') as out:
        out.write(header + '\n')
        for ce in sorted(ce_filtered, key=lambda e: (e[0], e[1])):
            if args.mode == 'sc':
                # filter_ce output tuples don't include cell_id; prepend it here
                row = cell_id + '\t' + '\t'.join(str(s) for s in ce)
            else:
                row = '\t'.join(str(s) for s in ce)
            out.write(row + '\n')


if __name__ == '__main__':
    main()
