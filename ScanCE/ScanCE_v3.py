#!/usr/bin/env python3
"""
ScanCE v3.0 — 综合版 Cryptic Exon Detector

直接基于 Scan_ce_hpc_longread_loose_multi_v2.py 扩展而来。

新增功能：
  --mode   sr  : 短读（Illumina bulk），使用 XS tag 推断链方向，计算 PSI
           lr  : 长读（ONT/PacBio bulk），与原 v2 逻辑完全一致
           sc  : 单细胞长读（PacBio CCS/MAS-seq），LR 逻辑 + cell_id + primary_only

  --ce_type single : 单外显子 CE（CE 区域内无 internal junction）
            multi  : 多外显子 CE（CE 区域内有 internal junction，v2 原有逻辑）
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
# 参数解析
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
                        help='Min PSI threshold (sr mode only).')
    parser.add_argument('--min_junction_reads', type=int, default=2,
                        help='Min reads for internal junction (multi mode). '
                             'Lower to 1 for sparse single-cell data.')
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
    annotation_ref  = config.get('sorted GENCODE annotation', 'annotation')
    annotation_ref2 = config.get('sorted NCBI annotation', 'annotation')
    return annotation_ref, annotation_ref2


# ─────────────────────────────────────────────────────────────────────────────
# Junction 提取（v2 原版 + SR 扩展）
# ─────────────────────────────────────────────────────────────────────────────

def find_introns(read_iterator, mode='lr', stranded='no'):
    """
    从比对 reads 中提取剪接 junction。

    lr/sc : 与 v2 完全一致 —— ts tag 推断链方向，大 deletion(>=30bp) 视为 intron，
            相邻 deletion 校正 junction 边界。
    sr    : XS tag（或 library strandedness）推断链方向，仅处理 N-CIGAR。

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

        # ── 短读：XS tag ──────────────────────────────────────────────────
        if mode == 'sr':
            strand = None
            if stranded == 'no':
                try:
                    strand = r.get_tag('XS')   # '+' or '-'
                except KeyError:
                    continue
            elif stranded == 'fr-firststrand':
                strand = ('+' if r.is_reverse else '-') if r.is_read1 \
                    else ('-' if r.is_reverse else '+')
            elif stranded == 'fr-secondstrand':
                strand = ('-' if r.is_reverse else '+') if r.is_read1 \
                    else ('+' if r.is_reverse else '-')
            if strand is None:
                continue

            for tag, nt in r.cigartuples:
                if tag in match:
                    base_position += nt
                elif tag == BAM_CREF_SKIP:
                    junc_start = base_position
                    base_position += nt
                    junc_end = base_position
                    introns[(junc_start, junc_end, strand)] += 1
                    reads[(junc_start, junc_end, strand)].append(r.query_name)
                elif tag == 2:
                    base_position += nt

        # ── 长读/单细胞：ts tag（v2 原版逻辑）────────────────────────────
        else:
            for i, (tag, nt) in enumerate(r.cigartuples):
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
                    # v2 Bug-fix: 相邻 deletion 校正
                    if r.cigartuples[i - 1][0] == 2:
                        junc_start -= r.cigartuples[i - 1][1]
                    if r.cigartuples[i + 1][0] == 2:
                        junc_end += r.cigartuples[i + 1][1]
                    introns[(junc_start, junc_end, strand)] += 1
                    reads[(junc_start, junc_end, strand)].append(r.query_name)

    return introns, reads


# ─────────────────────────────────────────────────────────────────────────────
# CE 检测（v2 原版逻辑 + mode/ce_type 分支）
# ─────────────────────────────────────────────────────────────────────────────

def ce_caller(bamfile, referencename, referencename2, chrm,
              mode='lr', stranded='no', mapq=0,
              primary_only=False, min_junc=2, ce_type='multi',
              stringency='loose'):
    """
    对单条染色体进行 CE 检测。核心算法与 v2 一致，新增：
      - SR 模式（PSI、a_count）
      - single CE 类型（无 internal junction）
      - primary_only 过滤（SC/LR）
    """
    known_splices_D = defaultdict(list)
    known_splices_A = defaultdict(list)

    gtf  = pysam.TabixFile(referencename, parser=pysam.asGTF())
    gtf2 = pysam.TabixFile(referencename2, parser=pysam.asGTF())
    db = gffutils.FeatureDB(referencename + '.db')

    # ── 读取 reads（新增 primary_only 过滤）─────────────────────────────
    all_reads = []
    for r in bamfile.fetch(chrm):
        if r.mapping_quality < mapq:
            continue
        if primary_only and (r.is_secondary or r.is_supplementary):
            continue
        all_reads.append(r)

    introns, reads = find_introns(all_reads, mode=mode, stranded=stranded)

    # ── 与注释比对（v2 原版逻辑不变）─────────────────────────────────────
    for intron in introns:
        intron_start    = intron[0] - 1
        intron_end      = intron[1] + 2
        intron_witnesses = introns[intron]
        try:
            intersection = list(gtf.fetch(chrm, intron_start, intron_end))
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
            # loose: junction end anywhere within exon; strict: exactly at exon start
            acceptor_range = range(region_start, region_start + 1) if stringency == 'strict' \
                else range(region_start, region_end + 1)
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

    # ── CE candidate 配对（v2 原版逻辑 + 分支）───────────────────────────
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

                # 检查 CE 候选区域无已知外显子覆盖
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

                # ── internal junction 收集（v2 Bug-fix：用 list 避免 key 冲突）
                junction_inside = []
                for intron in introns:
                    if intron[0] > i[1] and intron[1] < x[0]:
                        junction_inside.append(
                            (intron[0], intron[1], introns[intron]))

                # ce_type 分支
                if ce_type == 'single':
                    # 单外显子：要求无 internal junction
                    if len(junction_inside) != 0:
                        continue
                else:
                    # 多外显子（v2 原始行为）：要求有 internal junction 且满足最低支持数
                    junction_inside = [(s, e, c) for s, e, c in junction_inside
                                       if c > min_junc]
                    if len(junction_inside) == 0:
                        continue

                # ── AO 计算：SR vs LR/SC ─────────────────────────────────
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

                else:
                    # LR/SC：同一条 read 必须跨越两个外侧 junction（v2 原版 ao 计算）
                    list_sameread = list(
                        set(reads[(i[0], i[1] - 1, i[4])]).intersection(
                            reads[(x[0], x[1] - 1, x[4])]))
                    ao = len(list_sameread)
                    if ao == 0:
                        continue

                    if ce_type == 'single':
                        ce.append({
                            'chrom': chrms, 'D': ref_intron_start, 'A': ref_intron_end,
                            'ce_start': i[1], 'ce_end': x[0],
                            'ao': ao, 'ao1': i[2], 'ao2': x[2],
                            'strand': strand, 'gene_id': i[5], 'gene_name': i[6],
                        })
                    else:
                        # v2 Bug-fix：ce_start_2/ce_end_2 顺序已修正
                        for (js, je, jc) in junction_inside:
                            ce.append({
                                'chrom': chrms, 'D': ref_intron_start, 'A': ref_intron_end,
                                'ce_start_1': i[1],  'ce_end_1': js,
                                'ce_start_2': je,    'ce_end_2': x[0],
                                'ao': ao, 'ao1': i[2], 'ao2': jc, 'ao3': x[2],
                                'strand': strand, 'gene_id': i[5], 'gene_name': i[6],
                            })

    return ce


# ─────────────────────────────────────────────────────────────────────────────
# 过滤与去重
# ─────────────────────────────────────────────────────────────────────────────

def filter_ce(ces, ao_min=1, psi_min=0.0, mode='lr', ce_type='multi'):
    ce_filtered = set()
    for ce in ces:
        if ce['ao1'] < ao_min or ce['ao2'] < ao_min:
            continue
        if mode == 'sr' and ce.get('psi', 0.0) < psi_min:
            continue

        if ce_type == 'single':
            if mode == 'sr':
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start'], ce['ce_end'],
                    ce['ao1'], ce['ao2'], ce['ao3'], ce['a_count'], ce['psi'],
                    ce['strand'], ce['gene_id'], ce['gene_name'],
                ))
            else:
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start'], ce['ce_end'],
                    ce['ao'], ce['ao1'], ce['ao2'],
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
            else:
                ce_filtered.add((
                    ce['chrom'], ce['D'], ce['A'],
                    ce['ce_start_1'], ce['ce_end_1'],
                    ce['ce_start_2'], ce['ce_end_2'],
                    ce['ao'], ce['ao1'], ce['ao2'], ce['ao3'],
                    ce['strand'], ce['gene_id'], ce['gene_name'],
                ))

    return ce_filtered


# ─────────────────────────────────────────────────────────────────────────────
# 输出表头
# ─────────────────────────────────────────────────────────────────────────────

HEADERS = {
    ('sr', 'single'): 'chrom\tD\tA\tce_start\tce_end\tao1\tao2\tao3\ta_count\tPSI\tstrand\tgene_id\tgene_name',
    ('sr', 'multi'):  'chrom\tD\tA\tce_start_1\tce_end_1\tce_start_2\tce_end_2\tao1\tao2\tao3\tstrand\tgene_id\tgene_name',
    ('lr', 'single'): 'chrom\tD\tA\tce_start\tce_end\tao\tao1\tao2\tstrand\tgene_id\tgene_name',
    ('lr', 'multi'):  'chrom\tD\tA\tce_start_1\tce_end_1\tce_start_2\tce_end_2\tao\tao1\tao2\tao3\tstrand\tgene_id\tgene_name',
    ('sc', 'single'): 'cell_id\tchrom\tD\tA\tce_start\tce_end\tao\tao1\tao2\tstrand\tgene_id\tgene_name',
    ('sc', 'multi'):  'cell_id\tchrom\tD\tA\tce_start_1\tce_end_1\tce_start_2\tce_end_2\tao\tao1\tao2\tao3\tstrand\tgene_id\tgene_name',
}


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    annotation_ref, annotation_ref2 = config_getter('config.ini')

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

    # SC 模式：推断 cell_id
    cell_id = None
    if args.mode == 'sc':
        cell_id = args.cell_id or \
            os.path.basename(bamfile.filename.decode('UTF-8')).split('.')[0]
        print('Cell ID         : {}'.format(cell_id))

    print('Mode            : {}'.format(args.mode))
    print('CE type         : {}'.format(args.ce_type))
    print('MAPQ >=         : {}'.format(args.mapq))
    print('ao_min          : {}'.format(args.ao))
    if args.mode == 'sr':
        print('Stranded        : {} | PSI >= {}'.format(args.stranded, args.psi))
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
            bamfile, annotation_ref, annotation_ref2, chrm,
            mode=args.mode,
            stranded=args.stranded,
            mapq=args.mapq,
            primary_only=args.primary_only,
            min_junc=args.min_junction_reads,
            ce_type=args.ce_type,
            stringency=args.stringency,
        )
        ce_total.extend(ce)

    filter_mode = 'lr' if args.mode == 'sc' else args.mode
    ce_filtered = filter_ce(ce_total,
                            ao_min=args.ao,
                            psi_min=args.psi,
                            mode=filter_mode,
                            ce_type=args.ce_type)

    # 输出文件名
    if args.output:
        out_file_name = args.output
    else:
        prefix = cell_id if args.mode == 'sc' \
            else os.path.basename(bamfile.filename.decode('UTF-8')).split('.')[0]
        out_file_name = '{}.{}.{}.ce'.format(prefix, args.mode, args.ce_type)

    print('Finished. {} events → {}'.format(len(ce_filtered), out_file_name))

    header = HEADERS[(args.mode, args.ce_type)]

    with open(out_file_name, 'w') as out:
        out.write(header + '\n')
        for ce in sorted(ce_filtered, key=lambda e: (e[0], e[1])):
            row = '\t'.join(str(s) for s in ce)
            if args.mode == 'sc':
                row = cell_id + '\t' + row
            out.write(row + '\n')


if __name__ == '__main__':
    main()
