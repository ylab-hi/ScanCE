from ctypes.wintypes import PSIZE
import sys
import os
import argparse
from numpy import sort
import pysam
from collections import Counter, defaultdict
import gffutils
from configparser import ConfigParser

__version__ = 'v1.0'
def parse_args():
    parser = argparse.ArgumentParser(
        description="%(prog)s ",
        epilog="ScanCE v1.0: detecting cryptic exon splicing events using RNA-Seq data",
    )
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input BAM/CRAM file (if no index is found, an index file is created).",
        required=True,
    )
    parser.add_argument(
        "-r1",
        "--reference-annotations1",
        action="store",
        dest="annotation_ref1",
        default=None,
        required=True
    )
    parser.add_argument(
        "-r2",
        "--reference-annotations2",
        action="store",
        dest="annotation_ref2",
        default=None,
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="out",
        help="Output filename. (default: $bam_file_name.ce)",
        default=None,
    )

    parser.add_argument(
        '-m',
        '--mapq',
        action='store',
        dest='mapq',
        type=int,
        help="Consider only reads with MAPQ >= cutoff (default: %(default)s)",
        default=50)

    parser.add_argument(
        "-a",
        "--ao",
        action="store",
        dest="ao_min",
        type=int,
        help="AO cutoff (default: %(default)s)",
        default=1,
    )
    parser.add_argument(
        "-p",
        "--PSI",
        action="store",
        dest="PSI_min",
        type=float,
        help="PSI cutoff (default: %(default)s)",
        default=0,
    )
    

    parser.add_argument(
        "-s",
        "--stranded",
        action="store",
        dest="stranded",
        help="Determines how read strand is inferred. Options are 'no', 'fr-firststrand', 'fr-secondstrand'.  If 'no' for unstranded reads, the XS tag will be used. Otherwise, strand will be inferred based on library prep. See https://chipster.csc.fi/manual/library-type-summary.html for details. (default: %(default)s",
        type=str,
        choices=['no', 'fr-firststrand','fr-secondstrand'],
        default='no',
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {}".format(__version__)
    )
    args = parser.parse_args()
    return args

def find_introns(read_iterator,stranded):
    """
    Parameters
    ----------
    read_iterator : iterator of reads from a pysam.AlignmentFile.
        Expected that the iterator will be an entire chromosome. See exitron_caller
    Returns
    -------
    introns -- counter of (intron_start, intron_end, strand)
    reads -- dictionary of reads that support the junctions in introns.  This
        will be used in later filtering steps.
    """
    BAM_CREF_SKIP = 3

    introns = Counter()
    reads = defaultdict(list)

    # only M/=/X (0/7/8) and D (2) are related to genome position
    match = {0, 7, 8}
    for r in read_iterator:
        base_position = r.pos
        read_position = 0
        # if cigarstring is * (r.cigartuples == None), unmatched, continue
        if r.cigartuples == None:
            continue
        # iterate through cigar string looking for N
        for i, (tag, nt) in enumerate(r.cigartuples):
            # if (0, X), keep track of base_position.
            # if (3, X), which corresponds to N,
            # look at match before and after
            if tag in match or (tag == 2 and nt < 30):
                base_position += nt
                read_position += nt
            elif tag == BAM_CREF_SKIP or (tag == 2 and nt >= 30):
                junc_start = base_position
                base_position += nt
                junc_end = base_position
                try:
                    if r.get_tag('ts') == '+':
                        strand = '-' if r.is_reverse else '+'
                    elif r.get_tag('ts') == '-':
                        strand = '+' if r.is_reverse else '-'
                except KeyError:
                    # ts tag is not present, thus we don't know strand, continue
                    continue
                if r.cigartuples[i - 1][0] == 2:
                    junc_start -= r.cigartuples[i - 1][1]
                if r.cigartuples[i + 1][0] == 2:
                    junc_end += r.cigartuples[i + 1][1]
                introns[(junc_start, junc_end, strand)] += 1
                reads[(junc_start, junc_end, strand)].append(r.query_name)

    return introns, reads 

def ce_caller(bamfile,referencename1, referencename2, chrm, stranded = 'no', mapq = 50):

    known_splices_D = defaultdict(list)
    known_splices_A = defaultdict(list)
    try:
        gtf = pysam.TabixFile(referencename1, parser=pysam.asGTF())
    except OSError:
        print('Building tabix index.')
        pysam.tabix_index(referencename1, preset='gtf')
    try:
        gtf2 = pysam.TabixFile(referencename2, parser=pysam.asGTF())
    except OSError:
        print('Building tabix index.')
        pysam.tabix_index(referencename2, preset='gtf')

    # Prepage gffutils database
    try:
        db = gffutils.FeatureDB(referencename1 + '.db')
    except ValueError:
        print('Preparing annotation database...')
        db = gffutils.create_db(referencename1,
                                dbfn=referencename1 + '.db',
                                disable_infer_transcripts=True,
                                disable_infer_genes=True)
        db = gffutils.FeatureDB(referencename1 + '.db')

    introns, reads = find_introns(
       (read for read in bamfile.fetch(chrm) if read .mapping_quality >= mapq),  
        stranded
        )

    for intron in introns:
        intron_start = intron[0] - 1 
        intron_end = intron[1] + 2 
        #-1 and +2 so that we can capture the ends and beginning of adjacent transcripts
        intron_witnesses = introns[intron]
        intersection = list(gtf.fetch(chrm, intron_start, intron_end))
        for feature in intersection:      
            #if feature.strand != intron[2]: continue # gtf.fetch is strand agnostic
            #if intron[2] =='+':
            region_type = feature.feature
            if region_type == 'exon':
                region_start = feature.start
                region_end = feature.end
                transcript_id = feature.transcript_id 
                gene_name = feature.gene_name
                gene_id = feature.gene_id
                strand = feature.strand
                if intron_start+1 in range(region_end, region_end + 1):
                    db_exons_D =db.children(db[transcript_id], featuretype = 'exon', order_by = 'start')
                    a = list(db_exons_D)

                    for i, exon in enumerate(a):
                        if exon.start > intron_start+2:
                            ref_intron_end_D = exon.start
                            ref_intron_start_D = a[i-1].end
                            if ref_intron_end_D > intron_end-2+1:
                                known_splices_D[(chrm, ref_intron_start_D,ref_intron_end_D,strand)].append((intron_start+1, intron_end-2+1, intron_witnesses, transcript_id, strand, gene_id, gene_name))
                                break
                            if ref_intron_end_D <= intron_end-2+1:
                                break
                # intron matches a known acceptor
                if intron_end-2 in range(region_start, region_end + 1):
                    db_exons_A = db.children(db[transcript_id], featuretype = 'exon', order_by = 'end')
                    b = list(db_exons_A)[::-1]
                    for i, exon in enumerate(b):
                        if exon.end < intron_end-2+1:
                            ref_intron_start_A = exon.end
                            ref_intron_end_A = b[i-1].start
                            if ref_intron_start_A < intron_start+1:
                                known_splices_A[(chrm, ref_intron_start_A,ref_intron_end_A,strand)].append((intron_start+1, intron_end-2+1, intron_witnesses, transcript_id, strand, gene_id, gene_name))
                                break
                            if ref_intron_start_A >= intron_start+1:
                                break
   
   
    ce=[]
    known_exons = set()
    chrm_gencode= ['chr1','chr2','chr3','chr4', 'chr5',
         'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11','chr12', 'chr13', 'chr14', 'chr15',
         'chr16','chr17', 'chr18', 'chr19', 'chr20',
         'chr21', 'chr22', 'chrX','chrY', 'chrM']
    chrm_ncbi = ['NC_000001.11','NC_000002.12','NC_000003.12','NC_000004.12','NC_000005.10',
            'NC_000006.12','NC_000007.14','NC_000008.11','NC_000009.12','NC_000010.11',
            'NC_000011.10','NC_000012.12','NC_000013.11','NC_000014.9','NC_000015.10',
            'NC_000016.10','NC_000017.11','NC_000018.10','NC_000019.10','NC_000020.11',
            'NC_000021.9','NC_000022.11','NC_000023.11','NC_000024.10','NC_012920.1']
    dic_gencodetoncbi=dict(zip(chrm_gencode,chrm_ncbi))
    for (chrms, ref_intron_start,ref_intron_end,strand) in known_splices_D:
        for i in known_splices_D[(chrms, ref_intron_start,ref_intron_end,strand)]:
            for x in known_splices_A[(chrms, ref_intron_start,ref_intron_end,strand)]:
                if x[0] - i[1] > 20 and i[3]==x[3]:
                    intersection_ce = list(gtf.fetch(chrms, i[1], x[0]))
                    overlap_exon=[]
                    for feature in intersection_ce:      
                        region_type = feature.feature
                        if region_type == 'exon':
                             overlap_exon.append(feature.gene_id)  
                    intersection_ce2 = list(gtf2.fetch(dic_gencodetoncbi[chrms], i[1], x[0]))
                    for feature in intersection_ce2:      
                        region_type = feature.feature
                        if region_type == 'exon':
                             overlap_exon.append(feature.contig)
                    if len(overlap_exon)==0:
                        junction_inside=set()
                        a_count = 0
                        for intron in introns:
                            if intron[0] > i[1] and  intron[1] < x[0] and introns[intron]>2:
                                junction_inside.add((intron[0], intron[1]))
                        if len(junction_inside)==0:
                            a_area = bamfile.fetch(chrms, start = i[1], stop = x[0])
                            for r in a_area:
                                if r.mapping_quality >= 50:
                                    if r.reference_start > i[1] and r.reference_end < x[0]:
                                        a_count = a_count+1
                            
                            list_sameread=list(set(reads[(i[0],i[1]-1,i[4])]).intersection(reads[(x[0],x[1]-1,x[4])]))
                            ao = len(list_sameread)
                            if ao > 0:
                                ce.append({'chrom':chrms,
                                        "D":ref_intron_start,
                                        "A":ref_intron_end,
                                        'start':str(i[1])+"+"+str(i[0]),
                                        'end':str(x[0])+"+"+str(x[1]),
                                        'ao':ao,
                                        'ao1':i[2],
                                        'ao2':x[2],
                                        'ao3':introns[ref_intron_start,ref_intron_end-1],
                                        'a_count':a_count,
                                        'transcript_id':i[3],
                                        'strand':i[4],
                                        'gene_id':i[5],
                                        'gene_name':i[6]})



    return ce


def filter_ce(ces, ao_min=1, PSI_min=0):

    ce_filtered=set()
   
    for ce in ces:
        
        chrm = ce['chrom']
        ao1 = ce['ao1'] 
        ao2 = ce['ao2'] 
        ao3 = ce['ao3']
        ao = ce['ao']
        a_count = ce['a_count']
        ce_start = ce['start']
        ce_end = ce['end']
        #transcript_id=ce['transcript_id']
        strand = ce['strand']
        intron_start = ce['D']
        intron_end = ce['A']
        gene_id = ce['gene_id']
        gene_name = ce['gene_name']
        PSI = (a_count + ao1 + ao2)/ (a_count + ao1 + ao2 +ao3)
        if ao1 < ao_min or ao2 < ao_min: 
            continue
        if PSI < PSI_min: 
            continue
        ce_filtered.add((chrm,intron_start,intron_end,ce_start,ce_end,ao,ao1,ao2,ao3,a_count,PSI,strand,gene_id,gene_name))
   
    return ce_filtered

    
def main():
    args = parse_args()
    try:
        bamfile = pysam.AlignmentFile(args.input, 'rb', require_index = True)
    except FileNotFoundError:
        try:
            print(' bam index file is needed')
        except FileNotFoundError:
            print(f'There is a problem opening bam file at: {args.input}')
    
    chrms = ['chr1','chr2','chr3','chr4', 'chr5',
         'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11','chr12', 'chr13', 'chr14', 'chr15',
         'chr16','chr17', 'chr18', 'chr19', 'chr20',
         'chr21', 'chr22', 'chrX','chrY', 'chrM']

    ce_total=[]
    for chrm in chrms:
        print(f'Finding cryptic exon in {chrm}')
        sys.stdout.flush()
        ce  = ce_caller(bamfile,
                    args.annotation_ref1,
                    args.annotation_ref2,
                    chrm,
                    args.stranded,
                    args.mapq)
        ce_total.extend(ce)
    print("beforefilter")
    ce_filtered=filter_ce(ce_total, ao_min=1, PSI_min=0)
        
       

    out_file_name = args.out
    if not out_file_name:
        prefix = os.path.splitext(os.path.basename(bamfile.filename.decode('UTF-8')))[0]
        out_file_name = f'{prefix}.LR.loosece'
    print(f'Finished ce calling and filtering. Printing to {out_file_name}')
 
    with open(out_file_name, 'w') as out:
        out.write('chrom\tD\tA\tstart\tend\tao\tao1\tao2\tao3\ta_count\tPSI\tstrand\tgene_id\tgene_name\n')  
        ce_list = list(ce_filtered)

        ce_result=sorted(ce_list, key=lambda element: (element[0], element[1])) 

        for ce in ce_result:
            out.write('\t'.join(str(s) for s in ce) + '\n')
            

if __name__ == '__main__':

    main()