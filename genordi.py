#!/usr/bin/env python3

import sys
import logging
import argparse
import unittest

__author__ = 'Qiyun Zhu'
__license__ = 'BSD-3-Clause'
__version__ = '0.0.1-dev'
__email__ = 'qiyunzhu@gmail.com'

usage = """%(prog)s -i INPUT_ALIGN -g GENE_TABLE -o OUTPUT_PROFILE [options]"""

description = """examples:
  %(prog)s -i read_map.txt -g gene_table.txt -o profile.txt
"""

epilog = """Genordi: Match read alignments and annotated genes in an ordinal \
scale on genome."""


def parse_args():
    """Command-line interface.
    """
    parser = argparse.ArgumentParser(
        usage=usage, description=description, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    arg = parser.add_argument
    arg('-i', '--input', type=argparse.FileType('r'), default=sys.stdin,
        help='input read alignment, default: stdin')
    arg('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        help=('output profile, format: gene <tab> count, default: stdout'))
    arg('-g', '--genetab', type=argparse.FileType('rt'),
        help=('table of gene coordinates on genome sequences'))
    arg('-f', '--format', choices=['auto', 'b6o', 'sam', 'map'],
        default='auto',
        help=('format of read alignment: "auto": automatic determination '
              '(default), "b6o": BLAST tabular format (-outfmt 6), "sam": SAM '
              'format, "map": read-to-gene(s) map'))
    arg('-t', '--threshold', type=float, default=0.8,
        help=('minimum ratio of overlap length vs. alignment length to '
              'qualify for a match, default: 0.8'))
    arg('-a', '--ambiguity', choices=['uniq', 'all', 'norm'], default='norm',
        help=('how to treat reads that occur multiple times: "uniq": drop '
              'non-unique reads; "all": count each occurrence once; "norm": '
              'count each occurrence 1/k times (k = number of occurrences) '
              '(default)'))
    arg('-n', '--lines', type=int, default=1000000,
        help=('number of lines per chunck for parsing; higher value improves '
              'time efficiency but consumes more memory, default: 1 million'))
    arg('-m', '--outmap', action='store_true',
        help=('generate a read-to-gene(s) map instead of a gene-to-count '
              'profile'))
    arg('-p', '--prefix', action='store_true',
        help=('prefix nucleotide ID to gene ID in the format of "nucl_gene"'))
    for arg in parser._actions:
        arg.metavar = ''
    return parser.parse_args()


def main():
    """Main workflow.
    """
    # config logger
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s %(message)s')

    # parser arguments
    args = parse_args()

    # whether return read map or gene counts (latter saves compute)
    ismap = args.outmap or args.format == 'map' or args.ambiguity != 'all'

    # read gene coordinates
    if args.genetab:
        logging.info('Gene table parsing started.')
        genetab = read_gene_table(args.genetab)

        # sort gene coordinates per nucleotide
        for gene, queue in genetab.items():
            genetab[gene] = sorted(queue, key=lambda x: x[0])
        logging.info('Gene table parsing completed.')

    # gene profile to be constructed
    profile = {}

    # identifiers, alignment lengths and coordinates of reads
    # note: for a read, "id" is the identifier (a free text), "idx" is the
    # index (an incremental integer for all reads of the current chunck)
    rids, lenmap, locmap = [], {}, {}

    # read input read map
    logging.info('Input read alignment parsing started.')

    # line counter
    ii = 0

    # associate reads with genes for current chunck
    def _match_reads_genes():
        nonlocal rids, lenmap, locmap
        readmap = {}
        for nucl, loci in locmap.items():

            # merge and sort coordinates
            # question is to merge an unsorted list into a sorted one
            # Python's built-in "timesort" algorithm is efficient at this
            try:
                queue = sorted(genetab[nucl] + loci, key=lambda x: x[0])

            # it's possible that no gene was annotated on the nucleotide
            except KeyError:
                continue

            # map reads to genes
            res = match_read_gene(queue, lenmap[nucl], args.threshold, ismap)

            # prefix
            pref = nucl if args.prefix else None

            # merge into master read map (of current chunck)
            if args.outmap:
                add_match_to_readmap(readmap, res, pref)

            # merge into master profile (of entire sample)
            else:
                add_match_to_profile(profile, res, ismap, pref)

        # write read map
        if args.outmap:
            logging.info('Mapped reads: {}.'.format(len(readmap)))
            for i, rid in enumerate(rids):
                try:
                    args.output.write('{}\t{}\n'.format(rid, ','.join(sorted(
                        readmap[i]))))
                except KeyError:
                    pass

        # free memory
        lenmap, locmap = {}, {}
        logging.info('Parsed {} lines.'.format(ii))

    # parser for read-to-gene(s) map
    def _parse_map_line(line, *args):
        nonlocal profile, rids
        rid, genes = line.rstrip('\r\n').split('\t')
        rix = len(rids)
        rids.append(rid)
        for gene in genes.split(','):
            profile.setdefault(gene, []).append(rix)

    # determine parser
    def _assign_parser(fmt):
        """Assign parser function based on format code.
        """
        if fmt == 'map':  # read-to-gene(s) map
            return _parse_map_line
        if fmt == 'b6o':  # BLAST format
            return parse_b6o_line
        elif args.format == 'sam':  # SAM format
            return parse_sam_line
        else:
            logging.error('Invalid format code: {}.'.format(args.format))
            sys.exit(1)

    # determine alignment format
    if args.format == 'auto':  # auto-determine
        line = args.input.readline()
        try:
            args.format = infer_align_format(line)
            logging.info('Alignment format: {}.'.format(args.format))
        except ValueError:
            logging.error('Alignment format cannot be determined.')
            sys.exit(1)
        if args.format == 'map':
            ismap = True
        parser = _assign_parser(args.format)
        parser(line, rids, lenmap, locmap)
        ii += 1
    else:
        parser = _assign_parser(args.format)

    # parse read map in chuncks
    for line in args.input:
        parser(line, rids, lenmap, locmap)
        ii += 1
        if args.lines and ii % args.lines == 0:
            _match_reads_genes()
    _match_reads_genes()
    logging.info('Input read alignment parsing completed.')

    # read map is already written
    if args.outmap:
        logging.info('Output read-to-gene(s) map written.')
        return

    # convert read maps into counts
    if ismap:
        readmap_to_profile(profile, rids, args.ambiguity == 'norm')
        logging.info('Multiple-occurrence reads processed.')

    # write profile
    for gene, n in sorted(profile.items()):
        args.output.write('{}\t{}\n'.format(gene, n))
    logging.info('Output gene profile written.')


def read_gene_table(f):
    """Read coordinates of genes on genomes.

    Parameters
    ----------
    f : file handle
        Gene table file.

    Returns
    -------
    dict of list of tuple of (int, bool, bool, str)
        Flattened list of gene coordinates per nucleotide.
            Coordinate (nt).
            Whether start (True) or end (False).
            Whether gene (True) or read (False).
            Identifier of gene.

    See Also
    --------
    map_read_gene

    Notes
    -----
    This data structure is central to this algorithm. Starting and ending
    coordinates of each gene are separated and flattened into a sorted list.
    which enables only one round of list traversal for the entire set of genes
    plus reads.
    """
    res = {}
    nucl = None
    for line in f:
        line = line.rstrip('\r\n')

        # ">" or "#" indicates nucleotide name
        if line.startswith(('>', '#')):
            nucl = line[1:].strip()

            # double ">" or "#" indicates genome name
            if not nucl.startswith(('>', '#')):
                res[nucl] = []
        else:
            x = line.split('\t')
            idx = x[0]

            # start and end are based on genome, not gene itself
            start, end = sorted([int(x[1]), int(x[2])])
            res[nucl].extend((
                (start, True, True, idx),
                (end,  False, True, idx)))
    return res


def infer_align_format(line):
    """Guess the format of an alignment file based on first line.

    Parameters
    ----------
    line : str
        First line of alignment.

    Returns
    -------
    str
        Alignment file format (map, b6o or sam).

    Raises
    ------
    ValueError
        Format cannot be determined.

    See Also
    --------
    parse_b6o_line
    parse_sam_line
    """
    if line.split()[0] == '@HD':
        return 'sam'
    row = line.rstrip('\r\n').split('\t')
    if len(row) == 2:
        return 'map'
    if len(row) == 12:
        if all(row[i].isdigit() for i in range(3, 10)):
            return 'b6o'
    if len(row) >= 11:
        if all(row[i].isdigit() for i in (1, 3, 4)):
            return 'sam'


def parse_b6o_line(line, rids, lenmap, locmap):
    """Parse a line in a BLAST tabular file (b6o).

    Parameters
    ----------
    line : str
        Line to parse.
    rids : list
        Read identifiers.
    lenmap : dict of dict
        Query to alignment length map.
    locmap : dict of list
        Coordinates of features.

    See Also
    --------
    read_gene_table

    Notes
    -----
    BLAST tabular format:
        qseqid sseqid pident length mismatch gapopen qstart qend sstart send
        evalue bitscore

    .. _BLAST manual:
        https://www.ncbi.nlm.nih.gov/books/NBK279684/
    """
    x = line.rstrip('\r\n').split('\t')
    qseqid, sseqid, length = x[0], x[1], int(x[3])
    sstart, send = sorted([int(x[8]), int(x[9])])
    idx = len(rids)
    rids.append(qseqid)
    lenmap.setdefault(sseqid, {})[idx] = length
    locmap.setdefault(sseqid, []).extend((
        (sstart, True, False, idx),
        (send,  False, False, idx)))


def parse_sam_line(line, rids, lenmap, locmap):
    """Parse a line in a SAM format (sam).

    Parameters
    ----------
    line : str
        Line to parse.
    rids : list
        Read identifiers.
    lenmap : dict of dict
        Query to alignment length map.
    locmap : dict of list
        Coordinates of features.

    See Also
    --------
    read_gene_table

    Notes
    -----
    SAM format:
        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL,
        TAGS

    .. _Wikipedia:
        https://en.wikipedia.org/wiki/SAM_(file_format)
    .. _SAM format specification:
        https://samtools.github.io/hts-specs/SAMv1.pdf
    .. _Bowtie2 manual:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output
    """
    # skip header
    if line.startswith('@'):
        return
    x = line.rstrip('\r\n').split('\t')
    qname, rname = x[0], x[2]  # query and subject identifiers

    # skip unmapped
    if rname == '*':
        return
    pos = int(x[3])  # leftmost mapping position

    # parse CIGAR string
    length, offset = cigar_to_lens(x[5])

    # append strand to read Id if not already
    if not qname.endswith(('/1', '/2')):
        flag = int(x[1])
        if flag & (1 << 6):  # forward strand: bit 64
            qname += '/1'
        elif flag & (1 << 7):  # reverse strand: bit 128
            qname += '/2'

    # store read Id
    idx = len(rids)
    rids.append(qname)

    # update maps
    lenmap.setdefault(rname, {})[idx] = length
    locmap.setdefault(rname, []).extend((
        (pos, True, False, idx),
        (pos + offset - 1,  False, False, idx)))


def cigar_to_lens(cigar):
    """Extract lengths from a CIGAR string.

    Parameters
    ----------
    cigar : str
        CIGAR string.

    Returns
    -------
    int, int
        Alignment length.
        Offset in subject sequence.

    Raises
    ------
    ValueError
        CIGAR string is missing.
    """
    if cigar in ('', '*'):
        raise ValueError('Missing CIGAR string.')
    align, offset = 0, 0
    ops = 'DHIMNPSX='
    n = ''  # current step size
    for c in cigar:
        if c in ops:
            if c in 'M=X':
                align += int(n)
            if c in 'MDN=X':
                offset += int(n)
            n = ''
        else:
            n += c
    return align, offset


def match_read_gene(queue, lens, th=0.8, ismap=False):
    """Associate reads with genes based on a sorted queue of coordinates.

    Parameters
    ----------
    queue : list of tuple
        Sorted list of elements.
        (loc, is_start, is_gene, id)
    lens : dict
        Read-to-alignment length map.
    th : float, optional
        Threshold for read/gene overlapping fraction.
    ismap : bool, optional
        Return read-to-gene map (True) instead of per-gene counts (False).

    Returns
    -------
    dict
        Per-gene counts or read-to-gene map.

    See Also
    --------
    read_gene_table

    Notes
    -----
    This algorithm is the core of the program. It uses a flattened, sorted
    list to store starting and ending coordinates of both genes and reads.
    Only one round of traversal (O(n)) of this list is needed to accurately
    find all gene-read matches.
    """
    match = {}  # read/gene match
    genes = {}  # current genes
    reads = {}  # current reads

    def _add_to_match(rid, gid):
        if ismap:
            match.setdefault(rid, set()).add(gid)
        else:
            match[gid] = match.get(gid, 0) + 1

    for loc, is_start, is_gene, id_ in queue:
        if is_gene:

            # when a gene starts, added to current genes
            if is_start:
                genes[id_] = loc

            # when a gene ends,
            else:

                # check current reads
                for rid, rloc in reads.items():

                    # add to match if read/gene overlap is long enough
                    if loc - max(genes[id_], rloc) + 1 >= lens[rid] * th:
                        _add_to_match(rid, id_)

                # remove it from current genes
                del(genes[id_])

        # the same for reads
        else:
            if is_start:
                reads[id_] = loc
            else:
                for gid, gloc in genes.items():
                    if loc - max(reads[id_], gloc) + 1 >= lens[id_] * th:
                        _add_to_match(id_, gid)
                del(reads[id_])
    return match


def add_match_to_readmap(readmap, match, nucl=None):
    """Merge current read-gene matches to master read map.

    Parameters
    ----------
    readmap : dict
        Master read map.
    match : dict
        Current read map.
    nucl : str, optional
        Prefix nucleotide Id to gene Ids.
    """
    for rid, genes in match.items():
        if nucl:
            genes = {'{}_{}'.format(nucl, x) for x in genes}
        readmap.setdefault(rid, set()).update(genes)


def add_match_to_profile(profile, match, ismap=True, nucl=None):
    """Merge current read-gene matches to master profile.

    Parameters
    ----------
    profile : dict
        Master gene profile.
    match : dict
        Read-gene matches.
    ismap : bool, optional
        Whether matches are a read-to-gene(s) map or simple counts.
    nucl : str, optional
        Prefix nucleotide Id to gene Ids.

    See Also
    --------
    match_read_gene
    """
    # prefix gene Id with nucleotide Id
    def prefix(gene):
        return '{}_{}'.format(nucl, gene) if nucl else gene

    # read-to-gene(s) map
    if ismap:
        for ridx, genes in match.items():
            for gene in genes:
                profile.setdefault(prefix(gene), []).append(ridx)

    # simple counts
    else:
        for gene, count in match.items():
            gene = prefix(gene)
            profile[gene] = profile.get(gene, 0) + count


def readmap_to_profile(profile, rids, normalize=True):
    """Convert read-to-gene map into gene profile.

    Parameters
    ----------
    profile : dict
        Master gene profile.
    rids : list
        Read identifiers.
    normalize : bool, optional
        For a read occurring k times, normalize by k (True), or drop (False).
    """
    # count occurrences of reads
    freqs = {}
    for read in rids:
        freqs[read] = freqs.get(read, 0) + 1

    # treat multiple occurrences
    todel = []
    for gene, ridxx in profile.items():
        n = 0
        for ridx in ridxx:
            k = freqs[rids[ridx]]
            if k > 1:
                if normalize:
                    n += 1 / k
            else:
                n += 1
        n = round(n)
        if n > 0:
            profile[gene] = n
        else:
            todel.append(gene)

    # delete zero values
    for gene in todel:
        del profile[gene]


class Tests(unittest.TestCase):
    def setUp(self):
        self.sam = (
            '@HD	VN:1.0	SO:unsorted\n'
            'S1	77	NC_123456	26	0	100M	*	0	0	*	*\n'
            'S1	141	NC_123456	151	0	80M	*	0	0	*	*\n'
            'S2	0	NC_789012	186	0	50M5I20M5D20M	*	0	0	*	*\n'
            'S2	16	*	0	0	*	*	0	0	*	*\n')

    def test_read_gene_table(self):
        tbl = ('## GCF_000123456',
               '# NC_123456',
               '1	5	384',
               '2	410	933',
               '# NC_789012',
               '1	912	638',
               '2	529	75')
        obs = read_gene_table(tbl)
        exp = {'NC_123456': [
            (5,   True, True, '1'), (384, False, True, '1'),
            (410, True, True, '2'), (933, False, True, '2')],
               'NC_789012': [
            (638, True, True, '1'), (912, False, True, '1'),
            (75,  True, True, '2'), (529, False, True, '2')]}
        self.assertDictEqual(obs, exp)

    def test_parse_b6o_line(self):
        rids, lenmap, locmap = [], {}, {}
        b6o = ('S1/1	NC_123456	100	100	0	0	1	100	225	324	1.2e-30	345',
               'S1/2	NC_123456	95	98	2	1	2	99	608	708	3.4e-20	270')

        # first line
        parse_b6o_line(b6o[0], rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1'])
        self.assertDictEqual(lenmap, {'NC_123456': {0: 100}})
        self.assertDictEqual(locmap, {'NC_123456': [
            (225, True, False, 0), (324, False, False, 0)]})

        # second line
        parse_b6o_line(b6o[1], rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1', 'S1/2'])
        self.assertDictEqual(lenmap, {'NC_123456': {0: 100, 1: 98}})
        self.assertDictEqual(locmap, {'NC_123456': [
            (225, True, False, 0), (324, False, False, 0),
            (608, True, False, 1), (708, False, False, 1)]})

    def test_parse_sam_line(self):
        sam = self.sam.splitlines()
        rids, lenmap, locmap = [], {}, {}

        # header
        parse_sam_line(sam[0], rids, lenmap, locmap)
        self.assertFalse(len(rids) + len(lenmap) + len(locmap))

        # normal, fully-aligned, forward strand
        parse_sam_line(sam[1], rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1'])
        self.assertDictEqual(lenmap, {'NC_123456': {0: 100}})
        self.assertDictEqual(locmap, {'NC_123456': [
            (26, True, False, 0), (125, False, False, 0)]})

        # shortened, reverse strand
        parse_sam_line(sam[2], rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1', 'S1/2'])
        self.assertDictEqual(lenmap, {'NC_123456': {
            0: 100, 1: 80}})
        self.assertDictEqual(locmap, {'NC_123456': [
            (26,  True, False, 0), (125, False, False, 0),
            (151, True, False, 1), (230, False, False, 1)]})

        # not perfectly aligned, unpaired
        parse_sam_line(sam[3], rids, lenmap, locmap)
        self.assertListEqual(rids, ['S1/1', 'S1/2', 'S2'])
        self.assertEqual(lenmap['NC_789012'][2], 90)
        self.assertTupleEqual(locmap['NC_789012'][0], (186,  True, False, 2))
        self.assertTupleEqual(locmap['NC_789012'][1], (280, False, False, 2))

        # not aligned
        parse_sam_line(sam[4], rids, lenmap, locmap)
        self.assertEqual(len(rids), 3)
        self.assertEqual(len(lenmap['NC_789012']), 1)
        self.assertEqual(len(locmap['NC_789012']), 2)

    def test_cigar_to_lens(self):
        self.assertTupleEqual(cigar_to_lens('150M'), (150, 150))
        self.assertTupleEqual(cigar_to_lens('3M1I3M1D5M'), (11, 12))

        with self.assertRaises(ValueError) as context:
            cigar_to_lens('*')
        msg = 'Missing CIGAR string.'
        self.assertEqual(str(context.exception), msg)

    def test_match_read_gene(self):
        # illustration of map

        #  reads:             ---------r1---------
        #                           ---------r2---------
        #                               ---------r3---------
        #                                 ---------r4---------
        #                                         ---------r5---------
        #                                         ---------r6---------
        # genome:  1 ====>>>>>>>>>>>g1>>>>>>>>>>>>===>>>>>>>>>>>>>>g2>> 50

        #  reads:             ---------r7---------
        #                          ---------r8---------
        #                                           ------r9------
        # genome: 51 >>>>>>>>>>>===>>>>>>>>>>>>>>g3>>>>>>>>>>>>>======= 100
        # --------------------------------------------------

        # gene table
        genes = [('g1',  5, 29),
                 ('g2', 33, 61),
                 ('g3', 65, 94)]
        # read map
        reads = [('r1', 10, 29),
                 ('r2', 16, 35),
                 ('r3', 20, 39),
                 ('r4', 22, 41),
                 ('r5', 30, 49),
                 ('r6', 30, 49),  # identical
                 ('r7', 60, 79),
                 ('r8', 65, 84),
                 ('r9', 82, 95)]  # shorter

        # length map
        lens = {'r{}'.format(i): 20 for i in range(1, 10)}

        # flatten lists
        genes = [x for id_, start, end in genes for x in
                 ((start, True, True, id_),
                  (end,  False, True, id_))]
        reads = [x for id_, start, end in reads for x in
                 ((start, True, False, id_),
                  (end,  False, False, id_))]

        queue = sorted(genes + reads, key=lambda x: x[0])

        # default (threshold = 80%)
        obs = match_read_gene(queue, lens)
        exp = {'g1': 1, 'g2': 2, 'g3': 1}
        self.assertDictEqual(obs, exp)

        # return read map instead of counts
        obs = match_read_gene(queue, lens, th=0.8, ismap=True)
        exp = {'r1': {'g1'},
               'r5': {'g2'},
               'r6': {'g2'},
               'r8': {'g3'}}
        self.assertDictEqual(obs, exp)

        # threashold = 50%
        obs = match_read_gene(queue, lens, th=0.5, ismap=True)
        exp = {'r1': {'g1'},
               'r2': {'g1'},
               'r3': {'g1'},
               'r5': {'g2'},
               'r6': {'g2'},
               'r7': {'g3'},
               'r8': {'g3'},
               'r9': {'g3'}}
        self.assertDictEqual(obs, exp)

    def test_readmap_to_profile(self):
        # no ambiguity
        rids = ['R1', 'R2', 'R3', 'R4']
        profile = {'G1': [0, 1, 2], 'G2': [1, 3]}
        readmap_to_profile(profile, rids, True)
        self.assertDictEqual(profile, {'G1': 3, 'G2': 2})

        # R1 occurs 3 times; R2 occurs 2 times
        rids = ['R1', 'R2', 'R3', 'R1', 'R1', 'R2']

        # drop R1 and R2
        profile = {'G1': [0, 1], 'G2': [1, 2, 3]}
        readmap_to_profile(profile, rids, False)
        self.assertDictEqual(profile, {'G2': 1})

        # normalize R1 by 3, R2 by 2
        profile = {'G1': [0, 1], 'G2': [1, 2, 3]}
        readmap_to_profile(profile, rids, True)
        self.assertDictEqual(profile, {'G1': 1, 'G2': 2})


if __name__ == "__main__":
    main()
