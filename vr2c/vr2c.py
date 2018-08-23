import sys
import argparse
import logging

from tqdm import tqdm
import pandas as pd
import pysam

from vr2c.plot import plot_alignment

logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


def get_args():
    parser = argparse.ArgumentParser(
        description='plot r2c and c2g in one coordinate')
    parser.add_argument(
        '-c', '--contigs-to-genome', type=str, required=True,
        help='input contig-to-genome alignment BAM file'
    )
    parser.add_argument(
        '-r', '--reads-to-contigs', type=str, required=True,
        help='input read-to-contig alignment BAM file'
    )
    parser.add_argument(
        '-t', '--contig-id', type=str, required=True,
        help=('contig id, e.g. A1.J70460')
    )
    parser.add_argument(
        '-s', '--seqname', type=str, required=True,
        help=('e.g. chr12. This is required because one contig could be '
              'aligned to two chromosomes, e.g. hardclipping')
    )
    parser.add_argument(
        '-p', '--positions', type=int, nargs='+', required=True,
        help=('plot vertical lines for intended positions. e.g. cleavage site')
    )
    parser.add_argument(
        '--figure-width', type=int, default=16,
        help=('figure width, default to 16 based on experience')
    )
    parser.add_argument(
        '--seqname-beg-end', type=int, nargs=2, default=None,
        help=('if specified, the beginning and ending of the plot '
              'on the chromosome are enforced, e.g. [25357088, 25357993]')
    )
    parser.add_argument(
        '-o', '--output', type=str, default=None,
        help=('default to be <contig_id>.png')
    )
    return parser.parse_args()


def fetch_contig(c2g_bam, contig_id, seqname):
    for contig in tqdm(c2g_bam):
        if (
                contig.query_name == contig_id and
                contig.reference_name == seqname
        ):
            return contig


def log_contig_info(contig):
    for attr in [
            'query_name', 'is_reverse', 'cigarstring',
            'reference_start', 'reference_end',
    ]:
        logging.info('contig.{0}: {1}'.format(attr, getattr(contig, attr)))


def extract_read_info(contig, read):
    if contig.is_reverse:
        contig_len = contig.infer_query_length(always=True)
        read_info = [
            contig_len - read.reference_end,
            contig_len - read.reference_start,
            not read.is_reverse,
            f'rev({read.cigarstring})',
            tuple(reversed(read.cigartuples))
        ]
    else:
        contig_len = contig.infer_query_length(always=True)
        read_info = [
            read.reference_start,
            read.reference_end,
            read.is_reverse,
            f'{read.cigarstring}',
            read.cigartuples
        ]
    return read_info


def collect_reads(contig, r2c_bam, max_num=100):
    reads = r2c_bam.fetch(contig.query_name)
    reads_info = []
    for k, rd in enumerate(reads):
        if rd.is_unmapped:
            continue

        # TODO: add mechnism for selecting particular types of reads
        rd_info = extract_read_info(contig, rd)
        reads_info.append(rd_info)

        if len(reads_info) == max_num:
            break

    df_reads = pd.DataFrame(
        reads_info, columns=[
            'reference_start',
            'reference_end',
            'is_reverse',
            'cigarstring',
            'cigartuples'
        ])
    df_reads = df_reads.sort_values('reference_start').reset_index(drop=True)
    return df_reads


def gen_output(contig, clvs):
    return '{0}_{1}_{2}.png'.format(
        contig.reference_name,
        contig.query_name,
        '-'.join([str(_) for _ in sorted(clvs)]),
    )


def main():
    args = get_args()

    c2g_bam_file = args.contigs_to_genome
    r2c_bam_file = args.reads_to_contigs
    logging.info('c2g_bam_file: {0}'.format(c2g_bam_file))
    logging.info('r2g_bam_file: {0}'.format(r2c_bam_file))

    c2g_bam = pysam.AlignmentFile(c2g_bam_file)

    contig = fetch_contig(c2g_bam, args.contig_id, args.seqname)
    if contig is None:
        logging.info(f'The specified contig id: "{args.contig_id}" is NOT found in {args.contigs_to_genome}')
        sys.exit(1)

    log_contig_info(contig)

    r2c_bam = pysam.AlignmentFile(r2c_bam_file)
    df_reads = collect_reads(contig, r2c_bam)
    logging.info('collected {0} reads aligned to for {1} '
                 '(use --plot-all-reads to plot all reads)'.format(
                     df_reads.shape[0], contig.query_name))

    output = args.output
    if output is None:
        output = gen_output(contig, args.positions)

    plot_alignment(
        contig, df_reads, args.positions, output,
        args.figure_width,
        args.seqname_beg_end,
    )


if __name__ == "__main__":
    main()
