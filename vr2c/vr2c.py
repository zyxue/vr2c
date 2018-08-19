import sys
import argparse
import logging

from tqdm import tqdm
import pandas as pd
import pysam
import matplotlib.pyplot as plt

from kleat.evidence import bridge
from kleat.visaln import visaln


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
        '-l', '--clvs', type=int, nargs='+', required=True,
        help=('the predicted cleavage site')
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
        '--plot-all-reads', action='store_true',
        help=('By default, only bridge reads are collected and plotted. '
              'If this argument is specified, all reads will be plotted '
              'the output figure could be large if many reads are aligned '
              'to the given contig. Unmapped reads are always ignored')
    )
    parser.add_argument(
        '-o', '--output', type=str, default=None,
        help=('default to be <contig_id>.png')
    )
    return parser.parse_args()


def fetch_contig(c2g_bam, contig_id, seqname):
    for contig in tqdm(c2g_bam):
        if (contig.query_name == contig_id and
            contig.reference_name == seqname):
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


def collect_reads(contig, r2c_bam, plot_all_reads, max_num=100):
    reads = r2c_bam.fetch(contig.query_name)
    reads_info = []
    for k, rd in enumerate(reads):
        if rd.is_unmapped:
            continue

        if plot_all_reads or bridge.is_a_bridge_read(rd):
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


def calc_figure_height(df_reads):
    # +1 for the contig
    num_reads = df_reads.shape[0] + 1
    height_per_read = 0.3       # arbitray

    fig_height = num_reads * height_per_read
    return fig_height


def prepare_fig_axes(num_skips, df_xlims, fig_width, fig_height):
    nrows = 1
    ncols = num_skips + 1
    width_ratios = df_xlims.span.values

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(fig_width, fig_height), sharey=True,
        gridspec_kw={
            'width_ratios': width_ratios,
            'hspace': 0.5
        }
    )

    if ncols == 1:
        axes = [axes]
    else:
        axes = axes.ravel()

    return fig, axes


def plot_a_read(ax, contig, read_row, cy):
    cx = visaln.get_abs_start(read_row)
    cx = visaln.convert_contig2genome_coord(cx, contig)
    cy = cy
    visaln.draw_alignment(ax, cx, cy, read_row)


def plot_clv(ax, clv, ylim):
    ax.plot([clv, clv], ylim, color='black')
    ax.set_ylim(ylim)


def format_ax(ax, ax_idx):
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)

    if ax_idx > 0:
        ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if ax_idx == 0:
        ax.set_ylabel('# reads')
        ax.tick_params(top=False, right=False)
    if ax_idx > 0:
        ax.tick_params(top=False, right=False, left=False)

    ax.yaxis.grid()
    ax.xaxis.grid()


def add_figure_title(figure, contig, clvs):
    rev = 'reverse' if contig.is_reverse else 'forward'
    title = '{0}, {1}, {2}, ({3}), clv@{4}'.format(
        contig.query_name,
        contig.cigarstring,
        contig.reference_name,
        rev,
        ','.join([str(_) for _ in sorted(clvs)]),
    )
    figure.suptitle(title, fontsize=15, y=1.02, verticalalignment='bottom')


def plot_alignment(contig, df_reads, clvs, output, fig_width, seqname_beg_end=None):
    num_skips = visaln.calc_num_skips(contig)
    df_xlims = visaln.calc_xlim_pairs(contig, clvs)

    if seqname_beg_end is not None:
        min_idx = df_xlims.xmin.idxmin()
        max_idx = df_xlims.xmax.idxmax()
        df_xlims.loc[min_idx, 'xmin'] = seqname_beg_end[0]
        df_xlims.loc[max_idx, 'xmax'] = seqname_beg_end[1]

    fig_height = calc_figure_height(df_reads)
    fig, axes = prepare_fig_axes(
        num_skips, df_xlims, fig_width, fig_height)

    for idx, row in df_xlims.iterrows():
        ax = axes[idx]
        for k, read_row in df_reads.iterrows():
            cy = k + 1
            plot_a_read(ax, contig, read_row, cy)

        cx = visaln.get_abs_start(contig)
        visaln.draw_alignment(ax, cx, cy=0, aln=contig)

        for _clv in clvs:
            if row.xmin <= _clv <= row.xmax:
                ylim = [-1, df_reads.shape[0] + 1]
                plot_clv(ax, _clv, ylim)

        ax.set_xlim(row.xmin, row.xmax)
        format_ax(ax, idx)

    add_figure_title(fig, contig, clvs)

    # useful for zoom in a particular section, left here for later use
    # xlim = ax.get_xlim()
    # ax.set_xlim(xlim[0] + 2450, xlim[1] + 50)
    plt.savefig(output, bbox_inches='tight')


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
        logging.info('The specified contig id: {0} is NOT found in {1}'.format(args.contigs_to_genome))
        sys.exit(1)

    log_contig_info(contig)

    r2c_bam = pysam.AlignmentFile(r2c_bam_file)
    df_reads = collect_reads(contig, r2c_bam, args.plot_all_reads)
    logging.info('collected {0} reads aligned to for {1} '
                 '(use --plot-all-reads to plot all reads)'.format(
                     df_reads.shape[0], contig.query_name))

    output = args.output
    if output is None:
        output = gen_output(contig, args.clvs)

    plot_alignment(contig, df_reads, args.clvs, output,
                   args.figure_width,
                   args.seqname_beg_end
    )


if __name__ == "__main__":
    main()
