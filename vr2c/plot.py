import logging

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from vr2c import visaln         # noqa: E402


logger = logging.getLogger(__name__)


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


def plot_alignment(contig, df_reads, clvs, output, fig_width,
                   seqname_beg_end=None):
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
