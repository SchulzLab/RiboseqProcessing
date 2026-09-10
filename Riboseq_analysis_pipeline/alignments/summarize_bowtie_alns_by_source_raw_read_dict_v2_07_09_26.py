import sys
import os
import os.path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip


def assign_trans_IDs_to_source(ref_fastas):
    """
    Assign the IDs in reference fasta files to the source

    Return a dict with the sources as keys and the transcript IDs as values

    """
    ref_fastas = ref_fastas.split(',')
    reference_dict = {}
    for fasta in ref_fastas:
        ref_seqs = SeqIO.parse(fasta, 'fasta')
        source = os.path.basename(fasta)
        source = source.rsplit('.', 1)[0]
        for seq in ref_seqs:
            reference_dict[seq.id] = source
    return reference_dict


def plot_summary(summary_mapping_stats_df, outdir, raw=False, figsize=None,
                 inches_per_sample=0.45, min_label_pct=2.0):
    """
        plot summary statistics of the mapping to Ignolia like reference

        barplot colored by reference type with the sample on the x-axis

        figsize, font sizes and label density scale with the number of
        samples; pass figsize=(w, h) to override the automatic size
    """

    def add_labels(x, y, bottom):
        for i in range(len(x)):
            # CHANGED: skip segments too thin to hold a readable label
            if y.iloc[i] < min_label_pct or not annotate:
                continue
            if bottom is not None:
                # CHANGED: / instead of //, fontsize added
                plt.text(i, y.iloc[i] / 2 + bottom.iloc[i], round(y.iloc[i], 2),
                         ha='center', va='center', fontsize=label_fs)
                # Placing text at half the bar height
            else:
                plt.text(i, y.iloc[i] / 2, round(y.iloc[i], 2),
                         ha='center', va='center', fontsize=label_fs)

    summary_mapping_stats_df_t = summary_mapping_stats_df.transpose()
    summary_mapping_stats_df_t = summary_mapping_stats_df_t * 100
    summary_mapping_stats_df_t = summary_mapping_stats_df_t[summary_mapping_stats_df_t.index.str.contains(
        'percentage')]
    summary_mapping_stats_df_t = summary_mapping_stats_df_t.sort_index()

    # rename columns for better readability
    summary_mapping_stats_df_t = summary_mapping_stats_df_t.rename(columns={'Ens_Gencode_lncRNA_ncRNA': 'ncRNA',
                                                                            'Homo_sapiens.GRCh38.dna.chromosome.MT': 'MT',
                                                                            'MANE.GRCh38.v0.95.select_ensembl_rna': 'MANE_mRNA',
                                                                            'hg38-tRNAs': 'tRNA',
                                                                            'rRNA_ref_NCBI_Ens': 'rRNA'})
    # Make one pd.Series per source
    mRNA_df = summary_mapping_stats_df_t['MANE_mRNA']
    rRNA_df = summary_mapping_stats_df_t['rRNA']
    ncRNA_df = summary_mapping_stats_df_t['ncRNA']
    mt_df = summary_mapping_stats_df_t['MT']
    tRNA_df = summary_mapping_stats_df_t['tRNA']

    sample_names = ['_'.join(x[2:9]) if 'RR' in x else '_'.join(
        x[2:8]) for x in mRNA_df.index.str.split('_')]

    x = np.arange(len(sample_names))  # Numeric x-axis positions

    # CHANGED: scale the canvas and the type with the number of samples
    n = len(sample_names)
    if figsize is None:
        max_label_len = max([len(s) for s in sample_names] + [10])
        width = float(np.clip(2.5 + inches_per_sample * n, 8, 60))
        height = 6 + 0.07 * max_label_len   # room for the rotated ticklabels
        figsize = (width, height)
    # font sizes follow the width of one bar in points, not the sample count:
    # while the figure is still growing each bar keeps the same physical width,
    # so the text only has to shrink once the width hits its 60 inch cap
    bar_pts = 72.0 * figsize[0] * 0.75 / max(n, 1)
    label_fs = float(np.clip(bar_pts / 3.2, 4, 9))    # '12.3' is ~3.2 pts/size
    # rotated 90, so ~1 em wide
    tick_fs = float(np.clip(bar_pts * 0.8, 4, 10))
    annotate = bar_pts / 3.2 >= 4.5    # below this nothing is legible anyway
    bar_width = 0.8 if n > 3 else 0.5

    plt.figure(figsize=figsize)
    plt.bar(sample_names, ncRNA_df, width=bar_width,
            color='#40A3CD', label='ncRNA')
    plt.bar(x, tRNA_df, bottom=ncRNA_df, width=bar_width,
            color='#019c91', label='tRNA')
    plt.bar(x, rRNA_df, bottom=ncRNA_df + tRNA_df, width=bar_width,
            color='#b79165', label='rRNA')
    plt.bar(x, mt_df, bottom=ncRNA_df + tRNA_df + rRNA_df, width=bar_width,
            color='orange', label='MT')
    plt.bar(x, mRNA_df, bottom=ncRNA_df + tRNA_df +
            rRNA_df + mt_df, width=bar_width, color='#c73832', label='mRNA')

    # add percentage labels
    add_labels(x, mRNA_df, ncRNA_df + tRNA_df + rRNA_df + mt_df)
    add_labels(x, mt_df, ncRNA_df + tRNA_df + rRNA_df)
    add_labels(x, rRNA_df, ncRNA_df + tRNA_df)
    add_labels(x, tRNA_df, ncRNA_df)
    add_labels(x, ncRNA_df, None)

    if raw:
        unmapped_filtered_df = summary_mapping_stats_df_t['filtered_unmapped']
        plt.bar(x, unmapped_filtered_df, bottom=ncRNA_df +
                tRNA_df + rRNA_df + mt_df + mRNA_df, width=bar_width,
                color='#FDDA0D', label='filtered_unmapped')
        add_labels(x, unmapped_filtered_df, ncRNA_df +
                   tRNA_df + rRNA_df + mt_df + mRNA_df)

    plt.ylabel('Percentage')
    # CHANGED: scaled tick fontsize, xlim so bars are not squashed to the edge
    plt.xticks(x, sample_names, rotation=90, fontsize=tick_fs)
    plt.xlim(-0.6, n - 0.4)
    plt.legend(title='RNA Type', loc='upper left', bbox_to_anchor=(1, 1))
    plt.title('Stacked Mapping Percentages per Sample')

    if raw:
        plt.savefig(os.path.join(
            outdir, 'barplot_mapping_stats_per_sample_raw.png'), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(
            outdir, 'barplot_mapping_stats_per_sample.png'), bbox_inches='tight')
    plt.close()


def summarize_mapping_statistics(reference_dict, idx_stats, out_dir, raw_read_path='', raw=False):
    """
    Obtain mapping percentages of references from samtools idxstats.

    Each mutlifasta file given as input (list with commaseparated fasta files) is a different source.
    Summarize the alignments by source from mapping stats of single transcripts.

    """
    print('raw', raw)
    if raw:
        def count_reads(raw_read_path):
            open_func = gzip.open if raw_read_path.endswith('.gz') else open
            with open_func(raw_read_path, 'rt') as f:
                nr_lines = sum(1 for _ in f)
                nr_reads = nr_lines / 4
                assert (nr_lines % 4 == 0)
                print(nr_reads)
                return nr_reads

        raw_read_counts_dict = {}

        for f in os.listdir(raw_read_path):
            if f.endswith('.fastq') or f.endswith('.fastq.gz') or f.endswith('.fq') or f.endswith('.fq.gz'):
                sample = f.split('.')[0]
                print('sample_name_dict', sample)
                raw_read_counts_dict[sample] = count_reads(
                    os.path.join(raw_read_path, f))

        print(raw_read_counts_dict)

    summary_mapping_stats_df = None
    for idx_stat in idx_stats:
        idx_stat_df = pd.read_csv(idx_stat,
                                  sep='\t',
                                  index_col=False,
                                  header=None,
                                  names=['tID', 'length', 'mapped', 'unmapped'])
        idx_name = os.path.basename(idx_stat).split('.')[0]
        print('idx_name', idx_name)
        # skip the last row that is just a star and two dots
        idx_stat_df = idx_stat_df.iloc[:-1, :]
        # map the tids to their source
        idx_stat_df['source'] = idx_stat_df['tID'].map(reference_dict)
        if raw:
            nr_unmapped = raw_read_counts_dict[idx_name] - \
                sum(idx_stat_df['mapped'])
            add_unmapped_row = {'tID': 'unmapped', 'length': 0,
                                'mapped': nr_unmapped, 'unmapped': 0, 'source': 'filtered_unmapped'}
            idx_stat_df = pd.concat(
                [idx_stat_df, pd.DataFrame([add_unmapped_row])], ignore_index=True)

        # get the number of mapping reads per tID
        idx_stat_df_by_source = idx_stat_df[[
            'source', 'mapped']].groupby('source').agg('sum')

        # calculate the mapping percentage per source of all mapping reads
        # if raw equals to true than all the unmapped and filtered reads are counted as well
        idx_stat_df_by_source[f'{idx_name}_percentage'] = idx_stat_df_by_source['mapped'] / \
            sum(idx_stat_df_by_source['mapped'])

        if summary_mapping_stats_df is not None:
            summary_mapping_stats_df[f'{idx_name}_percentage'] = idx_stat_df_by_source[f'{idx_name}_percentage']
            summary_mapping_stats_df[f'{idx_name}_mapped'] = idx_stat_df_by_source['mapped']
        else:
            summary_mapping_stats_df = idx_stat_df_by_source.rename(
                columns={'mapped': f'{idx_name}_mapped'})

        # CHANGED: removed the per-sample pie chart that was written here

    plot_summary(summary_mapping_stats_df, out_dir, raw)

    return summary_mapping_stats_df


def main():
    ref_fastas = sys.argv[1]
    idx_dir = sys.argv[2]
    raw_read_path = sys.argv[3]

    print(raw_read_path)

    idx_stats = []

    for filename in os.listdir(idx_dir):
        if filename.endswith('idxstats.out'):
            idx_stats.append(os.path.join(idx_dir, filename))

    out_dir = os.path.dirname(idx_stats[0])

    reference_dict = assign_trans_IDs_to_source(ref_fastas)

    summarized_map_stats_df = summarize_mapping_statistics(
        reference_dict, idx_stats, out_dir=out_dir)
    summarized_map_stats_df = summarize_mapping_statistics(
        reference_dict, idx_stats, out_dir, raw_read_path, raw=True)

    summarized_map_stats_df.to_csv(f'{out_dir}/summarized_map_stats_df.csv')

    sys.stdout.flush()
    sys.stderr.flush()


if __name__ == '__main__':
    main()
