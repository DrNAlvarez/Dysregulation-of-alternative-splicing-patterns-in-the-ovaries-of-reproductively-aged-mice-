# module for Figure 1
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import viridis
import matplotlib.gridspec as gridspec

# Create the figure
fig = plt.figure(figsize=(8, 11))

def mean_and_sem_numpy(data):
    """Calculate the mean and the Standard Error of the Mean (SEM) using NumPy only.
    
    Parameters:
    - data (array-like): The data to compute the statistics on.
    
    Returns:
    - tuple: mean, SEM
    """
    mean = np.mean(data)
    sem = np.std(data, ddof=1) / np.sqrt(len(data))  # ddof=1 for sample standard deviation
    return mean, sem

def seq_output(file_paths, table=True, ax=None, base_size=12):
    # Data storage for plotting
    total_counts = []
    pass_counts = []
    fail_counts = []
    sample_labels = []

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    
    # Color mapping for viridis
    colors = viridis(np.linspace(0, 1, len(file_paths)))
    
    #sample_names = ['6 weeks', '6 weeks', '6 weeks']
    #group_names = ['PMSG','PMSG','PMSG']
    
    
    sample_names = ['6 weeks', '6 weeks', '6 weeks',
                    '6 weeks', '6 weeks', '6 weeks',
                    '14 months', '14 months', '14 months',
                    '14 months', '14 months', '14 months']
    group_names = ['PMSG','PMSG','PMSG', 'PMSG+hCG','PMSG+hCG','PMSG+hCG', 'PMSG','PMSG','PMSG', 'PMSG+hCG','PMSG+hCG','PMSG+hCG']

    # Loop through each file to collect data
    for i, file_path in enumerate(file_paths):
        # Load the sequence summary data
        df = pd.read_csv(file_path, delimiter='\t')
        
        # Calculate pass and fail counts
        pass_count = df['passes_filtering'].sum()
        fail_count = len(df) - pass_count
        total_count = pass_count + fail_count
        
        # Store data
        total_counts.append(total_count)
        pass_counts.append(pass_count)
        fail_counts.append(fail_count)
        sample_labels.append(sample_names[i])

        
    fail_counts = np.log(fail_counts)
    pass_counts = np.log(pass_counts)
    total_counts = np.log(total_counts)

    mean_total, sem_total = mean_and_sem_numpy(np.exp(total_counts))
    summary_df = pd.DataFrame({
        'Mean Total': [mean_total],
        'SEM Total': [sem_total]
    })

    #print(summary_df)
    
    # Save summary statistics to CSV if table=True
    if table:
        summary_df.to_csv('figure_files/Fig1A_summary_statistics.csv', index=False)
    # Dynamically scale font sizes
    fontsize = scale_fontsize(fig, base_size)
    # Plot passed and failed counts as stacked bars
    ax.bar(np.arange(len(file_paths)), total_counts, label='Total (Passed + Failed)', color=colors, alpha=0.3)
    ax.bar(np.arange(len(file_paths)), pass_counts, label='Passed', color=colors, alpha=0.7)
    
    # Add total counts as text above bars
    for i, (pass_count, fail_count) in enumerate(zip(pass_counts, fail_counts)):
        total_count = np.exp(pass_count) + np.exp(fail_count)
        ax.text(i, 1 + pass_count, str(np.round(total_count, decimals=0)), ha='center', rotation=90,size=fontsize*.75)
    
    ax.set_ylim(0, 30)
    ax.set_xticks(np.arange(len(file_paths)))
    ax.set_xticklabels(sample_names, ha='center')
    ax.tick_params(axis='x', which='major', labelsize=fontsize*.75,rotation=90)
    
    # Adding the secondary x-axis for group labels
    num_samples_per_group = 3
    label_positions = [i + 1 for i in range(0, len(sample_names), num_samples_per_group)]
    tick_mark_positions = [i for i in range(0, len(sample_names), num_samples_per_group)] + [i + 2 for i in range(0, len(sample_names), num_samples_per_group)]
        
    secaxA = ax.secondary_xaxis('bottom')

    secaxA.set_xticks(np.arange(len(file_paths)))

    secaxA.set_xticklabels(group_names,ha='center')
    secaxA.tick_params(axis='x', which='major', labelsize=fontsize*.75,rotation=90)
    secaxA.spines['bottom'].set_position(('outward', 45))
    secaxA.set_xlim(0, 11)
    
    
    ax.set_ylabel('RNA Molecules (Log10)', fontsize=fontsize)
    ax.set_title('Total RNA Molecules by Sequencing Run', fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)

def plot_seq_length(file_paths, table=True, ax=None, base_size=12):
    # Color mapping for viridis
    colors = viridis(np.linspace(0, 1, len(file_paths)))
    
    plt.rcParams['font.family'] = 'sans-serif'

    # Prioritize 'Arial' within the 'sans-serif' font list
    plt.rcParams['font.sans-serif'] = ['Arial'] + plt.rcParams['font.sans-serif']

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    
    # Dynamically scale font sizes
    fontsize = scale_fontsize(fig, base_size)

    # Panel B - Line frequency distribution using the previously defined logic
    # Loop through each file and plot the read length distribution on ax1
    #sample_names = ['6 weeks\nPMSG','6 weeks\nPMSG','6 weeks\nPMSG']
    sample_names = ['6 weeks\nPMSG','6 weeks\nPMSG','6 weeks\nPMSG',
                    '6 weeks\nPMSG+hCG','6 weeks\nPMSG+hCG','6 weeks\nPMSG+hCG',
                   '14 months\nPMSG','14 months\nPMSG','14 months\nPMSG',
                   '14 months\nPMSG+hCG','14 months\nPMSG+hCG','14 months\nPMSG+hCG']
    read_lengths = []
    for i, file_path in enumerate(file_paths):
        # Load the sequence summary data
        df = pd.read_csv(file_path, delimiter='\t')
    
        # Filter the data to only include reads that passed filtering
        filtered_data = df[df['sequence_length_template'] <= 10000]
        filtered_data = filtered_data[filtered_data['passes_filtering'] == 1]
        
        # Extract the sequence lengths of the reads that passed
        read_lengths_passed = filtered_data['sequence_length_template']
        read_lengths.append(read_lengths_passed)
        
        # Calculate the histogram
        hist, bin_edges = np.histogram(read_lengths_passed, bins=50)
        
        # Calculate the center of each bin
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # Plot as a line plot on ax
        ax.plot(bin_centers, hist, label=sample_names[i], color=colors[i], alpha=0.9)
    
    # Calculate mean and SEM for read lengths
    all_read_lengths = np.concatenate(read_lengths)
    mean_length, sem_length = mean_and_sem_numpy(all_read_lengths)
    summary_df = pd.DataFrame({
        'Mean Length': [mean_length],
        'SEM Length': [sem_length]
    })

    #print(summary_df)
    
    # Save summary statistics to CSV if table=True
    if table:
        summary_df.to_csv('figure_files/Fig1B_length_summary_statistics.csv', index=False)
    
    
    ax.set_title('Distribution of RNA Lengths', fontsize=fontsize)
    ax.set_xlabel('RNA Length (nt)', fontsize=fontsize)
    ax.set_ylabel('Number of Bases (Log10)', fontsize=fontsize)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.tick_params(axis='both', which='major', labelsize=fontsize*0.75)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3, frameon=False, fontsize=fontsize*.75)

    return read_lengths
    
# Load and annotate the transcript data
def annotate_transcripts(gtf_file_path: str,
                         ensembl_ids_tsv: str,
                         translation_tsv: str) -> pd.DataFrame:
    gtf = pd.read_csv(gtf_file_path, sep='\t', low_memory=False, header=0, dtype=str)
    tx = gtf[gtf['feature'] == 'transcript'][['transcript_id', 'gene_id']].drop_duplicates()

    ensembl = pd.read_csv(ensembl_ids_tsv, sep='\t', low_memory=False, dtype=str)
    ensembl = ensembl.rename(columns={
        'Transcript stable ID': 'transcript_id',
        'Gene stable ID': 'gene_id',
        'Transcript type': 'transcript_biotype',
        'Gene type': 'gene_biotype'
    })[['transcript_id', 'gene_id', 'transcript_biotype', 'gene_biotype']]

    df = tx.merge(ensembl, on=['transcript_id', 'gene_id'], how='left')
    df['transcript_biotype'] = df['transcript_biotype'].fillna('unannotated')
    df['gene_biotype'] = df['gene_biotype'].fillna('unannotated')

    translation = pd.read_csv(translation_tsv, sep='\t', low_memory=False, dtype=str, usecols=['transcript_id'])
    coding_set = set(translation['transcript_id'])

    df['bioseq2seq_biotype'] = np.where(
        df['transcript_id'].isin(coding_set),
        'coding',
        'noncoding'
    )

    return df[['transcript_id', 'gene_id', 'transcript_biotype', 'gene_biotype', 'bioseq2seq_biotype']]

def scale_fontsize(fig, base_size=12):
    width, height = fig.get_size_inches()
    return base_size * min(width / 8, height / 11)


def collapse_categories(proportions, cutoff):
    large = proportions[proportions >= cutoff]
    small = proportions[proportions < cutoff]
    collapsed = large.copy()
    if not small.empty:
        collapsed['Other'] = small.sum()
    return collapsed

def clean_label(x):
    if x == 'Other':
        return 'Other'
    parts = x.split('_')
    cleaned = []
    for part in parts:
        if 'RNA' in part:
            cleaned.append(part)  # leave casing untouched
        else:
            cleaned.append(part.capitalize())
    return ' '.join(cleaned)

def plot_biotype(annotation_df, table=True, ax1=None, ax2=None, base_size=12, cutoff=0.01):

    fontsize = scale_fontsize(fig, base_size)

    # Transcript biotype collapsing and cleanup
    transcript_counts = annotation_df['transcript_biotype'].value_counts()
    transcript_props = transcript_counts / transcript_counts.sum()
    transcript_collapsed = collapse_categories(transcript_props, cutoff)
    transcript_counts_collapsed = (transcript_collapsed * transcript_counts.sum()).astype(int)
    transcript_labels_clean = [clean_label(x) for x in transcript_collapsed.index]

    # Bioseq2seq collapsing and cleanup
    bioseq_counts = annotation_df['bioseq2seq_biotype'].value_counts()
    bioseq_props = bioseq_counts / bioseq_counts.sum()
    bioseq_collapsed = collapse_categories(bioseq_props, cutoff)
    bioseq_counts_collapsed = (bioseq_collapsed * bioseq_counts.sum()).astype(int)
    bioseq_labels_clean = [clean_label(x) for x in bioseq_collapsed.index]

    # Percentages
    transcript_percentages = [f"{(v * 100):.1f}%" for v in transcript_collapsed]
    bioseq_percentages = [f"{(v * 100):.1f}%" for v in bioseq_collapsed]

    # Output DataFrames
    bio_df = pd.DataFrame({
        'Biotype': transcript_labels_clean,
        'Count': transcript_counts_collapsed.values,
        'Percentage': transcript_percentages
    })

    bio2_df = pd.DataFrame({
        'Biotype': bioseq_labels_clean,
        'Count': bioseq_counts_collapsed.values,
        'Percentage': bioseq_percentages
    })

    if table:
        bio_df.to_csv('figure_files/Fig1C_biotype_percentages_cutoff.csv', index=False)
        bio2_df.to_csv('figure_files/Fig1C_bioseq2seq_percentages_cutoff.csv', index=False)

    # Colors
    colors1 = viridis(np.linspace(0, 1, len(transcript_counts_collapsed)))
    colors2 = viridis(np.linspace(0, 1, len(bioseq_counts_collapsed)))

    # Pie plots
    ax1.pie(transcript_counts_collapsed, labels=None, startangle=90, counterclock=False, colors=colors1, wedgeprops={'width': 0.4})
    ax1.legend(loc='upper center', bbox_to_anchor=(1, 0),
               labels=[f"{x}: {y}" for x, y in zip(transcript_labels_clean, transcript_percentages)],
               prop={'size': fontsize * 0.75}, frameon=False)
    ax1.set_title(f"Transcript Biotype\n(Total Transcripts: {transcript_counts.sum()})", size=fontsize)

    ax2.pie(bioseq_counts_collapsed, labels=None, startangle=90, counterclock=False, colors=colors2, wedgeprops={'width': 0.4})
    ax2.legend(loc='upper center', bbox_to_anchor=(1, 0),
               labels=[f"{x}: {y}" for x, y in zip(bioseq_labels_clean, bioseq_percentages)],
               prop={'size': fontsize * 0.75}, frameon=False)
    ax2.set_title(f"Bioseq2Seq\n(Total Transcripts: {bioseq_counts.sum()})", size=fontsize)

    for ax in (ax1, ax2):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')


