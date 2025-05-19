import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlinescd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import matplotlib.patches as patches
from collections import Counter
from scipy.signal import find_peaks
from tqdm.notebook import tqdm
from multiprocessing import Pool
import re
import os
from functools import partial
from scipy.stats import ks_2samp
from scipy.stats import fisher_exact
import scipy.stats as stats

def create_splicetype_plot(df, ax, title, splicetype, direction,sample_name,xlab=True,ylab=True,export=False,base_size=12):
    """
    Plot the fold change with credible intervals, adjusting for specified splicetype.

    Parameters:
    df (pd.DataFrame): DataFrame containing the differential expression analysis results.
    ax: The axis object to plot on.
    title (str): Title for the plot.
    splicetype (str): The type of splicing event to consider ('ce', 'a5', 'a3', 'i').
    direction (str): The direction (i.e. Inclusion or Exclusion): I or E
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    
    # Dynamically scale font sizes
    fontsize = scale_fontsize(fig, base_size)
    # Ensure splicetype is one of the expected values
    if splicetype not in ['ce', 'a5', 'a3', 'i','ap','ae']:
        raise ValueError("splicetype must be one of: 'ce', 'a5', 'a3', 'i'")

     # Ensure splicetype is one of the expected values
    if direction not in ['I', 'E']:
        raise ValueError("direction must be one of: 'I', 'E'")

    # Check if splicetype column exists
    if splicetype not in df.columns:
        raise ValueError(f"{splicetype} column is missing from the DataFrame")

    # Duplicate rows based on the splicetype value because the code has the number of events by transcript row. if transcripts 
    #sig different then all splicing events are sig different. So 1 transcript can have multi sig events
    if splicetype in ['ce', 'a5', 'a3', 'i']:
        df = df.reindex(df.index.repeat(df[splicetype+direction])).reset_index(drop=True)
    if splicetype in ['ap','ae']:
        df = df[df[splicetype]>0]
    # Continue with existing plotting code
    df['log2_fold_change_mean'] = np.log2(df['fold_change_mean'])

    # Define conditions for upregulated and downregulated transcripts
    upregulated_condition = (df['log2_fold_change_mean'] > 0) & (df['fold_change_hdi_low'] > 1) & (df['p_diff_greater_than_zero'] > 0.95)
    
    downregulated_condition = (df['log2_fold_change_mean'] < 0) & (df['fold_change_hdi_high'] < 1) & (df['p_diff_less_than_zero'] > 0.95)

    df['color'] = 'grey'  # Default color for transcripts not strongly up or downregulated
    df.loc[upregulated_condition, 'color'] = 'red'
    df.loc[downregulated_condition, 'color'] = 'blue'
    df = df.sort_values(by='log2_fold_change_mean', ascending=True)

    # Create a new DataFrame to store the results
    results_df = pd.DataFrame({
        'Splicetype': [splicetype],
        'Upregulated': [len(df[df['color'] == 'red'][splicetype])],
        'Downregulated': [len(df[df['color'] == 'blue'][splicetype])],
        'Total': [len(df[splicetype])]
    })
    
    # Save the results to a CSV file
    results_df.to_csv("figure_files/"+sample_name+"_"+title+"_"+direction+'_splicing_results.csv', index=False)
    #print(splicetype+" up: "+str(len(df[df['color'] == 'red'][splicetype])),splicetype+" down: "+str(len(df[df['color'] == 'blue'][splicetype])),splicetype+" total: "+str(len(df[splicetype])))
    if export:
        df.to_csv("figure_files/"+sample_name+"_"+splicetype+'_'+direction+'.csv',index=True)
    
    # Plotting
    ax.errorbar(x=df['log2_fold_change_mean'], y=np.arange(df.shape[0]),
                xerr=[abs(df['log2_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])),
                      abs(np.log2(df['fold_change_hdi_high']) - df['log2_fold_change_mean'])],
                fmt='.', ecolor=df['color'], color='none', elinewidth=1, alpha=0.1)

    for category, color in zip(["Up regulated", "Down regulated", "Not affected"], ["red", "blue", "grey"]):
        condition = df['color'] == color
        
        ax.scatter(x=df[condition]['log2_fold_change_mean'], y=np.arange(df.shape[0])[condition],
                   color=color, s=0.5, alpha=0.5)

    # Custom legend handles
    legend_handles = [mlines.Line2D([], [], color='red', marker='o', linestyle='None', markersize=6, label='Up regulated'),
                      mlines.Line2D([], [], color='blue', marker='o', linestyle='None', markersize=6, label='Down regulated'),
                      mlines.Line2D([], [], color='grey', marker='o', linestyle='None', markersize=6, label="Not affected"),
                      mlines.Line2D([], [], color='grey', marker=None, linestyle='-', markersize=6, label="Credible Interval")]

    #min_max = [abs(df['log_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])).max(),
    #                  abs(np.log2(df['fold_change_hdi_high']) - df['log_fold_change_mean']).max()]
    #min_max = [abs(df['log_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])).max(),
    #                  abs(np.log2(df['fold_change_hdi_high']) + df['log_fold_change_mean']).max()]
    #max = np.log2(df.iloc[-1]['fold_change_hdi_high'])
    #print(min_max)
    ax.set_xlim(-40, 40)
    ax.tick_params(axis='y',labelsize=fontsize*.75)
    ax.tick_params(axis='x',labelsize=fontsize*.75)
    # Set the axis limits directly if needed
    #ax.set_ylim([-10, df.shape[0]+100)

    # Set specific ticks
    # This example sets the minimum and maximum x-ticks explicitly
    #ax.set_yticks([0, df.shape[0]])  # Adjust this to set your desired y-ticks


    # Set y-axis limits explicitly to include 0 to 1281
    ax.set_ylim(bottom=-10, top=df.shape[0]+10)

    # Get current y-axis ticks
    current_ticks = ax.get_yticks()

    # Create a new ticks array that ensures the last tick is 1281, without any beyond
    new_ticks = current_ticks[current_ticks <= df.shape[0]]
    if new_ticks[-1] != df.shape[0]:
        new_ticks = np.append(new_ticks, df.shape[0])  # Append 1281 if not the last tick
    new_ticks = new_ticks[1:]

    # Your target end tick
    end_tick = df.shape[0]

    # Check and remove any tick within % units of the end_tick
    new_ticks = [tick for tick in new_ticks if abs(tick - end_tick) > .1*end_tick]

    # Append the end_tick if it's not already the last tick
    if not new_ticks or new_ticks[-1] != end_tick:
        new_ticks.append(end_tick)

    ax.set_yticks(new_ticks)  # Apply the modified ticks


    if xlab:
        ax.set_xlabel(r"${log_{2}({PS}_{FC})}$",fontsize=fontsize*.75)
    if ylab:
        ax.set_ylabel('Events',fontsize=fontsize*.75)
    if title:
        ax.set_title(title, fontsize=fontsize*.75)
    #ax.legend(handles=legend_handles, fontsize=6)
    ax.axvline(x=0, linestyle='--', color='k', alpha=0.7)
    # Define conditions for upregulated and downregulated transcripts
    upregulated_condition = (df['log2_fold_change_mean'] > 0) & (df['fold_change_hdi_low'] > 1) & (df['p_diff_greater_than_zero'] > 0.95)
    
    downregulated_condition = (df['log2_fold_change_mean'] < 0) & (df['fold_change_hdi_high'] < 1) & (df['p_diff_less_than_zero'] > 0.95)

    df['color'] = 'grey'  # Default color for transcripts not strongly up or downregulated
    df.loc[upregulated_condition, 'color'] = 'red'
    df.loc[downregulated_condition, 'color'] = 'blue'
    df = df.sort_values(by='log2_fold_change_mean', ascending=True)

    # Create a new DataFrame to store the results
    results_df = pd.DataFrame({
        'Splicetype': [splicetype],
        'Upregulated': [len(df[df['color'] == 'red'][splicetype])],
        'Downregulated': [len(df[df['color'] == 'blue'][splicetype])],
        'Total': [len(df[splicetype])]
    })
    
    # Save the results to a CSV file
    #results_df.to_csv(sample_name+"_"+title+"_"+direction+'_splicing_results.csv', index=False)
    #print(splicetype+" up: "+str(len(df[df['color'] == 'red'][splicetype])),splicetype+" down: "+str(len(df[df['color'] == 'blue'][splicetype])),splicetype+" total: "+str(len(df[splicetype])))
    #if export:
    #    df.to_csv(title+"_"+splicetype+'_'+direction+'.csv',index=True)
    
    # Plotting
    ax.errorbar(x=df['log2_fold_change_mean'], y=np.arange(df.shape[0]),
                xerr=[abs(df['log2_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])),
                      abs(np.log2(df['fold_change_hdi_high']) - df['log2_fold_change_mean'])],
                fmt='.', ecolor=df['color'], color='none', elinewidth=1, alpha=0.1)

    for category, color in zip(["Up regulated", "Down regulated", "Not affected"], ["red", "blue", "grey"]):
        condition = df['color'] == color
        
        ax.scatter(x=df[condition]['log2_fold_change_mean'], y=np.arange(df.shape[0])[condition],
                   color=color, s=0.5, alpha=0.5)

    # Custom legend handles
    legend_handles = [mlines.Line2D([], [], color='red', marker='o', linestyle='None', markersize=6, label='Up regulated'),
                      mlines.Line2D([], [], color='blue', marker='o', linestyle='None', markersize=6, label='Down regulated'),
                      mlines.Line2D([], [], color='grey', marker='o', linestyle='None', markersize=6, label="Not affected"),
                      mlines.Line2D([], [], color='grey', marker=None, linestyle='-', markersize=6, label="Credible Interval")]

    #min_max = [abs(df['log_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])).max(),
    #                  abs(np.log2(df['fold_change_hdi_high']) - df['log_fold_change_mean']).max()]
    #min_max = [abs(df['log_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])).max(),
    #                  abs(np.log2(df['fold_change_hdi_high']) + df['log_fold_change_mean']).max()]
    #max = np.log2(df.iloc[-1]['fold_change_hdi_high'])
    #print(min_max)
    ax.set_xlim(-40, 40)
    ax.tick_params(axis='y',labelsize=fontsize*.75)
    ax.tick_params(axis='x',labelsize=fontsize*.75)
    # Set the axis limits directly if needed
    #ax.set_ylim([-10, df.shape[0]+100)

    # Set specific ticks
    # This example sets the minimum and maximum x-ticks explicitly
    #ax.set_yticks([0, df.shape[0]])  # Adjust this to set your desired y-ticks


    # Set y-axis limits explicitly to include 0 to 1281
    ax.set_ylim(bottom=-10, top=df.shape[0]+10)

    # Get current y-axis ticks
    current_ticks = ax.get_yticks()

    # Create a new ticks array that ensures the last tick is 1281, without any beyond
    new_ticks = current_ticks[current_ticks <= df.shape[0]]
    if new_ticks[-1] != df.shape[0]:
        new_ticks = np.append(new_ticks, df.shape[0])  # Append 1281 if not the last tick
    new_ticks = new_ticks[1:]

    # Your target end tick
    end_tick = df.shape[0]

    # Check and remove any tick within % units of the end_tick
    new_ticks = [tick for tick in new_ticks if abs(tick - end_tick) > .1*end_tick]

    # Append the end_tick if it's not already the last tick
    if not new_ticks or new_ticks[-1] != end_tick:
        new_ticks.append(end_tick)

    ax.set_yticks(new_ticks)  # Apply the modified ticks


    if xlab:
        ax.set_xlabel(r"${log_{2}({PS}_{FC})}$",fontsize=fontsize*.75)
    if ylab:
        ax.set_ylabel('Events',fontsize=fontsize*.75)
    if title:
        ax.set_title(title, fontsize=fontsize*.75)
    #ax.legend(handles=legend_handles, fontsize=6)
    ax.axvline(x=0, linestyle='--', color='k', alpha=0.7)

def create_dist_plot(df_avsb, ax, title, splicetype, csv_file_path, test_result_path, xlab=True,ylab=True):
    """
    Plot the fold change with credible intervals, adjusting for specified splicetype.

    Parameters:
    df (pd.DataFrame): DataFrame containing the differential expression analysis results.
    ax: The axis object to plot on.
    title (str): Title for the plot.
    splicetype (str): The type of splicing event to consider ('ce', 'a5', 'a3', 'i').
    """

    if splicetype in ['ap','ae']:
        df_avsb = df_avsb[(df_avsb[splicetype] >= 0) & (df_avsb[splicetype] <= 1)]

    # Apply log2 transformation to the 'fold_change_mean' column
    df_avsb['log2_fold_change'] = np.log2(df_avsb['fold_change_mean'])

    # Define conditions for significantly affected transcripts
    upregulated_avsb = (df_avsb['log2_fold_change'] > 0) & (df_avsb['fold_change_hdi_low'] > 1) & (df_avsb['p_diff_greater_than_zero'] > 0.95)
    #downregulated_avsb = (df_avsb['log2_fold_change'] < 0) & (df_avsb['fold_change_hdi_high'] < 1) & (df_avsb['p_diff_greater_than_zero'] <= 0.05)
    downregulated_avsb = (df_avsb['log2_fold_change'] < 0) & (df_avsb['fold_change_hdi_high'] < 1) & (df_avsb['p_diff_less_than_zero'] > 0.95)

    # Define conditions for not significantly affected transcripts
    not_affected_avsb = ~(upregulated_avsb | downregulated_avsb)
    affected_avsb = (upregulated_avsb | downregulated_avsb)
    
    # Find indices for transcripts significantly affected in AvsB, CvsD, and both
    up_avsb = df_avsb.index[upregulated_avsb]
    down_avsb = df_avsb.index[downregulated_avsb]
    affected_only_in_avsb = df_avsb.index[upregulated_avsb | downregulated_avsb]

    # Find indices for transcripts not significantly affected in both
    not_affected_in_both = df_avsb.index[not_affected_avsb]
    #print(not_affected_in_both)

     # Adjust the plot_ecdf function to accept the ax parameter
    def plot_ecdf(data, label, color, ax):
        #sns.ecdfplot(data=data, ax=ax, label=label, color=color)
        sns.kdeplot(data=data, ax=ax, label=label, color=color,fill=True)
    
    # Ensure ax is defined
    if ax is None:
        fig, ax = plt.subplots()

    # Plot ECDF for all transcripts and specific groups
    plot_ecdf(df_avsb.loc[not_affected_in_both, splicetype], 'All Transcripts', 'lightgrey', ax)
    plot_ecdf(df_avsb.loc[affected_only_in_avsb, splicetype], 'Significant', 'skyblue', ax)
    #plot_ecdf(df_avsb.loc[up_avsb, splicetype], 'Significant', 'red', ax)
    #plot_ecdf(df_avsb.loc[down_avsb, splicetype], 'Significant', 'blue', ax)

    
    
    # Adjust the plot aesthetics
    ax.tick_params(axis='y',labelsize=6*.75)
    ax.tick_params(axis='x',labelsize=6*.75)
    ax.set_xlabel(xlab,fontsize=6)
    ax.set_ylabel(ylab,fontsize=6)
    ax.set_title(title,fontsize=6)
    ax.set_xlim([-0.5,1.5])
    #ax.legend()
    #ax.set_aspect('equal', 'box')
    if ax is None:
        plt.tight_layout()
        plt.show()

    # Create a DataFrame to save the data for statistical tests
    data_for_stats = pd.DataFrame({
        'splicetype': df_avsb[splicetype],
        'group': np.select(
            [affected_avsb, not_affected_avsb],
            ['affected', 'not_affected'],
            default='not_affected'
        )
    })

    # Save the DataFrame to a CSV file
    data_for_stats.to_csv(csv_file_path, index=False)

    # Extract the data for each group
    not_affected_data = data_for_stats[data_for_stats['group'] == 'not_affected']['splicetype']
    affected_data = data_for_stats[data_for_stats['group'] == 'affected']['splicetype']
    #upregulated_data = data_for_stats[data_for_stats['group'] == 'upregulated']['splicetype']
    #downregulated_data = data_for_stats[data_for_stats['group'] == 'downregulated']['splicetype']

    # Perform Kolmogorov-Smirnov tests
    ks_stat_affected, p_value_affected = ks_2samp(not_affected_data, affected_data)
    #ks_stat_down, p_value_down = ks_2samp(not_affected_data, downregulated_data)

    # Create a DataFrame for the test results
    test_results = pd.DataFrame({
        'Comparison': ['Affected vs Not Affected'],
        'KS Statistic': [ks_stat_affected],
        'p-value': [p_value_affected]
    })

    # Save the test results to a CSV file
    test_results.to_csv(test_result_path, index=False)

    return data_for_stats, test_results

def load_csv_as_dataframe(filename):
    return pd.read_csv(filename, index_col=0)

def scale_fontsize(fig, base_size=12):
    width, height = fig.get_size_inches()
    return base_size * min(width / 8, height / 11)

def fractional_distances(keys):
    # Convert dict_keys to a list if not already a list
    key_list = list(keys)
    # Find the minimum and maximum values
    min_key = min(key_list)
    max_key = max(key_list)
    # Calculate fractional distances
    fractional_distances = [(key - min_key) / (max_key - min_key) for key in key_list]
    return fractional_distances

def identify_alternative_promoters(gtf_data, gene_id, promoter_region_size=1):
    """
    Correctly identify alternative promoters for a given gene, ensuring unique association of transcripts with promoters.

    :param gtf_data: DataFrame containing GTF data.
    :param gene_id: The ID of the gene.
    :param promoter_region_size: The size of the promoter region to consider around the TSS.
    :return: Dictionary of alternative promoters and their associated transcripts.
    """
    # Filter data for the given gene and exons
    gene_exon_data = gtf_data[
        (gtf_data['gene_id'] == gene_id) &
        (gtf_data['feature'] == 'exon')
    ].copy()

    # Ensure 'strand' column is present
    if 'strand' not in gene_exon_data.columns:
        raise ValueError("The GTF data must contain a 'strand' column.")

    # Initialize dictionaries to store promoters and fractional distances
    unique_promoters = {}
    fract_dict = {}

    # Group data by transcript_id
    grouped = gene_exon_data.groupby('transcript_id')

    # Iterate over each transcript to find the TSS based on strand
    for transcript_id, group in grouped:
        strand = group['strand'].iloc[0]

        if strand == '+':
            # For positive strand, TSS is the minimum 'start' position
            tss = group['start'].min()
            promoter_region = (tss - promoter_region_size, tss + promoter_region_size)
        elif strand == '-':
            # For negative strand, TSS is the maximum 'end' position
            tss = group['end'].max()
            promoter_region = (tss - promoter_region_size, tss + promoter_region_size)
        else:
            raise ValueError(f"Unknown strand '{strand}' for transcript '{transcript_id}'.")

        # Check if the promoter region overlaps with any existing promoter
        found_overlap = False
        for unique_tss, data in unique_promoters.items():
            existing_promoter_region = data['promoter_region']
            # Check for overlap
            if not (promoter_region[1] < existing_promoter_region[0] or promoter_region[0] > existing_promoter_region[1]):
                # Overlaps with existing promoter
                data['transcripts'].add(transcript_id)
                found_overlap = True
                break

        # If no overlap, add as a new unique promoter
        if not found_overlap:
            unique_promoters[tss] = {
                'transcripts': set([transcript_id]),
                'strand': strand,
                'promoter_region': promoter_region
            }

    # Now, calculate fractional distances if there are multiple promoters
    if len(unique_promoters) > 1:
        # Extract TSS positions for fractional distance calculation
        tss_positions = list(unique_promoters.keys())
        min_tss = min(tss_positions)
        max_tss = max(tss_positions)

        # Calculate fractional distances based on strand
        fractional_distances = {}
        for tss, data in unique_promoters.items():
            strand = data['strand']
            if strand == '+':
                fractional_distance = (tss - min_tss) / (max_tss - min_tss)
            else:
                fractional_distance = (max_tss - tss) / (max_tss - min_tss)
            fractional_distances[tss] = fractional_distance

        # Assign fractional distances to transcripts
        for tss, data in unique_promoters.items():
            for transcript_id in data['transcripts']:
                fract_dict[transcript_id] = fractional_distances[tss]
    else:
        # Only one promoter, assign 0 fractional distance to all transcripts
        for data in unique_promoters.values():
            for transcript_id in data['transcripts']:
                fract_dict[transcript_id] = 0.0

    # Prepare the final promoters dictionary
    promoters_dict = {}
    for tss, data in unique_promoters.items():
        promoters_dict[tss] = data['transcripts']

    return promoters_dict, fract_dict

def identify_alternative_ends(gtf_data, gene_id, end_region_size=1):
    """
    Correctly identify alternative ends for a given gene, ensuring unique association of transcripts with ends.

    :param gtf_data: DataFrame containing GTF data.
    :param gene_id: The ID of the gene.
    :param end_region_size: The size of the end region to consider around the TES.
    :return: Dictionary of alternative ends and their associated transcripts.
    """
    # Filter data for the given gene and exons
    gene_exon_data = gtf_data[
        (gtf_data['gene_id'] == gene_id) &
        (gtf_data['feature'] == 'exon')
    ].copy()

    # Ensure 'strand' column is present
    if 'strand' not in gene_exon_data.columns:
        raise ValueError("The GTF data must contain a 'strand' column.")

    # Initialize dictionaries to store ends and fractional distances
    unique_ends = {}
    fract_dict = {}

    # Group data by transcript_id
    grouped = gene_exon_data.groupby('transcript_id')

    # Iterate over each transcript to find the TES based on strand
    for transcript_id, group in grouped:
        strand = group['strand'].iloc[0]

        if strand == '+':
            # For positive strand, TES is the maximum 'end' position
            tes = group['end'].max()
            end_region = (tes - end_region_size, tes + end_region_size)
        elif strand == '-':
            # For negative strand, TES is the minimum 'start' position
            tes = group['start'].min()
            end_region = (tes - end_region_size, tes + end_region_size)
        else:
            raise ValueError(f"Unknown strand '{strand}' for transcript '{transcript_id}'.")

        # Check if the end region overlaps with any existing end
        found_overlap = False
        for unique_tes, data in unique_ends.items():
            existing_end_region = data['end_region']
            # Check for overlap
            if not (end_region[1] < existing_end_region[0] or end_region[0] > existing_end_region[1]):
                # Overlaps with existing end
                data['transcripts'].add(transcript_id)
                found_overlap = True
                break

        # If no overlap, add as a new unique end
        if not found_overlap:
            unique_ends[tes] = {
                'transcripts': set([transcript_id]),
                'strand': strand,
                'end_region': end_region
            }

    # Now, calculate fractional distances if there are multiple ends
    if len(unique_ends) > 1:
        # Extract TES positions for fractional distance calculation
        tes_positions = list(unique_ends.keys())
        min_tes = min(tes_positions)
        max_tes = max(tes_positions)

        # Calculate fractional distances based on strand
        fractional_distances = {}
        for tes, data in unique_ends.items():
            strand = data['strand']
            if strand == '+':
                fractional_distance = (tes - min_tes) / (max_tes - min_tes)
            else:
                fractional_distance = (max_tes - tes) / (max_tes - min_tes)
            fractional_distances[tes] = fractional_distance

        # Assign fractional distances to transcripts
        for tes, data in unique_ends.items():
            for transcript_id in data['transcripts']:
                fract_dict[transcript_id] = fractional_distances[tes]
    else:
        # Only one end, assign 0 fractional distance to all transcripts
        for data in unique_ends.values():
            for transcript_id in data['transcripts']:
                fract_dict[transcript_id] = 0.0

    # Prepare the final ends dictionary
    ends_dict = {}
    for tes, data in unique_ends.items():
        ends_dict[tes] = data['transcripts']

    return ends_dict, fract_dict

def plot_alternative_promoters(transcript_data, alt_promoters, max_transcripts=30):
    """
    Plot the splicing events for the transcripts of a gene, showing transcripts grouped by alternative promoter usage.

    :param transcript_data: A Series object where each entry is a list of [start, end] positions for exons in a transcript.
    :param alt_promoters: A dictionary with alternative promoter start positions as keys and list of transcript IDs as values.
    :param max_transcripts: Maximum number of transcripts to plot.
    """
    # Sorting transcripts by promoter groups
    sorted_transcript_data = []
    promoter_group_indices = {}  # To mark the start and end indices of each promoter group
    current_index = 0

    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in transcript_data:
                sorted_transcript_data.append((transcript_id, transcript_data[transcript_id]))
                current_index += 1

    # Limiting the number of transcripts to plot for clarity
    if len(sorted_transcript_data) > max_transcripts:
        sorted_transcript_data = sorted_transcript_data[:max_transcripts]

    fig, ax = plt.subplots(figsize=(12, len(sorted_transcript_data) * 0.6))

    # Generate a unique color for each alternative promoter group using viridis colormap
    promoter_colors = plt.get_cmap('viridis', len(alt_promoters))
    #promoter_colors = plt.cm.get_cmap('viridis', len(alt_promoters))
    promoter_color_map = {start_pos: promoter_colors(i) for i, start_pos in enumerate(alt_promoters.keys())}

    # Plotting each transcript
    for i, (transcript_id, exons) in enumerate(sorted_transcript_data):
        #print(exons)
        promoter_color = 'black'  # Default color
        for start_pos, transcripts in alt_promoters.items():
            if transcript_id in transcripts:
                promoter_color = promoter_color_map[start_pos]
                break

        # Drawing the transcript line
        transcript_start = min([exon[0] for exon in exons])
        transcript_end = max([exon[1] for exon in exons])
        ax.plot([transcript_start+5, transcript_end-5], [i, i], color=promoter_color, lw=2)

        # Drawing the exons as rectangles
        for exon in exons:
            exon_start, exon_end = exon
            #print(exon)
            rect = patches.Rectangle((exon_start, i-.1), exon_end-exon_start, 0.2, color=promoter_color)
            ax.add_patch(rect)

    # Drawing lines to mark promoter groups
    for start_pos, (start_index, end_index) in promoter_group_indices.items():
        ax.axhline(y=start_index - 0.5, color='black', linestyle='--')
        ax.axhline(y=end_index + 0.5, color='black', linestyle='--')

    ax.set_yticks(np.arange(len(sorted_transcript_data)))
    ax.set_yticklabels([transcript_id for transcript_id, _ in sorted_transcript_data])
    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Transcripts')
    ax.set_xlim(min(transcript_data.explode().min())-100, max(transcript_data.explode().max())+100)  # Set x-axis limits
    ax.set_title('AS Events')

    # Creating a legend for the alternative promoters
    #ax.legend([patches.Patch(color=promoter_colors(0))], ['Alternative Promoter Group'])

    plt.tight_layout()

    return fig, ax

def create_global_exon_index_for_promoter_group(gtf_data, transcripts):
    """
    Create a global index of all exons for a given promoter group of transcripts, including identification of
    constitutive and alternative exons.

    :param gtf_data: DataFrame containing GTF data.
    :param transcripts: List of transcripts in the promoter group.
    :return: DataFrame with global indexed exons and splicing events for each transcript.
    """
    # Filter exons belonging to the transcripts in the promoter group
    promoter_exons = gtf_data[gtf_data['transcript_id'].isin(transcripts) & (gtf_data['feature'] == 'exon')]

    # Creating a unique identifier for each exon based on start and end positions
    promoter_exons['exon_identifier'] = promoter_exons['start'].astype(str) + '-' + promoter_exons['end'].astype(str)

    # Creating a global index for exons
    unique_exons = promoter_exons['exon_identifier'].drop_duplicates().reset_index(drop=True)
    exon_global_index = {exon: idx + 1 for idx, exon in enumerate(unique_exons)}

    # Mapping transcripts to the global exon index
    promoter_exons['global_exon_index'] = promoter_exons['exon_identifier'].map(exon_global_index)

    # Constitutive exons check
    exon_freq = promoter_exons['exon_identifier'].value_counts()
    constitutive_exons = exon_freq[exon_freq == len(transcripts)].index.tolist()

    # Identify splicing events
    for index, row in promoter_exons.iterrows():
        exon_id = row['exon_identifier']

        if exon_id in constitutive_exons:
            promoter_exons.at[index, 'exon_type'] = 'Constitutive'
        else:
            promoter_exons.at[index, 'exon_type'] = 'Alternative'

    return promoter_exons[['transcript_id', 'global_exon_index', 'start', 'end', 'exon_type']]


def normalize_exon_identifier(row, variation=2):
    """
    Normalize exon identifier based on a small variation in start and end positions.

    :param row: DataFrame row containing exon information.
    :param variation: Allowed variation in start and end positions.
    :return: Normalized exon identifier.
    """
    normalized_start = row['start'] // variation * variation
    normalized_end = row['end'] // variation * variation
    #return f"{normalized_start}-{normalized_end}"
    return normalized_start, normalized_end

def get_exon_ranges(gtf_data, alt_promoters):
    # Filter exons belonging to the transcripts in the promoter group
    promoter_exons = gtf_data[
        gtf_data['transcript_id'].isin(alt_promoters) &
        (gtf_data['feature'] == 'exon')
    ].copy()

    # Ensure 'strand' column is present
    if 'strand' not in promoter_exons.columns:
        raise ValueError("The GTF data must contain a 'strand' column.")

    # Get the strand of the gene (assuming all exons have the same strand)
    strand = promoter_exons['strand'].iloc[0]

    # Sort exons based on strand
    if strand == '+':
        promoter_exons = promoter_exons.sort_values(by=["start", "end"])
    else:
        promoter_exons = promoter_exons.sort_values(by=["end", "start"], ascending=[False, False])

    # Create exon identifiers
    promoter_exons['exon_identifier_start'] = promoter_exons['start']
    promoter_exons['exon_identifier_end'] = promoter_exons['end']
    promoter_exons['exon_identifier'] = promoter_exons.apply(lambda row: (row['start'], row['end']), axis=1)
    exon_ranges = promoter_exons['exon_identifier'].unique()
    return exon_ranges, promoter_exons


def find_alternative_5_prime_splice_sites(coverage, intronic_regions, peaks, strand):
    """
    Identify potential alternative 5' splice sites in the coverage data,
    accounting for the gene's strand.

    :param coverage: List of coverage values.
    :param intronic_regions: List of tuples representing intronic regions.
    :param peaks: List of peak positions.
    :param strand: '+' or '-' indicating the gene's strand.
    :return: List of positions of alternative 5' splice sites.
    """
    alt_5_prime_splice_sites = []

    if strand == '+':
        # For positive strand, 5' end is at the beginning
        last_peak = max(peaks)
        iterator = range(1, len(coverage))
        compare = lambda i: coverage[i] < coverage[i - 1]
        position_check = lambda i: i <= last_peak
    elif strand == '-':
        # For negative strand, 5' end is at the end
        last_peak = min(peaks)
        iterator = range(len(coverage) - 2, -1, -1)
        compare = lambda i: coverage[i] < coverage[i + 1]
        position_check = lambda i: i >= last_peak
    else:
        raise ValueError("Strand must be '+' or '-'.")

    for i in iterator:
        # Check if current position is within an intron or beyond the last peak
        if any(start <= i <= end for start, end in intronic_regions) or not position_check(i):
            continue

        # Check for a step (decrease in coverage)
        if compare(i):
            alt_5_prime_splice_sites.append(i)

    return alt_5_prime_splice_sites

def find_alternative_3_prime_splice_sites(coverage, intronic_regions, peaks, strand):
    """
    Identify potential alternative 3' splice sites in the coverage data,
    accounting for the gene's strand.

    :param coverage: List of coverage values.
    :param intronic_regions: List of tuples representing intronic regions.
    :param peaks: List of peak positions.
    :param strand: '+' or '-' indicating the gene's strand.
    :return: List of positions of alternative 3' splice sites.
    """
    alt_3_prime_splice_sites = []

    if strand == '+':
        # For positive strand, iterate backwards from each peak
        for peak in peaks:
            for i in range(peak - 1, -1, -1):
                if any(start <= i <= end for start, end in intronic_regions):
                    continue

                if coverage[i] < coverage[i + 1]:
                    alt_3_prime_splice_sites.append(i)
                    break
    elif strand == '-':
        # For negative strand, iterate forwards from each peak
        for peak in peaks:
            for i in range(peak + 1, len(coverage)):
                if any(start <= i <= end for start, end in intronic_regions):
                    continue

                if coverage[i] < coverage[i - 1]:
                    alt_3_prime_splice_sites.append(i)
                    break
    else:
        raise ValueError("Strand must be '+' or '-'.")

    return alt_3_prime_splice_sites


def find_alternative_end(coverage, intronic_regions, peaks, strand):
    """
    Identify alternative transcription end sites (TES) in the coverage data,
    accounting for the gene's strand.

    :param coverage: List of coverage values.
    :param intronic_regions: List of tuples representing intronic regions.
    :param peaks: List of peak positions.
    :param strand: '+' or '-' indicating the gene's strand.
    :return: List of positions of alternative ends.
    """
    alt_end = []
    coverage = coverage + [0]*20 + [1]*20  # Padding
    ir = intronic_regions + [(len(coverage) + 1, len(coverage) + 10)]
    pad_peak = len(coverage) - 10

    if strand == '+':
        last_peak = max(peaks)
        iterator = range(last_peak + 1, pad_peak)
        compare = lambda i: coverage[i] < coverage[i - 1]
    elif strand == '-':
        last_peak = min(peaks)
        iterator = range(last_peak - 1, 0, -1)
        compare = lambda i: coverage[i] < coverage[i + 1]
    else:
        raise ValueError("Strand must be '+' or '-'.")

    for i in iterator:
        # Check if current position is within an intron or beyond the padding peak
        if any(start <= i <= end for start, end in intronic_regions) or (strand == '+' and i > pad_peak) or (strand == '-' and i < 0):
            continue

        # Check for a step (decrease in coverage)
        if compare(i):
            alt_end.append(i)

    return alt_end

def plot_exons_with_introns_and_alt_splice_sites(exon_ranges, intron_ranges, alt_5splice_sites,alt_3splice_sites,alt_ends,exon_ranges_below_one):
    """
    Plot the exon structure with shaded intronic regions and alternative 5' splice sites.
    """
    # Normalize ranges for plotting
    min_start = min([start for start, _ in exon_ranges])
    normalized_exon_ranges = [(start - min_start, end - min_start) for start, end in exon_ranges]

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plotting exons
    for i, (start, end) in enumerate(normalized_exon_ranges):
        exon_rect = patches.Rectangle((start, i - 0.1), end - start, 0.2, color='skyblue', edgecolor='black')
        ax.add_patch(exon_rect)
        ax.text(start - 100, i, str(i + 1), verticalalignment='center', fontsize=8)  # Label each exon numerically

    # Highlighting introns
    for start, end in intron_ranges:
        normalized_start = start - min_start
        normalized_end = end - min_start
        intron_rect = patches.Rectangle((normalized_start, -1), normalized_end - normalized_start, len(exon_ranges) + 1, color='lightgreen', alpha=0.3)
        ax.add_patch(intron_rect)

    # Highlighting alternative 5' splice sites
    for site in alt_5splice_sites:
        normalized_site = site
        ax.axvline(x=normalized_site, color='orange', linestyle='--',lw=0.5)

    # Highlighting alternative 3' splice sites
    for site in alt_3splice_sites:
        normalized_site = site
        ax.axvline(x=normalized_site, color='red', linestyle='--',lw=0.5)

    # Highlighting alternative ends
    #use the key values from the identify_alternative_ends def
    for site in alt_ends:
        normalized_site = site - min_start
        ax.axvline(x=normalized_site, color='black', linestyle='--',lw=0.5)

     # Highlighting cassette exons
    for start, end in exon_ranges_below_one:
        normalized_start = start - min_start
        normalized_end = end - min_start
        ce_rect = patches.Rectangle((normalized_start, 0), normalized_end - normalized_start,  len(exon_ranges) + 1, color='red', alpha=0.3)
        ax.add_patch(ce_rect)


    # Setting plot parameters
    ax.set_yticks(range(len(exon_ranges)))
    ax.set_yticklabels(range(1, len(exon_ranges) + 1))
    ax.set_ylabel('Exon Number')
    ax.set_xlabel('Genomic Position (Normalized)')
    ax.set_title('Exon Arrangement with AS Events Highlight')
    ax.set_xlim(-200, max([end for _, end in normalized_exon_ranges]) + 200)

    plt.tight_layout()
    plt.show()

def get_cassette_exons(promoter_exons):
    # Step 1: Determine the span of each transcript
    transcript_spans = {}
    for _, row in promoter_exons.iterrows():
        transcript_id = row['transcript_id']
        start = row['exon_identifier_start']
        end = row['exon_identifier_end']

        if transcript_id not in transcript_spans:
            transcript_spans[transcript_id] = [start, end]
        else:
            transcript_spans[transcript_id][0] = min(transcript_spans[transcript_id][0], start)
            transcript_spans[transcript_id][1] = max(transcript_spans[transcript_id][1], end)

    # Step 2: Count the number of transcripts each exon is part of
    exon_transcript_count = {}
    for _, row in promoter_exons.iterrows():
        exon_range = (row['exon_identifier_start'], row['exon_identifier_end'])
        exon_transcript_count[exon_range] = sum(
            exon_range[0] >= span[0] and exon_range[1] <= span[1]
            for span in transcript_spans.values()
        )

    min_start_t = min(promoter_exons['exon_identifier_start'])
    max_end_t = max(promoter_exons['exon_identifier_end'])

    coverage_t = [0] * (max_end_t - min_start_t + 1)
    for _, row in promoter_exons.iterrows():
        start = row['exon_identifier_start']
        end = row['exon_identifier_end']
        exon_range = (start, end)
        exon_count = exon_transcript_count.get(exon_range, 0)

        for i in range(start - min_start_t, end - min_start_t + 1):
            coverage_t[i] += 1 / exon_count if exon_count else 0

    # Find peaks in the coverage
    peaks, _ = find_peaks(coverage_t, height=0)  # Adjust parameters as needed for peak detection
    valleys, _ = find_peaks([-p for p in coverage_t])

    # Identify peaks with coverage less than 1 and return their exon ranges
    peaks_below_one = [peak for peak in peaks if round(coverage_t[peak],3) < 1]

    # Find peaks with coverage less than 1 that are flanked by valleys
    exon_ranges_below_one = []
    for peak in peaks:
        if round(coverage_t[peak], 3) < 1:
            # Finding the closest valleys before and after the peak
            previous_valleys = [v for v in valleys if v < peak]
            next_valleys = [v for v in valleys if v > peak]

            if previous_valleys and next_valleys:
                previous_valley = previous_valleys[-1]+min_start_t
                next_valley = next_valleys[0]+min_start_t
                #print(previous_valley,next_valley)
                # Convert peak position back to exon range
                peak_position = peak + min_start_t
                for exon_range in exon_transcript_count.keys():
                    # Filter out the first and last exons
                    if exon_range[0] != min_start_t and exon_range[1] != max_end_t:
                    #print(exon_range)
                        if exon_range[0] > previous_valley and exon_range[1] < next_valley:
                            exon_ranges_below_one.append(exon_range)
                            break  # Break the loop once the corresponding exon range is found

    return exon_ranges_below_one

def ranges_overlap(range1, range2):
    return range1[0] <= range2[1] and range1[1] >= range2[0]

def get_intronic_region(promoter_exons):
    """
        The output from this function converts the python 0 based coords back to
        to genomic based  coords
    """

    # Step 1: Determine the span of each transcript
    transcript_spans = {}
    for _, row in promoter_exons.iterrows():
        transcript_id = row['transcript_id']
        start = row['exon_identifier_start']
        end = row['exon_identifier_end']

        if transcript_id not in transcript_spans:
            transcript_spans[transcript_id] = [start, end]
        else:
            transcript_spans[transcript_id][0] = min(transcript_spans[transcript_id][0], start)
            transcript_spans[transcript_id][1] = max(transcript_spans[transcript_id][1], end)

    # Step 2: Count the number of transcripts each exon is part of
    exon_transcript_count = {}
    for _, row in promoter_exons.iterrows():
        exon_range = (row['exon_identifier_start'], row['exon_identifier_end'])
        exon_transcript_count[exon_range] = sum(
            exon_range[0] >= span[0] and exon_range[1] <= span[1]
            for span in transcript_spans.values()
        )
    #Step 3: Create  a list of intron start and stop locations
    min_start_t = min(promoter_exons['exon_identifier_start'])
    max_end_t = max(promoter_exons['exon_identifier_end'])

    coverage_t = [0] * (max_end_t - min_start_t + 1)
    for _, row in promoter_exons.iterrows():
        start = row['exon_identifier_start']
        end = row['exon_identifier_end']
        exon_range = (start, end)
        exon_count = exon_transcript_count.get(exon_range, 0)
        #0 based to 1 based conversion
        for i in range(start - min_start_t + 1, end - min_start_t + 1):
            coverage_t[i] += 1 / exon_count if exon_count else 0


    # Find peaks in the coverage
    peaks, _ = find_peaks(coverage_t, height=0)  # Adjust parameters as needed for peak detection
    valleys, _ = find_peaks([-p for p in coverage_t])

    # Identify peaks with coverage less than 1 and return their exon ranges
    valleys_above_0 = [valley for valley in valleys if round(coverage_t[valley],3) >= 0]

    # Find valleys with coverage greater than or equal to 0 that are flanked by peaks
    intron_ranges = []

    for valley in valleys:
        start = valley
        end = valley

        # Scan backward
        while start > 0 and coverage_t[start - 1] <= coverage_t[start]:
            start -= 1

        # Scan forward
        while end < len(coverage_t) - 1 and coverage_t[end + 1] <= coverage_t[end]:
            end += 1

        # Add the start and end points to the list
        intron_ranges.append((start, end))

    return intron_ranges

def norm_coverage(promoter_exons):
    # Step 1: Determine the span of each transcript
    transcript_spans = {}
    for _, row in promoter_exons.iterrows():
        transcript_id = row['transcript_id']
        start = row['start']
        end = row['end']
        #normed range
        #start = row['exon_identifier_start']
        #end = row['exon_identifier_end']

        if transcript_id not in transcript_spans:
            transcript_spans[transcript_id] = [start, end]
        else:
            transcript_spans[transcript_id][0] = min(transcript_spans[transcript_id][0], start)
            transcript_spans[transcript_id][1] = max(transcript_spans[transcript_id][1], end)

    # Step 2: Count the number of transcripts each exon is part of
    exon_transcript_count = {}
    for _, row in promoter_exons.iterrows():
        exon_range = (row['start'], row['end'])
        #norm ranged
        #exon_range = (row['exon_identifier_start'], row['exon_identifier_end'])
        exon_transcript_count[exon_range] = sum(
            exon_range[0] >= span[0] and exon_range[1] <= span[1]
            for span in transcript_spans.values()
        )

    #min_start_t = min(promoter_exons['exon_identifier_start'])
    #max_end_t = max(promoter_exons['exon_identifier_end'])
    min_start_t = min(promoter_exons['start'])
    max_end_t = max(promoter_exons['end'])

    coverage_t = [0] * (max_end_t - min_start_t+1)
    for _, row in promoter_exons.iterrows():
        #start = row['exon_identifier_start']
        #end = row['exon_identifier_end']
        start = row['start']
        end = row['end']
        exon_range = (start, end)
        exon_count = exon_transcript_count.get(exon_range, 0)

        for i in range(start - min_start_t, end - min_start_t+1 ):
            coverage_t[i] += 1 / exon_count if exon_count else 0

    # Find peaks in the coverage
    # Pad the coverage array for peak detection
    coverage_t_pad = [0] + coverage_t + [0]
    peaks, _ = find_peaks(coverage_t_pad, height=0)  # Adjust parameters as needed for peak detection
    valleys, _ = find_peaks([-p for p in coverage_t])
    # Adjust the peak positions to account for the padding
    adjusted_peaks = [peak - 1 for peak in peaks]
    return coverage_t, adjusted_peaks, valleys



















def plot_alternative_events(
    gtf_data, 
    gene_id, 
    alt_promoters, 
    alt5sites, 
    alt3sites, 
    alt_ends, 
    exon_ranges_below_one, 
    nintrons, 
    min_start, 
    ncoverage, 
    strand, 
    ax=None, 
    max_transcripts=30, 
    fscale=10,
    hide_xticks=True,
    hide_yticks=False
):
    """
    Plot the splicing events for the transcripts of a gene, showing transcripts grouped by alternative promoter usage.

    Parameters:
        gtf_data (DataFrame): DataFrame containing GTF data.
        gene_id (str/int): The gene ID.
        alt_promoters (dict): Dictionary with alternative promoter start positions as keys and list of transcript IDs as values.
        alt5sites (list): List of alternative 5' splice sites.
        alt3sites (list): List of alternative 3' splice sites.
        alt_ends (list): List of alternative end splice sites.
        exon_ranges_below_one (list): List of exon ranges to highlight.
        nintrons (list): List of intronic regions.
        min_start (int): Minimum start position for normalization.
        ncoverage (list/array): Coverage data.
        strand (str): Strand information ('+' or '-').
        ax (matplotlib.axes.Axes, optional): Axes to plot on. If None, creates a new figure and axes.
        max_transcripts (int, optional): Maximum number of transcripts to plot.
        fscale (float, optional): Scaling factor for font sizes.

    Returns:
        matplotlib.axes.Axes: The axes with the plot.
    """
    # If no Axes provided, create a new figure and Axes
    created_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        created_ax = True

    # Extract transcript data
    transcript_data = gtf_data[
        (gtf_data['gene_id'] == gene_id) & (gtf_data['feature'] == 'exon')
    ].groupby('transcript_id').apply(lambda x: x[['start', 'end']].values.tolist())

    # Sorting transcripts by promoter groups
    sorted_transcript_data = []
    promoter_group_indices = {}  # To mark the start and end indices of each promoter group
    current_index = 0

    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in transcript_data:
                sorted_transcript_data.append((transcript_id, transcript_data[transcript_id]))
                current_index += 1

    # Limiting the number of transcripts to plot for clarity
    if len(sorted_transcript_data) > max_transcripts:
        sorted_transcript_data = sorted_transcript_data[:max_transcripts]

    # Calculate font scale
    font_scale_ax = fscale  # Adjust as needed or make dynamic based on ax size

    # Generate a unique color for each alternative promoter group using viridis colormap
    promoter_colors = plt.get_cmap('viridis', len(alt_promoters))
    promoter_color_map = {start_pos: promoter_colors(i) for i, start_pos in enumerate(alt_promoters.keys())}

    # Highlighting cassette exons
    for start, end in exon_ranges_below_one:
        normalized_start = start
        normalized_end = end
        ce_rect = patches.Rectangle(
            (normalized_start, 0), 
            normalized_end - normalized_start, 
            len(sorted_transcript_data)-1, 
            color='red', 
            alpha=0.1, 
            lw=0.5*0.75
        )
        ax.add_patch(ce_rect)
    
    # Highlighting introns
    for start, end in nintrons:
        normalized_start = start
        normalized_end = end
        g_start = start + min_start
        g_end = end + min_start
        ce_rect = patches.Rectangle(
            (g_start, 0), 
            g_end - g_start, 
            len(sorted_transcript_data)-1, 
            color='green', 
            alpha=0.1, 
            lw=0.5*0.75
        )
        ax.add_patch(ce_rect)
        if sum(ncoverage[normalized_start:normalized_end]) / len(ncoverage[normalized_start:normalized_end]) > 0:
            ce_rect = patches.Rectangle(
                (g_start, 0), 
                g_end - g_start, 
                len(sorted_transcript_data)-1, 
                color='green', 
                alpha=0.6, 
                lw=0.5*0.75
            )
            ax.add_patch(ce_rect)

    # Plotting each transcript
    for i, (transcript_id, exons) in enumerate(sorted_transcript_data):
        promoter_color = 'black'  # Default color
        for start_pos, transcripts in alt_promoters.items():
            if transcript_id in transcripts:
                promoter_color = promoter_color_map[start_pos]
                break

        # Drawing the transcript line
        transcript_start = min([exon[0] for exon in exons])
        transcript_end = max([exon[1] for exon in exons])
        ax.plot(
            [transcript_start + 5, transcript_end - 5], 
            [i, i], 
            color=promoter_color, 
            lw=1*0.75
        )

        # Drawing the exons as rectangles
        for exon in exons:
            exon_start, exon_end = exon
            rect = patches.Rectangle(
                (exon_start, i - 0.1), 
                exon_end - exon_start, 
                0.2, 
                color=promoter_color, 
                lw=0.5*0.75
            )
            ax.add_patch(rect)

    # Drawing lines to mark promoter groups
    for start_pos, (start_index, end_index) in promoter_group_indices.items():
        ax.axhline(y=start_index - 0.5, color='black', linestyle='--', lw=0.5*0.75)
        ax.axhline(y=end_index + 0.5, color='black', linestyle='--', lw=0.5*0.75)

    # Highlighting alternative 5' splice sites
    for site in alt5sites:
        normalized_site = site
        if site not in alt_ends:
            ax.axvline(x=normalized_site, color='orange', linestyle='-.', lw=0.5*0.75)
    
    # Highlighting alternative 3' splice sites
    for site in alt3sites:
        normalized_site = site
        if site not in alt_promoters:
            ax.axvline(x=normalized_site, color='red', linestyle=':', lw=0.5*0.75)
    
    # Highlighting alternative ends
    if len(alt_ends) >1:
        for site in alt_ends:
            normalized_site = site
            ax.axvline(x=normalized_site, color='black', linestyle='--', lw=0.5*0.75)

    # Setting labels and titles
    ax.set_yticks(np.arange(len(sorted_transcript_data)))
    ax.set_yticklabels([transcript_id for transcript_id, _ in sorted_transcript_data], fontsize=font_scale_ax)
    ax.set_xlabel('Genomic Position', fontsize=font_scale_ax * 1.5)
    ax.set_ylabel('Transcripts', fontsize=font_scale_ax * 1.5)
    ax.set_xlim(
        min(transcript_data.explode().min()) - 100, 
        max(transcript_data.explode().max()) + 100
    )  # Set x-axis limits
    ax.set_title(f'{gene_id} AS Events', fontsize=font_scale_ax*1.5)
    
    # Remove x-axis tick labels if requested
    if hide_xticks:
        ax.tick_params(axis='x', labelbottom=False)
        ax.set_xticklabels([])
    else:
        ax.tick_params(axis='x', labelsize=font_scale_ax)
        
    if hide_yticks:
        ax.tick_params(axis='y', labelleft=False)
        ax.set_yticklabels([])
    else:
        ax.tick_params(axis='y', labelsize=font_scale_ax)
        
    ax.ticklabel_format(style='plain', axis='x', useOffset=False)
    # Adjust font sizes for better readability
    #for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    #    label.set_fontsize(font_scale_ax)

    # Return the Axes object
    return ax

# Function to calculate scaling factor based on subplot size
def calculate_scale_factor(ax,fig):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    scale = min(width, height) * 25  # Adjust the multiplier as needed for appropriate scaling
    return scale

def adjust_fontsize(ax, scale_factor):
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(scale_factor)


def plot_cov_splice_sites(gene_id, ncoverage, alt5sites, alt3sites, alt_ends, exon_ranges_below_one,
                      nintrons, npeaks, nvalleys, min_start, fig, gs_pos=[0,0],fscale=10):
    """
    Plots splice site data with highlights for alternative splice sites, exon ranges, and introns.

    :param ncoverage: Normalized coverage data.
    :param alt5sites: Alternative 5' splice sites.
    :param alt3sites: Alternative 3' splice sites.
    :param alt_ends: Alternative ends.
    :param exon_ranges_below_one: Exon ranges with coverage below one.
    :param nintrons: Normalized intron ranges.
    :param npeaks: Normalized peak positions.
    :param nvalleys: Normalized valley positions.
    :param min_start: Minimum start position for normalization.
    """
    #fig, ax = plt.subplots(figsize=(12, 8))


    # Determine the width span of gs_pos
    width_span = gs_pos.get_geometry()[3] - gs_pos.get_geometry()[2]

    # Calculate the number of columns proportionally
    base_columns = 4
    new_columns = base_columns * width_span // 4  # Adjust proportionally based on a base span of 4

    nested_gs = gridspec.GridSpecFromSubplotSpec(8, new_columns+1, subplot_spec=gs_pos)

    ax = fig.add_subplot(nested_gs[0:7, 0:new_columns])
    #ax1s3 = fig.add_subplot(nested_gs[4:7, new_columns])

    # Adjust font size for ax1s1
    font_scale_ax1 = calculate_scale_factor(ax,fig) / fscale  # Adjust the division factor as needed
    adjust_fontsize(ax, font_scale_ax1*2)

    ax.plot(ncoverage,lw=0.5*.75)
    # Highlighting alternative 5' splice sites
    for site in alt5sites:
        normalized_site = site - min_start
        ax.axvline(x=normalized_site, color='orange', linestyle='--', lw=0.5*.75)

    # Highlighting alternative 3' splice sites
    for site in alt3sites:
        normalized_site = site - min_start
        ax.axvline(x=normalized_site, color='red', linestyle=':', lw=0.5*.75)

    # Highlighting alternative ends
    for site in alt_ends:
        normalized_site = site - min_start
        ax.axvline(x=normalized_site, color='black', linestyle='-.', lw=0.5*.75)

    # Highlighting cassette exons
    for start, end in exon_ranges_below_one:
        normalized_start = start - min_start
        normalized_end = end - min_start
        ce_rect = patches.Rectangle((normalized_start, 0), normalized_end - normalized_start, 1, color='red', alpha=0.1,lw=0.5*.75)
        ax.add_patch(ce_rect)

    # Highlighting introns
    for start, end in nintrons:
        normalized_start = start
        normalized_end = end
        ce_rect = patches.Rectangle((normalized_start, 0), normalized_end - normalized_start, 1, color='green', alpha=0.1,lw=0.5*.75)
        ax.add_patch(ce_rect)

    # Highlighting retained introns
    for start, end in nintrons:
        normalized_start = start
        normalized_end = end
        if sum(ncoverage[normalized_start:normalized_end]) / len(ncoverage[normalized_start:normalized_end]) > 0:
            ce_rect = patches.Rectangle((normalized_start, 0), normalized_end - normalized_start, 1, color='green', alpha=0.6,lw=0.5*.75)
            ax.add_patch(ce_rect)

    # Plotting peaks and valleys
    ax.plot(npeaks, [ncoverage[p] for p in npeaks], "p", label='Peaks',ms=0.5*.75)
    ax.plot(nvalleys, [ncoverage[v] for v in nvalleys], "v", label='Valleys',ms=0.5*.75)
    ax.set_xlabel('Genomic Position',fontsize=font_scale_ax1*3)
    ax.set_ylabel('Normalized Coverage',fontsize=font_scale_ax1*3)
    ax.set_title(gene_id+' Coverage',fontsize=font_scale_ax1*3)
    return ax

def gtf_alternative_events(gtf_data, gene_id, alt_promoters, alt5sites, alt3sites, alt_ends, exon_ranges_below_one, nintrons, min_start, ncoverage):
    """
    Label the splicing events for the transcripts of a gene, showing transcripts grouped by alternative promoter usage.

    :param gtf_data: DataFrame containing GTF data.
    :param gene_id: The ID of the gene to process.
    :param alt_promoters: A dictionary with alternative promoter start positions as keys and lists of transcript IDs as values.
    :param alt5sites: List of alternative 5' splice sites.
    :param alt3sites: List of alternative 3' splice sites.
    :param alt_ends: List of alternative ends.
    :param exon_ranges_below_one: List of exon ranges with coverage below one.
    :param nintrons: List of intronic regions.
    :param min_start: Minimum start position for normalization.
    :param ncoverage: Normalized coverage data.
    :param max_transcripts: Maximum number of transcripts to process.
    :return: A DataFrame slice of the GTF data with labeled exons for the specific gene_id.
    """
    # Filter GTF data for the specific gene_id and exons
    transcript_data = gtf_data[(gtf_data['gene_id'] == gene_id) & (gtf_data['feature'] == 'exon')]

    # Initialize a dictionary to store labeled exons
    labeled_exons = {}

    promoter_group_indices = {}  # To mark the start and end indices of each promoter group
    current_index = 0

    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in transcript_data:
                sorted_transcript_data.append((transcript_id, transcript_data[transcript_id]))
                current_index += 1

    alt_starts = [c for c in alt_promoters]

    # Drawing lines to mark promoter groups
    for site, (start_index, end_index) in promoter_group_indices.items():
        labeled_exons[site] = labeled_exons.get(site, set())
        labeled_exons[site].add("_alt_pro")
        #ax.axhline(y=start_index - 0.5, color='black', linestyle='--',lw=0.5*.75)
        #ax.axhline(y=end_index + 0.5, color='black', linestyle='--',lw=0.5*.75)

    # Highlighting alternative 5' splice sites
    for site in alt5sites:
        if site not in alt_ends:
            labeled_exons[site] = labeled_exons.get(site, set())
            labeled_exons[site].add("_alt5")

    # Highlighting alternative 3' splice sites
    for site in alt3sites:
        if site not in alt_promoters:
            labeled_exons[site] = labeled_exons.get(site, set())
            labeled_exons[site].add("_alt3")

    # Highlighting alternative ends
    for site in alt_ends:
        labeled_exons[site] = labeled_exons.get(site, set())
        labeled_exons[site].add("_alt_end")

    # Highlighting cassette exons
    for start, end in exon_ranges_below_one:
        for pos in range(start, end + 1):
            labeled_exons[pos] = labeled_exons.get(pos, set())
            labeled_exons[pos].add("_ce")

    # Highlighting retained introns
    for start, end in nintrons:
        for pos in range(start, end + 1):
            labeled_exons[pos] = labeled_exons.get(pos, set())
            labeled_exons[pos].add("_ir")

    # Assign the labeled exon names back to the GTF data
    for idx, row in transcript_data.iterrows():
        exon_start = row['start']
        exon_end = row['end']
        exon_labels = set()

        # Aggregate the labels for the exon
        for pos in range(exon_start, exon_end + 1):
            if pos in labeled_exons:
                exon_labels.update(labeled_exons[pos])

        # Combine the labels and assign to the feature column
        exon_label = "exon" + "".join(sorted(exon_labels))
        transcript_data.at[idx, 'feature'] = exon_label

    return transcript_data


def identify_alternative_ends_per_group(exon_data, end_region_size=1):
    """
    Identify alternative ends for a given set of exons (e.g., promoter group), ensuring unique association of transcripts with ends.

    :param exon_data: DataFrame containing exon data for the transcripts in the promoter group.
    :param end_region_size: The size of the end region to consider around the TES.
    :return: Dictionary of alternative ends and their associated transcripts, and fractional distances.
    """
    # Ensure 'strand' column is present
    if 'strand' not in exon_data.columns:
        raise ValueError("The exon data must contain a 'strand' column.")

    # Initialize dictionaries to store ends and fractional distances
    unique_ends = {}
    fract_dict = {}

    # Group data by transcript_id
    grouped = exon_data.groupby('transcript_id')

    # Iterate over each transcript to find the TES based on strand
    for transcript_id, group in grouped:
        strand = group['strand'].iloc[0]

        if strand == '+':
            # For positive strand, TES is the maximum 'end' position
            tes = group['end'].max()
            end_region = (tes - end_region_size, tes + end_region_size)
        elif strand == '-':
            # For negative strand, TES is the minimum 'start' position
            tes = group['start'].min()
            end_region = (tes - end_region_size, tes + end_region_size)
        else:
            raise ValueError(f"Unknown strand '{strand}' for transcript '{transcript_id}'.")

        # Check if the end region overlaps with any existing end
        found_overlap = False
        for unique_tes, data in unique_ends.items():
            existing_end_region = data['end_region']
            # Check for overlap
            if not (end_region[1] < existing_end_region[0] or end_region[0] > existing_end_region[1]):
                # Overlaps with existing end
                data['transcripts'].add(transcript_id)
                found_overlap = True
                break

        # If no overlap, add as a new unique end
        if not found_overlap:
            unique_ends[tes] = {
                'transcripts': set([transcript_id]),
                'strand': strand,
                'end_region': end_region
            }

    # Now, calculate fractional distances if there are multiple ends
    if len(unique_ends) > 1:
        # Extract TES positions for fractional distance calculation
        tes_positions = list(unique_ends.keys())
        min_tes = min(tes_positions)
        max_tes = max(tes_positions)

        # Calculate fractional distances based on strand
        fractional_distances = {}
        for tes, data in unique_ends.items():
            strand = data['strand']
            if strand == '+':
                fractional_distance = (tes - min_tes) / (max_tes - min_tes)
            else:
                fractional_distance = (max_tes - tes) / (max_tes - min_tes)
            fractional_distances[tes] = fractional_distance

        # Assign fractional distances to transcripts
        for tes, data in unique_ends.items():
            for transcript_id in data['transcripts']:
                fract_dict[transcript_id] = fractional_distances[tes]
    else:
        # Only one end, assign 0 fractional distance to all transcripts
        for data in unique_ends.values():
            for transcript_id in data['transcripts']:
                fract_dict[transcript_id] = 0.0

    # Prepare the final ends dictionary
    ends_dict = {}
    for tes, data in unique_ends.items():
        ends_dict[tes] = data['transcripts']

    return ends_dict, fract_dict

def plot_promoter_group(
    promoter_group_data, 
    promoter_group_id_to_plot, 
    gtf_data, 
    gene_id, 
    ax_ae=None, 
    ax_cov=None, 
    asplot=True, 
    covplot=False, 
    save_path=None,
    fscale1=10,
    fscale2=10,
    hide_xticks=True,
    hide_yticks=False
):
    """
    Plots the alternative events and coverage for a specified promoter group.

    Parameters:
        promoter_group_data (list): List of dictionaries containing promoter group data.
        promoter_group_id_to_plot (str/int): The ID of the promoter group to plot.
        gtf_data (DataFrame): DataFrame containing GTF data.
        gene_id (str/int): The gene ID.
        ax_ae (matplotlib.axes.Axes, optional): Axes for the alternative events plot.
        ax_cov (matplotlib.axes.Axes, optional): Axes for the coverage plot.
        asplot (bool, optional): Whether to plot alternative events. Defaults to True.
        covplot (bool, optional): Whether to plot coverage. Defaults to False.
        save_path (str, optional): Path to save the figure. If None, the figure will be displayed.
        fscale1 (float, optional): Scaling factor for the alternative events plot.
        fscale2 (float, optional): Scaling factor for the coverage plot.
    """
    # Set up font
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial'] + plt.rcParams['font.sans-serif']

    # Find the promoter group data for the specified promoter_group_id
    data = next((d for d in promoter_group_data if d['promoter_group_id'] == promoter_group_id_to_plot), None)
    if data is None:
        print(f"Promoter group ID {promoter_group_id_to_plot} not found.")
        return

    # Extract data
    promoter_group_id = data['promoter_group_id']
    transcripts = data['transcripts']
    alt5sites = data['alt5sites']
    alt3sites = data['alt3sites']
    alt_ends = data['alt_ends']
    exon_ranges_below_one = data['exon_ranges_below_one']
    intronic_regions = data['intronic_regions']
    min_start = data['min_start']
    ncoverage = data['ncoverage']
    npeaks = data['npeaks']
    nvalleys = data['nvalleys']
    promoter_exons = data['promoter_exons']
    alternative_promoters = data['alternative_promoters']
    nintrons = data.get('nintrons', [])
    gene_strand = data['gene_strand']

    # Extract the start position for the alternative promoters
    start_pos = list(alternative_promoters.keys())[0]

    # Initialize figure if axes are not provided
    created_fig = False
    if (asplot and ax_ae is None) or (covplot and ax_cov is None):
        fig = plt.figure(figsize=(10, 8))
        created_fig = True
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.4)
        if asplot and ax_ae is None:
            ax_ae = fig.add_subplot(gs[0, 0])
        if covplot and ax_cov is None:
            ax_cov = fig.add_subplot(gs[1, 0])

    # Plot Alternative Events
    if asplot and ax_ae is not None:
        plot_alternative_events(
            gtf_data, 
            gene_id, 
            {start_pos: transcripts}, 
            alt5sites, 
            alt3sites, 
            alt_ends,
            exon_ranges_below_one, 
            nintrons, 
            min_start, 
            ncoverage, 
            gene_strand, 
            ax=ax_ae,
            fscale=fscale1,
            hide_xticks=False,
            hide_yticks=False
        )
        ax_ae.set_title(f"{gene_id} Alternative Events for Promoter Group {promoter_group_id}",fontsize=fscale1*1.5)

    # Plot Coverage
    if covplot and ax_cov is not None:
        plot_cov_splice_sites(
            gene_id, 
            ncoverage, 
            alt5sites, 
            alt3sites, 
            alt_ends, 
            exon_ranges_below_one,
            nintrons, 
            npeaks, 
            nvalleys, 
            min_start, 
            ax=ax_cov,
            fscale=fscale2
        )
        #ax_cov.set_title(f"Coverage Plot for Promoter Group {promoter_group_id}")

    # Save or show the figure
    if save_path and created_fig:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    elif created_fig:
        plt.show()

    # Close the figure to free memory if it was created within the function
    if created_fig:
        plt.close(fig)




def parse_gtf_attributes(attribute_string):
    """
    Parses the attribute field of a GTF entry into a dictionary.
    """
    attributes = {}
    # Split the string by semicolons
    for attribute in attribute_string.strip().split(';'):
        if attribute.strip() == '':
            continue
        key_value = attribute.strip().split(' ', 1)
        if len(key_value) == 2:
            key = key_value[0]
            value = key_value[1].strip('"')
            attributes[key] = value
    return attributes

def get_splicing_events_V5(gtf_data, gene_id, gtf=True):
    alternative_promoters, alt_promoter_fract_dist = identify_alternative_promoters(gtf_data, gene_id)
    gene_strand = gtf_data[gtf_data['gene_id'] == gene_id]['strand'].iloc[0]

    # Initialize variables to collect results across promoter groups
    transcript_data = {}
    labeled_gtf_data = pd.DataFrame()
    promoter_group_ids = {start_pos: idx for idx, start_pos in enumerate(sorted(alternative_promoters.keys()))}

    # Initialize list to collect data for plotting
    promoter_group_data = []
    # Iterate over each promoter group
    for start_pos, transcripts in alternative_promoters.items():
        promoter_group_id = promoter_group_ids[start_pos]

        # Get exons for transcripts in this promoter group
        exon_ranges, promoter_exons = get_exon_ranges(gtf_data, transcripts)

        # Initialize variables for labeling
        alt5sites = []
        alt3sites = []
        alt_ends = []
        exon_ranges_below_one = []
        intronic_regions = []
        min_start = None
        ncoverage = None
        npeaks = []
        nvalleys = []
        nintrons = []
        promoter_alternative_ends = {}
        promoter_alt_end_fract_dist = {}

        # Handle cases based on the structure of the promoter_exons
        if len(promoter_exons['transcript_id'].unique()) == len(promoter_exons['exon_identifier'].unique()):
            # Single-exon transcripts
            for _, row in promoter_exons.iterrows():
                transcript_id = row['transcript_id']
                if transcript_id not in transcript_data:
                    transcript_data[transcript_id] = {
                        'gene_id': gene_id,
                        'ap': 0, 'ce': 0, 'cet': 0, 'a5': 0,
                        'a5t': 0, 'a3': 0, 'a3t': 0,
                        'ae': 0, 'aet': 0, 'i': 0, 'it': 0,
                        'promoter_group': promoter_group_id
                    }
        elif len(promoter_exons['transcript_id'].unique()) == 1 and len(promoter_exons['exon_identifier'].unique()) > 1:
            # Single transcript with multiple exons
            for _, row in promoter_exons.iterrows():
                transcript_id = row['transcript_id']
                if transcript_id not in transcript_data:
                    transcript_data[transcript_id] = {
                        'gene_id': gene_id,
                        'ap': 0, 'ce': 0, 'cet': 0, 'a5': 0,
                        'a5t': 0, 'a3': 0, 'a3t': 0,
                        'ae': 0, 'aet': 0, 'i': 0, 'it': 0,
                        'promoter_group': promoter_group_id
                    }
        else:
            # Multi-exon, multi-isoform transcripts in this promoter group
            if len(promoter_exons['transcript_id'].unique()) > 1 and len(promoter_exons['exon_identifier'].unique()) > len(promoter_exons['transcript_id'].unique()):
                min_start = promoter_exons['exon_identifier_start'].min()
                max_end = promoter_exons['exon_identifier_end'].max()

                # Normalize coverage for this promoter group
                ncoverage, npeaks, nvalleys = norm_coverage(promoter_exons)
                nintrons = get_intronic_region(promoter_exons)

                # Find alternative splice sites
                alternative_5splice_sites = find_alternative_5_prime_splice_sites(ncoverage, nintrons, npeaks, gene_strand)
                alternative_3splice_sites = find_alternative_3_prime_splice_sites(ncoverage, nintrons, npeaks, gene_strand)

                # Identify cassette exons
                exon_ranges_below_one = get_cassette_exons(promoter_exons)

                # Identify alternative ends for this promoter group
                promoter_alternative_ends, promoter_alt_end_fract_dist = identify_alternative_ends_per_group(promoter_exons)
                alt_ends = list(promoter_alternative_ends.keys())

                # Adjust positions based on min_start
                alt5sites = [site + min_start - 1 for site in alternative_5splice_sites]
                alt3sites = [site + min_start + 1 for site in alternative_3splice_sites]
                intronic_regions = [(start + min_start, end + min_start - 1) for start, end in nintrons]

                # Collect data needed for plotting
                promoter_group_data.append({
                    'promoter_group_id': promoter_group_id,
                    'transcripts': transcripts,
                    'alt5sites': alt5sites,
                    'alt3sites': alt3sites,
                    'alt_ends': alt_ends,
                    'exon_ranges_below_one': exon_ranges_below_one,
                    'intronic_regions': intronic_regions,
                    'min_start': min_start,
                    'ncoverage': ncoverage,
                    'npeaks': npeaks,
                    'nvalleys': nvalleys,
                    'promoter_exons': promoter_exons,
                    'alternative_promoters': {start_pos: transcripts},
                    'alternative_ends': promoter_alternative_ends,
                    'gene_strand': gene_strand,
                    'nintrons': nintrons
                })
                # Iterate through the promoter_exons and compute values
                for _, row in promoter_exons.iterrows():
                    transcript_id = row['transcript_id']
                    exon_start = row['exon_identifier_start']
                    exon_end = row['exon_identifier_end']
                    exon_range = (exon_start, exon_end)

                    if transcript_id not in transcript_data:
                        transcript_data[transcript_id] = {
                            'gene_id': gene_id,
                            'ap': 0, 'ce': 0, 'cet': len(exon_ranges_below_one), 'a5': 0,
                            'a5t': len(alt5sites), 'a3': 0, 'a3t': len(alt3sites),
                            'ae': 0, 'aet': len(promoter_alternative_ends), 'i': 0, 'it': len(intronic_regions),
                            'promoter_group': promoter_group_id
                        }
                    # Compute and aggregate the values
                    transcript_data[transcript_id]['a5'] += exon_end in alt5sites
                    transcript_data[transcript_id]['a3'] += exon_start in alt3sites
                    transcript_data[transcript_id]['i'] += sum(
                        ranges_overlap(exon_range, intronic) for intronic in intronic_regions
                    )
                    transcript_data[transcript_id]['ce'] += exon_range in exon_ranges_below_one

                # Update 'ap' and 'ae' values if there are alternative promoters or ends
                if len(alternative_promoters) > 1:
                    for transcript_id in transcripts:
                        if transcript_id in alt_promoter_fract_dist:
                            transcript_data[transcript_id]['ap'] = alt_promoter_fract_dist[transcript_id]
                if len(promoter_alternative_ends) > 1:
                    for transcript_id in transcripts:
                        if transcript_id in promoter_alt_end_fract_dist:
                            transcript_data[transcript_id]['ae'] = promoter_alt_end_fract_dist[transcript_id]


        # Labeling exons for this promoter group
        if gtf:
            # Get exons and transcripts for this promoter group
            group_gtf_data = gtf_data[
                (gtf_data['gene_id'] == gene_id) &
                (gtf_data['transcript_id'].isin(transcripts)) &
                (gtf_data['feature'].isin(['exon', 'transcript']))
            ].copy()

            # Initialize a dictionary to store labeled exons
            labeled_exons = {}

            # Label alternative 5' splice sites
            if alt5sites:
                for site in alt5sites:
                    labeled_exons[site] = labeled_exons.get(site, set())
                    labeled_exons[site].add("_alt5")

            # Label alternative 3' splice sites
            if alt3sites:
                for site in alt3sites:
                    labeled_exons[site] = labeled_exons.get(site, set())
                    labeled_exons[site].add("_alt3")

            # Label cassette exons
            if exon_ranges_below_one:
                for start, end in exon_ranges_below_one:
                    for pos in range(start, end + 1):
                        labeled_exons[pos] = labeled_exons.get(pos, set())
                        labeled_exons[pos].add("_ce")

            # Label intronic regions
            if intronic_regions:
                for start, end in intronic_regions:
                    for pos in range(start, end + 1):
                        labeled_exons[pos] = labeled_exons.get(pos, set())
                        labeled_exons[pos].add("_ir")

            # Label alternative ends
            if len(promoter_alternative_ends) > 1:
                # Create a mapping of end positions to transcripts
                end_positions = {}
                for transcript_id in promoter_exons['transcript_id'].unique():
                    exons = promoter_exons[promoter_exons['transcript_id'] == transcript_id]
                    if not exons.empty:
                        strand = exons.iloc[0]['strand']
                        if strand == '+':
                            end_pos = exons['end'].max()
                        else:
                            end_pos = exons['start'].min()
                        if end_pos in promoter_alternative_ends:
                            end_positions.setdefault(end_pos, []).append(transcript_id)
                # Now label the last exons
                for end_pos, transcripts_in_end in end_positions.items():
                    for transcript_id in transcripts_in_end:
                        # Get exons for this transcript
                        exons = group_gtf_data[
                            (group_gtf_data['transcript_id'] == transcript_id) &
                            (group_gtf_data['feature'] == 'exon')
                        ]
                        if not exons.empty:
                            # Determine the last exon depending on the strand
                            strand = exons.iloc[0]['strand']
                            if strand == '+':
                                last_exon = exons.sort_values(by='start').iloc[-1]
                            else:
                                last_exon = exons.sort_values(by='end', ascending=False).iloc[-1]
                            # Label only the last exon with '_alt_end'
                            for pos in range(last_exon['start'], last_exon['end'] + 1):
                                labeled_exons[pos] = labeled_exons.get(pos, set())
                                labeled_exons[pos].add("_alt_end")

            # Assign labels back to group_gtf_data
            for idx, row in group_gtf_data.iterrows():
                if 'exon' in row['feature']:
                    exon_start = row['start']
                    exon_end = row['end']
                    exon_labels = set()
                    for pos in range(exon_start, exon_end + 1):
                        if pos in labeled_exons:
                            exon_labels.update(labeled_exons[pos])
                    if exon_labels:
                        exon_label = "exon" + "".join(sorted(exon_labels))
                        group_gtf_data.at[idx, 'feature'] = exon_label
                    else:
                        group_gtf_data.at[idx, 'feature'] = "exon"
                # Assign promoter group ID
                group_gtf_data.at[idx, 'promoter_group'] = promoter_group_id

            # Assign exon numbers within promoter group
            exons = group_gtf_data[group_gtf_data['feature'].str.contains('exon')]
            if not exons.empty:
                # Create a mapping of (transcript_id, exon_start, exon_end) to exon_number
                exon_number_mapping = {}
                for transcript_id, transcript_exons in exons.groupby('transcript_id'):
                    transcript_exons = transcript_exons.sort_values(by=['start', 'end']) if gene_strand == '+' else transcript_exons.sort_values(by=['end', 'start'], ascending=[False, False])
                    exon_numbers = range(1, len(transcript_exons) + 1)
                    for exon_num, (idx_exon, exon_row) in zip(exon_numbers, transcript_exons.iterrows()):
                        exon_number_mapping[(transcript_id, exon_row['start'], exon_row['end'])] = exon_num
                # Assign promoter group exon numbers
                exons_unique = exons.drop_duplicates(subset=['start', 'end'])
                exons_unique = exons_unique.sort_values(by=['start', 'end']) if gene_strand == '+' else exons_unique.sort_values(by=['end', 'start'], ascending=[False, False])
                promoter_group_exon_numbers = range(1, len(exons_unique) + 1)
                exons_unique['promoter_group_exon_number'] = promoter_group_exon_numbers
                exon_num_mapping = exons_unique.set_index(['start', 'end'])['promoter_group_exon_number'].to_dict()
                # Update group_gtf_data
                for idx, row in exons.iterrows():
                    exon_start = row['start']
                    exon_end = row['end']
                    promoter_group_exon_num = exon_num_mapping.get((exon_start, exon_end), None)
                    if promoter_group_exon_num is not None:
                        group_gtf_data.at[idx, 'promoter_group_exon_number'] = promoter_group_exon_num
                    # Get the original exon_number
                    transcript_id = row['transcript_id']
                    exon_num = exon_number_mapping.get((transcript_id, exon_start, exon_end), None)
                    if exon_num is not None:
                        group_gtf_data.at[idx, 'exon_number'] = exon_num
            else:
                group_gtf_data['promoter_group_exon_number'] = None

            # Update 'attribute' field to include 'promoter_group' and 'promoter_group_exon_number'
            for idx, row in group_gtf_data.iterrows():
                if 'exon' in row['feature']:
                    promoter_group_exon_number = row.get('promoter_group_exon_number', None)
                    promoter_group = row.get('promoter_group', None)
                    exon_number = row.get('exon_number', None)
                    # Parse the attribute field into a dictionary
                    attributes = parse_gtf_attributes(row['attribute'])
                    # Preserve the original exon_number
                    if exon_number is not None:
                        attributes['exon_number'] = str(exon_number)
                    if promoter_group_exon_number is not None:
                        attributes['promoter_group_exon_number'] = str(promoter_group_exon_number)
                    if promoter_group is not None:
                        attributes['promoter_group'] = str(promoter_group)
                    # Reconstruct the attribute string
                    attribute_string = "; ".join(f'{key} "{value}"' for key, value in attributes.items()) + ";"
                    group_gtf_data.at[idx, 'attribute'] = attribute_string
                elif row['feature'] == 'transcript':
                    promoter_group = promoter_group_id
                    # Parse the attribute field into a dictionary
                    attributes = parse_gtf_attributes(row['attribute'])
                    if promoter_group is not None:
                        attributes['promoter_group'] = str(promoter_group)
                    # Reconstruct the attribute string
                    attribute_string = "; ".join(f'{key} "{value}"' for key, value in attributes.items()) + ";"
                    group_gtf_data.at[idx, 'attribute'] = attribute_string

            # Append to labeled_gtf_data
            labeled_gtf_data = pd.concat([labeled_gtf_data, group_gtf_data], ignore_index=True)

    # After processing all promoter groups
    if not gtf:
        labeled_gtf_data = pd.DataFrame()  # Return an empty DataFrame if gtf is False

    return transcript_data, labeled_gtf_data, promoter_group_data

def calculate_proximal_distal(up_df, down_df,end_type='ae'):
    """
    Processes the input dataframes to calculate the proportions of 'Proximal' and 'Distal' based on the 'ae' column.

    Parameters:
        up_df (pd.DataFrame): DataFrame for the "up" direction.
        down_df (pd.DataFrame): DataFrame for the "down" direction.
        not_df (pd.DataFrame): DataFrame for the "not" direction.

    Returns:
        pd.DataFrame: A DataFrame containing the proportions of Proximal and Distal for each direction.
    """

    # Dataframes dictionary
    dataframes = {
       r"log_{2}(PS_{FC}) > 0": up_df,
        r"log_{2}(PS_{FC}) < 0": down_df
    }

    # Process dataframes to calculate proportions
    results = {}
    for direction, df in dataframes.items():
        if end_type == 'ae':
            if "ae" in df.columns:
                total = len(df)
                proximal_count = len(df[df['ae'] == 0])
                distal_count = len(df[df['ae'] > 0])
    
                results[direction] = {
                    "Proximal": proximal_count,
                    "Distal": distal_count
                }
        if end_type == 'ap':
            if "ap" in df.columns:
                total = len(df)
                proximal_count = len(df[df['ap'] == 0])
                distal_count = len(df[df['ap'] > 0])
    
                results[direction] = {
                    "Proximal": proximal_count,
                    "Distal": distal_count
                }            

    # Convert results to a DataFrame for visualization
    proportion_df = pd.DataFrame.from_dict(results, orient='index')
    proportion_df.reset_index(inplace=True)
    proportion_df.rename(columns={'index': 'Direction'}, inplace=True)
    proportion_df = proportion_df.set_index("Direction")
    return proportion_df

def custom_horizontal_bar_chart_with_fisher_ae_ap(affected_df, unaffected_df, end_type='ae', fscale=2, ax=None):
    """
    Create a custom bar chart using patches.Rectangle for fine control of colors and layout,
    and perform Fisher's exact test for statistical significance.

    Parameters:
    affected_df (pd.DataFrame): DataFrame with 'Last' and 'No Last' columns for the affected group.
    unaffected_df (pd.DataFrame): DataFrame with 'Last' and 'No Last' columns for the unaffected group.
    fscale: float (optional) - Font scale for customizing chart labels.
    ax (matplotlib.axes.Axes): Matplotlib Axes object for plotting.

    Returns:
    fisher_results_df: pd.DataFrame - Results of the Fisher's exact test for each category.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    # Combine proportions and sort by 'No Last' in affected
    affected_proportions = affected_df.div(affected_df.sum(axis=1), axis=0)
    unaffected_proportions = unaffected_df.div(unaffected_df.sum(axis=1), axis=0)
    sort_index = affected_df["Distal"].sort_values(ascending=False).index
    affected_proportions = affected_proportions.loc[sort_index]
    unaffected_proportions = unaffected_proportions.loc[sort_index]

    # Settings for layout
    categories = affected_proportions.index
    bar_height = 0.4
    y_positions = range(len(categories))

    # Fisher's exact test results
    fisher_results = []

    # Draw bars using patches.Rectangle
    for i, category in enumerate(categories):
        y = i  # y-position for the bar

        # Affected group
        ax.add_patch(patches.Rectangle(
            (0, y - bar_height / 2), -affected_proportions.loc[category, "Proximal"], bar_height,
            facecolor="blue", alpha=0.5, edgecolor="none", label="Affected Proximal" if i == 0 else None))
        ax.add_patch(patches.Rectangle(
            (0, y - bar_height / 2), affected_proportions.loc[category, "Distal"], bar_height,
            facecolor="blue", alpha=1, edgecolor="none", label="Affected Distal" if i == 0 else None))

        # Unaffected group
        ax.add_patch(patches.Rectangle(
            (0, y + bar_height / 2), -unaffected_proportions.loc[category, "Proximal"], bar_height,
            facecolor="orange", alpha=1, edgecolor="none", label="Unaffected Proximal" if i == 0 else None))
        ax.add_patch(patches.Rectangle(
            (0, y + bar_height / 2), unaffected_proportions.loc[category, "Distal"], bar_height,
            facecolor="orange", alpha=0.5, edgecolor="none", label="Unaffected Distal" if i == 0 else None))

        # Fisher's exact test
        contingency_table = [
            [affected_df.loc[category, "Proximal"], affected_df.loc[category, "Distal"]],
            [unaffected_df.loc[category, "Proximal"], unaffected_df.loc[category, "Distal"]]
        ]
        odds_ratio, p_value = fisher_exact(contingency_table)
        fisher_results.append({"Category": category, "Odds Ratio": odds_ratio, "p-value": p_value})

    # Calculate font scale
    font_scale_ax = fscale

    # Customize axes
    #ax.set_aspect('equal')  # Ensure equal aspect ratio
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_yticks([r + (bar_height / 2) for r in y_positions])
    #ax.set_yticklabels(categories, fontsize=font_scale_ax)

    # MODIFY HERE: set the labels to be the categories, but wrap in $ signs
    ax.set_yticklabels([r"${}$".format(cat) for cat in categories], fontsize=font_scale_ax)
    
    # Use axis-relative coordinates for text positioning
    if end_type == 'ae':
        ax.text(0.5, 1.4, "Proportion of Proximal/Distal \nAlternative End usage",
                ha="center", va="top", fontsize=font_scale_ax * 1.5, transform=ax.transAxes)
    if end_type == 'ap':
        ax.text(0.5, 1.4, "Proportion of Proximal/Distal \nAlternative Start usage",
                ha="center", va="top", fontsize=font_scale_ax * 1.5, transform=ax.transAxes)
    
    ax.text(0.25, 1.2, "Proximal", ha="center", va="top",
            fontsize=font_scale_ax * 1.2, transform=ax.transAxes)
    
    ax.text(0.75, 1.2, "Distal", ha="center", va="top",
            fontsize=font_scale_ax * 1.2, transform=ax.transAxes)
    
    ax.set_xlim(-1, 1)
    ax.set_ylim(-0.5, len(categories))
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    custom_ticks = [-1, -0.5, 0, 0.5, 1]
    custom_labels = ["1", "0.5", "0", "0.5", "1"]
    ax.set_xticks(custom_ticks)
    ax.set_xticklabels(custom_labels, fontsize=font_scale_ax)
    ax.invert_yaxis()

    # Add legend at the bottom of the plot
    ax.legend(
        loc="upper center",            # Position legend at the top-center of the bounding box
        bbox_to_anchor=(0.5, -0.10),  # Center the legend below the axis
        ncol=2,                       # Arrange legend items in two columns
        fontsize=font_scale_ax * 1.1,  # Adjust font size for clarity
        frameon=False
    )
    # Convert Fisher's exact test results to DataFrame
    fisher_results_df = pd.DataFrame(fisher_results)

    return fisher_results_df

def get_updown(asevents,direction):
    asevents['log2_fold_change'] = np.log2(asevents['fold_change_mean'])
    upregulated = (asevents['log2_fold_change'] > 0) & (asevents['fold_change_hdi_low'] > 1) & (asevents['p_diff_greater_than_zero'] > 0.95)
    downregulated = (asevents['log2_fold_change'] < 0) & (asevents['fold_change_hdi_high'] < 1) & (asevents['p_diff_less_than_zero'] > 0.95)

    if direction == 'up':
        affected = asevents.loc[asevents.index[upregulated]]
    if direction == 'down':
        affected = asevents.loc[asevents.index[downregulated]]
    if direction == 'not':
        affected = asevents.loc[asevents.index[~(upregulated | downregulated)]]
    if direction == 'both':
        affected = asevents.loc[asevents.index[(upregulated | downregulated)]]
    if direction == 'all':
        affected = asevents
    return affected

# Update the function to add the calculated columns to the original DataFrame
def calculate_and_add_splicing_events(df):
    df["ceI"] = df["ce"]
    df["ceE"] = df["cet"] - df["ce"]
    df["a5I"] = df["a5"]
    df["a5E"] = df["a5t"] - df["a5"]
    df["a3I"] = df["a3"]
    df["a3E"] = df["a3t"] - df["a3"]
    df["iI"] = df["i"]
    df["iE"] = df["it"] - df["i"]
    return df


def binomial_pvalue(up, down):
    """
    Two‑sided exact binomial test for imbalance between up‑ and down‑regulated counts.
    """
    n = up + down
    if n == 0:
        return 1.0
    # scipy ≤1.12: binomtest; older: binom_test
    try:
        return stats.binomtest(up, n, 0.5, alternative='two-sided').pvalue
    except AttributeError:
        return stats.binom_test(up, n, 0.5, alternative='two-sided')

def analyse_panel(path, label):
    """
    Load a *merged_I/E_splicing_results.csv* file, add grey counts,
    run a two‑sided binomial test for Up ≠ Down inside each splice class,
    and return a tidy DataFrame.
    """
    df = pd.read_csv(path, sep=',')
    df['Grey'] = df['Total'] - df['Upregulated'] - df['Downregulated']
    # per‑class binomial test
    df['Direction']   = df.apply(lambda r: 'Up>Down' if r.Upregulated > r.Downregulated
                                 else ('Down>Up' if r.Upregulated < r.Downregulated else 'Equal'), axis=1)
    df['p_two_sided'] = df.apply(lambda r: binomial_pvalue(r.Upregulated, r.Downregulated), axis=1)
    df.insert(0, 'Panel', label)
    return df[['Panel', 'Splicetype', 'Upregulated', 'Downregulated',
               'Grey', 'Total', 'Direction', 'p_two_sided']]


def aggregate_proximal_distal(**named_dfs):
    """
    Combine multiple 'proximal/distal' DataFrames into one tidy DataFrame.

    Parameters
    ----------
    **named_dfs : dict
        Keyword arguments where each key is the variable name (e.g. 
        'ac_ae_prox_dist_affected') and each value is the corresponding
        DataFrame with index ['log_{2}(PS_{FC}) > 0', 'log_{2}(PS_{FC}) < 0']
        and columns ['Proximal', 'Distal'].

    Returns
    -------
    pd.DataFrame
        Long-format table with columns:
        ['panel', 'end_type', 'status', 'direction', 'site', 'count']
    """
    records = []
    pat = re.compile(r'^(ac|bd)_(ae|ap)_prox_dist_(affected|not_affected)$')

    for name, df in named_dfs.items():
        m = pat.match(name)
        if m is None:
            raise ValueError(f"Name {name!r} does not match expected pattern.")
        panel, end_type, status = m.groups()

        # reshape: index -> column, stack Proximal/Distal into rows
        tidy = (df.reset_index()
                  .melt(id_vars='Direction', var_name='site', value_name='count')
                  .rename(columns={'Direction': 'direction'}))
        tidy['panel'] = panel
        tidy['end_type'] = end_type
        tidy['status'] = status
        records.append(tidy)

    return pd.concat(records, ignore_index=True)[
        ['panel', 'end_type', 'status', 'direction', 'site', 'count']
    ]

