import pandas as pd
import numpy as np
import re
from tqdm.notebook import tqdm
from multiprocessing import Pool
import pickle

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec

from matplotlib.patches import ConnectionPatch
from scipy.stats import fisher_exact

import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde, mannwhitneyu
from itertools import combinations



import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy.stats import chi2_contingency


def filter_and_count_features_orf(gtf_data, direction_data, keyword,diag=False):
    """
    Filters the GTF data based on the transcript IDs in the direction file and counts feature occurrences 
    based on the specified keyword ('orf', 'pfam', or 'stop').

    Parameters:
        gtf_data (pd.DataFrame): DataFrame containing GTF data.
        direction_data (pd.DataFrame): DataFrame containing transcript IDs.
        keyword (str): Either 'orf', 'pfam', or 'stop' to filter and count features.

    Returns:
        pd.DataFrame: A DataFrame containing counts of alternative splicing events.
    """
    # Filter GTF data by transcript IDs from the direction file
    filtered_gtf = gtf_data[gtf_data['transcript_id'].isin(direction_data['transcript_id'].values)]

    # Define keywords for alternative splicing events
    keywords = ['ce', 'alt3', 'alt5', 'alt_end', 'ir']

    if keyword.lower() == 'stop':
        # Special filtering for 'stop': exclude entries that contain 'orf'
        filtered_features = filtered_gtf[filtered_gtf['feature'].str.contains("stop", case=False)
            & ~filtered_gtf['feature'].str.contains("orf", case=False) #& ~filtered_gtf['feature'].str.contains("last", case=False)
             #& ~filtered_gtf['feature'].str.contains("codon", case=False)
        ]['feature'].value_counts()

        counts = {kw: filtered_features[filtered_features.index.str.contains(kw)].sum() for kw in keywords}

        if diag:
            return filtered_features
        else:
            return pd.DataFrame([counts])

    if keyword.lower() == 'start':
        # Special filtering for 'stop': exclude entries that contain 'orf'
        filtered_features = filtered_gtf[filtered_gtf['feature'].str.contains("start", case=False) 
            & ~filtered_gtf['feature'].str.contains("orf", case=False) & ~filtered_gtf['feature'].str.contains("codon", case=False)
        ]['feature'].value_counts()
        #filtered_gtf['feature']=="start_pfam"]['feature'].value_counts()
    
        counts = {kw: filtered_features[filtered_features.index.str.contains(kw)].sum() for kw in keywords}

        if diag:
            return filtered_features
        else:
            return pd.DataFrame([counts])
        
    if keyword.lower() == 'orf':
        # Special filtering for 'stop': exclude entries that contain 'orf'
        filtered_features = filtered_gtf[filtered_gtf['feature'].str.contains("orf", case=False) 
            #& ~filtered_gtf['feature'].str.contains("stop", case=False) #& ~filtered_gtf['feature'].str.contains("start", case=False)
        ]['feature'].value_counts()
        #filtered_gtf['feature']=="start_pfam"]['feature'].value_counts()
    
        counts = {kw: filtered_features[filtered_features.index.str.contains(kw)].sum() for kw in keywords}

        if diag:
            return filtered_features
        else:
            return pd.DataFrame([counts])



def perform_chi_squared_test(left_csv_path, right_csv_path, splicing_type):
    # Load the data from the CSV files
    left_data = pd.read_csv(left_csv_path)
    right_data = pd.read_csv(right_csv_path)
    
    # Extract relevant counts for the chi-squared test
    left_counts = left_data.set_index('Event Type')['Count']
    right_counts = right_data.set_index('Event Type')['Count']
    
    # Create contingency table for chi-squared test
    contingency_table = pd.DataFrame({
        'Left': [left_counts[splicing_type], left_counts.sum() - left_counts[splicing_type]],
        'Right': [right_counts[splicing_type], right_counts.sum() - right_counts[splicing_type]]
    })
    
    # Perform chi-squared test
    chi2, p, _, _ = chi2_contingency(contingency_table)
    
    return chi2, p

def perform_chi_squared_test_panel_a(left_csv_path, right_csv_path):
    # Load the data from the CSV files
    left_data = pd.read_csv(left_csv_path)
    right_data = pd.read_csv(right_csv_path)
    
    # Extract relevant counts for the chi-squared test
    left_counts = left_data.set_index('Event Type')['Count']
    right_counts = right_data.set_index('Event Type')['Count']
    
    # Create contingency table for chi-squared test
    contingency_table = pd.DataFrame({
        'Left': left_counts,
        'Right': right_counts
    })
    
    # Perform chi-squared test
    chi2, p, _, _ = chi2_contingency(contingency_table)
    
    return chi2, p
    


def load_csv_as_dataframe(filename):
    return pd.read_csv(filename, index_col=0)

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


def plot_alternative_events(
    gtf_data, 
    gene_id,
    promoter_group_id,
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
    Plot the splicing events for the transcripts of a gene, showing transcripts grouped by alternative promoter usage and overlaying ORF regions.

    Parameters:
        (See original function parameters)
    """
    # If no Axes provided, create a new figure and Axes
    created_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        created_ax = True

    # Extract transcript and ORF data
    transcript_data = gtf_data[
        (gtf_data['gene_id'] == gene_id) & (gtf_data['promoter_group'] == promoter_group_id) & (gtf_data['feature'].str.contains('exon',case=False))].groupby('transcript_id').apply(lambda x: x[['start', 'end']].values.tolist())

    orf_data = gtf_data[(gtf_data['gene_id'] == gene_id) & (gtf_data['promoter_group'] == promoter_group_id) & 
                        (gtf_data['feature'].str.contains('orf',case=False))].groupby('transcript_id',group_keys=True).apply(lambda x: x[['start', 'end']].values.tolist())
                         #& ~(gtf_data['feature'].str.contains('pf',case=False))].groupby('transcript_id').apply(lambda x: x[['start', 'end']].values.tolist())

    pfam_data = gtf_data[(gtf_data['gene_id'] == gene_id) & (gtf_data['promoter_group'] == promoter_group_id) & 
                        (gtf_data['feature'].str.contains('PF',case=False))].groupby('transcript_id',group_keys=True).apply(lambda x: x[['start', 'end']].values.tolist())

    start_data = gtf_data[(gtf_data['gene_id'] == gene_id) & (gtf_data['promoter_group'] == promoter_group_id) & (gtf_data['feature'] == 'start_codon')].groupby('transcript_id',group_keys=True).apply(lambda x: x[['start', 'end']].values.tolist())
    stop_data = gtf_data[(gtf_data['gene_id'] == gene_id) & (gtf_data['promoter_group'] == promoter_group_id) & (gtf_data['feature'] == 'stop_codon')].groupby('transcript_id',group_keys=True).apply(lambda x: x[['start', 'end']].values.tolist())

    pfam_start_data = gtf_data[(gtf_data['gene_id'] == gene_id) & (gtf_data['promoter_group'] == promoter_group_id) & (gtf_data['feature'] == 'start_pfam')].groupby('transcript_id',group_keys=True).apply(lambda x: x[['start', 'end']].values.tolist())
    pfam_stop_data = gtf_data[(gtf_data['gene_id'] == gene_id) & (gtf_data['promoter_group'] == promoter_group_id) & (gtf_data['feature'] == 'stop_pfam')].groupby('transcript_id',group_keys=True).apply(lambda x: x[['start', 'end']].values.tolist())
    
    
    # Sorting transcripts by promoter groups
    sorted_transcript_data = []
    sorted_orf_data = []
    sorted_pfam_data = []
    sorted_start_data = []
    sorted_stop_data = []
    sorted_pfam_start_data = []
    sorted_pfam_stop_data = []
    promoter_group_indices = {}
    current_index = 0

    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in transcript_data:
                sorted_transcript_data.append((transcript_id, transcript_data[transcript_id]))
                current_index += 1
                
    current_index = 0 #reset index
    
    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in orf_data:
                sorted_orf_data.append((transcript_id, orf_data[transcript_id]))
                current_index += 1

    current_index = 0 #reset index

    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in pfam_data:
                sorted_pfam_data.append((transcript_id, pfam_data[transcript_id]))
                current_index += 1

    current_index = 0 #reset index
    
    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in start_data:
                sorted_start_data.append((transcript_id, start_data[transcript_id]))
                current_index += 1

    current_index = 0 #reset index
    
    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in stop_data:
                sorted_stop_data.append((transcript_id, stop_data[transcript_id]))
                current_index += 1

    ##sort pfam_start/pfam_stop
    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in pfam_start_data:
                sorted_pfam_start_data.append((transcript_id, pfam_start_data[transcript_id]))
                current_index += 1

    current_index = 0 #reset index
    
    for start_pos, transcripts in alt_promoters.items():
        promoter_group_indices[start_pos] = [current_index, current_index + len(transcripts) - 1]
        for transcript_id in transcripts:
            if transcript_id in pfam_stop_data:
                sorted_pfam_stop_data.append((transcript_id, pfam_stop_data[transcript_id]))
                current_index += 1

    
    # Limiting the number of transcripts to plot for clarity
    if len(sorted_transcript_data) > max_transcripts:
        sorted_transcript_data = sorted_transcript_data[:max_transcripts]

    # Calculate font scale
    font_scale_ax = fscale

    # Generate a unique color for each alternative promoter group using viridis colormap
    promoter_colors = plt.get_cmap('viridis', len(alt_promoters))
    promoter_color_map = {start_pos: promoter_colors(i) for i, start_pos in enumerate(alt_promoters.keys())}

    # Highlighting cassette exons
    for start, end in exon_ranges_below_one:
        ce_rect = patches.Rectangle(
            (start, 0), end - start, len(sorted_transcript_data) - 1, color='red', alpha=0.1, lw=0.5 * 0.75
        )
        ax.add_patch(ce_rect)

   # Highlighting introns
    for start, end in nintrons:
        normalized_start = start
        normalized_end = end
        g_start = start + min_start
        g_end = end + min_start
        ce_rect = patches.Rectangle((g_start, 0), g_end - g_start, len(sorted_transcript_data)-1, color='green', alpha=0.1,lw=0.5*.75)
        ax.add_patch(ce_rect)
   
        if sum(ncoverage[normalized_start:normalized_end]) / len(ncoverage[normalized_start:normalized_end]) > 0:
          ce_rect = patches.Rectangle((g_start, 0), g_end - g_start, len(sorted_transcript_data)-1, color='green', alpha=0.6,lw=0.5*.75)
          ax.add_patch(ce_rect)

                
    # Plotting each transcript
    for i, (transcript_id, exons) in enumerate(sorted_transcript_data):
        # Default color for promoter
        promoter_color = 'black'
        for start_pos, transcripts in alt_promoters.items():
            if transcript_id in transcripts:
                promoter_color = promoter_color_map[start_pos]
                break
    
        # Drawing the transcript line with a lower zorder to place it in the background
        transcript_start = min([exon[0] for exon in exons])
        transcript_end = max([exon[1] for exon in exons])
        ax.plot([transcript_start + 5, transcript_end - 5], [i, i], color=promoter_color, lw=1 * .75, zorder=1)
    
        # Drawing the exons as rectangles with a higher zorder to place them in the foreground
        for exon in exons:
            exon_start, exon_end = exon
            rect = patches.Rectangle((exon_start, i - .1), exon_end - exon_start, 0.2, color=promoter_color, lw=0.5 * .75, zorder=2)
            ax.add_patch(rect)
    
        # Drawing the ORFs as rectangles with the highest zorder to make sure they are in front
        for trans_id_orf, orf_data in sorted_orf_data:
            if trans_id_orf == transcript_id:
                for orf in orf_data:
                    orf_start, orf_end = orf
                    orf_rect = patches.Rectangle((orf_start, i - .1), orf_end - orf_start, 0.2, color='red', lw=0.5 * .75, zorder=3)
                    ax.add_patch(orf_rect)

        # Drawing the pfam as rectangles with the highest zorder to make sure they are in front
        for trans_id_pfam, pfam_data in sorted_pfam_data:
            if trans_id_pfam == transcript_id:
                for pfam in pfam_data:
                    pfam_start, pfam_end = pfam
                    pfam_rect = patches.Rectangle((pfam_start, i - .1), pfam_end - pfam_start, 0.2, color='blue', lw=0.5 * .75, zorder=4)
                    ax.add_patch(pfam_rect)

    
    

    
    # Drawing lines to mark promoter groups
    for start_pos, (start_index, end_index) in promoter_group_indices.items():
        ax.axhline(y=start_index - 0.5, color='black', linestyle='--', lw=0.5 * 0.75)
        ax.axhline(y=end_index + 0.5, color='black', linestyle='--', lw=0.5 * 0.75)

    # Highlighting splice sites
    for site in alt5sites:
        if site not in alt_ends:
            ax.axvline(x=site, color='orange', linestyle='-.', lw=0.5 * 0.75)
    for site in alt3sites:
        if site not in alt_promoters:
            ax.axvline(x=site, color='red', linestyle=':', lw=0.5 * 0.75)
    if len(alt_ends) > 1: 
        for site in alt_ends:
            ax.axvline(x=site, color='black', linestyle='--', lw=0.5 * 0.75)

    # Setting labels and titles
    
    ax.set_yticks(np.arange(len(sorted_transcript_data)))
    ax.set_yticklabels([transcript_id for transcript_id, _ in sorted_transcript_data], fontsize=font_scale_ax)
    ax.set_xlabel('Genomic Position', fontsize=font_scale_ax * 1.5)
    ax.set_ylabel('Transcripts', fontsize=font_scale_ax * 1.5)
    ax.set_title(f'{gene_id} AS Events', fontsize=font_scale_ax * 1.5)

    # Hide ticks if specified
    if hide_xticks:
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    else:
        ax.tick_params(axis='x', labelsize=font_scale_ax)
        ax.xaxis.get_offset_text().set_fontsize(font_scale_ax)
        
    if hide_yticks:
        ax.tick_params(axis='y', labelleft=False)

    #ax.ticklabel_format(style='plain', axis='x', useOffset=False)

    return ax

import matplotlib.pyplot as plt
from matplotlib import gridspec

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
    hide_xtick=True,
    hide_ytick=False
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
            promoter_group_id,
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
            hide_xticks=hide_xtick,
            hide_yticks=hide_ytick

            
        )
        ax_ae.set_title(f"{gene_id} Alternative Events for Promoter Group {promoter_group_id} ({gene_strand} Strand)",fontsize=fscale1*1.5)

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


def find_all_promoter_groups_by_gene_id(promoter_group_data, gene_id_to_find):
    """
    Searches for a specified gene_id within the promoter_group_data list and returns all matching elements.

    Parameters:
        promoter_group_data (list): List of dictionaries, each containing promoter group data.
        gene_id_to_find (str): The gene ID to search for.

    Returns:
        list: A list of dictionaries from promoter_group_data that contain the specified gene ID.
              Returns an empty list if no matches are found.
    """
    matches = []
    for element in promoter_group_data:
        # Check if 'promoter_exons' and 'gene_id' exist within the element and if there's a match
        unique_gene_ids = element['promoter_exons']['gene_id'].unique()
        if gene_id_to_find in unique_gene_ids:
            matches.append(element)
    
    return matches




def custom_pie_v4(asevents, title, events='splice', event_name=None, event_column=None, table_name=None, event_order=None, innerauto=False, outerauto=True,
                 innerman=False, outerman=False, legend=False, inner_pad=1.0, outer_pad=1.4, fscale=5, axvis=True, ax=None,
                 table=True):
    """
    Draws a custom pie chart with a hole in the center, customizable label positioning, adjustable padding for labels,
    and optional legend. Adjusts label padding to prevent overlap of long text labels with the circle and adds
    leader lines if labels overlap.

    Parameters:
    - asevents: DataFrame containing the events data. For 'splice' events, it should include the following columns:
        - 'ce'
        - 'alt3'
        - 'alt5'
        - 'alt_end'
        - 'ir'
        - 'exon_stop'
    - title: Title of the plot.
    - events: Type of events to plot ('splice' or 'trans').
    - table_name: Name of the CSV file to save the table. If None, defaults based on event type.
    - event_order: List of event names to control the plotting order.
    - innerauto: Automatically places inner labels horizontally if True.
    - outerauto: Automatically places outer labels horizontally if True.
    - innerman: Manually positions inner labels if True.
    - outerman: Manually positions outer labels if True.
    - legend: Adds a legend to the plot if True.
    - inner_pad: Padding factor for inner labels. Default is 1.0 (no padding).
    - outer_pad: Padding factor for outer labels. Default is 1.4 (increased to prevent overlaps).
    - axvis: Adjusts axis visibility. True shows the axis; False hides it.
    - ax: Optional matplotlib Axes object. If not provided, a new figure and axes will be created.
    - table: Saves the table of events to a CSV file if True.
    """
    
    if events == 'splice':
        if event_name and event_column:
            event_names = event_name
            event_columns = event_column
        else:
            # Define the event names and corresponding columns based on the new input
            #event_names = ['ce', 'alt3', 'alt5', 'alt_end', 'ir', 'exon_stop']
            event_names = ['Alternative Exon', 'Alternative 3′ End', 'Alternative 5′ End', 'Alternative Last Exon', 'Intron Retention', 'Constitutive Exon']
            event_columns = ['ce', 'alt3', 'alt5', 'alt_end', 'ir', 'exon_stop']
        
        # Verify that all required columns are present
        missing_cols = [col for col in event_columns if col not in asevents.columns]
        if missing_cols:
            raise ValueError(f"Missing columns in asevents DataFrame: {missing_cols}")
        
        # Calculate event counts by summing each column
        events_counts = asevents[event_columns].sum()
        
        # Calculate percentages
        sizes = events_counts.values
        total = sizes.sum()
        percentages = (sizes / total) * 100
        
        # Create ASsum and AStotal dictionaries with event names
        ASsum = dict(zip(event_names, percentages))
        AStotal = dict(zip(event_names, sizes))
    
    elif events == 'trans':
        # Count unique trans_id for each gene_id
        gene_trans_counts = asevents.groupby('gene_id')['transcript_id'].nunique()
        events_counts = pd.Series({
            "Single Transcript Gene": (gene_trans_counts == 1).sum(),
            "Multi Transcript Gene": (gene_trans_counts > 1).sum()
        })
    
        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100
    
        ASsum = dict(zip(events_counts.index, percentages))
        AStotal = dict(zip(events_counts.index, sizes))
    
    else:
        raise ValueError(f"Unsupported event type: {events}. Supported types are 'splice' and 'trans'.")
    
    # Remove zero percentage values from ASsum and corresponding counts from AStotal
    ASsum_nonzero = {k: v for k, v in ASsum.items() if v > 0}
    AStotal_nonzero = {k: AStotal[k] for k in ASsum_nonzero.keys()}
    
    # If event_order is provided, reorder the dictionaries
    if event_order:
        # Filter out events not in ASsum_nonzero
        event_order = [event for event in event_order if event in ASsum_nonzero]
        ASsum_nonzero = {k: ASsum_nonzero[k] for k in event_order}
        AStotal_nonzero = {k: AStotal_nonzero[k] for k in event_order}
    else:
        # Ensure consistent order for plotting
        ASsum_nonzero = dict(sorted(ASsum_nonzero.items()))
        AStotal_nonzero = dict(sorted(AStotal_nonzero.items()))
    
    # Combine ASsum_nonzero and AStotal_nonzero into a single DataFrame
    ASsum_df = pd.DataFrame({
        'Event Type': list(ASsum_nonzero.keys()),
        'Percentage': list(ASsum_nonzero.values()),
        'Count': list(AStotal_nonzero.values())
    })
    
    # Save ASsum_df to CSV if table=True
    if table:
        if table_name is None:
            table_name = 'splice_events.csv' if events == 'splice' else 'trans_events.csv'
        ASsum_df.to_csv(table_name, index=False)
    
    # Check if ax is None, create a new figure and axes if necessary
    create_new_fig = ax is None
    if create_new_fig:
        fig, ax = plt.subplots(figsize=(10, 10))  # Increased figure size for better label accommodation
        ax.set_aspect('equal')
    
    golden_ratio = (1 + np.sqrt(5)) / 2
    hole_size_relative = 1 / golden_ratio
    outer_radius = 1
    inner_radius = hole_size_relative
    
    total_percentage = sum(ASsum_nonzero.values())
    start_angle = 0
    colors = plt.cm.viridis(np.linspace(0, 1, len(ASsum_nonzero)))
    patches = []  # For legend
    
    # Lists to hold label positions
    label_infos = []
    
    for (label, value), color in zip(ASsum_nonzero.items(), colors):
        slice_angle = (value / total_percentage) * 360
        end_angle = start_angle + slice_angle
    
        mid_angle = np.deg2rad((start_angle + end_angle) / 2)
        is_right_side = np.cos(mid_angle) >= 0  # True if the label is on the right side of the pie
        oh_align = 'left' if is_right_side else 'right'  # Adjusting horizontal alignment based on label position
        ih_align = 'center'  # Center align inner labels
    
        # Draw the pie slice using ax.fill
        theta = np.linspace(start_angle, end_angle, 100)
        x = np.cos(np.deg2rad(theta)) * outer_radius
        y = np.sin(np.deg2rad(theta)) * outer_radius
        wedge = ax.fill(np.concatenate(([0], x, [0])), np.concatenate(([0], y, [0])), color=color)
        if isinstance(wedge, list) and len(wedge) > 0:
            patches.append((wedge[0], label))  # Store patch and label for legend
    
        # Calculate initial label position
        label_radius = outer_radius * outer_pad
        x_label = label_radius * np.cos(mid_angle)
        y_label = label_radius * np.sin(mid_angle)
        label_text = f'{label}: {round(value, 1):.1f}%'
        
        label_infos.append({
            'label': label_text,
            'x': x_label,
            'y': y_label,
            'angle': mid_angle,
            'is_right': is_right_side
        })
    
        start_angle = end_angle

    # Calculate font scale
    font_scale_ax = fscale
    
    # Adjust label positions to prevent overlapping
    def adjust_labels(label_infos, vertical_spacing=0.05):
        """
        Adjust label positions to prevent overlapping by shifting labels vertically.
        """
        # Separate labels into left and right
        left_labels = [label for label in label_infos if not label['is_right']]
        right_labels = [label for label in label_infos if label['is_right']]
        
        def sort_and_adjust(labels, side):
            # Sort labels by y position
            labels_sorted = sorted(labels, key=lambda x: x['y'], reverse=True)
            for i in range(1, len(labels_sorted)):
                if labels_sorted[i]['y'] + vertical_spacing > labels_sorted[i-1]['y']:
                    labels_sorted[i]['y'] = labels_sorted[i-1]['y'] - vertical_spacing
            return labels_sorted
        
        left_labels = sort_and_adjust(left_labels, 'left')
        right_labels = sort_and_adjust(right_labels, 'right')
        
        return left_labels + right_labels
    
    adjusted_labels = adjust_labels(label_infos)
    
    # Draw labels and leader lines
    for label_info in adjusted_labels:
        label = label_info['label']
        x = label_info['x']
        y = label_info['y']
        angle = label_info['angle']
        is_right = label_info['is_right']
        
        # Determine horizontal alignment
        ha = 'left' if is_right else 'right'
        
        # Place the label with a bounding box for better readability
        ax.text(
            x,
            y,
            label,
            ha=ha,
            va='center',
            fontsize=font_scale_ax,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.7)
        )
        
        # Draw a leader line from the pie slice to the label
        # Start point: edge of the outer circle
        line_start_x = outer_radius * np.cos(angle)
        line_start_y = outer_radius * np.sin(angle)
        # End point: near the label position (95% towards the label)
        leader_line_fraction = 0.95
        line_end_x = x * leader_line_fraction
        line_end_y = y * leader_line_fraction
        
        # Connection points for the line
        con = ConnectionPatch(
            (line_start_x, line_start_y),
            (line_end_x, line_end_y),
            "data", "data",
            arrowstyle="-",
            linewidth=0.8,
            color='gray'
        )
        ax.add_artist(con)
    
    # Adjust axes limits based on outer_pad to ensure all labels and leader lines are visible
    max_label_distance = max([np.sqrt(label['x']**2 + label['y']**2) for label in adjusted_labels]) + 0.2  # Additional margin
    ax.set_xlim(-max_label_distance, max_label_distance)
    ax.set_ylim(-max_label_distance, max_label_distance)
    
    # Draw the inner circle to create a donut chart effect
    inner_circle = plt.Circle((0, 0), inner_radius, color='white')
    ax.add_artist(inner_circle)
    ax.set_aspect('equal', adjustable='box')

    
    
    # Add legend if requested and patches are present
    if legend and patches:
        try:
            # Extract patches and labels
            legend_patches, legend_labels = zip(*patches)
            ax.legend(
                legend_patches,
                legend_labels,
                loc='upper center',
                bbox_to_anchor=(0.5, -0.05),
                fancybox=True,
                shadow=True,
                ncol=1,
                frameon=False,
                fontsize=font_scale_ax
            )
        except ValueError:
            # Handle cases where patches might be empty
            print("No patches available to create a legend.")
    
    # Adjust axis visibility
    if not axvis:
        ax.axis('off')
    
    # Set the title
    if title:
        total_count = int(sum(AStotal_nonzero.values()))
        #title_text = f"{title}\n (Total Stop Codons: {total_count})"
        title_text = f"{title}: {total_count})"
        ax.set_title(title_text, fontsize=font_scale_ax*1.5)
    
    # Save the figure if a new figure was created
    if create_new_fig:
        #plt.tight_layout()
        plt.show()

def levenshtein_ratio(s1, s2):
    """
    Calculate the Levenshtein similarity ratio between two strings.
    """
    n, m = len(s1), len(s2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        dp[i][0] = i
    for j in range(m + 1):
        dp[0][j] = j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost)

    distance = dp[n][m]
    return (n + m - distance) / (n + m)


def plot_similarity_heatmap(data, gene_id, promoter_group, fscale = 2, ax=None):
    """
    Generates a heatmap with squares and circles overlaid to represent
    Levenshtein similarity of AA sequences for a specified gene_id and promoter group.

    Parameters:
    - data (pd.DataFrame): The input dataframe containing the AA sequences and metadata.
    - gene_id (str): The gene_id to filter data by.
    - promoter_group (int): The promoter group to filter data by.
    - ax (matplotlib.axes._axes.Axes): The axes to plot on. If None, a new figure and axes are created.

    Returns:
    - None: Displays the plot.
    """
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import Normalize

    # Filter data for the specified gene_id and promoter group
    filtered_data = data[(data['gene_id'] == gene_id) & (data['promoter_group'] == promoter_group)]
    if filtered_data.empty:
        print("No data available for the specified gene_id and promoter group.")
        return

    # Extract AA sequences and transcript_id
    aa_sequences = filtered_data['AA sequence'].tolist()
    labels = filtered_data['transcript_id']

    # Calculate Levenshtein similarity matrix
    n = len(aa_sequences)
    similarity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1):
            similarity_matrix[i, j] = levenshtein_ratio(aa_sequences[i], aa_sequences[j])
            similarity_matrix[j, i] = similarity_matrix[i, j]  # Symmetric matrix

    # Create a new figure and axes if no axes are provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_aspect('equal')  # Ensure equal aspect ratio

    # Create a normalization instance for color scaling
    norm = Normalize(vmin=similarity_matrix.min(), vmax=1)

    # Draw the grid of squares and overlay circles
    for i in range(len(similarity_matrix)):
        for j in range(i + 1):  # Only lower triangle and diagonal
            similarity = similarity_matrix[i, j]

            # Draw the square with a border
            square = patches.Rectangle(
                (j - 0.01, i + 0.01), 1, 1, fill=False, edgecolor="black", linewidth=0.5
            )
            ax.add_patch(square)

            # Overlay the circle
            circle = patches.Circle(
                (j + 0.5, i + 0.5),  # Centered in the square
                radius=0.4,          # Size of the circle
                color=plt.cm.viridis(norm(similarity)),  # Use viridis color map
                ec=None  # No border for the circle
            )
            ax.add_patch(circle)

    # Calculate font scale
    font_scale_ax = fscale
    
    # Configure the axis with transcript_id labels
    ax.set_xlim(0, len(labels))
    ax.set_ylim(0, len(labels))
    ax.set_xticks(np.arange(len(labels)) + 0.5)
    ax.set_yticks(np.arange(len(labels)) + 0.5)
    ax.set_xticklabels(labels, rotation=90,fontsize=font_scale_ax)
    ax.set_yticklabels(labels,fontsize=font_scale_ax)
    ax.invert_xaxis()
    
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.04)
    cbar.ax.tick_params(labelsize=font_scale_ax*1.2)  # Adjust tick font size
    cbar.set_label('Levenshtein Similarity',fontsize=font_scale_ax*1.5)
    ax.set_title(f'AA Sequence Similarity (Gene: {gene_id}, Promoter Group: {promoter_group})', fontsize=font_scale_ax*1.5)
    # Add the colorbar if this is a new figure
    if ax is None:
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.04,)
        cbar.set_label('Levenshtein Similarity')

    # Add a title only if this is a new figure
    #if ax is None:
    #    plt.title(f'AA Sequence Similarity (Gene: {gene_id}, Promoter Group: {promoter_group})', fontsize=16)


def custom_horizontal_bar_chart_with_fisher(affected_df, unaffected_df, fscale=2, ax=None):
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
    sort_index = affected_df["No Last"].sort_values(ascending=False).index
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
            (0, y - bar_height / 2), -affected_proportions.loc[category, "Last"], bar_height,
            facecolor="blue", alpha=0.5, edgecolor="none", label="Affected Last" if i == 0 else None))
        ax.add_patch(patches.Rectangle(
            (0, y - bar_height / 2), affected_proportions.loc[category, "No Last"], bar_height,
            facecolor="blue", alpha=1, edgecolor="none", label="Affected Not Last" if i == 0 else None))

        # Unaffected group
        ax.add_patch(patches.Rectangle(
            (0, y + bar_height / 2), -unaffected_proportions.loc[category, "Last"], bar_height,
            facecolor="orange", alpha=1, edgecolor="none", label="Unaffected Last" if i == 0 else None))
        ax.add_patch(patches.Rectangle(
            (0, y + bar_height / 2), unaffected_proportions.loc[category, "No Last"], bar_height,
            facecolor="orange", alpha=0.5, edgecolor="none", label="Unaffected Not Last" if i == 0 else None))

        # Fisher's exact test
        contingency_table = [
            [affected_df.loc[category, "Last"], affected_df.loc[category, "No Last"]],
            [unaffected_df.loc[category, "Last"], unaffected_df.loc[category, "No Last"]]
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
    ax.set_yticklabels(categories, fontsize=font_scale_ax)
    
    # Use axis-relative coordinates for text positioning
    ax.text(0.5, 1.4, "Proportion of Stop Codon Exon Location\n by Splicing Type",
            ha="center", va="top", fontsize=font_scale_ax * 1.5, transform=ax.transAxes)
    
    ax.text(0.25, 1.2, "Last Exon", ha="center", va="top",
            fontsize=font_scale_ax * 1.2, transform=ax.transAxes)
    
    ax.text(0.75, 1.2, "Not Last Exon", ha="center", va="top",
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






def analyze_direction(gtf_data, direction, keywords=None, event_names=None):
    """
    Analyzes the GTF data for the given direction and computes the keyword breakdown.

    Parameters:
    - gtf_data: pd.DataFrame - GTF data containing transcript and feature information.
    - direction: pd.DataFrame - DataFrame with a column 'transcript_id' to filter GTF data.
    - keywords: list (optional) - List of keywords to analyze. Defaults to ['ce', 'alt3', 'alt5', 'alt_end', 'ir'].
    - event_names: list (optional) - Full names for keywords. Must match the length and order of `keywords`.

    Returns:
    - result_df: pd.DataFrame - DataFrame with counts of occurrences for each keyword.
    """
    if keywords is None:
        keywords = ['ce', 'alt3', 'alt5', 'alt_end', 'ir']
    if event_names is None:
        event_names = ['Alternative Exon', 'Alternative 3′ End', 'Alternative 5′ End', 
                       'Alternative Last Exon', 'Intron Retention']
    
    # Ensure the length of keywords and event_names matches
    if len(keywords) != len(event_names):
        raise ValueError("The length of `keywords` and `event_names` must be the same.")
    
    # Filter GTF data based on the direction
    filt_data = gtf_data[gtf_data['transcript_id'].isin(direction['transcript_id'].values)]
    
    # Count features containing "stop" but not "orf"
    feature_counts = filt_data[
        filt_data['feature'].str.contains("stop", case=False) &
        ~filt_data['feature'].str.contains("orf", case=False)
    ]['feature'].value_counts()
    
    # Filter counts for rows containing "stop"
    stop_counts = feature_counts[feature_counts.index.str.contains('stop')]
    
    # Initialize a dictionary to store keyword breakdown
    keyword_breakdown = {kw: {'last': 0, 'no_last': 0} for kw in keywords}
    
    # Loop through the keywords and count occurrences
    for kw in keywords:
        keyword_last = stop_counts[stop_counts.index.str.contains(kw) & stop_counts.index.str.contains('last')].sum()
        keyword_no_last = stop_counts[stop_counts.index.str.contains(kw) & ~stop_counts.index.str.contains('last')].sum()
        keyword_breakdown[kw]['last'] = keyword_last
        keyword_breakdown[kw]['no_last'] = keyword_no_last
    
    # Convert the breakdown to a DataFrame
    result_df = pd.DataFrame(keyword_breakdown).T
    result_df.columns = ['Last', 'No Last']
    
    # Rename the index with full event names
    result_df.index = event_names
    
    return result_df

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

def perform_mannwhitneyu_test(values1, values2, label1="Group1", label2="Group2"):
    """
    Performs a Mann-Whitney U test between two sets of numeric values.
    Returns a small DataFrame with the statistic and p-value.
    """
    values1_arr = np.array(values1.dropna(), dtype=float)
    values2_arr = np.array(values2.dropna(), dtype=float)
    
    if len(values1_arr) == 0 or len(values2_arr) == 0:
        return pd.DataFrame({
            "Comparison": [f"{label1} vs. {label2}"],
            "U-statistic": [np.nan],
            "p-value": [np.nan]
        })
    
    u_stat, p_val = mannwhitneyu(values1_arr, values2_arr, alternative='two-sided')
    
    return pd.DataFrame({
        "Comparison": [f"{label1} vs. {label2}"],
        "U-statistic": [u_stat],
        "p-value": [p_val]
    })



def find_last_feature(feature_string, match_list):
    """
    Determines which feature from the match list comes last in the feature string.

    Args:
        feature_string (str): The input string containing features separated by '_'.
        match_list (list): A list of features to search for in the string.

    Returns:
        str: The feature from the match list that comes last, or None if no match is found.
    """
    # Split the string by '_'
    split_parts = feature_string.split('_')
    
    # Track the last position and feature
    last_feature = None
    last_position = -1

    # Iterate over the match list to find the last occurrence
    for feature in match_list:
        if feature in split_parts:
            position = split_parts.index(feature)
            if position > last_position:
                last_position = position
                last_feature = feature

    return last_feature



def add_multi_trans_column_to_data(splice_abundance_data, gtf_data):
    """
    This function takes already loaded splice abundance data (TSV) and GTF data,
    identifies multi-isoform genes in the GTF data, and adds a new column 'multi_trans'
    to the splice abundance data.

    Parameters:
        splice_abundance_data (pd.DataFrame): DataFrame containing splice abundance data.
        gtf_data (pd.DataFrame): DataFrame containing GTF data.

    Returns:
        pd.DataFrame: Updated splice abundance DataFrame with the 'multi_trans' column added.
    """
    # Identify multi-isoform genes by grouping by gene_id and counting unique transcript IDs
    multi_isoform_genes = (
        gtf_data.groupby("gene_id")["transcript_id"].nunique()
    )
    multi_isoform_genes = multi_isoform_genes[multi_isoform_genes > 1].index.tolist()

    # Add 'multi_trans' column to splice abundance data based on the gene_id column
    splice_abundance_data['multi_trans'] = splice_abundance_data['gene_id'].apply(
        lambda x: 1 if x in multi_isoform_genes else 0
    )

    return splice_abundance_data



