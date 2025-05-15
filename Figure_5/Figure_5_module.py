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

from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2_contingency

def plot_protein_variants(ax, img_path1, img_path2,title="", label1="", label2="",fscale=5):
    import matplotlib.image as mpimg
    from matplotlib.patches import Rectangle

    # Load images
    img1 = mpimg.imread(img_path1)
    img2 = mpimg.imread(img_path2)

    # Clear axis and hide it
    ax.axis('off')

    # Compute aspect ratio
    img_aspect = img1.shape[1] / img1.shape[0]
    img_display_height = 0.9
    img_display_width = img_display_height * img_aspect

    # Set limits to accommodate both images side-by-side
    ax.set_xlim(0, 2 * img_display_width + 0.2)
    ax.set_ylim(0, img_display_height + 0.1)

    # Plot first image
    left_x = 0.05
    left_y = 0.05
    ax.imshow(img1, extent=[left_x, left_x + img_display_width, left_y, left_y + img_display_height])
    ax.add_patch(Rectangle((left_x, left_y), img_display_width, img_display_height, linewidth=0.5, edgecolor='black', facecolor='none'))

    # Plot second image
    right_x = left_x + img_display_width + 0.1
    right_y = 0.05
    ax.imshow(img2, extent=[right_x, right_x + img_display_width, right_y, right_y + img_display_height])
    ax.add_patch(Rectangle((right_x, right_y), img_display_width, img_display_height, linewidth=0.5, edgecolor='black', facecolor='none'))

    # Add text labels
    ax.set_title(title,
                 fontsize=fscale*1.5)
    ax.text(left_x + img_display_width / 2, left_y - 0.05, label1, ha='center', va='top', fontsize=8)
    ax.text(right_x + img_display_width / 2, right_y - 0.05, label2, ha='center', va='top', fontsize=8)


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
    
def perform_chi_squared_test(left_csv_path, right_csv_path, splicing_type):
    # Load the data from the CSV files
    left_data = pd.read_csv(left_csv_path)
    right_data = pd.read_csv(right_csv_path)
    
    # Extract relevant counts for the chi-squared test
    left_counts = left_data.set_index('Event Type')['Count']
    right_counts = right_data.set_index('Event Type')['Count']
    
    # Use .get() to default missing splice types to zero
    left_splice = left_counts.get(splicing_type, 0)
    right_splice = right_counts.get(splicing_type, 0)
    
    # Create contingency table for chi-squared test
    contingency_table = pd.DataFrame({
        'Left': [left_splice, left_counts.sum() - left_splice],
        'Right': [right_splice, right_counts.sum() - right_splice]
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
        ~filt_data['feature'].str.contains("pfam", case=False)
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


def process_and_filter_pfam_features(gtf_data, ac_direction, diag=False, filter_option=None):
    """
    Processes a GTF file to extract and filter PFAM features, then filters based on transcript IDs from ac_direction.
    Ensures multiple PFAM IDs per transcript are retained.
    Adds binary columns for specific splicing event keywords and optionally filters and sums based on 'start' or 'stop'.

    Parameters:
        gtf_data (pd.DataFrame): DataFrame containing GTF data.
        ac_direction (pd.DataFrame): DataFrame containing transcript IDs to filter.
        filter_option (str, optional): If 'start' or 'stop', filters the DataFrame and returns column sums.

    Returns:
        pd.DataFrame: Filtered PFAM features with added keyword indicators or aggregated sums if 'start' or 'stop' is specified.
    """
    
    # Filter rows where the feature contains "PFAM" (case-sensitive)
    pfam_filtered = gtf_data[gtf_data["feature"].str.contains("PFAM", case=True, na=False)].copy()

    # Extract PFAM ID and exon feature
    pfam_filtered[["pfam_id", "exon_feature"]] = pfam_filtered["feature"].str.extract(r'(PFAM_PF\d+\.\d+)_(.*)')

    # Remove "PFAM_" prefix from pfam_id
    pfam_filtered["pfam_id"] = pfam_filtered["pfam_id"].str.replace("PFAM_", "", regex=False)

    # Drop unwanted exon features and ensure we're working with a copy
    pfam_filtered = pfam_filtered[~pfam_filtered["exon_feature"].isin(["exon_start", "exon_stop", "exon"])].copy()

    # Filter PFAM features based on transcript IDs in ac_direction (keep all PFAM IDs)
    filtered_pfam = pfam_filtered[pfam_filtered["transcript_id"].isin(ac_direction["transcript_id"])].copy()

    # Ensure no unexpected dropping of duplicate PFAM IDs per transcript
    filtered_pfam = filtered_pfam.sort_values(by=["transcript_id", "pfam_id"]).reset_index(drop=True)

    # Define the keyword columns to be added
    keywords = ['stop', 'start', 'ce', 'alt3', 'alt5', 'alt_end', 'ir']

    # Add new columns based on keyword presence in exon_feature
    for kw in keywords:
        filtered_pfam.loc[:, kw] = filtered_pfam['exon_feature'].apply(lambda x: 1 if kw in x else 0)

    # Select only required columns
    filtered_pfam = filtered_pfam[["gene_id", "transcript_id", "pfam_id", "exon_feature"] + keywords]

    if diag:
        return filtered_pfam

    # If 'start' or 'stop' is specified, filter and return summed counts of specific columns
    if filter_option in ['start', 'stop']:
        filtered_subset = filtered_pfam[filtered_pfam[filter_option] == 1]
        summed_values = filtered_subset[['ce', 'alt3', 'alt5', 'alt_end', 'ir']].sum().to_frame().T
        return summed_values

    return filtered_pfam[['ce', 'alt3', 'alt5', 'alt_end', 'ir']].sum().to_frame().T

def process_and_filter_pfam_features(gtf_data, ac_direction, diag=False, filter_option=None):
    """
    Processes a GTF file to extract and filter PFAM features, then filters based on transcript IDs from ac_direction.
    Ensures multiple PFAM IDs per transcript are retained.
    Adds binary columns for specific splicing event keywords and optionally filters and sums based on 'start' or 'stop'.

    Parameters:
        gtf_data (pd.DataFrame): DataFrame containing GTF data.
        ac_direction (pd.DataFrame): DataFrame containing transcript IDs to filter.
        filter_option (str, optional): If 'start' or 'stop', filters the DataFrame and returns column sums.

    Returns:
        pd.DataFrame: Filtered PFAM features with added keyword indicators or aggregated sums if 'start' or 'stop' is specified.
    """
    
    # Filter rows where the feature contains "PFAM" (case-sensitive)
    pfam_filtered = gtf_data[gtf_data["feature"].str.contains("PFAM", case=True, na=False)].copy()

    # Extract PFAM ID and exon feature
    pfam_filtered[["pfam_id", "exon_feature"]] = pfam_filtered["feature"].str.extract(r'(PFAM_PF\d+\.\d+)_(.*)')

    # Remove "PFAM_" prefix from pfam_id
    pfam_filtered["pfam_id"] = pfam_filtered["pfam_id"].str.replace("PFAM_", "", regex=False)

    # Drop unwanted exon features and ensure we're working with a copy
    pfam_filtered = pfam_filtered[~pfam_filtered["exon_feature"].isin(["exon_start", "exon_stop", "exon"])].copy()

    # Filter PFAM features based on transcript IDs in ac_direction (keep all PFAM IDs)
    filtered_pfam = pfam_filtered[pfam_filtered["transcript_id"].isin(ac_direction["transcript_id"])].copy()

    # Ensure no unexpected dropping of duplicate PFAM IDs per transcript
    filtered_pfam = filtered_pfam.sort_values(by=["transcript_id", "pfam_id"]).reset_index(drop=True)

    # Define the keyword columns to be added
    keywords = ['stop', 'start', 'ce', 'alt3', 'alt5', 'alt_end', 'ir']

    # Add new columns based on keyword presence in exon_feature
    for kw in keywords:
        filtered_pfam.loc[:, kw] = filtered_pfam['exon_feature'].apply(lambda x: 1 if kw in x else 0)

    # Select only required columns
    filtered_pfam = filtered_pfam[["gene_id", "transcript_id", "pfam_id", "exon_feature"] + keywords]

    if diag:
        return filtered_pfam

    # If 'start' or 'stop' is specified, filter and return summed counts of specific columns
    if filter_option in ['start', 'stop']:
        filtered_subset = filtered_pfam[filtered_pfam[filter_option] == 1]
        summed_values = filtered_subset[['ce', 'alt3', 'alt5', 'alt_end', 'ir']].sum().to_frame().T
        return summed_values

    return filtered_pfam[['ce', 'alt3', 'alt5', 'alt_end', 'ir']].sum().to_frame().T


# Define the function to update ac_pfam_hits with target_name and description
def annotate_pfam_hits(ac_pfam_hits_df, pfam_hits_df):
    """
    Annotates ac_pfam_hits dataframe with target_name and description from pfam_hits dataframe.

    Parameters:
    ac_pfam_hits_df (pd.DataFrame): DataFrame containing pfam_id to be matched.
    pfam_hits_df (pd.DataFrame): DataFrame containing target_accession, target_name, and description.

    Returns:
    pd.DataFrame: Updated ac_pfam_hits_df with target_name and description columns.
    """
    # Drop unnecessary columns from pfam_hits
    pfam_hits_df = pfam_hits_df.iloc[:, 2:]

    # Drop the first column from ac_pfam_hits if it exists
    ac_pfam_hits_df = ac_pfam_hits_df.iloc[:, 1:]

    # Create a dictionary for quick lookup of target_name and description
    pfam_lookup = {row['target_accession']: (row['target_name'], row['description']) 
                   for _, row in pfam_hits_df.iterrows()}

    # Function to add target_name and description row by row
    def add_pfam_annotations(row):
        target_name, description = pfam_lookup.get(row['pfam_id'], (None, None))
        row['target_name'] = target_name
        row['description'] = description
        return row

    # Apply the function row-wise without modifying the index or removing duplicates
    return ac_pfam_hits_df.apply(add_pfam_annotations, axis=1)


def plot_clustered_heatmap_v2(
    dataframe,
    event_name=None,
    event_column=None,
    event_order=None,
    title=None,
    fscale=5,
    ax_main=None,
    ax_colorbar=None,
    show_dendrogram=True,
    export_tsv=False,
    tsv_filename="heatmap_data_normalized.tsv"
):
    """
    Generates a clustered heatmap of normalized feature counts per unique target_name,
    with an optional dendrogram showing the hierarchical clustering of rows.
    Optionally exports the clustered, normalized DataFrame to a TSV file.

    Parameters:
    -----------
    dataframe (pd.DataFrame):
        Input DataFrame with 'target_name' and feature columns.
    event_name (List[str]):
        Labels for feature columns, e.g. ['Alternative Exon', ...].
    event_column (List[str]):
        Corresponding column names in the DataFrame.
    event_order (List[str]):
        Desired order of labels on the x-axis (matches event_name).
    title (str):
        Title for the heatmap.
    fscale (float):
        Font and line scale factor.
    ax_main (plt.Axes):
        Optional axis for the heatmap.
    ax_colorbar (plt.Axes):
        Optional axis for the colorbar.
    show_dendrogram (bool):
        Whether to display the dendrogram on the left side.
    export_tsv (bool):
        If True, exports the final clustered normalized DataFrame to a TSV file.
    tsv_filename (str):
        TSV filename to use for exporting the data.

    Returns:
    --------
    fig : matplotlib.figure.Figure
    ax_dendro : matplotlib.axes._axes.Axes or None
    ax_main : matplotlib.axes._axes.Axes
    ax_colorbar : matplotlib.axes._axes.Axes
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy.cluster.hierarchy as sch
    from matplotlib import gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    lw = fscale / 10

    # Default feature columns if not provided
    if event_column is None:
        event_column = ["ce", "alt3", "alt5", "alt_end"]
    if event_name is None:
        event_name = event_column
    if event_order is None:
        event_order = event_name

    # Map column names to event names and reorder columns as desired
    column_label_map = dict(zip(event_column, event_name))
    ordered_columns = [col for label in event_order for col, lbl in column_label_map.items() if lbl == label]

    # Prepare and normalize data
    heatmap_data = dataframe.groupby("target_name")[event_column].sum()
    heatmap_data_normalized = heatmap_data.div(heatmap_data.sum(axis=1), axis=0).fillna(0)
    heatmap_data_normalized = heatmap_data_normalized[ordered_columns]

    # Perform hierarchical clustering on rows
    linkage = sch.linkage(heatmap_data_normalized, method='ward')
    dendro_order = sch.leaves_list(linkage)
    heatmap_data_clustered = heatmap_data_normalized.iloc[dendro_order]

    # Optionally export the clustered normalized data to TSV
    if export_tsv:
        heatmap_data_clustered.to_csv(tsv_filename, sep='\t', index=True)

    # Set up figure and axes
    if ax_main is None:
        if show_dendrogram:
            # Create a new figure with three columns: dendrogram, heatmap, and colorbar
            fig = plt.figure(figsize=(12, 12))
            gs = gridspec.GridSpec(
                nrows=1, ncols=3, width_ratios=[0.2, 0.7, 0.1], wspace=0.05
            )
            ax_dendro = fig.add_subplot(gs[0])
            ax_main = fig.add_subplot(gs[1])
            ax_colorbar = fig.add_subplot(gs[2])
        else:
            fig = plt.figure(figsize=(10, 12))
            gs = gridspec.GridSpec(
                nrows=1, ncols=2, width_ratios=[0.9, 0.05], wspace=0.05
            )
            ax_main = fig.add_subplot(gs[0])
            ax_colorbar = fig.add_subplot(gs[1])
            ax_dendro = None
    else:
        fig = ax_main.figure
        if show_dendrogram:
            # Create an inset dendrogram axis when a main axis is provided
            ax_dendro = inset_axes(
                ax_main,
                width="15%",
                height="100%",
                loc='upper left',
                bbox_to_anchor=(-0.2, 0, 1, 1),
                bbox_transform=ax_main.transAxes,
                borderpad=0
            )
        else:
            ax_dendro = None
        if ax_colorbar is None:
            ax_colorbar = inset_axes(
                ax_main,
                width="3%",
                height="15%",
                loc="lower left",
                bbox_to_anchor=(1.02, 0, 1, 1),
                bbox_transform=ax_main.transAxes,
                borderpad=0
            )

    # Plot dendrogram if required
    if show_dendrogram and ax_dendro is not None:
        sch.dendrogram(
            linkage,
            orientation='left',
            no_labels=True,
            ax=ax_dendro,
            color_threshold=0
        )
        ax_dendro.invert_yaxis()

        ax_dendro.set_xticks([])
        ax_dendro.set_yticks([])
        ax_dendro.axis('off')
        # Set dendrogram line styles after drawing
        plt.setp(ax_dendro.collections, linewidth=0.5, linestyle='-', color='black')

    # Plot heatmap
    sns.heatmap(
        heatmap_data_clustered,
        cmap="coolwarm",
        linewidths=0,
        cbar=True,
        cbar_ax=ax_colorbar,
        ax=ax_main,
        annot=False
    )

    # Set x-axis labels and remove y ticks and labels
    ax_main.set_xticklabels(event_order, rotation=45, ha='right', fontsize=fscale)
    ax_main.set_xlabel(None)
    ax_main.set_ylabel("Pfam Domains", fontsize=fscale*1.2,labelpad=20)

    # Remove y ticks and y tick labels
    ax_main.set_yticks([])
    ax_main.set_yticklabels([])

    if title:
        ax_main.set_title(title, fontsize=fscale * 1.5)

    ax_main.tick_params(labelsize=fscale)
    ax_colorbar.tick_params(labelsize=fscale, width=lw)

    for ax in [ax_main, ax_colorbar]:
        if ax is not None:
            for spine in ax.spines.values():
                spine.set_linewidth(lw)

    # Apply tight_layout only if ax_main is the first axis (i.e., new figure)
    if ax_main is fig.axes[0]:
        fig.tight_layout()

    return fig, ax_dendro, ax_main, ax_colorbar


def plot_alternative_events_v2(
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
    hide_yticks=False,
    transcript_space=False
):
    created_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        created_ax = True

    # 1) transcript_data: only true exons, no "orf_" labels, deduped & sorted
    transcript_data = (
        gtf_data
        .loc[
            (gtf_data['gene_id'] == gene_id) &
            (gtf_data['promoter_group'] == promoter_group_id) &
            (gtf_data['feature'].str.contains('exon', case=False)) &
            (~gtf_data['feature'].str.contains('orf',  case=False))
        ]
        .groupby('transcript_id', group_keys=False)
        .apply(lambda df: sorted({tuple(x) for x in df[['start','end']].values},
                                 key=lambda x: x[0]))
    ).to_dict()

    # 2) orf_data: all ORF‐exon segments, deduped & sorted
    orf_data = (
        gtf_data
        .loc[
            (gtf_data['gene_id'] == gene_id) &
            (gtf_data['promoter_group'] == promoter_group_id) &
            (gtf_data['feature'].str.contains('orf', case=False))
        ]
        .groupby('transcript_id', group_keys=False)
        .apply(lambda df: sorted({tuple(x) for x in df[['start','end']].values},
                                 key=lambda x: x[0]))
    ).to_dict()

    # 3) pfam, start/stop codons, pfam_start/stop
    pfam_data = (
        gtf_data
        .loc[
            (gtf_data['gene_id'] == gene_id) &
            (gtf_data['promoter_group'] == promoter_group_id) &
            (gtf_data['feature'].str.contains('PF', case=False))
        ]
        .groupby('transcript_id', group_keys=False)
        .apply(lambda df: sorted({tuple(x) for x in df[['start','end']].values},
                                 key=lambda x: x[0]))
    ).to_dict()

    start_data = (
        gtf_data
        .loc[
            (gtf_data['gene_id'] == gene_id) &
            (gtf_data['promoter_group'] == promoter_group_id) &
            (gtf_data['feature'] == 'start_codon')
        ]
        .groupby('transcript_id', group_keys=False)
        .apply(lambda df: sorted({tuple(x) for x in df[['start','end']].values},
                                 key=lambda x: x[0]))
    ).to_dict()

    stop_data = (
        gtf_data
        .loc[
            (gtf_data['gene_id'] == gene_id) &
            (gtf_data['promoter_group'] == promoter_group_id) &
            (gtf_data['feature'] == 'stop_codon')
        ]
        .groupby('transcript_id', group_keys=False)
        .apply(lambda df: sorted({tuple(x) for x in df[['start','end']].values},
                                 key=lambda x: x[0]))
    ).to_dict()

    pfam_start_data = (
        gtf_data
        .loc[
            (gtf_data['gene_id'] == gene_id) &
            (gtf_data['promoter_group'] == promoter_group_id) &
            (gtf_data['feature'] == 'start_pfam')
        ]
        .groupby('transcript_id', group_keys=False)
        .apply(lambda df: sorted({tuple(x) for x in df[['start','end']].values},
                                 key=lambda x: x[0]))
    ).to_dict()

    pfam_stop_data = (
        gtf_data
        .loc[
            (gtf_data['gene_id'] == gene_id) &
            (gtf_data['promoter_group'] == promoter_group_id) &
            (gtf_data['feature'] == 'stop_pfam')
        ]
        .groupby('transcript_id', group_keys=False)
        .apply(lambda df: sorted({tuple(x) for x in df[['start','end']].values},
                                 key=lambda x: x[0]))
    ).to_dict()

    # 4) sort transcripts by promoter group & limit
    sorted_transcript_data = []
    for p_start, txs in alt_promoters.items():
        for tx in txs:
            if tx in transcript_data:
                sorted_transcript_data.append((tx, transcript_data[tx]))
    sorted_transcript_data = sorted_transcript_data[:max_transcripts]
   
    # 5) (optional) collapse to transcript‐space
    if transcript_space:
        # build per‐tx exon → t‐coord map
        tx_exon_map = {}
        for tx, exons in sorted_transcript_data:
            cum = 0
            mapped = []
            for g0, g1 in sorted(exons, key=lambda x: x[0]):
                L = g1 - g0
                mapped.append((cum, cum + L, g0, g1))
                cum += L
            tx_exon_map[tx] = mapped

        def g2t(pos, mapping):
            for t0, t1, g0, g1 in mapping:
                if g0 <= pos <= g1:
                    return t0 + (pos - g0)
            return None

        # rebuild exon lists in t‐space
        sorted_transcript_data = [
            (tx, [(t0, t1) for t0, t1, _, _ in tx_exon_map[tx]])
            for tx, _ in sorted_transcript_data
        ]
        

        def map_feats(raw_dict):
            out = []
            for tx, feats in raw_dict.items():
                if tx not in tx_exon_map:
                    continue
                vals = []
                for g0, g1 in feats:
                    t0, t1 = g2t(g0, tx_exon_map[tx]), g2t(g1, tx_exon_map[tx])
                    if t0 is not None and t1 is not None:
                        vals.append((t0, t1))
                if vals:
                    out.append((tx, vals))
            return out
            
        def fix_transcript_space(transcripts,strand):
            """
            For each (tx_id, exons) in transcripts, if strand=='-',
            flips each interval around the transcript length and re-sorts.
            
            transcripts: list of (tx_id, [(start, end), ...])
            strand: '+' or '-'
            
            Returns a new list of (tx_id, fixed_exons).
            """
            fixed_list = []
            for tx_id, exons in transcripts:
                # always sort input exons by start first
                exons = sorted(exons, key=lambda iv: iv[0])
                if strand == '-':
                    # total transcript length = max end
                    L = max(end for _, end in exons)
                    # flip each interval, then sort by new start
                    flipped = sorted(
                        [(L - end, L - start) for start, end in exons],
                        key=lambda iv: iv[0]
                    )
                    fixed_list.append((tx_id, flipped))
                else:
                    # plus-strand: just keep sorted order
                    fixed_list.append((tx_id, exons))
            return fixed_list



        def fix_transcript_space_with_features(exon_list, feature_list, strand):
            """
            exon_list:   [(tx_id, [(s,e),…]), …]
            feature_list: [(tx_id, [(s,e),…]), …]
            strand:      '+' or '-'
            
            Returns two lists:
              fixed_exons, fixed_features
            each in the same [(tx_id, intervals), …] form.
            """
            # build a lookup for your features
            feat_map = {tx: feats for tx, feats in feature_list}
            
            fixed_exons   = []
            fixed_features = []
            
            for tx, exons in exon_list:
                # sort input just in case
                exons = sorted(exons, key=lambda iv: iv[0])
                feats = sorted(feat_map.get(tx, []), key=lambda iv: iv[0])
                
                if strand == '-':
                    L = max(e for _, e in exons)
                    # flip exons
                    fx = sorted([(L - end, L - start) for start, end in exons],
                                key=lambda iv: iv[0])
                    # flip features
                    ff = sorted([(L - end, L - start) for start, end in feats],
                                key=lambda iv: iv[0])
                else:
                    fx = exons
                    ff = feats
                
                fixed_exons.append((tx, fx))
                fixed_features.append((tx, ff))
            
            return fixed_exons, fixed_features
        
        sorted_transcript_data,_ = fix_transcript_space_with_features(sorted_transcript_data,map_feats(orf_data),strand)
        _,sorted_orf_data        = fix_transcript_space_with_features(sorted_transcript_data,map_feats(orf_data),strand)
        _,sorted_pfam_data       = fix_transcript_space_with_features(sorted_transcript_data,map_feats(pfam_data),strand)
        _,sorted_start_data      = fix_transcript_space_with_features(sorted_transcript_data,map_feats(start_data),strand)
        _,sorted_stop_data       = fix_transcript_space_with_features(sorted_transcript_data,map_feats(stop_data),strand)
        _,sorted_pfam_start_data = fix_transcript_space_with_features(sorted_transcript_data,map_feats(pfam_start_data),strand)
        _,sorted_pfam_stop_data  = fix_transcript_space_with_features(sorted_transcript_data,map_feats(pfam_stop_data),strand)
        

    else:
        # genomic‐space: just list()
        sorted_orf_data        = list(orf_data.items())
        sorted_pfam_data       = list(pfam_data.items())
        sorted_start_data      = list(start_data.items())
        sorted_stop_data       = list(stop_data.items())
        sorted_pfam_start_data = list(pfam_start_data.items())
        sorted_pfam_stop_data  = list(pfam_stop_data.items())

    # 6) plotting (unchanged)
    promoter_colors = plt.get_cmap("viridis", len(alt_promoters))
    pcm = {p: promoter_colors(i) for i, p in enumerate(alt_promoters)}

    for i, (tx, exons) in enumerate(sorted_transcript_data):
        col = next((pcm[p] for p,txs in alt_promoters.items() if tx in txs), "black")

        # backbone
        tmin, tmax = min(e[0] for e in exons), max(e[1] for e in exons)
        ax.plot([tmin, tmax], [i, i], color=col, lw=1, zorder=1)

        # exon boxes
        for t0, t1 in exons:
            ax.add_patch(patches.Rectangle((t0, i-0.1), t1-t0, 0.2, color=col, zorder=2))

        # ORF
        for tx2, feats in sorted_orf_data:
            if tx2==tx:
                for t0,t1 in feats:
                    ax.add_patch(patches.Rectangle((t0, i-0.1), t1-t0, 0.2,
                                                  color="red", zorder=3))

        # PFAM
        for tx2, feats in sorted_pfam_data:
            if tx2==tx:
                for t0,t1 in feats:
                    ax.add_patch(patches.Rectangle((t0, i-0.1), t1-t0, 0.2,
                                                  color="blue", zorder=4))

    ax.set_yticks(np.arange(len(sorted_transcript_data)))
    ax.set_yticklabels([tx for tx,_ in sorted_transcript_data], fontsize=fscale)
    ax.set_xlabel(
        "Transcript coordinate (nt)" if transcript_space else "Genomic position",
        fontsize=fscale*1.5
    )
    ax.tick_params(axis='x', labelsize=fscale)
    ax.set_ylabel("Transcripts", fontsize=fscale*1.5)
    ax.set_title(f"{gene_id} AS events ({'transcript' if transcript_space else 'genomic'} space)",
                 fontsize=fscale*1.5)

    if hide_xticks:
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
    if hide_yticks:
        ax.tick_params(axis="y", labelleft=False)

    if created_ax:
        plt.show()
        plt.close(ax.figure)

    return ax


def plot_promoter_group_v2(
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
    hide_ytick=False,
    transcript_space=False
):
    data = next((d for d in promoter_group_data
                 if d['promoter_group_id']==promoter_group_id_to_plot), None)
    if data is None:
        raise ValueError(f"Promoter group {promoter_group_id_to_plot} not found")

    # unpack
    transcripts           = data['transcripts']
    alt5sites             = data['alt5sites']
    alt3sites             = data['alt3sites']
    alt_ends              = data['alt_ends']
    exon_ranges_below_one = data['exon_ranges_below_one']
    nintrons              = data.get('nintrons', [])
    min_start             = data['min_start']
    ncoverage             = data['ncoverage']
    gene_strand           = data['gene_strand']
    alternative_promoters = data['alternative_promoters']

    created = False
    if (asplot and ax_ae is None) or (covplot and ax_cov is None):
        fig = plt.figure(figsize=(10, 8))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.4)
        if asplot: ax_ae = fig.add_subplot(gs[0])
        if covplot: ax_cov = fig.add_subplot(gs[1])
        created = True

    if asplot:
        plot_alternative_events_v2(
            gtf_data, gene_id, promoter_group_id_to_plot,
            {list(alternative_promoters.keys())[0]: transcripts},
            alt5sites, alt3sites, alt_ends,
            exon_ranges_below_one, nintrons, min_start,
            ncoverage, gene_strand,
            ax=ax_ae,
            fscale=fscale1,
            hide_xticks=hide_xtick,
            hide_yticks=hide_ytick,
            transcript_space=transcript_space
        )
        ax_ae.set_title(
            f"{gene_id} Alternative Events for PG {promoter_group_id_to_plot}",
            fontsize=fscale1*1.5
        )

    if created and save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    elif created:
        plt.show()

    return ax_ae, ax_cov


