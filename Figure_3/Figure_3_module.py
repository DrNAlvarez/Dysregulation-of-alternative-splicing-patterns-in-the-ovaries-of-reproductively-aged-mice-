import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlinescd
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
from scipy.stats import fisher_exact
import pandas as pd
from scipy.stats import chi2_contingency


def custom_pie(asevents, title,events='splice', innerauto=False, outerauto=False, innerman=False, outerman=False, legend=False,
               inner_pad=1.0, outer_pad=1.2, axvis=True, ax=None):
    """
    Draws a custom pie chart with a hole in the center, customizable label positioning, adjustable padding for labels, and optional legend.
    Adjusts label padding to prevent overlap of long text labels with the circle.

    Parameters:
    - value_dict: A dictionary where keys are labels and values are percentages of each slice.
    - innerauto: Automatically places inner labels horizontally if True.
    - outerauto: Automatically places outer labels horizontally if True.
    - innerman: Manually positions inner labels if True.
    - outerman: Manually positions outer labels if True.
    - legend: Adds a legend to the plot if True.
    - inner_pad: Padding factor for inner labels. Default is 1.0 (no padding).
    - outer_pad: Padding factor for outer labels. Default is 1.2 (slightly outside the outer circle).
    - axvis: Adjusts axis visibility. True shows the axis; False hides it.
    - ax: Optional matplotlib Axes object. If not provided, a new figure and axes will be created.
    """
    if events == 'splice':
        events_counts = pd.concat([asevents[['ce', 'a5', 'a3', 'i']].sum(),
                                               ((asevents[['ap']] > 0) & (asevents[['ap']] < 1)).sum(),
                                   ((asevents[['ae']] > 0) & (asevents[['ae']] < 1)).sum()],axis=0)
        sizes = events_counts.values
        ASsum=dict(zip(['Alternative Exon','Alternative 5′ End','Alternative 3′ End','Intron Retention',
                           'Alternative First Exon','Alternative Last Exon'],[round(a,1) for a in (sizes/sizes.sum())*100]))
    if events == 'trans':
        # Count unique trans_id for each gene_id
        gene_trans_counts = asevents.groupby('gene_id')['transcript_id'].nunique()
        events_counts = pd.Series({"single":(gene_trans_counts == 1).sum(),"multi": (gene_trans_counts > 1).sum()})
        sizes = events_counts.values
        ASsum=dict(zip(['Single Transcript Gene','Multi Transcript Gene'],[round(a,1) for a in (sizes/sizes.sum())*100]))

    # Check if ax is None, create a new figure and axes if necessary
    create_new_fig = ax is None
    if create_new_fig:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')

    golden_ratio = (1 + np.sqrt(5)) / 2
    hole_size_relative = 1 / golden_ratio
    outer_radius = 1
    inner_radius = hole_size_relative

    total = sum(ASsum.values())
    start_angle = 0
    colors = plt.cm.viridis(np.linspace(0, 1, len(ASsum)))
    patches = []  # For legend

    for (label, value), color in zip(ASsum.items(), colors):
        slice_angle = (value / total) * 360
        end_angle = start_angle + slice_angle

        mid_angle = np.deg2rad((start_angle + end_angle) / 2)
        is_right_side = np.cos(mid_angle) > 0  # True if the label is on the right side of the pie
        oh_align = 'left' if is_right_side else 'right'  # Adjusting horizontal alignment based on label position
        ih_align = 'center' if is_right_side else 'center'  # Adjusting horizontal alignment based on label position

        start_rad = np.deg2rad(start_angle)
        end_rad = np.deg2rad(end_angle)
        rads = np.linspace(start_rad, end_rad, 100)

        x = np.cos(rads)
        y = np.sin(rads)

        patch = ax.fill(np.append(x, 0), np.append(y, 0), color=color)
        patches.append((patch[0], label))  # Store patch and label for legend

        # Outer labels with adjusted alignment
        if outerauto and not legend:
            ax.text(outer_radius * np.cos(mid_angle) * outer_pad, outer_radius * np.sin(mid_angle) * outer_pad, label, horizontalalignment=oh_align, verticalalignment='center', fontsize=6)

        # Inner labels with adjusted alignment
        if innerauto:
            ax.text(inner_radius * np.cos(mid_angle) * inner_pad, inner_radius * np.sin(mid_angle) * inner_pad, f'{value}%', horizontalalignment=ih_align, verticalalignment='center', fontsize=6)

        start_angle = end_angle

    inner_circle = plt.Circle((0, 0), inner_radius, color='white')
    ax.add_artist(inner_circle)
    ax.set_aspect('equal', adjustable='box')

    if legend:
        legend_patches, legend_labels = zip(*patches)
        ax.legend(legend_patches, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=1, frameon=False)

    if not axvis:
        ax.axis('off')

    if title and events == 'trans':
        title = title+"\n (Total Genes: "+str(sizes.sum())+")"
        ax.set_title(title, fontsize=8)

    if title and events == 'splice':
        title = title+"\n (Total Events: "+str(sizes.sum())+")"
        ax.set_title(title, fontsize=8)


def custom_pie_v2(asevents, title,events='splice',table_name=None, innerauto=False, outerauto=False, innerman=False, outerman=False, legend=False,
               inner_pad=1.0, outer_pad=1.2, axvis=True, ax=None, table=True):
    """
    Draws a custom pie chart with a hole in the center, customizable label positioning, adjustable padding for labels, and optional legend.
    Adjusts label padding to prevent overlap of long text labels with the circle.

    Parameters:
    - value_dict: A dictionary where keys are labels and values are percentages of each slice.
    - innerauto: Automatically places inner labels horizontally if True.
    - outerauto: Automatically places outer labels horizontally if True.
    - innerman: Manually positions inner labels if True.
    - outerman: Manually positions outer labels if True.
    - legend: Adds a legend to the plot if True.
    - inner_pad: Padding factor for inner labels. Default is 1.0 (no padding).
    - outer_pad: Padding factor for outer labels. Default is 1.2 (slightly outside the outer circle).
    - axvis: Adjusts axis visibility. True shows the axis; False hides it.
    - ax: Optional matplotlib Axes object. If not provided, a new figure and axes will be created.
    """

    # Apply log2 transformation to the 'fold_change_mean' column
    #asevents['log2_fold_change'] = np.log2(asevents['fold_change_mean'])

    # Define conditions for significantly affected transcripts
    #upregulated = (asevents['log2_fold_change'] > 0) & (asevents['fold_change_hdi_low'] > 1) & (asevents['p_diff_greater_than_zero'] > 0.95)
    #downregulated = (asevents['log2_fold_change'] < 0) & (asevents['fold_change_hdi_high'] < 1) & (asevents['p_diff_greater_than_zero'] <= 0.05)

    # Find indices for transcripts significantly affected in AvsB, CvsD, and both
    #asevents_up = asevents.loc[asevents.index[upregulated]
    #asevents_down = asevents.loc[asevents.index[downregulated]
    #not_affected = asevents.loc[~(upregulated | downregulated)]

    if events == 'splice':
        events_counts = pd.concat([asevents[['ce', 'a5', 'a3', 'i']].sum(),
                                               ((asevents[['ap']] > 0) & (asevents[['ap']] < 1)).sum(),
                                   ((asevents[['ae']] > 0) & (asevents[['ae']] < 1)).sum()],axis=0)
        sizes = events_counts.values
        ASsum=dict(zip(['Alternative Exon','Alternative 5′ End','Alternative 3′ End','Intron Retention',
                           'Alternative First Exon','Alternative Last Exon'],[round(a,1) for a in (sizes/sizes.sum())*100]))

        AStotal=dict(zip(['Alternative Exon','Alternative 5′ End','Alternative 3′ End','Intron Retention',
                           'Alternative First Exon','Alternative Last Exon'],[a for a in sizes]))
    if events == 'trans':

        # Count unique trans_id for each gene_id
        gene_trans_counts = asevents.groupby('gene_id')['transcript_id'].nunique()
        events_counts = pd.Series({"single":(gene_trans_counts == 1).sum(),"multi": (gene_trans_counts > 1).sum()})
        sizes = events_counts.values
        ASsum=dict(zip(['Single Transcript Gene','Multi Transcript Gene'],[round(a,1) for a in (sizes/sizes.sum())*100]))
        AStotal=dict(zip(['Alternative Exon','Alternative 5′ End','Alternative 3′ End','Intron Retention',
                           'Alternative First Exon','Alternative Last Exon'],[a for a in sizes]))

    # Calculate percentages
    #percentages = [round((value / sum(sizes)) * 100, 1) for value in sizes]

   # Combine ASsum and AStotal into a single DataFrame
    ASsum_df = pd.DataFrame({
        'Event Type': ASsum.keys(),
        'Percentage': ASsum.values(),
        'Count': AStotal.values()
    })

    # Save ASsum_df to CSV if table=True
    if table:
        if table_name is None:
            table_name = 'splice_events.csv' if events == 'splice' else 'trans_events.csv'
        ASsum_df.to_csv(table_name, index=False)
    
    # Check if ax is None, create a new figure and axes if necessary
    create_new_fig = ax is None
    if create_new_fig:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')

    golden_ratio = (1 + np.sqrt(5)) / 2
    hole_size_relative = 1 / golden_ratio
    outer_radius = 1
    inner_radius = hole_size_relative

    total = sum(ASsum.values())
    start_angle = 0
    colors = plt.cm.viridis(np.linspace(0, 1, len(ASsum)))
    patches = []  # For legend

    for (label, value), color in zip(ASsum.items(), colors):
        slice_angle = (value / total) * 360
        end_angle = start_angle + slice_angle

        mid_angle = np.deg2rad((start_angle + end_angle) / 2)
        is_right_side = np.cos(mid_angle) > 0  # True if the label is on the right side of the pie
        oh_align = 'left' if is_right_side else 'right'  # Adjusting horizontal alignment based on label position
        ih_align = 'center' if is_right_side else 'center'  # Adjusting horizontal alignment based on label position

        start_rad = np.deg2rad(start_angle)
        end_rad = np.deg2rad(end_angle)
        rads = np.linspace(start_rad, end_rad, 100)

        x = np.cos(rads)
        y = np.sin(rads)

        patch = ax.fill(np.append(x, 0), np.append(y, 0), color=color)
        patches.append((patch[0], label))  # Store patch and label for legend

        # Outer labels with adjusted alignment
        if outerauto and not legend:
            ax.text(outer_radius * np.cos(mid_angle) * outer_pad, outer_radius * np.sin(mid_angle) * outer_pad, label, horizontalalignment=oh_align, verticalalignment='center', fontsize=4)

        # Inner labels with adjusted alignment
        if innerauto:
            ax.text(inner_radius * np.cos(mid_angle) * inner_pad, inner_radius * np.sin(mid_angle) * inner_pad, f'{round(value,0):.0f}%', horizontalalignment=ih_align, verticalalignment='center', fontsize=3)

        start_angle = end_angle

    inner_circle = plt.Circle((0, 0), inner_radius, color='white')
    ax.add_artist(inner_circle)
    ax.set_aspect('equal', adjustable='box')

    if legend:
        legend_patches, legend_labels = zip(*patches)
        ax.legend(legend_patches, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=1, frameon=False)

    if not axvis:
        ax.axis('off')

    if title and events == 'trans':
        
        title = title+"\n (Total Genes: "+str(sizes.sum())+")"
        ax.set_title(title, fontsize=8)

    if title and events == 'splice':
        
        title = title+"\n (Total Events: "+str(sizes.sum())+")"
        ax.set_title(title, fontsize=8)

def custom_pie_v3(asevents, title, events='splice', table_name=None, event_order=None, innerauto=False, outerauto=False, innerman=False, outerman=False, legend=False,
                  inner_pad=1.0, outer_pad=1.2, axvis=True, ax=None, table=True):
    """
    Draws a custom pie chart with a hole in the center, customizable label positioning, adjustable padding for labels, and optional legend.
    Adjusts label padding to prevent overlap of long text labels with the circle.

    Parameters:
    - asevents: DataFrame containing the events data.
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
    - outer_pad: Padding factor for outer labels. Default is 1.2 (slightly outside the outer circle).
    - axvis: Adjusts axis visibility. True shows the axis; False hides it.
    - ax: Optional matplotlib Axes object. If not provided, a new figure and axes will be created.
    - table: Saves the table of events to a CSV file if True.
    """

    if events == 'splice_total':
        # Define default event names and their corresponding columns
        event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
                       'Alternative First Exon', 'Alternative Last Exon']
        event_columns = ['cet', 'a5t', 'a3t', 'it', 'ap', 'ae']
        #event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
        #               'Alternative Last Exon']
        #event_columns = ['ce', 'a5', 'a3', 'i', 'ae']

        # Calculate event counts
        events_counts = pd.Series({
            'ce': asevents['cet'].sum(),
            'a5': asevents['a5t'].sum(),
            'a3': asevents['a3t'].sum(),
            'i': asevents['it'].sum(),
            'ap': ((asevents['ap'] >= 0) & (asevents['ap'] <= 1)).sum(),
            'ae': ((asevents['ae'] >= 0) & (asevents['ae'] <= 1)).sum()
        })

        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100

        # Create ASsum and AStotal dictionaries with event names
        ASsum = dict(zip(event_names, percentages))
        AStotal = dict(zip(event_names, sizes))

    if events == 'splice_I':
        # Define default event names and their corresponding columns
        event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
                       'Alternative First Exon', 'Alternative Last Exon']
        event_columns = ['ceI', 'a5I', 'a3I', 'iI', 'ap', 'ae']
        #event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
        #               'Alternative Last Exon']
        #event_columns = ['ce', 'a5', 'a3', 'i', 'ae']

        # Calculate event counts
        events_counts = pd.Series({
            'ce': asevents['ceI'].sum(),
            'a5': asevents['a5I'].sum(),
            'a3': asevents['a3I'].sum(),
            'i': asevents['iI'].sum(),
            'ap': ((asevents['ap'] > 0) & (asevents['ap'] <= 1)).sum(),
            'ae': ((asevents['ae'] > 0) & (asevents['ae'] <= 1)).sum()
        })

        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100

        # Create ASsum and AStotal dictionaries with event names
        ASsum = dict(zip(event_names, percentages))
        AStotal = dict(zip(event_names, sizes))
      
    if events == 'splice_E':
        # Define default event names and their corresponding columns
        event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
                       'Alternative First Exon', 'Alternative Last Exon']
        event_columns = ['ceE', 'a5E', 'a3E', 'iE', 'ap', 'ae']
        #event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
        #               'Alternative Last Exon']
        #event_columns = ['ce', 'a5', 'a3', 'i', 'ae']

        # Calculate event counts
        events_counts = pd.Series({
            'ce': asevents['ceE'].sum(),
            'a5': asevents['a5E'].sum(),
            'a3': asevents['a3E'].sum(),
            'i': asevents['iE'].sum(),
            'ap': (asevents['ap'] == 0).sum(),
            'ae': (asevents['ae'] == 0).sum()
        })

        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100

        # Create ASsum and AStotal dictionaries with event names
        ASsum = dict(zip(event_names, percentages))
        AStotal = dict(zip(event_names, sizes))
    
    elif events == 'trans':
           # Reduce to unique gene_id and multi_trans status
        gene_mult_trans = asevents[['gene_id', 'multi_trans']].drop_duplicates()
        
        # Determine multi-transcript status for each gene_id
        # A gene is multi-transcript if it has any transcript marked as multi_trans = 1
        gene_mult_trans = gene_mult_trans.groupby('gene_id')['multi_trans'].max()
    
        # Calculate single and multi-transcript gene counts
        events_counts = pd.Series({
            "Single Transcript Gene": (gene_mult_trans == 0).sum(),
            "Multi Transcript Gene": (gene_mult_trans == 1).sum()
        })
        
        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100

        ASsum = dict(zip(events_counts.index, percentages))
        AStotal = dict(zip(events_counts.index, sizes))

 
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
        fig, ax = plt.subplots()
        ax.set_aspect('equal')

    golden_ratio = (1 + np.sqrt(5)) / 2
    hole_size_relative = 1 / golden_ratio
    outer_radius = 1
    inner_radius = hole_size_relative

    total = sum(ASsum_nonzero.values())
    start_angle = 0
    colors = plt.cm.viridis(np.linspace(0, 1, len(ASsum_nonzero)))
    patches = []  # For legend

    for (label, value), color in zip(ASsum_nonzero.items(), colors):
        slice_angle = (value / total) * 360
        end_angle = start_angle + slice_angle

        mid_angle = np.deg2rad((start_angle + end_angle) / 2)
        is_right_side = np.cos(mid_angle) > 0  # True if the label is on the right side of the pie
        oh_align = 'left' if is_right_side else 'right'  # Adjusting horizontal alignment based on label position
        ih_align = 'center'  # Center align inner labels

        start_rad = np.deg2rad(start_angle)
        end_rad = np.deg2rad(end_angle)
        rads = np.linspace(start_rad, end_rad, 100)

        x = np.cos(rads)
        y = np.sin(rads)

        patch = ax.fill(np.append(x, 0), np.append(y, 0), color=color)
        patches.append((patch[0], label))  # Store patch and label for legend

        # Outer labels with adjusted alignment
        if outerauto and not legend:
            ax.text(outer_radius * np.cos(mid_angle) * outer_pad, outer_radius * np.sin(mid_angle) * outer_pad,
                    label, horizontalalignment=oh_align, verticalalignment='center', fontsize=4)

        # Inner labels with adjusted alignment
        if innerauto:
            ax.text(inner_radius * np.cos(mid_angle) * inner_pad, inner_radius * np.sin(mid_angle) * inner_pad,
                    f'{round(value, 1):.1f}%', horizontalalignment=ih_align, verticalalignment='center', fontsize=3)

        start_angle = end_angle

    inner_circle = plt.Circle((0, 0), inner_radius, color='white')
    ax.add_artist(inner_circle)
    ax.set_aspect('equal', adjustable='box')

    if legend:
        legend_patches, legend_labels = zip(*patches)
        ax.legend(legend_patches, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=1, frameon=False)

    if not axvis:
        ax.axis('off')

    if title:
        total_count = int(sum(AStotal_nonzero.values()))
        title_text = f"{title}\n (Total {'Events' if 'splice' in events else 'Genes'}: {total_count})"
        ax.set_title(title_text, fontsize=8)


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

def create_splicetype_plot(df, ax, title, splicetype,xlab=True,ylab=True,export=False):
    """
    Plot the fold change with credible intervals, adjusting for specified splicetype.

    Parameters:
    df (pd.DataFrame): DataFrame containing the differential expression analysis results.
    ax: The axis object to plot on.
    title (str): Title for the plot.
    splicetype (str): The type of splicing event to consider ('ce', 'a5', 'a3', 'i').
    """
    # Ensure splicetype is one of the expected values
    if splicetype not in ['ce', 'a5', 'a3', 'i','ap','ae']:
        raise ValueError("splicetype must be one of: 'ce', 'a5', 'a3', 'i'")

    # Check if splicetype column exists
    if splicetype not in df.columns:
        raise ValueError(f"{splicetype} column is missing from the DataFrame")

    # Duplicate rows based on the splicetype value
    if splicetype in ['ce', 'a5', 'a3', 'i']:
        df = df.reindex(df.index.repeat(df[splicetype])).reset_index(drop=True)

    if splicetype in ['ap','ae']:
        df = df[df[splicetype]>0]
    # Continue with existing plotting code
    df['log_fold_change_mean'] = np.log2(df['fold_change_mean'])

    # Define conditions for upregulated and downregulated transcripts
    upregulated_condition = (df['log_fold_change_mean'] > 0) & (df['fold_change_hdi_low'] > 1) & (df['p_diff_greater_than_zero'] > 0.95)
    downregulated_condition = (df['log_fold_change_mean'] < 0) & (df['fold_change_hdi_high'] < 1) & (df['p_diff_greater_than_zero'] < 0.05)

    df['color'] = 'grey'  # Default color for transcripts not strongly up or downregulated
    df.loc[upregulated_condition, 'color'] = 'red'
    df.loc[downregulated_condition, 'color'] = 'blue'
    df = df.sort_values(by='log_fold_change_mean', ascending=True)

    if export:
        df.to_csv(splicetype+'.csv',index=True)
    
    # Plotting
    ax.errorbar(x=df['log_fold_change_mean'], y=np.arange(df.shape[0]),
                xerr=[abs(df['log_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])),
                      abs(np.log2(df['fold_change_hdi_high']) - df['log_fold_change_mean'])],
                fmt='.', ecolor=df['color'], color='none', elinewidth=1, alpha=0.1)

    for category, color in zip(["Up regulated", "Down regulated", "Not affected"], ["red", "blue", "grey"]):
        condition = df['color'] == color
        ax.scatter(x=df[condition]['log_fold_change_mean'], y=np.arange(df.shape[0])[condition],
                   color=color, s=0.5, alpha=0.5)

    # Custom legend handles
    legend_handles = [mlines.Line2D([], [], color='red', marker='o', linestyle='None', markersize=6, label='Up regulated'),
                      mlines.Line2D([], [], color='blue', marker='o', linestyle='None', markersize=6, label='Down regulated'),
                      mlines.Line2D([], [], color='grey', marker='o', linestyle='None', markersize=6, label="Not affected"),
                      mlines.Line2D([], [], color='grey', marker='none', linestyle='-', markersize=6, label="Credible Interval")]

    #min_max = [abs(df['log_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])).max(),
    #                  abs(np.log2(df['fold_change_hdi_high']) - df['log_fold_change_mean']).max()]
    #min_max = [abs(df['log_fold_change_mean'] - np.log2(df['fold_change_hdi_low'])).max(),
    #                  abs(np.log2(df['fold_change_hdi_high']) + df['log_fold_change_mean']).max()]
    #max = np.log2(df.iloc[-1]['fold_change_hdi_high'])
    #print(min_max)
    ax.set_xlim(-2, 2)
    ax.tick_params(axis='y',labelsize=6)
    ax.tick_params(axis='x',labelsize=6)
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

    # Check and remove any tick within 50 units of the end_tick
    new_ticks = [tick for tick in new_ticks if abs(tick - end_tick) > 50]

    # Append the end_tick if it's not already the last tick
    if not new_ticks or new_ticks[-1] != end_tick:
        new_ticks.append(end_tick)

    ax.set_yticks(new_ticks)  # Apply the modified ticks


    if xlab:
        ax.set_xlabel('Fold Change (log2)',fontsize=8)
    if ylab:
        ax.set_ylabel('Transcripts',fontsize=8)
    if title:
        ax.set_title(title, fontsize=8)
    #ax.legend(handles=legend_handles, fontsize=6)
    ax.axvline(x=0, linestyle='--', color='k', alpha=0.7)

def create_dist_plot(df_avsb, ax, title, splicetype,xlab=True,ylab=True):
    """
    Plot the fold change with credible intervals, adjusting for specified splicetype.

    Parameters:
    df (pd.DataFrame): DataFrame containing the differential expression analysis results.
    ax: The axis object to plot on.
    title (str): Title for the plot.
    splicetype (str): The type of splicing event to consider ('ce', 'a5', 'a3', 'i').
    """

    if splicetype in ['ap','ae']:
        df_avsb = df_avsb[(df_avsb[splicetype] > 0) & (df_avsb[splicetype] < 1)]

    # Apply log2 transformation to the 'fold_change_mean' column
    df_avsb['log2_fold_change'] = np.log2(df_avsb['fold_change_mean'])

    # Define conditions for significantly affected transcripts
    upregulated_avsb = (df_avsb['log2_fold_change'] > 0) & (df_avsb['fold_change_hdi_low'] > 1) & (df_avsb['p_diff_greater_than_zero'] > 0.95)
    downregulated_avsb = (df_avsb['log2_fold_change'] < 0) & (df_avsb['fold_change_hdi_high'] < 1) & (df_avsb['p_diff_greater_than_zero'] <= 0.05)

    # Define conditions for not significantly affected transcripts
    not_affected_avsb = ~(upregulated_avsb | downregulated_avsb)

    # Find indices for transcripts significantly affected in AvsB, CvsD, and both
    up_avsb = df_avsb.index[upregulated_avsb]
    down_avsb = df_avsb.index[downregulated_avsb]
    affected_only_in_avsb = df_avsb.index[upregulated_avsb | downregulated_avsb]

    # Find indices for transcripts not significantly affected in both
    not_affected_in_both = df_avsb.index[not_affected_avsb]
    #print(not_affected_in_both)

     # Adjust the plot_ecdf function to accept the ax parameter
    def plot_ecdf(data, label, color, ax):
        sns.ecdfplot(data=data, ax=ax, label=label, color=color)
        #sns.kdeplot(data=data, ax=ax, label=label, color=color)

    # Ensure ax is defined
    if ax is None:
        fig, ax = plt.subplots()

    # Plot ECDF for all transcripts and specific groups
    plot_ecdf(df_avsb.loc[not_affected_in_both, splicetype], 'All Transcripts', 'lightgrey', ax)
    #plot_ecdf(df_avsb.loc[affected_only_in_avsb, splicetype], 'Significant', 'green', ax)
    plot_ecdf(df_avsb.loc[up_avsb, splicetype], 'Significant', 'red', ax)
    plot_ecdf(df_avsb.loc[down_avsb, splicetype], 'Significant', 'blue', ax)


    # Adjust the plot aesthetics
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    #ax.legend()
    #ax.set_aspect('equal', 'box')
    if ax is None:
        plt.tight_layout()
        plt.show()

def load_csv_as_dataframe(filename):
    return pd.read_csv(filename, index_col=0)

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

import pandas as pd
import numpy as np

def save_event_table(asevents, events='splice_total', table_name=None, event_order=None, table=True):
    """
    Generates a table with event count and percentage data and saves it to a CSV file.

    Parameters:
    - asevents: DataFrame containing the events data.
    - events: Type of events to process. Options include:
        'splice_total', 'splice_I', 'splice_E', or 'trans'.
    - table_name: Name of the CSV file to save the table. If None, defaults based on event type.
    - event_order: List specifying the desired order of event types.
    - table: If True, the table is saved to a CSV file.
    
    Returns:
    - DataFrame with columns 'Event Type', 'Percentage', and 'Count'.
    """

    if events == 'splice_total':
        event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
                       'Alternative First Exon', 'Alternative Last Exon']
        events_counts = pd.Series({
            'ce': asevents['cet'].sum(),
            'a5': asevents['a5t'].sum(),
            'a3': asevents['a3t'].sum(),
            'i': asevents['it'].sum(),
            'ap': ((asevents['ap'] >= 0) & (asevents['ap'] <= 1)).sum(),
            'ae': ((asevents['ae'] >= 0) & (asevents['ae'] <= 1)).sum()
        })
        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100
        ASsum = dict(zip(event_names, percentages))
        AStotal = dict(zip(event_names, sizes))

    elif events == 'splice_I':
        event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
                       'Alternative First Exon', 'Alternative Last Exon']
        events_counts = pd.Series({
            'ce': asevents['ceI'].sum(),
            'a5': asevents['a5I'].sum(),
            'a3': asevents['a3I'].sum(),
            'i': asevents['iI'].sum(),
            'ap': ((asevents['ap'] > 0) & (asevents['ap'] <= 1)).sum(),
            'ae': ((asevents['ae'] > 0) & (asevents['ae'] <= 1)).sum()
        })
        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100
        ASsum = dict(zip(event_names, percentages))
        AStotal = dict(zip(event_names, sizes))

    elif events == 'splice_E':
        event_names = ['Alternative Exon', 'Alternative 5′ End', 'Alternative 3′ End', 'Intron Retention',
                       'Alternative First Exon', 'Alternative Last Exon']
        events_counts = pd.Series({
            'ce': asevents['ceE'].sum(),
            'a5': asevents['a5E'].sum(),
            'a3': asevents['a3E'].sum(),
            'i': asevents['iE'].sum(),
            'ap': (asevents['ap'] == 0).sum(),
            'ae': (asevents['ae'] == 0).sum()
        })
        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100
        ASsum = dict(zip(event_names, percentages))
        AStotal = dict(zip(event_names, sizes))

    elif events == 'trans':
        # For 'trans', reduce asevents to unique gene_id and calculate multi-transcript status
        gene_mult_trans = asevents[['gene_id', 'multi_trans']].drop_duplicates()
        gene_mult_trans = gene_mult_trans.groupby('gene_id')['multi_trans'].max()
        events_counts = pd.Series({
            "Single Transcript Gene": (gene_mult_trans == 0).sum(),
            "Multi Transcript Gene": (gene_mult_trans == 1).sum()
        })
        sizes = events_counts.values
        percentages = (sizes / sizes.sum()) * 100
        ASsum = dict(zip(events_counts.index, percentages))
        AStotal = dict(zip(events_counts.index, sizes))
    
    else:
        raise ValueError("Unsupported events type")

    # Remove events with zero percentage and the corresponding counts
    ASsum_nonzero = {k: v for k, v in ASsum.items() if v > 0}
    AStotal_nonzero = {k: AStotal[k] for k in ASsum_nonzero.keys()}

    # Reorder the events if event_order is provided
    if event_order:
        event_order = [event for event in event_order if event in ASsum_nonzero]
        ASsum_nonzero = {k: ASsum_nonzero[k] for k in event_order}
        AStotal_nonzero = {k: AStotal_nonzero[k] for k in event_order}
    else:
        ASsum_nonzero = dict(sorted(ASsum_nonzero.items()))
        AStotal_nonzero = dict(sorted(AStotal_nonzero.items()))

    # Combine the data into a DataFrame
    table_df = pd.DataFrame({
        'Event Type': list(ASsum_nonzero.keys()),
        'Percentage': list(ASsum_nonzero.values()),
        'Count': list(AStotal_nonzero.values())
    })

    table_df = table_df.set_index('Event Type')
    # Save the DataFrame to a CSV file if requested
    if table:
        if table_name is None:
            table_name = 'splice_events.csv' if 'splice' in events else 'trans_events.csv'
        table_df.to_csv(table_name, index=True)

    return table_df



def custom_horizontal_bar_chart_with_fisher_I_E(aff_I_df, unaff_I_df, aff_E_df, unaff_E_df, title=None, fscale=2, ax=None):
    """
    Create a custom horizontal bar chart comparing exon Inclusion (I) vs Exclusion (E)
    between affected and unaffected groups, and perform Fisher's exact test on each category.
    
    The function expects four DataFrames:
      - aff_I_df: Affected group inclusion counts
      - unaff_I_df: Unaffected group inclusion counts
      - aff_E_df: Affected group exclusion counts
      - unaff_E_df: Unaffected group exclusion counts
      
    Each DataFrame should have the same index (categories) and a single column of raw counts.
    
    In the plot:
      - For each category (row), two horizontal bars are shown:
         • The affected group is plotted at y - (bar_height/2) using blue patches:
             - Inclusion (I) is drawn as a bar to the left (negative direction, alpha=0.5)
             - Exclusion (E) is drawn as a bar to the right (positive direction, alpha=1)
         • The unaffected group is plotted at y + (bar_height/2) using orange patches:
             - Inclusion (I) is drawn as a bar to the left (negative, alpha=1)
             - Exclusion (E) is drawn as a bar to the right (positive, alpha=0.5)
      - The x-axis spans from -1 to 1 (representing proportions) with a vertical dashed line at 0.
      - Fisher’s exact test is performed for each category using the raw counts:
            contingency_table = [[aff_I, aff_E], [unaff_I, unaff_E]]
    
    Parameters:
      - aff_I_df (pd.DataFrame): Affected group inclusion counts.
      - unaff_I_df (pd.DataFrame): Unaffected group inclusion counts.
      - aff_E_df (pd.DataFrame): Affected group exclusion counts.
      - unaff_E_df (pd.DataFrame): Unaffected group exclusion counts.
      - fscale (float): Font scale for labels and title.
      - ax (matplotlib.axes.Axes): Axes object for plotting; if None, one is created.
    
    Returns:
      - fisher_results_df (pd.DataFrame): A DataFrame containing Fisher's exact test results
          for each category (including raw counts, odds ratios, and p-values).
    """
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    # Ensure all DataFrames have the same ordering of categories
    sort_index = aff_E_df.iloc[:, 0].sort_values(ascending=False).index
    aff_I_df = aff_I_df.loc[sort_index]
    aff_E_df = aff_E_df.loc[sort_index]
    unaff_I_df = unaff_I_df.loc[sort_index]
    unaff_E_df = unaff_E_df.loc[sort_index]
    
    categories = sort_index
    bar_height = 0.4
    y_positions = list(range(len(categories)))

    fisher_results = []
    
    def safe_float(x):
        """Convert a value to float; if conversion fails, return 0.0."""
        try:
            return float(x)
        except Exception:
            return 0.0
    
    # Loop over each category to plot bars and compute Fisher's exact test
    for i, category in enumerate(categories):

        aff_I = safe_float(aff_I_df.loc[category, "Count"])
        aff_E = safe_float(aff_E_df.loc[category, "Count"])
        unaff_I = safe_float(unaff_I_df.loc[category, "Count"])
        unaff_E = safe_float(unaff_E_df.loc[category, "Count"])

        aff_I_total = aff_I_df["Count"].sum()
        aff_E_total = aff_E_df["Count"].sum()
        unaff_I_total = unaff_I_df["Count"].sum()
        unaff_E_total = unaff_E_df["Count"].sum()
        
       
        aff_I_prop = aff_I / aff_I_total if aff_I_total != 0 else 0
        aff_E_prop = aff_E / aff_E_total if aff_E_total != 0 else 0
        unaff_I_prop = unaff_I / unaff_I_total if unaff_I_total != 0 else 0
        unaff_E_prop = unaff_E / unaff_E_total if unaff_E_total != 0 else 0
        
       
        # Determine y-positions for the two groups
        y_aff = i - bar_height / 2
        y_unaff = i + bar_height / 2
        
        # Plot affected group bars (blue)
        # Inclusion (I): plotted to the left (negative)
        ax.add_patch(patches.Rectangle(
            (0, y_aff), -aff_I_prop, bar_height,
            facecolor="blue", alpha=0.5, edgecolor="none",
            label="Affected Inclusion" if i == 0 else None))
        # Exclusion (E): plotted to the right (positive)
        ax.add_patch(patches.Rectangle(
            (0, y_aff), aff_E_prop, bar_height,
            facecolor="blue", alpha=1, edgecolor="none",
            label="Affected Exclusion" if i == 0 else None))
        
        # Plot unaffected group bars (orange)
        # Inclusion (I): plotted to the left (negative)
        ax.add_patch(patches.Rectangle(
            (0, y_unaff), -unaff_I_prop, bar_height,
            facecolor="orange", alpha=1, edgecolor="none",
            label="Unaffected Inclusion" if i == 0 else None))
        # Exclusion (E): plotted to the right (positive)
        ax.add_patch(patches.Rectangle(
            (0, y_unaff), unaff_E_prop, bar_height,
            facecolor="orange", alpha=0.5, edgecolor="none",
            label="Unaffected Exclusion" if i == 0 else None))
        
        # Perform Fisher's exact test separately for inclusion and exclusion per category,
        # and then a global test comparing inclusion vs exclusion across all categories.
        
        # Initialize lists to store results
        fisher_results_inclusion = []
        fisher_results_exclusion = []
        
        # Loop over each category and perform Fisher's exact test for inclusion and exclusion separately
        for category in aff_I_df.index:
            aff_I = aff_I_df.loc[category, "Count"]
            unaff_I = unaff_I_df.loc[category, "Count"]
            aff_E = aff_E_df.loc[category, "Count"]
            unaff_E = unaff_E_df.loc[category, "Count"]
            
            # Corrected Fisher test for inclusion (Affected vs Unaffected)
            contingency_inclusion = [[aff_I, aff_E], [unaff_I, unaff_E]]
            odds_ratio_inc, p_value_inc = fisher_exact(contingency_inclusion)
            fisher_results_inclusion.append({
                "Category": category,
                "Odds Ratio": odds_ratio_inc,
                "p-value": p_value_inc,
                "Test": "Inclusion (Affected vs Unaffected)"
            })
            
            # Corrected Fisher test for exclusion (Affected vs Unaffected)
            contingency_exclusion = [[aff_E, aff_I], [unaff_E, unaff_I]]
            odds_ratio_exc, p_value_exc = fisher_exact(contingency_exclusion)
            fisher_results_exclusion.append({
                "Category": category,
                "Odds Ratio": odds_ratio_exc,
                "p-value": p_value_exc,
                "Test": "Exclusion (Affected vs Unaffected)"
            })
        
        # Global test comparing inclusion vs exclusion across all categories
        total_aff_I = aff_I_df["Count"].sum()
        total_unaff_I = unaff_I_df["Count"].sum()
        total_aff_E = aff_E_df["Count"].sum()
        total_unaff_E = unaff_E_df["Count"].sum()
        
        contingency_global = [[total_aff_I, total_aff_E], [total_unaff_I, total_unaff_E]]
        odds_ratio_global, p_value_global = fisher_exact(contingency_global)
        
        global_test_result = pd.DataFrame([{
            "Category": "Global Inclusion vs Exclusion",
            "Odds Ratio": odds_ratio_global,
            "p-value": p_value_global,
            "Test": "Global Inclusion vs Exclusion"
        }])
        
        # Combine all results into a single DataFrame
        fisher_results_df = pd.DataFrame(fisher_results_inclusion + fisher_results_exclusion)
        fisher_results_df = pd.concat([fisher_results_df, global_test_result], ignore_index=True)


    font_scale_ax = fscale
    
    # Customize axes appearance
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_yticks([r for r in y_positions])
    ax.set_yticklabels(categories, fontsize=font_scale_ax)
    
    # Add title and custom text labels (using axis-relative coordinates)
    ax.text(0.5, 1.4,title,
            ha="center", va="top", fontsize=font_scale_ax * 1.5, transform=ax.transAxes)
    ax.text(0.25, 1.2, "Inclusion", ha="center", va="top",
            fontsize=font_scale_ax * 1.2, transform=ax.transAxes)
    ax.text(0.75, 1.2, "Exclusion", ha="center", va="top",
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
    
    # Add legend below the plot
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=2,
        fontsize=font_scale_ax * 1.1,
        frameon=False
    )
    
    #fisher_results_df = pd.DataFrame(fisher_results)
    return fisher_results_df
    
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

# Function to compute percentages within each group
def compute_groupwise_percentages(df, label):
    result = df.groupby(["Group", "Event Type"])["Count"].sum().reset_index()
    result["Total_in_Group"] = result.groupby("Group")["Count"].transform("sum")
    result["Percentage"] = 100 * result["Count"] / result["Total_in_Group"]
    result["Figure"] = label
    return result


