import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import numpy as np
import sys


SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 16
LARGE_SIZE = 20


plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=LARGE_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title

# plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', family='sans-serif')  # Specify a different font family




# def BinCollector_plot(input_file, plot_title):
#     df = pd.read_csv(input_file, sep='\t', header=0)
    
#     new_df = pd.DataFrame()
#     new_df = pd.concat([df[col].apply(lambda x: int(col) if x != 0 else 0) for col in df.columns], axis=1)

    
#     df1 = df.applymap(lambda x: 1 if x != 0 else x)
#     df1['RowSum'] = df1.sum(axis=1)
    
    
#     df['RowSum'] = new_df.sum(axis=1)
    
#     multi = [r * one for r, one in zip(df['RowSum'].tolist(), df1['RowSum'].tolist())]
    
#     df['RowSum'] = multi

    
#     # Reorder the DataFrame by the sum of row values
#     df = df.sort_values(['RowSum'], ascending=True)
#     # Remove the 'RowSum' column if desired
#     df = df.drop('RowSum', axis=1)

#     bars = [df[column].sum() for column in df]

#     # Create the figure and axes
#     fig, (ax_bar, ax_heatmap) = plt.subplots(2, 1, figsize=(8, 12), gridspec_kw={'height_ratios': [1, 6]}, sharex = True)

#     # Create the bar chart
#     x = range(len(bars))
#     bar_width = 0.9
#     ax_bar.bar([i + bar_width/2 for i in x], bars, width=bar_width)

#     # Remove the top and right spines from the bar chart
#     ax_bar.spines['top'].set_visible(False)
#     ax_bar.spines['right'].set_visible(False)

#     # Set the desired frequency of tick labels
#     tick_frequency = 10

#     # Get the number of columns in the DataFrame
#     num_columns = len(df.columns)

#     # Compute the step value to skip ticks
#     step = max(1, num_columns // tick_frequency)

#     # Set the x-axis tick positions
#     x_ticks = np.arange(0, num_columns, step)

#     # Set the x-axis tick labels based on the column names at the chosen positions
#     x_tick_labels = df.columns[x_ticks]

#     # Set the x-axis ticks and labels for the bar chart
#     ax_bar.set_xticks(x_ticks + bar_width / 2)
#     ax_bar.set_xticklabels(x_tick_labels, rotation=90)

#     # Remove the first tick (zero value) on the y-axis of the bar chart
#     ax_bar.set_yticks(ax_bar.get_yticks()[1:])
#     ax_bar.set_yticklabels(ax_bar.get_yticks())

#     # Create the heatmap
#     if len(df) < 100:
#         sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=True)
#     else:
#         sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=False)

#     # Add colorbar next to the heatmap
#     cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.1])  # Define the position and size of the colorbar
#     fig.colorbar(ax_heatmap.get_children()[0], cax=cbar_ax)  # Add the colorbar to the specified position

#     # Adjust the spacing between subplots
#     plt.subplots_adjust(hspace=0)

#     # Set the overall title
#     plt.suptitle(plot_title)

#     # Show the plot
#     plt.show()





# def BinCollector_plot(input_file, plot_title):
#     df = pd.read_csv(input_file, sep='\t', header=0)

#     new_df = pd.DataFrame()
#     new_df = pd.concat([df[col].apply(lambda x: int(col) if x != 0 else 0) for col in df.columns], axis=1)

#     df1 = df.applymap(lambda x: 1 if x != 0 else x)
#     df1['RowSum'] = df1.sum(axis=1)

#     df['RowSum'] = new_df.sum(axis=1)

#     multi = [r * one for r, one in zip(df['RowSum'].tolist(), df1['RowSum'].tolist())]

#     df['RowSum'] = multi

#     # Reorder the DataFrame by the sum of row values
#     df = df.sort_values(['RowSum'], ascending=True)
#     # Remove the 'RowSum' column if desired
#     df = df.drop('RowSum', axis=1)

#     bars = [df[column].sum() for column in df]

#     # Create the figure and axes
#     fig, (ax_bar, ax_heatmap) = plt.subplots(2, 1, figsize=(8, 12), gridspec_kw={'height_ratios': [1, 6]}, sharex=True)

#     # Set the desired frequency of tick labels
#     tick_frequency = 10

#     # Get the number of columns in the DataFrame
#     num_columns = len(df.columns)

#     # Compute the step value to skip ticks
#     step = max(1, num_columns // tick_frequency)

#     # Set the x-axis tick positions
#     x_ticks = np.arange(0, num_columns, step)

#     # Set the x-axis tick labels based on the column names at the chosen positions
#     x_tick_labels = [df.columns[i] for i in x_ticks]

#     # Create the bar chart
#     x = range(len(bars))
#     bar_width = 0.9
#     ax_bar.bar([i + bar_width/2 for i in x], bars, width=bar_width)

#     # Remove the top and right spines from the bar chart
#     ax_bar.spines['top'].set_visible(False)
#     ax_bar.spines['right'].set_visible(False)

#     # Set the x-axis ticks and labels for the bar chart
#     ax_bar.set_xticks([i + bar_width/2 for i in x_ticks])
#     ax_bar.set_xticklabels(x_tick_labels, rotation=90)

#     # Remove the first tick (zero value) on the y-axis of the bar chart
#     ax_bar.set_yticks(ax_bar.get_yticks()[1:])
#     ax_bar.set_yticklabels(ax_bar.get_yticks())

#     # Create the heatmap
#     if len(df) < 100:
#         sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=True)
#     else:
#         sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=False)

#     # Add colorbar next to the heatmap
#     cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.1])  # Define the position and size of the colorbar
#     fig.colorbar(ax_heatmap.get_children()[0], cax=cbar_ax)  # Add the colorbar to the specified position

#     # Adjust the spacing between subplots
#     plt.subplots_adjust(hspace=0)

#     # Set the overall title
#     plt.suptitle(plot_title)

#     # Show the plot
#     plt.show()






# def BinCollector_plot(input_file, plot_title, figure_folder, figure_file_name):
#     df = pd.read_csv(input_file, sep='\t', header=0)
    

#     # top200 = pd.read_csv(TPMPW_file, sep = '\t', header = 0)['Gene'].tolist()[:200]

#     # df = df[df.index.isin(top200)]


#     new_df = pd.DataFrame()
#     new_df = pd.concat([df[col].apply(lambda x: int(col) if x != 0 else 0) for col in df.columns], axis=1)

#     df1 = df.applymap(lambda x: 1 if x != 0 else x)
#     df1['RowSum'] = df1.sum(axis=1)

#     df['RowSum'] = new_df.sum(axis=1)

#     multi = [r * one for r, one in zip(df['RowSum'].tolist(), df1['RowSum'].tolist())]

#     df['RowSum'] = multi

#     # Reorder the DataFrame by the sum of row values
#     df = df.sort_values(['RowSum'], ascending=True)
#     # Remove the 'RowSum' column if desired
#     df = df.drop('RowSum', axis=1)

#     bars = [df[column].sum() for column in df]

#     # Create the figure and axes
#     fig, (ax_bar, ax_heatmap) = plt.subplots(2, 1, figsize=(8, 12), gridspec_kw={'height_ratios': [1, 6]}, sharex=True)

#     # Create the bar chart
#     x = range(len(bars))
#     bar_width = 0.9
#     ax_bar.bar([i + bar_width/2 for i in x], bars, width=bar_width)

#     # Remove the top and right spines from the bar chart
#     ax_bar.spines['top'].set_visible(False)
#     ax_bar.spines['right'].set_visible(False)

#     # Set the x-axis ticks and labels for the bar chart
#     x_ticks = range(len(df.columns))
#     ax_bar.set_xticks(x_ticks)
#     ax_bar.set_xticklabels(df.columns, rotation=90)

#     # Adjust the visibility of x-axis tick labels
#     tick_frequency = 10  # Show every 10th tick label
#     for i, label in enumerate(ax_bar.xaxis.get_ticklabels()):
#         if i % tick_frequency != 0:
#             label.set_visible(False)

#     # Remove the first tick (zero value) on the y-axis of the bar chart
#     ax_bar.set_yticks(ax_bar.get_yticks()[1:])
#     ax_bar.set_yticklabels(ax_bar.get_yticks())

#     # Create the heatmap
#     if len(df) < 100:
#         sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=True)
#     else:
#         sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=False)
#     ax_heatmap.tick_params(axis='y', labelsize=10)

#     # Add colorbar next to the heatmap
#     cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.1])  # Define the position and size of the colorbar
#     fig.colorbar(ax_heatmap.get_children()[0], cax=cbar_ax)  # Add the colorbar to the specified position

#     # Adjust the spacing between subplots
#     plt.subplots_adjust(hspace=0)

#     # Set the overall title
#     plt.suptitle(plot_title)

#     # Show the plot
#     plt.savefig(os.path.join(figure_folder, figure_file_name + '.pdf'), bbox_inches='tight')
#     plt.show()
#     plt.close()


def BinCollector_plot(input_file, plot_title, figure_folder, figure_file_name, genes_file):
    df = pd.read_csv(input_file, sep='\t', header=0)

    genes_df = pd.read_csv(genes_file, sep = '\t', names = ['gene'])
    genes_df = genes_df[genes_df.index <= 100]

    df['root'] = [i.split('.')[0] for i in df.index.values.tolist()]

    df = df[df['root'].isin(genes_df['gene'].tolist())]
    
    df = df.drop('root', axis = 1)

    row_values = [i for i in range(1, int((len(df.columns) + 2) / 2))]
    row_values = row_values + row_values[::-1]

    new_df = pd.DataFrame()
    for i in range(len(row_values)):
        new_df[i] = [-n * row_values[i] for n in df[f'{i + 1}'].tolist()]

    # Ensure that the indices of df and new_df match
    new_df.index = df.index

    df['RowSum2'] = new_df.sum(axis=1)

    df1 = df.applymap(lambda x: 1 if x != 0 else x)
    df['RowSum1'] = df1.sum(axis=1)


    # Reorder the DataFrame by the sum of row values
    df = df.sort_values(['RowSum2'], ascending=True)
    # Remove the 'RowSum' column if desired
    df = df.drop('RowSum1', axis=1)
    df = df.drop('RowSum2', axis=1)
    # if 'liaG.1' in df.index.values.tolist():
    #     print('lia yes')
    # else:
    #     print('lia no')
    
    # if 'pghB.2' in df.index.values.tolist():
    #     print('pghB yes')
    # else:
    #     print('pghB no')

    df = df[0:100]

    bars = [df[column].sum() for column in df]

    # Create the figure and axes
    fig, (ax_bar, ax_heatmap) = plt.subplots(2, 1, figsize=(8, 12), gridspec_kw={'height_ratios': [1, 6]}, sharex=True)

    # Create the bar chart
    x = range(len(bars))
    bar_width = 0.9
    ax_bar.bar([i + bar_width/2 for i in x], bars, width=bar_width)

    # Remove the top and right spines from the bar chart
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)

    # Set the x-axis ticks and labels for the bar chart
    x_ticks = range(len(df.columns))
    ax_bar.set_xticks(x_ticks)
    ax_bar.set_xticklabels(df.columns, rotation=90)

    # Adjust the visibility of x-axis tick labels
    tick_frequency = 10  # Show every 10th tick label
    for i, label in enumerate(ax_bar.xaxis.get_ticklabels()):
        if i % tick_frequency != 0:
            label.set_visible(False)

    # Remove the first tick (zero value) on the y-axis of the bar chart
    ax_bar.set_yticks(ax_bar.get_yticks()[1:])
    ax_bar.set_yticklabels(ax_bar.get_yticks())

    # # Create the heatmap
    # if len(df) < 100:
    #     sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=True)
    # else:
    #     sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=False)
    sns.heatmap(df, cmap='YlGnBu', cbar=False, ax=ax_heatmap, yticklabels=True)

    ax_heatmap.tick_params(axis='y', labelsize=10)

    # Add colorbar next to the heatmap
    cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.1])  # Define the position and size of the colorbar
    fig.colorbar(ax_heatmap.get_children()[0], cax=cbar_ax)  # Add the colorbar to the specified position

    # Adjust the spacing between subplots
    plt.subplots_adjust(hspace=0)

    # Set the overall title
    plt.suptitle(plot_title)

    # Show the plot
    plt.savefig(os.path.join(figure_folder, figure_file_name + '.pdf'), bbox_inches='tight')
    plt.show()
    plt.close()