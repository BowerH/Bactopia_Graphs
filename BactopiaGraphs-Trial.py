import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mpld3
from mpld3 import plugins
import json

# Directory is the root directory of the Bactopia run
# Save_directory is the directory where the graphs will be saved
directory = '/Users/hannah/Lab_Runs/Greg'
save_directory = '/Users/hannah/Lab_Runs/Greg'


#for the sourmash
# https://www.researchgate.net/figure/Krona-chart-of-the-bacteria-represented-by-16S-rRNA-gene-amplicon-based-bacterial_fig1_328672913


def process_prokka_files(directory):
    all_data = []
    try:
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file.endswith('.txt') and 'annotator/prokka' in root:
                    file_path = os.path.join(root, file)
                    data_dict = {}
                    try:
                        with open(file_path, 'r') as f:
                            for line in f:
                                if ':' in line:
                                    key, value = line.strip().split(':', 1)
                                    data_dict[key.strip()] = value.strip()
                        data_dict['organism'] = file.replace('.txt', '')
                        all_data.append(data_dict)
                    except IOError:
                        print(f"Error reading file: {file_path}")
        df = pd.DataFrame(all_data)
        df.to_csv('/Users/hannah/Lab_Runs/Greg/prokka_results.tsv', sep='\t', index=False)
        return df
    except Exception as e:
        print(f"Error in process_prokka_files: {str(e)}")
        return pd.DataFrame()


prokkadf = process_prokka_files(directory)


def process_snp_density_files(directory):
    try:
        snp_files = []
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file == 'snp_disty.txt' and 'fna_files' in root:
                    snp_files.append(os.path.join(root, file))

        if snp_files:
            first_df = pd.read_csv(snp_files[0], sep='\s+')
            column_names = first_df.columns.tolist()

            dfs = []
            for file in snp_files:
                try:
                    df = pd.read_csv(file, sep='\s+')
                    dfs.append(df)
                except pd.errors.EmptyDataError:
                    print(f"Empty file: {file}")
                except Exception as e:
                    print(f"Error reading file {file}: {str(e)}")

            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df.to_csv('/Users/hannah/Lab_Runs/Greg/combined_snp_density.tsv', sep='\t', index=False)
            return combined_df
        return None
    except Exception as e:
        print(f"Error in process_snp_density_files: {str(e)}")
        return None


snpdensity_df = process_snp_density_files(directory)


def process_sourmash_files(directory):
    try:
        dfs = []
        for root, _, files in os.walk(directory):
            for file in files:
                if 'sourmash' in file:
                    try:
                        df = pd.read_csv(os.path.join(root, file))
                        df.columns = [col.split('__')[-1] if '__' in col else col for col in df.columns]
                        df = df.apply(lambda x: x.str.split('__').str[-1] if x.dtype == 'object' else x)
                        dfs.append(df)
                    except Exception as e:
                        print(f"Error processing file {file}: {str(e)}")

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df.to_csv('/Users/hannah/Lab_Runs/Greg/sourmash_results.tsv', sep='\t', index=False)
            return combined_df
        return None
    except Exception as e:
        print(f"Error in process_sourmash_files: {str(e)}")
        return None


sourmash_df = process_sourmash_files(directory)


def process_samples_files(directory):
    try:
        dfs = []
        for root, _, files in os.walk(directory):
            for file in files:
                if file == 'GregD_samples.txt':
                    try:
                        df = pd.read_csv(os.path.join(root, file), delim_whitespace=True)
                        dfs.append(df)
                    except Exception as e:
                        print(f"Error reading file {file}: {str(e)}")

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df.to_csv('/Users/hannah/Lab_Runs/Greg/sample_info.tsv', sep='\t', index=False)
            return combined_df
        return None
    except Exception as e:
        print(f"Error in process_samples_files: {str(e)}")
        return None


sampleinfo_df = process_samples_files(directory)


def process_bactopia_summary(directory):
    try:
        summary_data = []
        for root, _, files in os.walk(directory):
            for file in files:
                if file == 'bactopia-summary.txt':
                    file_path = os.path.join(root, file)
                    try:
                        data = {
                            'total_samples': 0, 'passed': 0, 'gold': 0, 'silver': 0,
                            'bronze': 0, 'excluded': 0, 'failed_cutoff': 0, 'qc_failure': 0
                        }
                        with open(file_path, 'r') as f:
                            for line in f:
                                if 'Total Samples:' in line:
                                    data['total_samples'] = int(line.split(':')[-1].strip())
                                elif 'Passed:' in line and 'Gold:' not in line:
                                    data['passed'] = int(line.split(':')[-1].strip())
                                elif 'Gold:' in line:
                                    data['gold'] = int(line.split(':')[-1].strip())
                                elif 'Silver:' in line:
                                    data['silver'] = int(line.split(':')[-1].strip())
                                elif 'Bronze:' in line:
                                    data['bronze'] = int(line.split(':')[-1].strip())
                                elif 'Excluded:' in line:
                                    data['excluded'] = int(line.split(':')[-1].strip())
                                elif 'Failed Cutoff:' in line:
                                    data['failed_cutoff'] = int(line.split(':')[-1].strip())
                                elif 'QC Failure:' in line:
                                    data['qc_failure'] = int(line.split(':')[-1].strip())
                        summary_data.append(data)
                    except Exception as e:
                        print(f"Error processing file {file_path}: {str(e)}")

        if summary_data:
            df = pd.DataFrame(summary_data)
            df.to_csv('/Users/hannah/Lab_Runs/Greg/bactopia_summary.tsv', sep='\t', index=False)
            return df
        return None
    except Exception as e:
        print(f"Error in process_bactopia_summary: {str(e)}")
        return None


bactopiasummary_df = process_bactopia_summary(directory)


def process_sample_list(directory):
    try:
        dfs = []
        for root, _, files in os.walk(directory):
            for file in files:
                if 'Sample_List.txt' in file:
                    try:
                        df = pd.read_csv(os.path.join(root, file), header=None, names=['sample'])
                        dfs.append(df)
                    except Exception as e:
                        print(f"Error reading file {file}: {str(e)}")

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            combined_df.to_csv('/Users/hannah/Lab_Runs/Greg/sample_list.tsv', sep='\t', index=False)
            return combined_df
        return None
    except Exception as e:
        print(f"Error in process_sample_list: {str(e)}")
        return None


samplelist_df = process_sample_list(directory)


def find_and_process_tsv_files(directory):
    try:
        dfs = {}
        target_files = {
            'agrvate.tsv', 'spatyper.tsv', 'staphopiasccmec.tsv', 'amrfinderplus-proteins.tsv'
        }
        for root, _, files in os.walk(directory):
            for file in files:
                if file in target_files:
                    try:
                        file_path = os.path.join(root, file)
                        df = pd.read_csv(file_path, sep='\t')
                        dfs[file.replace('.tsv', '')] = df
                        output_name = f"{file.replace('.tsv', '')}_results.csv"
                        df.to_csv(output_name, index=False)
                    except Exception as e:
                        print(f"Error processing file {file}: {str(e)}")
        return dfs
    except Exception as e:
        print(f"Error in find_and_process_tsv_files: {str(e)}")
        return {}


b = find_and_process_tsv_files(directory)


def find_tsv_files(directory, target_files):
    try:
        found_files = {}
        for root, dirs, files in os.walk(directory):
            for target_file in target_files:
                if target_file in files:
                    found_files[target_file] = os.path.join(root, target_file)
        return found_files
    except Exception as e:
        print(f"Error in find_tsv_files: {str(e)}")
        return {}


tsv_filenames = ["agrvate.tsv", "staphopiasccmec.tsv", "spatyper.tsv"]
found_file_paths = find_tsv_files(directory, tsv_filenames)

bactopiasummaries = {}
for filename, file_path in found_file_paths.items():
    try:
        df_name = filename.replace('.tsv', '')
        bactopiasummaries[df_name] = pd.read_csv(file_path, sep='\t')
        print(f"File loaded into DataFrame: {df_name}")
    except Exception as e:
        print(f"Error loading file {filename}: {str(e)}")

not_found = set(tsv_filenames) - set(found_file_paths.keys())
if not_found:
    print("Files not found:", not_found)


def plot_agrvate(df, title):
    try:
        df['group'] = df['agr_group'].str.replace('gp', '')
        group_counts = df['group'].value_counts().sort_index()
        plt.figure(figsize=(8, 5))
        group_counts.plot(kind='bar', color='skyblue')
        plt.title(f'Number of Samples in Each Group - {title}')
        plt.xlabel('AGR Group')
        plt.ylabel('Number of Samples')
        plt.xticks(rotation=0)
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        filename = 'Agrvate.png'
        full_path = os.path.join(save_directory, filename)
        plt.savefig(full_path)

        print(f"Figure saved to: {full_path}")
    except Exception as e:
        print(f"Error in plot_agrvate: {str(e)}")


plot_agrvate(bactopiasummaries['agrvate'], 'Agrvate')


def plot_staphopiasccmec(df):
    try:
        true_counts = df.drop(columns='sample').apply(lambda x: x == True).sum()
        plt.figure(figsize=(15, 6))
        plt.grid(axis='y', linestyle='--', alpha=0.7, color='gray', linewidth=0.5)
        plt.bar(true_counts.index, true_counts.values, color='skyblue')
        plt.title("Staphopiasccmec True Counts per Column")
        plt.xlabel('Column')
        plt.ylabel("Count of 'True'")
        plt.xticks(rotation=0)
        plt.tight_layout()

        filename = 'Staphopiasccmec.png'
        full_path = os.path.join(save_directory, filename)
        plt.savefig(full_path)

        print(f"Figure saved to: {full_path}")
    except Exception as e:
        print(f"Error in plot_staphopiasccmec: {str(e)}")


plot_staphopiasccmec(bactopiasummaries['staphopiasccmec'])


def plot_spatyper(df):
    try:
        type_counts = df['Type'].value_counts()
        plt.figure(figsize=(12, 6))
        plt.bar(type_counts.index, type_counts.values, color='skyblue')
        plt.title('Frequency of Each Type- Spatyper')
        plt.xlabel('Type')
        plt.ylabel('Count')
        plt.xticks(rotation=90)
        plt.tight_layout()

        filename = 'Spatyper.png'
        full_path = os.path.join(save_directory, filename)
        plt.savefig(full_path)
        print(f"Figure saved to: {full_path}")

        return type_counts
    except Exception as e:
        print(f"Error in plot_spatyper: {str(e)}")
        return pd.Series()

plot_spatyper(bactopiasummaries['spatyper'])

import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import os

import base64

import json

def plot_snpdensity(df):
    try:
        def extract_below_diagonal(df):
            rows, cols = df.shape
            below_diagonal = []
            for i in range(rows):
                for j in range(cols):
                    if i > j:
                        value = df.iloc[i, j]
                        if value != 0:
                            below_diagonal.append({
                                'x_position': j,
                                'value': value
                            })
            return pd.DataFrame(below_diagonal)

        below_diag_df = extract_below_diagonal(df)

        # Create a Plotly figure
        fig = go.Figure()

        # Add histogram trace
        fig.add_trace(go.Histogram(
            x=below_diag_df['value'],
            xbins=dict(start=0, end=below_diag_df['value'].max(), size=1000),
            marker_color='skyblue',
            opacity=0.7,
            marker_line=dict(color='black', width=1)
        ))

        # Update layout
        fig.update_layout(
            title='Distribution of SNP Distances',
            xaxis_title='SNP Distance',
            yaxis_title='Number of Occurrences',
            bargap=0.1,
            template='plotly_white'
        )

        # Add grid
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')

        # Create TSV content
        tsv_content = below_diag_df.to_csv(index=False, sep='\t')
        tsv_b64 = base64.b64encode(tsv_content.encode()).decode()

        # Create download button HTML
        download_button = f'''
        <a href="data:text/tab-separated-values;base64,{tsv_b64}" download="snp_density_data.tsv">
            <button style="font-size: 16px; padding: 10px; margin-bottom: 20px;">
                Download Data as TSV
            </button>
        </a>
        '''

        # Save as interactive HTML with download functionality
        filename = 'SnpDensity_with_download.html'
        full_path = os.path.join(save_directory, filename)

        # Write HTML file
        with open(full_path, 'w') as f:
            f.write(download_button)  # Write the button first
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

        print(f"Interactive HTML with download functionality saved to: {full_path}")

    except Exception as e:
        print(f"An error occurred: {str(e)}")

'''
def plot_snpdensity(df):
    try:
        def extract_below_diagonal(df):
            rows, cols = df.shape
            below_diagonal = []
            for i in range(rows):
                for j in range(cols):
                    if i > j:
                        value = df.iloc[i, j]
                        if value != 0:
                            below_diagonal.append({
                                'x_position': j,
                                'value': value
                            })
            return pd.DataFrame(below_diagonal)

        below_diag_df = extract_below_diagonal(df)

        fig, ax = plt.subplots(figsize=(12, 6))
        n, bins, patches = ax.hist(below_diag_df['value'],
                                   bins=range(0, int(below_diag_df['value'].max()) + 1000, 1000),
                                   color='skyblue', alpha=0.7, edgecolor='black', density=False)

        ax.set_xlabel('SNP Distance')
        ax.set_ylabel('Number of Occurrences')
        ax.set_title('Distribution of SNP Distances')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()

        # Prepare data for download
        download_data = [{'bin_start': bins[i], 'bin_end': bins[i+1], 'count': n[i]} for i in range(len(n))]
        json_data = json.dumps(download_data)

        # Create tooltip for bars
        labels = [f"Distance: {bins[i]}-{bins[i + 1]}, Count: {n[i]}" for i in range(len(n))]
        for patch, label in zip(patches, labels):
            tooltip = mpld3.plugins.LineLabelTooltip(patch, label)
            mpld3.plugins.connect(fig, tooltip)



        # Save static image
        filename_png = 'SnpDensity.png'
        full_path_png = os.path.join(save_directory, filename_png)
        plt.savefig(full_path_png)
        print(f"Static figure saved to: {full_path_png}")

        # Save interactive HTML
        filename_html = 'SnpDensity_interactive.html'
        full_path_html = os.path.join(save_directory, filename_html)
        mpld3.save_html(fig, full_path_html)
        print(f"Interactive figure saved to: {full_path_html}")

        # Display the interactive plot
        mpld3.show()
        print(f"Figure saved to: {full_path_html}")
    except Exception as e:
        print(f"Error in plot_snpdensity: {str(e)}")
'''
# Call the function
plot_snpdensity(snpdensity_df)

def plot_prokka_results(df):
    try:
        df['Contigs'] = pd.to_numeric(df['contigs'])
        df['MBasePairs'] = pd.to_numeric(df['bases'])
        df['CDS'] = pd.to_numeric(df['CDS'])

        if len(df) > 300:
            fig, axes = plt.subplots(3, 1, figsize=(6, 18))  # Stacked vertically
            vert_orientation = False  # Horizontal box plots for large datasets
        else:
            fig, axes = plt.subplots(1, 3, figsize=(12, 6))  # Side by side
            vert_orientation = True  # Vertical box plots for smaller datasets

        fig.suptitle('Prokka Results Distribution', fontsize=14, y=0.95)
        data_columns = ['Contigs', 'MBasePairs', 'CDS']

        for ax, column in zip(axes, data_columns):
            # Create boxplot with dynamic orientation
            bp = ax.boxplot(df[column], vert=vert_orientation, patch_artist=True,
                            boxprops=dict(facecolor='skyblue', alpha=0.6))

            # Add scatter points for smaller datasets
            if len(df) < 200 and vert_orientation:
                ax.scatter(np.ones(len(df[column])), df[column], color='darkblue', alpha=0.5, zorder=2)

            if vert_orientation:  # Labels for vertical orientation
                ax.set_title(column)
                ax.set_xticks([])
                ax.set_xlim(0.5, 1.5)
                if column == 'MBasePairs':
                    ax.set_ylabel('Mbasepairs')
                else:
                    ax.set_ylabel(column)
            else:  # Labels for horizontal orientation
                ax.set_title(column)
                ax.set_yticks([])
                ax.set_ylim(0.5, 1.5)
                if column == 'MBasePairs':
                    ax.set_xlabel('Mbasepairs')
                else:
                    ax.set_xlabel(column)

        plt.tight_layout()
        plt.subplots_adjust(top=0.85)

        filename = 'ProkkaResults.png'
        full_path = os.path.join(save_directory, filename)
        plt.savefig(full_path)

        print(f"Figure saved to: {full_path}")
    except Exception as e:
        print(f"Error in plot_prokka_results: {str(e)}")

plot_prokka_results(prokkadf)
