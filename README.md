# Bactopia Results Processor & Visualizer

This repository contains a comprehensive Python script that parses, processes, and visualizes results from a **Bactopia** analysis pipeline. It automatically generates summary files and plots based on key outputs like SNP distances, Prokka annotations, sample information, sourmash, and more.

## 📁 Directory Structure

Set your directory paths to the root of your Bactopia run and where you’d like to save output plots:

```python
directory = '/path/to/bactopia/output'
save_directory = '/path/to/save/figures'
```

## 🧩 Features

✅ Parsing and Processing

Prokka Annotation Files → Generates ```prokka_results.tsv```

SNP Density Matrices → Combines into ```combined_snp_density.tsv```

Sourmash Files → Merged into ```sourmash_results.tsv```

Sample Info (GregD_samples.txt) → Saved as ```sample_info.tsv```

Sample List Files → Saved as ```sample_list.tsv```

Bactopia Summary Files → Merged into ```bactopia_summary.tsv```

Additional TSVs (agrvate, spatyper, sccmec, amrfinder) → Saved and visualized

## 📊 Visualization

AGR Groups (agrvate.tsv) → Bar plot

SCCmec Elements (staphopiasccmec.tsv) → Bar plot of presence/absence

Spatyper Types → Frequency bar chart

SNP Distance Matrix → Interactive histogram + downloadable TSV

Prokka Results → Box plots of contigs, CDS, and base pairs

## 🌐 Interactive HTML

SNP Density plot includes a download button to export binned data as TSV directly from the figure.

## 🚀 Usage

Run the main script:

```
python bactopia_summary_visualizer.py
```

## 📂 Output Files

All processed tables and plots will be saved in the save_directory. Example outputs include:

prokka_results.tsv

combined_snp_density.tsv

Agrvate.png

SnpDensity_with_download.html

ProkkaResults.png
