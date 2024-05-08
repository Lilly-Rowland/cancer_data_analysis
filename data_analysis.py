import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set the color scheme
colors = ['#FAC6D2', '#282F44', '#645070', '#9E6471', '#98A8D7']

# Total protein count
total_protein_counts = {}
total_chrom_counts = {**{str(key): 0 for key in range(1, 23)}, 'X': 0, 'Y': 0}
base_changes = {}
total_cases = 0
total_subs = 0

def graph_chroms(data_dict, cancer_type, total):
    """
    Generate a bar graph showing the percentages of chromosomal mutations for a given cancer type.

    Parameters:
        data_dict (dict): Dictionary containing chromosome counts.
        cancer_type (str): Name of the cancer type.
        total (int): Total number of cases for the given cancer type.
    """
    plt.figure(facecolor=colors[1])
    percentages = {chrom: (count / total) * 100 for chrom, count in data_dict.items()}
    keys = list(percentages.keys())
    values = list(percentages.values())
    plt.bar(keys, values, color=colors[2])
    std_error = np.sqrt(np.array(values) * (100 - np.array(values)) / total)
    plt.errorbar(keys, values, yerr=std_error, fmt='none', capsize=2, capthick=1, ecolor=colors[1])
    plt.xlabel('Chromosome', color=colors[0])
    plt.ylabel('Percentages (%)', color=colors[0])
    plt.title(("Chromosomal Mutation Percentages in {} Cancer".format(cancer_type)).upper(), color=colors[0])
    plt.xticks(rotation=45, color=colors[0])
    plt.yticks(color=colors[0])
    plt.savefig("plots/{}_chrom_plot.png".format(cancer_type))
    plt.clf()

def graph_substitution_types(data_dict, total):
    """
    Generate a bar graph showing the percentages of substitution types.

    Parameters:
        data_dict (dict): Dictionary containing substitution counts.
        total (int): Total number of cases.
    """
    plt.figure(facecolor=colors[1])
    percentages = {sub: (count / total) * 100 for sub, count in data_dict.items()}
    keys = list(percentages.keys())
    values = list(percentages.values())
    plt.bar(keys, values, color=colors[2])
    std_error = np.sqrt(np.array(values) * (100 - np.array(values)) / total)
    plt.errorbar(keys, values, yerr=std_error, fmt='none', capsize=2, capthick=1, ecolor=colors[1])
    plt.xlabel('Substitution', color=colors[0])
    plt.ylabel('Percentages (%)', color=colors[0])
    plt.xticks(rotation=45, color=colors[0])
    plt.yticks(color=colors[0])
    plt.savefig("plots/substitution_plot.png")
    plt.clf()

def plot_transition_transversion_pie(transition_vs_transversions, cancer_type):
    """
    Generate a pie chart showing the frequency of transition vs. transversion mutations.

    Parameters:
        transition_vs_transversions (dict): Dictionary containing counts of transitions and transversions.
        cancer_type (str): Name of the cancer type.
    """
    plt.figure(facecolor=colors[1])
    legend_labels = ['Transitions', 'Transversions']
    plt.pie(transition_vs_transversions.values(), labels=legend_labels, autopct='%1.1f%%', startangle=140, colors=[colors[2], colors[3]], textprops={'color': colors[0]})
    plt.title('Transition vs Transversion Frequency in {} cancer'.format(cancer_type).upper(), color=colors[0])
    plt.axis('equal')
    plt.savefig("plots/{}_transitions_transversions_plot.png".format(cancer_type))
    plt.clf()

def graph_protein_mutations(protein_counts, num_cases, show_percents=False, cancer_type="aggregate", num_to_show=20):
    """
    Generate a bar graph showing the percentages of protein mutations for a given cancer type.

    Parameters:
        protein_counts (dict): Dictionary containing protein mutation counts.
        num_cases (int): Total number of cases for the given cancer type.
        show_percents (bool): Flag to display percentages on the graph (default is False).
        cancer_type (str): Name of the cancer type (default is "aggregate").
        num_to_show (int): Number of top proteins to display on the graph (default is 20).
    """
    plt.figure(facecolor=colors[1])
    top_proteins = sorted(protein_counts.items(), key=lambda x: x[1], reverse=True)[:num_to_show]
    keys = [item[0] for item in top_proteins]
    counts = [item[1] for item in top_proteins]
    percentages = [(count / num_cases) * 100 for count in counts]
    bars = plt.bar(keys, percentages, color=colors[2])
    plt.title(("Protein Mutation Percentages in {} Cancer".format(cancer_type)).upper(), color=colors[0])
    plt.xlabel("Protein", color=colors[0])
    plt.xticks(rotation=45, ha='right', color=colors[0])
    plt.yticks(color=colors[0])
    if show_percents:
        for bar, percentage in zip(bars, percentages):
            plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{percentage:.1f}%", ha='center', va='bottom', fontsize=8, color=colors[1])
    plt.tight_layout()
    plt.savefig("plots/{}_protein_mutation_plot.png".format(cancer_type))
    plt.clf()

def plot_mutation_types(mutations, num_cases, cancer_type):
    """
    Generate a bar graph showing the percentages of mutation types for a given cancer type.

    Parameters:
        mutations (dict): Dictionary containing mutation type counts.
        num_cases (int): Total number of cases for the given cancer type.
        cancer_type (str): Name of the cancer type.
    """
    plt.figure(facecolor=colors[1])
    top_mutations = sorted(mutations.items(), key=lambda x: x[1], reverse=True)[:10]
    keys = [item[0] for item in top_mutations]
    counts = [item[1] for item in top_mutations]
    percentages = [(count / num_cases) * 100 for count in counts]
    bars = plt.bar(keys, percentages, color=colors[2])
    plt.title(("Mutation Types in {} Cancer".format(cancer_type)).upper(), color=colors[0])
    plt.xlabel("Mutation Type", color=colors[0])
    plt.xticks(rotation=45, ha='right', color=colors[0])
    plt.yticks(color=colors[0])
    plt.tight_layout()
    plt.savefig("plots/{}_mutation_type_plot.png".format(cancer_type))
    plt.clf()

def loop_files(folder_path):
    """
    Loop through all files in a folder and process each file.

    Parameters:
        folder_path (str): Path to the folder containing the files.
    """
    for file_name in os.listdir(folder_path):
        if not file_name.endswith("tsv"):
            continue
        cancer_type = file_name[:-4]
        protein_counts = {}
        mutations = {}
        chrom_counts = total_chrom_counts = {**{str(key): 0 for key in range(1, 23)}, 'X': 0, 'Y': 0}
        transition_vs_transversions = {"num_transitions": 0, "num_transversions": 0}
        file_path = os.path.join(folder_path, file_name)
        df = pd.read_csv(file_path, sep='\t', encoding='utf-8')
        total_cohort_cases = 0
        for index, row in df.iterrows():
            if pd.isna(row['protein_change']) or pd.isna(row['dna_change']):
                continue
            protein = row['protein_change'].split()[0]
            dna_change = row['dna_change']
            chrom = dna_change.split(":")[0][3:]
            mutation = row['consequence']
            affected_cases = row['num_cohort_ssm_affected_cases']
            if protein_counts.get(protein):
                protein_counts[protein] += affected_cases
            else:
                protein_counts[protein] = affected_cases
            if total_protein_counts.get(protein):
                protein_counts[protein] += affected_cases
            else:
                total_protein_counts[protein] = affected_cases
            if mutations.get(mutation):
                mutations[mutation] += affected_cases
            else:
                mutations[mutation] = affected_cases
            chrom_counts[chrom] += affected_cases
            total_chrom_counts[chrom] += affected_cases
            if dna_change[-2] == '>':
                base_change = dna_change[-3:]
                global total_subs
                total_subs += 1
                if base_changes.get(base_change):
                    base_changes[base_change] += affected_cases
                else:
                    base_changes[base_change] = affected_cases
                old_base = base_change[0]
                new_base = base_change[2]
                if (old_base in ('A', 'G') and new_base in ('T', 'C')
                    or old_base in ('T', 'C') and new_base in ('A', 'G')):
                    transition_vs_transversions["num_transversions"] += affected_cases
                else:
                    transition_vs_transversions["num_transitions"] += affected_cases
            total_cohort_cases += affected_cases
            global total_cases
            total_cases += affected_cases
        show_percents = False
        if cancer_type == "aggregate":
            show_percents = True
        graph_protein_mutations(protein_counts, total_cohort_cases, cancer_type=cancer_type, show_percents=show_percents)
        plt.close()
        graph_chroms(chrom_counts, cancer_type, total_cohort_cases)
        plt.close()
        plot_transition_transversion_pie(transition_vs_transversions, cancer_type)
        plt.close()
        plot_mutation_types(mutations, total_cohort_cases, cancer_type)
        plt.close()
    
def main():
    """
    Main function to execute the program.
    """
    folder_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"data")
    loop_files(folder_path)

if __name__ == "__main__":
    main()
