File Structure
gui.py: Runs python GUI to present data
data_analysis.py: Runs data analysis on tsv files for different types of cancer
data/: Directory containing the cancer data files in TSV format.
plots/: Directory where the generated plots will be saved.
requirements.txt: File listing the Python dependencies required for the project.
README.md: This file, providing an overview of the project.

Usage
Ensure you have the required cancer data files in the data directory. These files should be in TSV format, containing information about protein changes, DNA changes, and affected cases.
Run the data_analysis with data_analysis.py (not necessary if plots are already downloads in plots/)
Run the GUI with gui.py