# HIV-1 Protease Docking Analysis

This project provides a Python-based tool for analyzing protein-ligand docking results, specifically for HIV-1 protease inhibitors. It calculates binding energies, visualizes protein-ligand interactions, and offers interpretative analysis of the results.

## Features

- Automated download of PDB structures
- Analysis of protein-ligand interactions
- Calculation and visualization of binding energies
- Generation of interaction heatmaps
- Interpretative output for easy understanding of results

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/hiv1-protease-docking-analysis.git
   cd hiv1-protease-docking-analysis
   ```

2. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```

## Usage

1. Ensure you have a CSV file named `hiv1_protease_docking_results.csv` in the project directory. This file should contain your docking results with at least two columns: 'pose_id' and 'binding_energy'.

2. Run the analysis script:
   ```
   python hiv1_protease_docking_analysis.py
   ```

3. The script will generate two output files:
   - `hiv1_protease_interaction_heatmap.png`: A heatmap showing the interactions between ligand atoms and the protein.
   - `hiv1_protease_binding_energy_distribution.png`: A histogram showing the distribution of binding energies.

4. Check the console output for additional interpretative information about the docking results.

## Output Interpretation

- The interaction heatmap shows the number of interactions between each ligand atom and the protein. Darker colors indicate more interactions.
- The binding energy distribution plot shows the range of binding energies across all poses. Lower energies indicate stronger binding.
- The console output provides additional context, including the best binding pose and average number of interactions.

