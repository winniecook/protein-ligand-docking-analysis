# hiv1_protease_docking_analysis.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBList
from scipy.spatial.distance import cdist
import os
import seaborn as sns
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

def download_pdb(pdb_id, file_type="pdb"):
    pdbl = PDBList()
    filename = pdbl.retrieve_pdb_file(pdb_id, pdir=".", file_format=file_type)
    return filename

def read_docking_results(filename):
    return pd.read_csv(filename)

def parse_pdb(filename):
    parser = PDBParser()
    structure = parser.get_structure("molecule", filename)
    return structure

def extract_ligand(structure, ligand_id):
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_id:
                    return residue
    return None

def calculate_interactions(ligand_coords, protein_coords, cutoff=3.5):
    distances = cdist(ligand_coords, protein_coords)
    interactions = np.sum(distances < cutoff, axis=1)
    return interactions

def plot_interaction_heatmap(interactions, ligand_atoms):
    plt.figure(figsize=(12, 10))
    heatmap_data = pd.DataFrame(interactions, index=ligand_atoms, columns=['Interactions'])
    ax = sns.heatmap(heatmap_data, cmap='YlOrRd', annot=True, fmt='d', cbar_kws={'label': 'Number of Interactions'})
    plt.xlabel("HIV-1 Protease")
    plt.ylabel("Ligand Atoms")
    plt.title("HIV-1 Protease-Inhibitor Interactions")
    plt.tight_layout()
    plt.savefig("hiv1_protease_interaction_heatmap.png")
    plt.close()

def plot_binding_energies(energies):
    plt.figure(figsize=(10, 6))
    sns.histplot(energies, bins=20, kde=True)
    plt.xlabel("Binding Energy (kcal/mol)")
    plt.ylabel("Number of Poses")
    plt.title("Distribution of Binding Energies for HIV-1 Protease Inhibitor")
    
    # Add mean and best energy lines
    mean_energy = np.mean(energies)
    best_energy = np.min(energies)
    plt.axvline(mean_energy, color='r', linestyle='--', label=f'Mean Energy: {mean_energy:.2f} kcal/mol')
    plt.axvline(best_energy, color='g', linestyle='--', label=f'Best Energy: {best_energy:.2f} kcal/mol')
    
    plt.legend()
    plt.savefig("hiv1_protease_binding_energy_distribution.png")
    plt.close()

def main():
    # Download PDB file
    pdb_id = "1HVR"
    pdb_filename = download_pdb(pdb_id)
    
    # Rename the downloaded file
    os.rename(pdb_filename, f"{pdb_id}.pdb")
    pdb_filename = f"{pdb_id}.pdb"

    docking_results = read_docking_results("hiv1_protease_docking_results.csv")
    structure = parse_pdb(pdb_filename)
    
    # Extract ligand (XK2 is the ligand ID for 1HVR)
    ligand = extract_ligand(structure, "XK2")
    
    if ligand is None:
        print("Ligand not found in the PDB file.")
        return

    protein_coords = np.array([atom.coord for atom in structure.get_atoms() if atom.element != 'H' and atom.get_parent() != ligand])
    ligand_coords = np.array([atom.coord for atom in ligand.get_atoms() if atom.element != 'H'])
    ligand_atoms = [atom.name for atom in ligand.get_atoms() if atom.element != 'H']

    interactions = calculate_interactions(ligand_coords, protein_coords)

    plot_binding_energies(docking_results['binding_energy'])
    plot_interaction_heatmap(interactions, ligand_atoms)

    # Find the best binding pose
    best_pose = docking_results.loc[docking_results['binding_energy'].idxmin()]
    print("Best binding pose:")
    print(f"Pose ID: {best_pose['pose_id']}")
    print(f"Binding Energy: {best_pose['binding_energy']} kcal/mol")

    # Calculate average number of interactions
    avg_interactions = np.mean(interactions)
    print(f"Average number of interactions per ligand atom: {avg_interactions:.2f}")
    
    # Print additional interpretative information
    print(f"\nInterpretation:")
    print(f"1. The binding energy range is from {docking_results['binding_energy'].min():.2f} to {docking_results['binding_energy'].max():.2f} kcal/mol.")
    print(f"2. Lower binding energies indicate stronger binding. The best pose has a binding energy of {best_pose['binding_energy']} kcal/mol.")
    print(f"3. Ligand atoms with more than {avg_interactions:.0f} interactions are considered to have strong interactions with the protein.")
    print("4. Check the heatmap to identify which ligand atoms have the strongest interactions with the protein.")

if __name__ == "__main__":
    main()
