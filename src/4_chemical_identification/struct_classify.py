import pandas as pd
from rdkit import Chem
import os

# File paths (you may need to update these based on your local system)
negative_file = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\rdkit\positive_matches_with_fingerprints.csv'
output_dir = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\rdkit'

# Load the data
df = pd.read_csv(negative_file)

# Define the updated priority order based on abundance and detectability by LCMS
priority_order = [
    "Lipids", "Carbohydrates", "Phenolics", "Amino Acids", "Glycosides",
    "Terpenoids", "Proteins", "Alkaloids", "Oxylipins", "Phytochelatins", 
    "Volatile Organic Compounds (VOCs)"
]

# Updated SMARTS patterns for matching
biomolecule_smarts = {
    "Carbohydrates": [
        "[C@H]([OH])[C@H](O)[C@H](O)[C@H](O)[C@H](O)[CH2OH]",  # Linear sugar (e.g., glucose, fructose)
        "[OX2H][CX4H]1[OX2][CX4H][CX4H][CX4H][CX4H]1[OX2H]",  # Cyclic sugar (e.g., monosaccharides)
        "[C@H]([OH])[C@H](O)[C@H](O)[C@H](O)[C@H](O)[OX2H]"  # Another variation for sugars
    ],
    "Lipids": [
        "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)O",  # Fatty acids with a long hydrocarbon chain and carboxyl group
        "C(C(CO)O)O",  # Glycerolipid core
        "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=O)[O]",  # Long-chain hydrocarbon with carboxyl
        "[C](=O)[O]"  # Simple fatty acid carboxyl group
    ],
    "Proteins": [
        "[NX3][CX4H]([*])[CX3](=[OX1])[NX3H]",  # Peptide bond pattern
        "[NX3][CX4H]([*])[CX3](=[OX1])[OX2H]",  # General amino acid backbone
        "[#6](-[#7])-[#6](=[#8])-[#7]",  # General peptide backbone pattern
        "[N][C][C](=O)[N]"  # Short peptide segment pattern
    ],
    "Amino Acids": [
        "[NX3H2][CX4H]([*])[CX3](=[OX1])[OX2H]",  # General amino acid
        "[#7][#6](-[*])[#6](=[#8])[#8]",  # General backbone of amino acids
        "[NH2][CH](*)[C](=O)[OH]",  # Common amino acid backbone pattern
        "N[C@@H](C(=O)O)"  # Alanine or similar simple amino acids
    ],
    "Phenolics": [
        "[OX2H][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1",  # Phenolic hydroxyl and benzene ring
        "[OX2H][c]1[c][c][c][c][c]1",  # Another phenol ring variation
        "[OX2H][c]1[c][c]([OX2H])[c][c][c]1"  # Hydroxylated benzene (phenol)
    ],
    "Terpenoids": [
        "[CH3][C](=[CH2])[CH2][CH2][C](=[CH2])[CH3]",  # Isoprene unit (building block of terpenoids)
        "[C](=[CH2])[C](=[CH2])[C](=[CH2])[C](=[CH2])",  # Polyisoprene chain
        "[C@H]1[C@H]([C@H]([C@H](C1)C)[CH2])[CH2]"  # Cyclic terpenoid core
    ],
    "Alkaloids": [
        "[$([NX3H2]),$([NX3H](C)(C)),$([NX4H2+](C)(C)(C))]",  # General alkaloid structure
        "[#7]",  # Any nitrogen-containing structure (alkaloids often contain nitrogen)
        "[#6]~[#7]~[#6]~[#7]",  # Carbon-nitrogen backbone common in alkaloids
        "C1CC[NH+](C1)C2=CC=CC=C2"  # Heterocyclic alkaloid example
    ],
    "Glycosides": [
        "[C,O][C,O][C,O][C,O][C,O][OX2R0,OX2R1]",  # Sugar moiety with potential glycosidic bond
        "[#6]1[#8][#6]([#8])[#6]([#8])[#6]([#8])[#6]1[#8]",  # Cyclic sugar (glycoside linkage)
        "[C@H]1([C@H]([C@H](O[C@H]1O)[OH])[OH])[OH]"  # Glycoside structure
    ],
    "Volatile Organic Compounds (VOCs)": [
        "[#6][#6][#6]",  # Short carbon chain (volatile hydrocarbons)
        "[C,H][C,H][C,H]",  # General small hydrocarbon chains
        "[#6]~[#6]~[#8,#7]",  # Carbon chain with oxygen or nitrogen (volatile functional groups)
        "[CX3H1](=O)[CX3H2]"  # Ketones or aldehydes (common VOCs)
    ],
    "Phytochelatins": [
        "[NX3][CX4H]([CH2]S)[CX3](=[OX1])[NX3H]",  # Cysteine-containing peptide (phytochelatin core)
        "[NX3][CX4H]([*])[CX3](=[OX1])[NX3H]",  # General peptide bond pattern
        "[#7][#6](-[#6]-[#16])[#6](=[#8])[#7]"  # Phytochelatin with sulfur for metal binding
    ],
    "Oxylipins": [
        "[CH2][CH2][CH2]C=O",  # Carbonyl-containing lipid fragment (corrected pattern)
        "[CX4H2][CX3](=[OX1])[CX4H2][CX3H]=[CX3H]",  # Oxygenated fatty acid backbone (oxylipin)
        "[#6]~[#6]~[#6]~[#6]~[#6](=[#8])[#8]"  # Long chain fatty acid with carboxylic group
    ]
}

reviewed_smarts = {
    "Abscisic acid": "CC(=C(C(=O)O)C1=CC(=O)C(=CC1=O)C(C)C)C",
    "Brassinolids": "[#6H3]-[#6H](-[#6H3])-[#6@H](-[#6H3])-[#6@H](-[#6@@H](-[#6@@H](-[#6H3])-[#6@H]1-[#6H2]-[#6H2]-[#6@H]2-[#6@@H]3-[#6H2]-[#8]-[#6](-[#6@H]4-[#6H2]-[#6@@H](-[#6@@H](-[#6H2]-[#6@]-4(-[#6H3])-[#6@H]-3-[#6H2]-[#6H2]-[#6@]-1-2-[#6H3])-[#8H])-[#8H])=[#8])-[#8H])-[#8H]",
    "Ethylene": "C=C",
    "Salicylic Acid": "OC(=O)c1ccccc1O",
    "Jasmonates": "CC(=O)CCC=CC=CC(C)C(=O)O",
    "Auxins": "C1=CC=C2C(=C1)C(=O)NC2=O",
    "Gibberellins": "C1=CC2=C(C=C1)C3C(C2)C4C(C3)C5C(C4)C(C5)C(=O)O",
    "Cytokinins": "C1=NC2=C(N1)C(=NC=N2)N",
    "Strigolactones": "CC1=CC2=C(C=C1)C(=O)OC2=O",
    "Glycerolipids": "C(C(CO)O)O",
    "Phospholipids": "C(C(COP(=O)(O)O)O)O",
    "Galactolipids": "C1C(C(C(C(C(O1)O)O)O)O)CO",
    "Sphingolipids": "NC(CO)C(O)CCCCCCCCCCCCC=C",
    "Superoxide anion": "[O-][O]",
    "Hydrogen peroxide": "OO",
    "Hydroxyl radical": "[OH]",
    "Singlet oxygen": "O=[O]",
    "Fructose": "C(C(C(C(C=O)O)O)O)O",
    "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
    "Proline": "C1CC(NC1)C(=O)O",
    "Lignin": [
        "[C]1=C[C]=C[C]=C1",  # General aromatic ring
        "[C]1=C([O])C=C([C])C=C1",  # Guaiacyl (G) unit
        "[C]1=C([O])C([O])=C([C])C([O])=C1",  # Syringyl (S) unit
        "[C]1=C([O])C=C([C])C=C1"  # p-Hydroxyphenyl (H) unit
    ],
    "Malic acid": "OC(CC(O)=O)C(O)=O"
}

import pandas as pd
from rdkit import Chem
import os

# Function to calculate a score based on SMARTS pattern complexity, match size, and generate confidence score
def get_smarts_score(smiles, smarts_list):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return 0, 0, 0  # Return 0 if the molecule cannot be parsed
    
    best_complexity = 0
    best_atoms_matched = 0
    confidence_score = 0  # Initialize confidence score
    
    # Ensure smarts_list is iterable (either list or single SMARTS pattern)
    if isinstance(smarts_list, str):
        smarts_list = [smarts_list]  # Wrap string into a list if not already a list
    
    # Iterate through the list of SMARTS patterns
    for smarts in smarts_list:
        substruct = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(substruct):
            match = mol.GetSubstructMatch(substruct)
            num_atoms_matched = len(match)  # Number of atoms matched in the molecule
            complexity_score = len(smarts)  # Complexity based on the length of the SMARTS pattern
            
            # Confidence score formula: complexity * atoms matched (adjust this as needed)
            current_confidence_score = complexity_score * num_atoms_matched
            
            # Keep track of the best match and highest confidence score
            if current_confidence_score > confidence_score:
                best_complexity = complexity_score
                best_atoms_matched = num_atoms_matched
                confidence_score = current_confidence_score
                
    return best_complexity, best_atoms_matched, confidence_score  # Return best scores

# Function to determine the best match based on the SMARTS complexity, atom count, and return a scaled confidence score
def best_match_smarts(smiles, biomolecule_smarts, reviewed_smarts, priority_order, group_weights=None):
    mol = Chem.MolFromSmiles(smiles)
    best_match = None
    best_score = 0
    best_atoms_matched = 0
    best_group = None
    best_confidence = 0  # Initialize best confidence score
    all_confidences = []  # Track all confidence scores for scaling
    
    if mol:
        # Check biomolecule_smarts (lists of patterns)
        for group, smarts_list in biomolecule_smarts.items():
            complexity_score, atoms_matched, confidence_score = get_smarts_score(smiles, smarts_list)
            group_weight = group_weights.get(group, 1) if group_weights else 1
            overall_score = (complexity_score * group_weight) + atoms_matched
            
            # Collect all confidence scores for scaling later
            all_confidences.append(confidence_score)
            
            if overall_score > best_score:
                best_score = overall_score
                best_group = group
                best_atoms_matched = atoms_matched
                best_confidence = confidence_score

        # Check reviewed_smarts (single SMARTS patterns, not a list)
        for group, smarts in reviewed_smarts.items():
            complexity_score, atoms_matched, confidence_score = get_smarts_score(smiles, smarts)
            group_weight = group_weights.get(group, 1) if group_weights else 1
            overall_score = (complexity_score * group_weight) + atoms_matched
            
            # Collect all confidence scores for scaling later
            all_confidences.append(confidence_score)
            
            if overall_score > best_score:
                best_score = overall_score
                best_group = group
                best_atoms_matched = atoms_matched
                best_confidence = confidence_score
    
    # Scale the confidence score between min and max confidence across all patterns
    if len(all_confidences) > 0:
        min_confidence = min(all_confidences)
        max_confidence = max(all_confidences)
        if max_confidence != min_confidence:
            # Continuous scaling (min-max normalization)
            scaled_confidence = (best_confidence - min_confidence) / (max_confidence - min_confidence)
        else:
            # If all confidence scores are the same, return 1 (best confidence)
            scaled_confidence = 1
    else:
        scaled_confidence = 0  # No match, confidence is 0
    
    # Apply the priority order if there are ties in scores (use first group with highest priority)
    if best_group in priority_order:
        return best_group, scaled_confidence  # Return best match and scaled confidence score
    else:
        return "Unclassified", 0  # No match, confidence is 0

# Optional: Define group-specific manual weights (higher weight gives priority)
group_weights = {
    "Lipids": 1.2,  # Example of higher weight for Lipids
    "Carbohydrates": 1.1  # Example of slightly higher weight for Carbohydrates
}

# Apply this to your dataframe: looping over each group and applying get_smarts_score
for group, smarts_list in biomolecule_smarts.items():
    df[group] = df['SMILES'].apply(lambda x: get_smarts_score(x, smarts_list)[0])  # Store the complexity score

# For reviewed compounds, use the same approach, but handle strings and lists appropriately
for group, smarts in reviewed_smarts.items():
    df[group] = df['SMILES'].apply(lambda x: get_smarts_score(x, smarts)[0])  # Handle both strings and lists

# Assign priority-based classification using the best match and store the scaled confidence score
df['Best_Classification'], df['Scaled_Confidence_Score'] = zip(*df['SMILES'].apply(
    lambda x: best_match_smarts(x, biomolecule_smarts, reviewed_smarts, priority_order, group_weights)
))

# Define output file path
output_file = os.path.join(output_dir, 'positive_matches_with_fingerprints_with_scaled_confidence.csv')

# Save the updated dataframe with scaled confidence scores
df.to_csv(output_file, index=False)

print(f"File saved at {output_file}")
