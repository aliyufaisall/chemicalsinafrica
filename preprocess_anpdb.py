import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors

# --- CONFIGURE PATHS ---
INPUT_CSV_PATH = "raw_anpdb_molecules.csv" 
OUTPUT_CSV_PATH = "processed_anpdb_molecules.csv"

def preprocess_data():
    print(f"Loading data from {INPUT_CSV_PATH}...")
    df = pd.read_csv(INPUT_CSV_PATH)
    
    # 1. Isolate only the required columns
    id_col = 'molecule_id' 
    smiles_col = 'smiles'
    
    if id_col not in df.columns or smiles_col not in df.columns:
        raise KeyError(f"Could not find '{id_col}' or '{smiles_col}' in your CSV columns: {list(df.columns)}")
        
    df = df[[id_col, smiles_col]].copy()
    print(f"Dropped all other columns. Kept: {id_col}, {smiles_col}")
    
    # Clean up empty rows
    df.dropna(subset=[smiles_col], inplace=True)
    
    # 2. Initialize feature lists for your exact requested properties
    mws = []
    logps = []
    tpsas = []
    h_donors = []
    h_acceptors = []
    rot_bonds = []
    ring_counts = []
    fingerprints = []
    valid_mask = []

    print("Calculating requested RDKit properties...")
    for index, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row[smiles_col]))
        
        if mol is not None:
            # 1. Molecular Weight
            mws.append(Descriptors.MolWt(mol))
            # 2. LogP (Lipophilicity)
            logps.append(Descriptors.MolLogP(mol))
            # 3. TPSA (Polarity)
            tpsas.append(Descriptors.TPSA(mol))
            # 4. H-Bond Donors
            h_donors.append(rdMolDescriptors.CalcNumHBD(mol))
            # 5. H-Bond Acceptors
            h_acceptors.append(rdMolDescriptors.CalcNumHBA(mol))
            # 6. Rotatable Bonds
            rot_bonds.append(rdMolDescriptors.CalcNumRotatableBonds(mol))
            # 7. Ring Count
            ring_counts.append(rdMolDescriptors.CalcNumRings(mol))
            # 8. Molecular Fingerprint (Morgan Fingerprint, radius 2, 2048 bits)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            fingerprints.append(fp.ToBitString())  # Storing as standard binary string representation
            
            valid_mask.append(True)
        else:
            mws.append(None)
            logps.append(None)
            tpsas.append(None)
            h_donors.append(None)
            h_acceptors.append(None)
            rot_bonds.append(None)
            ring_counts.append(None)
            fingerprints.append(None)
            valid_mask.append(False)
            print(f"Warning: Invalid SMILES skipped at index {index} (ID: {row[id_col]})")

    # Append new columns to DataFrame
    df['molecular_weight'] = mws
    df['log_p'] = logps
    df['tpsa'] = tpsas
    df['h_bond_donors'] = h_donors
    df['h_bond_acceptors'] = h_acceptors
    df['rotatable_bonds'] = rot_bonds
    df['ring_count'] = ring_counts
    df['fingerprint'] = fingerprints
    
    # Filter out rows that had bad SMILES
    df = df[valid_mask]

    # 3. Save to the final output CSV
    df.to_csv(OUTPUT_CSV_PATH, index=False)
    print(f"Data preparation complete! Saved clean dataset to: {OUTPUT_CSV_PATH}")

if __name__ == "__main__":
    preprocess_data()
