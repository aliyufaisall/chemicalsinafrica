from django.shortcuts import render
from .models import Molecule
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
from rdkit import DataStructs
# Create your views here.

def calculate_tanimoto(fp_str1, fp_str2):
    """Calculates Tanimoto similarity between two bit strings."""
    v1 = DataStructs.CreateFromBitString(str(fp_str1))
    v2 = DataStructs.CreateFromBitString(str(fp_str2))
    
    intersection = (v1 & v2).GetNumOnBits()
    union = (v1 | v2).GetNumOnBits()
    
    if union == 0:
        return 0.0
    return float(intersection) / union

def search_molecule(request):
    context = {}
    
    if request.method == "POST":
        query_smiles = request.POST.get("smiles", "").strip()
        
        # Parse query molecule
        query_mol = Chem.MolFromSmiles(query_smiles)
        if not query_mol:
            context["error"] = "Invalid SMILES string provided. Please try again."
            return render(request, "molecules/search.html", context)
            
        # Extract features from query molecule
        q_mw = Descriptors.MolWt(query_mol)
        q_logp = Descriptors.MolLogP(query_mol)
        q_tpsa = Descriptors.TPSA(query_mol)
        q_hbd = rdMolDescriptors.CalcNumHBD(query_mol)
        q_hba = rdMolDescriptors.CalcNumHBA(query_mol)
        q_rot = rdMolDescriptors.CalcNumRotatableBonds(query_mol)
        q_ring = rdMolDescriptors.CalcNumRings(query_mol)
        
        # Generate fingerprint bit-string for query
        q_fp_obj = AllChem.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)
        q_fp_str = q_fp_obj.ToBitString()
        
        # Load all candidate records from our database table
        db_molecules = Molecule.objects.all()
        ranked_results = []
        
        for db_mol in db_molecules:
            # 1. Structural Fingerprint Similarity (Tanimoto)
            tanimoto_sim = calculate_tanimoto(q_fp_str, db_mol.fingerprint)
            
            # 2. Continuous Physicochemical Scaling (closer to 1.0 is better)
            mw_sim = 1 / (1 + (abs(q_mw - db_mol.molecular_weight) / 100.0))
            logp_sim = 1 / (1 + abs(q_logp - db_mol.log_p))
            tpsa_sim = 1 / (1 + (abs(q_tpsa - db_mol.tpsa) / 50.0))
            
            # --- IMPLEMENTATION OF DECISION 2 ---
            # Continuous decay for count-based features instead of rigid binary true/false switches.
            # Small count variances stay close to 1.0; extreme count deltas drop to near 0.0.
            hbd_sim = 1 / (1 + abs(q_hbd - db_mol.h_bond_donors))
            hba_sim = 1 / (1 + abs(q_hba - db_mol.h_bond_acceptors))
            rot_sim = 1 / (1 + abs(q_rot - db_mol.rotatable_bonds))
            ring_sim = 1 / (1 + abs(q_ring - db_mol.ring_count))
            
            # --- IMPLEMENTATION OF DECISION 3 ---
            # Rebalanced structural weighting matrix. 
            # Structural fingerprint similarity is prioritized at 65%. 
            # Bulk properties scale at 20%, and binding site/functional mechanics map the remaining 15%.
            composite_score = (
                (tanimoto_sim * 0.65) +
                (mw_sim * 0.08) +
                (logp_sim * 0.08) +
                (tpsa_sim * 0.04) +
                (((hbd_sim + hba_sim + rot_sim + ring_sim) / 4.0) * 0.15)
            )
            
            percentage_match = round(composite_score * 100, 2)
            
            ranked_results.append({
                "molecule_id": db_mol.molecule_id,
                "smiles": db_mol.smiles,
                "molecular_weight": db_mol.molecular_weight,
                "log_p": db_mol.log_p,
                "percentage": percentage_match
            })
            
        # Sort values and pick the Top 3 entries
        ranked_results.sort(key=lambda x: x["percentage"], reverse=True)
        context["results"] = ranked_results[:3]
        context["query_smiles"] = query_smiles
        
    return render(request, "molecules/search.html", context)
