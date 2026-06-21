import streamlit as st
import pandas as pd
import requests
import json
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# -----------------------------------------------------------------------------
# DATABASE QUERY ENGINES (RCSB Search & Fetch API)
# -----------------------------------------------------------------------------

def search_pdb_by_name(query_text):
    """Programmatically queries RCSB text service to find relevant structural IDs."""
    if not query_text.strip():
        return []
    
    search_payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": f'"{query_text.strip()}"'
            }
        },
        "return_type": "entry"
    }
    
    try:
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        response = requests.post(url, json=search_payload, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if "result_set" in data:
                return [item["identifier"] for item in data["result_set"]]
        return []
    except Exception:
        return []

def fetch_sequences_from_pdb(pdb_id):
    """Fetches multi-FASTA format and maps candidate structural chains."""
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return None
        
        lines = response.text.split('\n')
        chains = {}
        current_header = None
        current_seq = []
        
        for line in lines:
            if line.startswith('>'):
                if current_header and current_seq:
                    chains[current_header] = "".join(current_seq)
                current_header = line.strip()
                current_seq = []
            elif line.strip():
                current_seq.append(line.strip())
        if current_header and current_seq:
            chains[current_header] = "".join(current_seq)
            
        return chains
    except Exception:
        return None

# -----------------------------------------------------------------------------
# PHYSICOCHEMICAL ANALYSIS LOOP (Biopython Integration)
# -----------------------------------------------------------------------------

def calculate_biopython_metrics(seq):
    """Computes exact empirical physical parameters via ProtParam package."""
    # Clean string input to capture canonical structural amino acids
    clean_seq = "".join([aa for aa in seq.upper() if aa in "ACDEFGHIKLMNPQRSTVWY"])
    if not clean_seq:
        return {"Isoelectric Point": 7.0, "GRAVY": 0.0, "Instability Index": 0.0, "Aromaticity": 0.0}
    
    analysis = ProteinAnalysis(clean_seq)
    try:
        pi = analysis.isoelectric_point()
        gravy = analysis.gravy()
        instability = analysis.instability_index()
        aromaticity = analysis.aromaticity()
    except Exception:
        pi, gravy, instability, aromaticity = 7.0, 0.0, 0.0, 0.0
        
    return {
        "Isoelectric Point": round(pi, 2),
        "GRAVY": round(gravy, 2),
        "Instability Index": round(instability, 2),
        "Aromaticity": round(aromaticity, 2)
    }

def score_for_choka(metrics):
    """
    Heuristic Evaluation Matrix for Human Choline Kinase alpha:
    - High pI indicates strong positive charge needed to mimic choline headgroup.
    - Balanced GRAVY preserves aqueous stability while ensuring target pocket entry.
    - Low Instability Index (< 40) validates physiological lifespan integrity.
    - Aromaticity tracks required framework for cation-pi interaction stacking.
    """
    score = 0
    
    # Target 1: Charge Distribution (Mimic choline headgroup)
    if metrics["Isoelectric Point"] >= 9.0:
        score += 6
    elif metrics["Isoelectric Point"] >= 7.5:
        score += 3
    else:
        score -= 4
        
    # Target 2: Hydrophobic Boundary Balance
    if -0.5 <= metrics["GRAVY"] <= 0.6:
        score += 5
        
    # Target 3: Structural Lifespan Integrity
    if metrics["Instability Index"] < 40:
        score += 4  # Classed as a stable chemical peptide architecture
    else:
        score -= 2  # Classed as biochemically unstable
        
    # Target 4: Core Aromatic Scaffolding 
    if metrics["Aromaticity"] >= 0.15:
        score += 5
        
    return round(score, 1)

# -----------------------------------------------------------------------------
# GRAPHICAL INTERACTION TERMINAL
# -----------------------------------------------------------------------------
st.set_page_config(page_title="pepTAC Automated Framework", layout="wide")

st.title("🧬 pepTAC Automated Sequence Designer")
st.markdown("Automated Protein Ingestion Search, Multi-Property Evaluation, and AlphaFold 3 Handoff")
st.write("---")

col_left, col_right = st.columns([1, 2])

with col_left:
    st.header("1. Target Discovery & Mining")
    
    search_mode = st.radio("Target Search Ingestion Mode:", ["Search by Protein Name", "Enter PDB ID Directly"])
    pdb_target = ""
    
    if search_mode == "Search by Protein Name":
        name_input = st.text_input("Enter Target Protein Name:", value="Choline Kinase alpha")
        if name_input:
            with st.spinner("Searching Protein Data Bank..."):
                discovered_ids = search_pdb_by_name(name_input)
            if discovered_ids:
                st.success(f"Found {len(discovered_ids)} matching structural files.")
                pdb_target = st.selectbox("Select Target PDB ID:", discovered_ids)
            else:
                st.warning("No entries found matching that query name string.")
    else:
        pdb_target = st.text_input("Enter Known PDB ID:", value="2CKO").strip()

    target_protein_seq = ""
    seed_peptide_seq = ""

    if pdb_target:
        with st.spinner("Harvesting Sequence Records..."):
            extracted_chains = fetch_sequences_from_pdb(pdb_target)
            
        if extracted_chains:
            st.success(f"Loaded Core Records for ID: {pdb_target.upper()}")
            
            proteins = {}
            peptides = {}
            for header, seq in extracted_chains.items():
                chain_info = header.split('|')[1] if '|' in header else header
                if len(seq) > 50:
                    proteins[chain_info] = seq
                elif len(seq) >= 5:
                    peptides[chain_info] = seq
            
            if proteins:
                selected_target_label = st.selectbox("Confirm Target Receptor Chain:", list(proteins.keys()))
                target_protein_seq = proteins[selected_target_label]
                st.caption(f"Receptor String Length: {len(target_protein_seq)} Amino Acids")
            
            if peptides:
                st.info(f"Detected {len(peptides)} potential co-crystallized peptide fragments.")
                selected_pep_label = st.selectbox("Select Extracted Seed Fragment:", list(peptides.keys()))
                seed_peptide_seq = peptides[selected_pep_label]
            else:
                st.warning("No short peptide fragments found. Initializing manual sequence anchor:")
                seed_peptide_seq = st.text_input("Manual Peptide Sequence Anchor:", value="RKVFLF").upper().strip()
                
            seed_peptide_seq = st.text_input("Confirm Sequence Mutation Blueprint:", value=seed_peptide_seq).upper().strip()
        else:
            st.error("Invalid database record reference hook.")

with col_right:
    st.header("2. Systematic Mutagenesis Terminal")
    
    if seed_peptide_seq and target_protein_seq:
        st.markdown("Calculating exact physicochemical metrics for point-mutation profiles...")
        
        variants_pool = []
        
        # Ingest Reference Sequence parameters
        base_metrics = calculate_biopython_metrics(seed_peptide_seq)
        base_score = score_for_choka(base_metrics)
        variants_pool.append({
            "Mutation Node": "Wildtype (Baseline)",
            "Sequence Matrix": seed_peptide_seq,
            **base_metrics,
            "Target Affinity Score": base_score
        })
        
        # Systematically execute amino acid substitutions
        search_residues = ['K', 'R', 'F', 'W', 'Y', 'L', 'A']
        for idx in range(len(seed_peptide_seq)):
            for aa in search_residues:
                if seed_peptide_seq[idx] != aa:
                    mut_list = list(seed_peptide_seq)
                    mut_list[idx] = aa
                    mut_seq = "".join(mut_list)
                    
                    metrics = calculate_biopython_metrics(mut_seq)
                    score = score_for_choka(metrics)
                    
                    variants_pool.append({
                        "Mutation Node": f"Position {idx+1} ➔ {aa}",
                        "Sequence Matrix": mut_seq,
                        **metrics,
                        "Target Affinity Score": score
                    })
                    
        df_variants = pd.DataFrame(variants_pool).drop_duplicates(subset=["Sequence Matrix"])
        df_variants = df_variants.sort_values(by="Target Affinity Score", ascending=False).reset_index(drop=True)
        
        st.dataframe(df_variants, use_container_width=True)
        
        # -----------------------------------------------------------------------------
        # ALPHAFOLD 3 INTEGRATION PARSER
        # -----------------------------------------------------------------------------
        st.write("---")
        st.header("3. AlphaFold 3 Integration Handoff")
        st.markdown("Export structural config blueprints to validate binding models on the AlphaFold Server platform.")
        
        optimal_seq = st.selectbox("Select Target Lead for Structural Prediction:", df_variants["Sequence Matrix"].tolist())
        
        af3_manifest = {
            "name": f"pepTAC_Design_Target_{pdb_target.upper()}",
            "modelModifiers": [],
            "sequences": [
                {"protein": {"sequence": target_protein_seq, "count": 1}},
                {"protein": {"sequence": optimal_seq, "count": 1}}
            ]
        }
        
        manifest_string = json.dumps(af3_manifest, indent=2)
        st.code(manifest_string, language="json")
        
        st.download_button(
            label="💾 Download AlphaFold 3 Manifest File",
            data=manifest_string,
            file_name=f"af3_config_{pdb_target.lower()}.json",
            mime="application/json"
        )
        
        st.markdown(
            """
            > 🚀 **How to process your run:**
            > 1. Download the `.json` file containing your chosen structural parameters.
            > 2. Navigate directly to the web submission terminal at [AlphaFold Server](https://alphafoldserver.com/).
            > 3. Click **Upload JSON** on your panel dashboard to load your target complex sequence without tedious copy-pasting.
            """
        )
    else:
        st.info("Awaiting protein database lookup coordinates.")
