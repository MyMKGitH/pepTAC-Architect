\# 🧬 pepTAC-Architect: Automated In Silico Peptide Design Platform



pepTAC-Architect is a streamlined computational biology application engineered to automate the pre-screening, mutation, and optimization of peptide candidates targeting therapeutic proteins, with a specific focus on \*\*Human Choline Kinase alpha (ChoKα)\*\*. 



By integrating direct text-based Protein Data Bank (PDB) queries, formal \*\*Biopython\*\* physicochemical profiling, and structured configuration generation for \*\*AlphaFold 3\*\*, this app serves as a high-throughput sequence optimization pipeline that filters out unviable candidates before intensive structural folding simulations.



\## 🚀 Features



\- \*\*Automated Target Ingestion:\*\* Search the RCSB Protein Data Bank natively using standard protein names or unique PDB IDs to instantly fetch structural chains and isolate co-crystallized peptide fragments.

\- \*\*Biopython Screening Engine:\*\* Computes exact empirical physical parameters including Isoelectric Point (pI), Hydropathy (GRAVY score), Instability Index, and Aromaticity.

\- \*\*Custom Heuristic Scorer:\*\* Automatically evaluates mutated candidate variants against a specific target binding grid (favoring positive charge and aromatic stacking required for Choline mimicry).

\- \*\*AlphaFold 3 Manifest Generator:\*\* Generates and exports structured `.json` system manifests perfectly formatted for direct upload into the AlphaFold 3 Server dashboard.



\## 🛠️ Installation \& Usage



1\. \*\*Clone the repository:\*\*

```bash

&#x20;  git clone \[https://github.com/MyMKGitH/pepTAC-Architect.git](https://github.com/MyMKGitH/pepTAC-Architect.git)

&#x20;  cd pepTAC-Architect



```

&#x20;2. \*\*Install dependencies:\*\*

```bash

&#x20;  pip install -r requirements.txt



```

&#x20;3. \*\*Run the local Streamlit server:\*\*

```bash

&#x20;  streamlit run app.py



```

\## 📊 Pipeline Workflow Overview

&#x20;1. \*\*Query \& Mine:\*\* Query target sequences via names or database tokens. Extract natural binder templates.

&#x20;2. \*\*Mutate \& Screen:\*\* Run high-throughput systematic single-point substitution matrices across 7 distinct amino acid properties.

&#x20;3. \*\*Export \& Fold:\*\* Download the programmatic JSON blueprint file and upload it directly into the AlphaFold 3 interface for final 3D complex alignment confirmation.

\## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.



Author

MK

