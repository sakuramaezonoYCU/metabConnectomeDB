import pandas as pd
import os

def generate_sub_definitions(cols, definitions_map, output_csv, output_md, title, desc):
    """
    Modular utility to map columns to their definitions, categories, and provenance,
    writing out clean CSV and formatted Markdown summaries.
    """
    data = []
    for col in sorted(cols):
        if col in definitions_map:
            db, category, definition = definitions_map[col]
        else:
            db = 'multiple'
            category = 'Interaction Evidence & Scores'
            definition = 'Consolidated metadata field.'
        data.append({
            'header': col,
            'database': db,
            'category': category,
            'definition': definition
        })
        
    df_out = pd.DataFrame(data)
    df_out.to_csv(output_csv, index=False)
    print(f"🎉 Successfully created column definitions at {output_csv}")
    
    # Save beautiful Markdown grouped by Category
    with open(output_md, 'w', encoding='utf-8') as f:
        f.write(f"# {title}\n\n")
        f.write(f"{desc}\n\n")
        
        categories = df_out['category'].unique()
        for cat in sorted(categories):
            f.write(f"## 📁 {cat}\n\n")
            f.write("| Header | Database | Definition |\n")
            f.write("| --- | --- | --- |\n")
            df_cat = df_out[df_out['category'] == cat]
            for _, r in df_cat.iterrows():
                f.write(f"| `{r['header']}` | **{r['database']}** | {r['definition']} |\n")
            f.write("\n")
            
    print(f"🎉 Successfully created markdown documentation at {output_md}")

def main():
    project_root = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB"
    
    # Core Database Inputs
    pair_csv = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
    metab_csv = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_with_HMDB_Info.csv')
    
    # Unique Metabolite Exploration Notebook Output Inputs
    tier1_csv = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_tier1.csv')
    tier2_csv = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_tier2.csv')
    
    # Metabolite-Target Pair Analysis Notebook Output Inputs
    cancer_csv = os.path.join(project_root, 'output', 'human_metab_target_pairs_cancer_annotated.csv')

    # Main Database Outputs
    main_csv = os.path.join(project_root, 'output', 'merged_database_col_definition.csv')
    main_md = os.path.join(project_root, 'output', 'merged_database_col_definition.md')
    
    # Notebook Output Definitions Outputs
    unique_csv = os.path.join(project_root, 'output', 'unique_metab_exploration_col_definition.csv')
    unique_md = os.path.join(project_root, 'output', 'unique_metab_exploration_col_definition.md')
    
    pair_analysis_csv = os.path.join(project_root, 'output', 'metab_target_pair_analysis_col_definition.csv')
    pair_analysis_md = os.path.join(project_root, 'output', 'metab_target_pair_analysis_col_definition.md')
    
    # Comprehensive Column Definitions Map
    # Maps exact column name -> (database provenance, functional category, description)
    definitions_map = {
        'Metabolite_Name': ('multiple', 'Metabolite Properties & Structure', 'The standardized common name of the metabolite.'),
        'Target': ('multiple', 'Target Protein Properties', 'The standardized HGNC gene symbol of the target protein (receptor, enzyme, or transporter).'),
        'HMDB_ID': ('multiple', 'Metabolite Properties & Structure', 'The unique identifier in the Human Metabolome Database (HMDB).'),
        'Sensor_Gene': ('MEBOCOST', 'Target Protein Properties', 'The sensing gene in MEBOCOST, representing the protein interacting with the metabolite.'),
        'Protein_Name': ('multiple', 'Target Protein Properties', 'Full descriptive name of the target protein.'),
        'Evidence': ('multiple', 'Interaction Evidence & Scores', 'Supporting evidence or category for the interaction.'),
        'Sensor_Type': ('MEBOCOST', 'Target Protein Properties', 'Functional type of the sensor gene (e.g. Receptor, Transporter, Enzyme).'),
        'database': ('multiple', 'Interaction Evidence & Scores', 'The source database(s) where this entry originates.'),
        'Evidence_Score': ('multiple', 'Interaction Evidence & Scores', 'Numerical score representing interaction confidence.'),
        'Text_Evidence': ('multiple', 'Interaction Evidence & Scores', 'Literature text snippets or PMIDs showing interaction evidence.'),
        'blood_concentration': ('MEBOCOST', 'Metabolite Localization & Physiology', 'Physiological concentration of the metabolite in human blood.'),
        'Secondary_HMDB_ID': ('HMDB', 'Metabolite Properties & Structure', 'Alternate HMDB IDs associated with the metabolite.'),
        'Kegg_ID': ('multiple', 'Metabolite Properties & Structure', 'Unique identifier in the Kyoto Encyclopedia of Genes and Genomes (KEGG) database.'),
        'Synonyms': ('multiple', 'Metabolite Properties & Structure', 'Alternative names and chemical synonyms for the metabolite.'),
        'Organ and components': ('MEBOCOST', 'Metabolite Localization & Physiology', 'Anatomical organs/components where the metabolite is localized.'),
        'Cell and elements': ('MEBOCOST', 'Metabolite Localization & Physiology', 'Cell types where the metabolite is localized.'),
        'Tissue and substructures': ('MEBOCOST', 'Metabolite Localization & Physiology', 'Tissue types and substructures where the metabolite is localized.'),
        'Biofluid and excreta': ('MEBOCOST', 'Metabolite Localization & Physiology', 'Human biofluids or excretions containing the metabolite.'),
        'Subcellular': ('MEBOCOST', 'Metabolite Localization & Physiology', 'Subcellular compartment or organelle localization.'),
        'Kingdom': ('HMDB', 'Metabolite Properties & Structure', 'Top-level taxonomic category for the chemical classification of the metabolite.'),
        'Super_Class': ('HMDB', 'Metabolite Properties & Structure', 'Broad chemical category of the metabolite (ClassyFire super_class level).'),
        'Sub_Class': ('HMDB', 'Metabolite Properties & Structure', 'Lower-level chemical subclass of the metabolite (ClassyFire sub_class level).'),
        'Class': ('multiple', 'Metabolite Properties & Structure', 'Standard chemical class of the metabolite.'),
        'associated_gene': ('MEBOCOST', 'Synthesis & Transport Genes', 'Metabolic genes associated with the synthesis or transport of the metabolite.'),
        'Software_Predicted': ('MEBOCOST', 'Interaction Evidence & Scores', 'Indicator of whether the metabolite\'s localization or interaction was computationally predicted.'),
        'BioLocation_Summary': ('MEBOCOST', 'Metabolite Localization & Physiology', 'Text summary outlining the biological locations of the metabolite.'),
        'Receptor_Gene_Symbol': ('multiple', 'Target Protein Properties', 'HGNC gene symbol of the receptor protein interacting with the metabolite.'),
        'Transporter': ('multiple', 'Target Protein Properties', 'HGNC gene symbol of the transporter protein importing/exporting the metabolite.'),
        'Enzyme': ('multiple', 'Target Protein Properties', 'HGNC gene symbol of the enzyme interacting with or processing the metabolite.'),
        'Task': ('multiple', 'Metabolic Modeling (scCellFie)', 'High-level signaling task or metabolic complex (e.g. in scCellFie or CellPhoneDBv5).'),
        'transmembrane': ('multiple', 'Target Protein Properties', 'Indicates if the target protein is a transmembrane protein.'),
        'peripheral': ('multiple', 'Target Protein Properties', 'Indicates if the target protein is a peripheral membrane protein.'),
        'secreted': ('multiple', 'Target Protein Properties', 'Indicates if the target protein is secreted into the extracellular space.'),
        'secreted_desc': ('MEBOCOST', 'Target Protein Properties', 'Descriptive text detailing the secretion properties of the target protein.'),
        'secreted_highlight': ('MEBOCOST', 'Target Protein Properties', 'Highlighted annotation specifying secreted status in MEBOCOST.'),
        'integrin': ('multiple', 'Target Protein Properties', 'Indicates if the target protein belongs to the integrin family.'),
        'other': ('MEBOCOST', 'Target Protein Properties', 'Indicates if the target protein has other functional classifications in MEBOCOST.'),
        'other_desc': ('MEBOCOST', 'Target Protein Properties', 'Functional description of other classifications in MEBOCOST.'),
        'pdb_structure': ('HMDB', 'Target Protein Properties', 'Protein Data Bank (PDB) structural identifier for the metabolite/target.'),
        'comments_complex': ('multiple', 'Target Protein Properties', 'Descriptive comments regarding the protein complex structure.'),
        'reactome_reaction': ('multiple', 'Metabolic Modeling (scCellFie)', 'Reactome reaction database identifier.'),
        'reactome_complex': ('multiple', 'Target Protein Properties', 'Reactome complex database identifier.'),
        'rhea_reaction': ('multiple', 'Metabolic Modeling (scCellFie)', 'Rhea reaction database identifier.'),
        'curator': ('multiple', 'Interaction Evidence & Scores', 'Database curator who verified or uploaded the entry.'),
        'CellphoneDB_ver': ('CellPhoneDBv5', 'Interaction Evidence & Scores', 'CellPhoneDBv5 database version supporting this interaction.'),
        'Uniprot': ('multiple', 'Target Protein Properties', 'UniProt accession identifier for the target protein.'),
        'Protein_Uniprot': ('multiple', 'Target Protein Properties', 'UniProt accession identifier for metabolic enzyme/transporter target proteins (historically UNIPROT_ID).'),
        'uniprot_1': ('multiple', 'Target Protein Properties', 'UniProt accession identifier for the first subunit/protein in a complex.'),
        'uniprot_2': ('multiple', 'Target Protein Properties', 'UniProt accession identifier for the second subunit/protein in a complex.'),
        'uniprot_3': ('multiple', 'Target Protein Properties', 'UniProt accession identifier for the third subunit/protein in a complex.'),
        'uniprot_4': ('multiple', 'Target Protein Properties', 'UniProt accession identifier for the fourth subunit/protein in a complex.'),
        'Class_source': ('multiple', 'Metabolite Properties & Structure', 'Source database for the chemical classification.'),
        'GEM_id': ('multiple', 'Metabolic Modeling (scCellFie)', 'Genome-scale Metabolic Model (GEM) identifier.'),
        'Cell_Compartment': ('multiple', 'Metabolite Localization & Physiology', 'Subcellular cellular compartment where the reaction occurs.'),
        'Synthetic_genes': ('multiple', 'Synthesis & Transport Genes', 'Metabolic genes responsible for synthesizing the metabolite.'),
        'Synthetic_genes_lv2': ('multiple', 'Synthesis & Transport Genes', 'Lower-confidence or level-2 metabolic genes responsible for synthesis.'),
        'Metaligand_transporter_genes_out': ('multiple', 'Synthesis & Transport Genes', 'Genes responsible for exporting or transporting the metabolite (historically Transporter_genes).'),
        'Metaligand_transporter_genes_in': ('multiple', 'Synthesis & Transport Genes', 'Genes responsible for importing the metabolite into cells (historically Transporter_genes_in).'),
        'precursor_transport': ('multiple', 'Synthesis & Transport Genes', 'Notes regarding transport of precursor molecules.'),
        'MRID': ('multiple', 'Interaction Evidence & Scores', 'Unique identifier in MRCLinkDB.'),
        'PubChem CID/SID': ('HMDB', 'Metabolite Properties & Structure', 'PubChem chemical compound or substance identifier.'),
        'Molecular Formula': ('HMDB', 'Metabolite Properties & Structure', 'Standard chemical molecular formula of the metabolite.'),
        'CHEMICAL_FORMULA': ('HMDB', 'Metabolite Properties & Structure', 'Standard chemical molecular formula of the metabolite.'),
        'Canonical SMILES': ('HMDB', 'Metabolite Properties & Structure', 'Standard Simplified Molecular-Input Line-Entry System (SMILES) representation.'),
        'SMILES': ('HMDB', 'Metabolite Properties & Structure', 'Standard Simplified Molecular-Input Line-Entry System (SMILES) representation.'),
        'Receptor_Gene_ID': ('multiple', 'Target Protein Properties', 'Entrez gene identifier for the receptor.'),
        'Receptor_Uniprot': ('multiple', 'Target Protein Properties', 'UniProt accession identifier for the receptor protein.'),
        'PMID': ('multiple', 'Interaction Evidence & Scores', 'PubMed identifier(s) for literature citations supporting the interaction.'),
        'Other.DB': ('multiple', 'Interaction Evidence & Scores', 'References or IDs in other chemical/biological databases.'),
        'HMDB_Protein_ID': ('multiple', 'Target Protein Properties', 'Identifier in other database consolidations (historically HAMDBP_ID).'),
        'Enzyme_Full_Name': ('multiple', 'Target Protein Properties', 'IUPAC name or common name of the enzyme (historically ENZYME_NAME).'),
        'Human_geneID': ('multiple', 'Target Protein Properties', 'Entrez Gene ID for the human target.'),
        'REACTIONS': ('multiple', 'Metabolic Modeling (scCellFie)', 'Biochemical reactions associated with the pair.'),
        'enzyme product/substrate': ('multiple', 'Synthesis & Transport Genes', 'Describes if the metabolite is a substrate or product of the target enzyme.'),
        'Interaction': ('MRCLinkDB', 'Clinical Oncology & TME Phenotypes', 'Direct clinical interaction type (e.g., regulation, activation, inhibition).'),
        'Cell type': ('MRCLinkDB', 'Clinical Oncology & TME Phenotypes', 'Clinical immune or somatic cell type mediating the interaction (T-cells, macrophages, etc.).'),
        'Experimental subject': ('MRCLinkDB', 'Clinical Oncology & TME Phenotypes', 'Experimental model or cell line used to verify the interaction.'),
        'Disease': ('MRCLinkDB', 'Clinical Oncology & TME Phenotypes', 'Explicit cancer/tumor clinical disease annotation.'),
        'Effect': ('MRCLinkDB', 'Clinical Oncology & TME Phenotypes', 'Functional effect of the metabolite on cells (e.g., pro-tumorigenic).'),
        'Effect_detail': ('MRCLinkDB', 'Clinical Oncology & TME Phenotypes', 'Detail regarding the functional effect on cells.'),
        'Target_Gene': ('multiple', 'Target Protein Properties', 'General gene symbol identifier (historically Gene_Name).'),
        'scCellFie_value': ('scCellFie', 'Metabolic Modeling (scCellFie)', 'Reaction pathway activity score or value.'),
        'System': ('scCellFie', 'Metabolic Modeling (scCellFie)', 'Top-level metabolic system category in scCellFie (e.g., Carbohydrate metabolism).'),
        'Subsystem': ('scCellFie', 'Metabolic Modeling (scCellFie)', 'Metabolic subsystem pathway in scCellFie (e.g., Glycolysis).'),
        'ensembl_id': ('scCellFie', 'Target Protein Properties', 'Ensembl gene identifier.'),
        'sccellfie_threshold': ('scCellFie', 'Metabolic Modeling (scCellFie)', 'Expression threshold value from scCellFie.'),
        'databases_count': ('multiple', 'Interaction Evidence & Scores', 'Number of integrated databases supporting this pair/metabolite.'),
        'HMDB_Name': ('HMDB', 'Metabolite Properties & Structure', 'Standardized metabolite name from HMDB.'),
        'INCHIKEY': ('HMDB', 'Metabolite Properties & Structure', 'Standard IUPAC International Chemical Identifier key.'),
        'AVERAGE_MASS': ('HMDB', 'Metabolite Properties & Structure', 'Average molecular mass of the metabolite.'),
        'MONO_MASS': ('HMDB', 'Metabolite Properties & Structure', 'Monoisotopic mass of the metabolite.'),
        'Pair_Confidence_Tier': ('multiple', 'Interaction Evidence & Scores', 'Confidence tier of the specific metabolite-target interaction pair.'),
        'Metabolite_Confidence_Tier': ('multiple', 'Interaction Evidence & Scores', 'Confidence tier of the metabolite itself.'),
        'Confidence_Tier': ('multiple', 'Interaction Evidence & Scores', 'Confidence tier category (Tier 1: high, Tier 2: medium, Tier 3: low) based on database overlap.'),
        'receptor': ('multiple', 'Target Protein Properties', 'Indicates if the target protein acts as a receptor in CellPhoneDBv5 or other source databases.'),
    }

    # 1. Main merged database column definitions
    main_cols = set()
    if os.path.exists(pair_csv):
        df_pair = pd.read_csv(pair_csv, nrows=1)
        main_cols.update(df_pair.columns)
    if os.path.exists(metab_csv):
        df_metab = pd.read_csv(metab_csv, nrows=1)
        main_cols.update(df_metab.columns)
        
    generate_sub_definitions(
        main_cols, definitions_map, main_csv, main_md,
        "MetabConnectomeDB Merged Database Column Definitions",
        "This file lists all unified database column headers in the consolidated database, their functional categories, their original source database (or `multiple` if shared/common), and their exact definitions."
    )
    
    # 2. Unique metabolite exploration notebook outputs definitions
    unique_cols = set()
    if os.path.exists(tier1_csv):
        df_t1 = pd.read_csv(tier1_csv, nrows=1)
        unique_cols.update(df_t1.columns)
    if os.path.exists(tier2_csv):
        df_t2 = pd.read_csv(tier2_csv, nrows=1)
        unique_cols.update(df_t2.columns)
        
    generate_sub_definitions(
        unique_cols, definitions_map, unique_csv, unique_md,
        "Unique Metabolite Exploration Column Definitions (Notebook Outputs)",
        "This file lists the column headers, categories, and definitions for the outputs of the unique metabolite exploration Jupyter notebook (`unique_metab_data_exploration.ipynb`), specifically `human_database_merge_unique_metab_tier1.csv` and `human_database_merge_unique_metab_tier2.csv`."
    )
    
    # 3. Metabolite-Target pair analysis notebook outputs definitions
    pair_analysis_cols = set()
    if os.path.exists(cancer_csv):
        df_cancer = pd.read_csv(cancer_csv, nrows=1)
        pair_analysis_cols.update(df_cancer.columns)
        
    generate_sub_definitions(
        pair_analysis_cols, definitions_map, pair_analysis_csv, pair_analysis_md,
        "Metabolite Target Pair Analysis Column Definitions (Notebook Outputs)",
        "This file lists the column headers, categories, and definitions for the outputs of the metabolite-target pair analysis Jupyter notebook (`metab_targetPair_analysis.ipynb`), specifically the clinical and cancer-annotated file `human_metab_target_pairs_cancer_annotated.csv`."
    )

if __name__ == '__main__':
    main()
