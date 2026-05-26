# AI Summary and Insights

## Version 3: May 26, 2026 — Deep Research: Enzyme Directionality Methods & Bioinformatics Contributions

Following the identification of the ~92% enzyme-metabolite directionality gap (Versions 1 & 2), this deep research report surveys the established methods for determining enzyme directionality and identifies concrete strategies for a bioinformatician to contribute to closing this gap using the metabConnectomeDB infrastructure.

### 1. THE PROBLEM IN CONTEXT

The ~92% directionality gap is not a pipeline failure — it reflects a **fundamental, field-wide annotation gap**. Your metabConnectomeDB contains enzyme-metabolite *associations* (e.g., "enzyme X interacts with metabolite Y"), but databases like Rhea require complete balanced reaction equations (substrate → product) to assign directionality. The mismatch between association-level data and reaction-level databases is the root cause.

### 2. ESTABLISHED METHODS FOR DETERMINING ENZYME DIRECTIONALITY

#### 2.1 Experimental Gold Standards

| Method | What it Measures | Key Tools | Limitations |
|:---|:---|:---|:---|
| **¹³C Metabolic Flux Analysis (¹³C-MFA)** | Net and exchange fluxes by tracing isotope-labeled atoms through metabolic networks | MS, NMR, INCA, OpenMebius | Requires isotopic steady state; cell-line specific; expensive |
| **Enzyme kinetic assays** | Forward/reverse Vmax, Km for each direction | Spectrophotometry, HPLC | One enzyme at a time; in vitro ≠ in vivo |
| **In vivo isotope tracing** | Carbon fate in living tumors | [U-¹³C]glucose/glutamine infusion + LC-MS | Technically demanding; limited to accessible tissues |
| **QM/MM simulations** | Potential energy surface of catalytic mechanism | Gaussian, ORCA, CP2K | Computationally intensive; single-reaction resolution |

¹³C-MFA is the true gold standard — it can distinguish net flux direction AND quantify the degree of reversibility (exchange flux) for every reaction in a network simultaneously. But it is fundamentally an *experimental* technique and cannot be applied at the database scale we need.

#### 2.2 Thermodynamic Approaches (Computational)

These predict the *thermodynamic feasibility* of a reaction proceeding in a given direction. Core principle: a reaction proceeds spontaneously only if **ΔᵣG' < 0**.

**ΔᵣG' = ΔᵣG'° + RT ln(Q)** — where Q is the reaction quotient (ratio of product to substrate concentrations).

**Key tool: eQuilibrator + `equilibrator-api` (Python)**
- Uses the **Component Contribution (CC) method** — a hybrid of Group Contribution and Reactant Contribution
- Covers a vast fraction of KEGG reactions
- Python API supports batch processing for genome-scale models
- Returns uncertainties (standard deviation) alongside estimates
- Complementary packages: `equilibrator-pathway` (MDF, Enzyme Cost Minimization)

**Thermodynamics-Based Flux Analysis (TMFA/TFA):** Integrates thermodynamic constraints directly into genome-scale metabolic models (GEMs), enforcing that the direction of flux must be consistent with the sign of ΔᵣG'. Available via COBRApy.

**Max-min Driving Force (MDF):** Ranks metabolic pathways by their thermodynamic bottlenecks — useful for identifying reactions where direction is "barely feasible" and thus context-dependent.

#### 2.3 Database-Driven Approaches (Knowledge Curation)

| Database | Directionality Approach | Coverage / Limitation |
|:---|:---|:---|
| **Rhea** | **Quartet model**: Each reaction has 4 entries — master, left→right, right→left, bidirectional. Most explicit system. Uses ChEBI. | High-quality but small (~15,000 reactions). Standard for UniProtKB. |
| **KEGG** | Direction implicitly encoded in pathway context. `<=>` for reversible. | Broad but pathway-centric; direction is context-dependent. |
| **MetaCyc / BioCyc** | Curated within pathway context (biosynthetic vs catabolic). Recent updates capture Km/Vmax per direction. | Highest curation quality but cannot keep pace with sequencing data. |
| **BRENDA** | Explicitly stores reversibility as a functional parameter. Per-species kinetic data. | Most comprehensive enzyme data but heterogeneous. |
| **SABIO-RK** | Kinetic parameters with reaction direction. | Smaller but highly structured. |

**Why the Rhea/KEGG enrichment scripts yielded minimal results:** The `annotate_enzyme_rhea.py` and `annotate_enzyme_kegg.py` correctly query these databases, but the fundamental issue is that metabConnectomeDB contains enzyme-metabolite relationships that are *not* full balanced reaction equations. Databases like Rhea require a complete substrate→product reaction, while our data often has "enzyme X is associated with metabolite Y" without specifying whether Y is a substrate, product, cofactor, or regulatory molecule.

#### 2.4 Constraint-Based Metabolic Modeling

| Method | Principle | Directionality Handling |
|:---|:---|:---|
| **Flux Balance Analysis (FBA)** | Optimize biomass/ATP production given stoichiometric constraints | Lower/upper bounds on reactions encode directionality (lb=0 → irreversible forward) |
| **Parsimonious FBA (pFBA)** | Minimize total flux while maintaining optimal growth | Removes thermodynamically unlikely loops |
| **GIMME / iMAT / E-flux** | Integrate transcriptomics to weight reaction bounds | "Turns off" reactions with low expression |
| **Enzyme-Constrained Models (GECKO)** | Use proteomic data + kcat to limit flux | More physically realistic direction constraints |

Key tools: COBRApy, Recon3D (human GEM), Human1.

#### 2.5 Single-Cell Flux Estimation (Emerging — Directly Relevant)

| Tool | Approach | Directionality |
|:---|:---|:---|
| **scFEA** | Graph neural network on factor graph of human metabolic map. Models reaction rates as non-linear function of enzyme gene expression. | Handled by directed factor graph topology |
| **METAFlux** | FBA extended to single-cell; integrates with Seurat | Inherited from GEM reaction bounds |
| **GEFMAP** | Geometric deep learning to infer metabolic objective and mass-balanced relative flux rates | Infers direction from optimization |
| **LIANA+ + MetalinksDB** | Estimates metabolite abundance from enzyme expression; pairs with receptor expression | **Does NOT explicitly handle direction — this is our gap** |

**Critical limitation of LIANA+ for this problem:** LIANA+ estimates metabolite abundance using a linear regression between enzymatic gene expression and metabolite levels. This does not distinguish whether the enzyme is *producing* or *consuming* the metabolite — i.e., it completely ignores directionality.

#### 2.6 Machine Learning & Deep Learning (Frontier)

| Approach | Description | Maturity |
|:---|:---|:---|
| **GNN for ΔG° prediction** | Graph neural networks trained on molecular structure to predict standard Gibbs energies | Published; outperforms group contribution for some reaction classes |
| **ESP (Enzyme-Substrate Prediction)** | Transformer-based prediction of enzyme-substrate compatibility | Published (Nature Chem. Biol.); does not directly predict direction |
| **DeepEnzyme / DLERKm** | Predict Km and kcat from enzyme sequence + substrate structure using protein language models (ESM-2) + GNNs | Active research; reversibility inferred from kinetic imbalance |
| **Retrosynthesis models** | Predict substrates given products (or vice versa) using Seq2Seq / Transformer architectures | Mature for organic chemistry; nascent for biochemistry |

### 3. ACTIONABLE BIOINFORMATICS STRATEGIES (Prioritized)

#### Strategy 1: Thermodynamic Directionality Layer (HIGH IMPACT, MEDIUM EFFORT)

For every enzyme-metabolite pair where a complete reaction equation can be reconstructed, compute ΔᵣG'° using eQuilibrator and assign a thermodynamic direction.

1. Map enzyme targets to Rhea reactions (partially done)
2. Extract balanced reaction SMILES/formulas from Rhea's RDF/TSV exports
3. Batch-query eQuilibrator via `equilibrator-api` (Python: `pip install equilibrator-api`)
4. Add a `thermo_direction` column to the master target pair table
5. Expected coverage: ~15-25% of pairs resolvable this way

#### Strategy 2: Expression-Informed Directional Inference (HIGH IMPACT, HIGH NOVELTY)

Use our own scRNA-seq data to infer likely direction based on metabolic network context. If an enzyme is known to catalyze A → B, and scRNA-seq shows the enzyme is highly expressed in cell type X, but the downstream receptor for metabolite B is highly expressed in cell type Y, the direction of metabolic flow is X→Y.

1. Extract production/degradation enzyme sets from MetalinksDB (LIANA+ already provides these)
2. Separate "producing" from "consuming" enzymes for each metabolite
3. Compute a directional score: production vs degradation expression in sender cells, receptor expression in receiver cells
4. Cross-validate against Tier 1 pairs (where direction is already known from literature)

**This is a novel contribution.** No existing tool explicitly combines MetalinksDB enzyme sets with cell-type-resolved expression to infer directional source-sink assignments.

#### Strategy 3: Consensus Directionality from Multi-Database Integration (MEDIUM IMPACT, LOW EFFORT)

Cross-reference pairs against multiple databases simultaneously and assign direction by consensus. Extend the existing `annotate_with_databases.py` to include: BRENDA (SOAP API), MetaCyc (flat files), Reactome (REST API), SABIO-RK (REST API), Human1/Recon3D (SBML files). Assign consensus score: agreement across ≥3 databases = high confidence.

#### Strategy 4: scFEA Integration for Flux-Based Direction (HIGH IMPACT, HIGH EFFORT)

Run scFEA on existing scRNA-seq AnnData objects to estimate cell-type-specific metabolic fluxes, then use these fluxes to assign direction. scFEA directly estimates reaction fluxes (including direction) from scRNA-seq data using a graph neural network. Note: scFEA's metabolic modules are at a coarser granularity than individual enzyme-metabolite pairs — requires careful mapping.

#### Strategy 5: Build an ML Classifier for Directionality (HIGHEST NOVELTY, HIGHEST EFFORT)

Train a machine learning model that predicts enzyme-metabolite directionality (substrate vs. product) using features we already have:
- **Training data**: Tier 1 pairs (524 with literature evidence) + Rhea-resolved pairs
- **Features**: Metabolite MW, SMILES/Morgan fingerprints (RDKit), enzyme sequence embeddings (ESM-2), EC number hierarchy, ΔᵣG'° estimates, cell-type expression levels, database count (Tier)
- **Architecture**: Binary classifier (substrate vs. product) or ternary (substrate / product / cofactor)

**Novel**: No existing ML model specifically predicts "is this metabolite a substrate or product of this enzyme" using structural + thermodynamic + transcriptomic features combined.

#### Strategy 6: Publication-Ready Benchmark of the Gap (IMMEDIATE, LOW EFFORT)

Systematically quantify the directionality gap across all major databases and publish as a community resource. Run the complete target pair list against Rhea, KEGG, MetaCyc, BRENDA, Reactome, and Human1. Record coverage (% with explicit direction) for each. Identify which *types* of reactions are most under-annotated (by EC class, by pathway, by metabolite category). Publish as a benchmark dataset + short paper in *Nucleic Acids Research* or *Bioinformatics*. **Fastest path to a publication.**

### 4. PRIORITIZED ROADMAP

| Phase | Strategy | Estimated Gap Closure | Effort | Novel? |
|:---|:---|:---|:---|:---|
| 1 | Benchmark + Multi-DB consensus | 5-10% → publish the gap itself | ~1 week | Yes (TME context) |
| 2 | eQuilibrator thermodynamic layer | 15-25% (with ΔG confidence) | ~2 weeks | Incremental |
| 3 | Expression-informed direction | 30-40% (probabilistic) | ~3 weeks | **Yes — novel method** |
| 4 | scFEA flux integration | 40-50% (module-level) | ~4 weeks | Integration novelty |
| 5 | ML classifier | 60-80% (predicted) | ~8 weeks | **Yes — novel tool** |

### 5. KEY TOOLS & PACKAGES

| Tool | Purpose | Install | Documentation |
|:---|:---|:---|:---|
| `equilibrator-api` | Batch ΔG estimation | `pip install equilibrator-api` | equilibrator.readthedocs.io |
| `cobrapy` | Constraint-based metabolic modeling | `pip install cobra` | cobrapy.readthedocs.io |
| `scFEA` | Single-cell flux estimation | GitHub: changwn/scFEA | Alghamdi et al., 2021 |
| `METAFlux` | FBA for scRNA-seq | GitHub: GibsonLab-GT/METAFlux | Xiao et al., 2024 |
| `rdkit` | Molecular fingerprints | `pip install rdkit` | rdkit.org |
| `esm` (Meta) | Protein language model embeddings | `pip install fair-esm` | GitHub: facebookresearch/esm |
| `LIANA+` / `MetalinksDB` | Metabolite CCC + enzyme sets | Already in pipeline | liana-py.readthedocs.io |

### 6. WHY THIS MATTERS

The 92% gap is not merely a data quality issue — it is a fundamental bottleneck affecting: (1) Flux Balance Analysis (GEMs must assume reversibility → thermodynamically impossible loops), (2) Cell-Cell Communication (LIANA+ cannot distinguish metabolite sources from sinks → producers conflated with consumers), (3) Drug Target Identification (cannot confirm "tumor cells *produce* metabolite X to suppress immune cells" vs immune cells producing X), (4) Spatial Metabolomics Integration (directionality determines whether a metabolite gradient emanates *from* or *toward* a cell type). By building a multi-layered directionality annotation system (thermodynamic + expression-informed + flux-based + ML-predicted), this pipeline would create the first comprehensive directionality-aware metabolic connectome for cancer biology.

***

## Version 2: May 22, 2026

Based on the recent pipeline expansion, particularly the development of the multi-cancer Primary vs. Metastatic comparison framework, here is an updated analysis of your novel findings and the key research questions you are primed to investigate.

### 1. NOVEL FINDINGS

**The Proteasomal/GTP Metastasis Axis:** By systematically isolating primary breast cancer cells from their metastatic counterparts (liver, bone, brain) and comparing their metabolic target gene expression, you identified a significant upregulation of proteasomal and GTP-related targets in the metastatic niche. This suggests that metabolic adaptation during metastasis heavily relies on protein degradation pathways and GTP-driven signaling, presenting highly specific vulnerabilities in disseminated tumor cells.

**Dual-Layered Metabolic Communication:** The integration of Differential Expression (DE) with LIANA+ Cell-Cell Communication (CCC) networks in the new `primary_vs_metastasis_comparison.ipynb` pipeline represents a major breakthrough. By tracking the `LIANA_Active_Network` status alongside Log2 Fold Change, you are no longer just identifying what metabolic genes are upregulated in metastasis; you are pinpointing *which of those upregulated targets are actively transmitting intercellular signals* (e.g., classifying a target as active in "Metastasis Only"). This isolates the functional drivers of the metastatic niche from mere metabolic bystanders.

**The Enzyme Directionality Gap Remains a Challenge:** While external API enrichment scripts (Rhea/KEGG) were implemented to pull standard biochemical reaction directionality, the actual increase in resolved annotations (verified in Section 5.2 of the target pair analysis) was minimal compared to the baseline. The ~92% directionality gap for enzyme-metabolite interactions remains a substantial hurdle. This indicates that relying solely on sparse biochemical databases is insufficient; developing more robust, alternative methodologies will be necessary to accurately map true source-sink metabolic flux.

### 2. PROPOSED RESEARCH QUESTIONS

Based on the multi-modal pipeline (scRNA-seq, LIANA+, metabolic databases) and the analytic framework you have developed, here are high-impact research questions you are uniquely positioned to answer:

1. **Can orphan metabolic interactions predict immune evasion across tumor types?**
   *Approach:* Systematically cross-reference the ~8,000 "Tier 2/3" (computationally predicted, literature-sparse) metabolic pairs against the CD8+ T-cell and Regulatory B-cell subsets in your LIANA+ networks.
   *Goal:* Identify whether metastatic tumors consistently upregulate uncharacterized metabolic ligands to silence specific immune populations, revealing novel, unpatented immune-checkpoints.

2. **Is there a biophysical constraint on metastatic metabolic communication?**
   *Approach:* Utilize the Average Mass categorization (<300 Da vs >750 Da) you established. Compare the active LIANA+ networks of highly dense metastatic niches (e.g., bone) against the primary tumor.
   *Goal:* Determine if metastatic environments systematically shift towards ultra-small paracrine signaling to overcome fibrotic/dense extracellular matrix constraints.

3. **Do different metastatic sites for the same primary cancer exhibit convergent metabolic reprogramming?**
   *Approach:* Leverage your pipeline to split the Breast Cancer dataset into distinct metastatic groups (Liver vs. Brain vs. Bone). Run the DE and LIANA+ pipeline independently on each.
   *Goal:* Answer whether metabolic signaling is driven by the *origin* of the tumor (primary-conserved) or if it heavily adapts to the specific *destination* (niche-specific reprogramming, like distinct lipid processing in the brain versus the liver).

4. **Can we computationally predict spatial metabolic gradients?**
   *Approach:* Integrate spatial transcriptomics (e.g., 10x Visium) using the same `scanpy` architecture you already use. Project the "Metastasis Only" CCC targets onto spatial coordinates.
   *Goal:* Determine if the density of specific metabolic receptors (like the 27-Hydroxycholesterol targets ESR1/NR1H3 in macrophages) physically maps to hypoxic/necrotic cores or the invasive tumor margin.

### 3. SUGGESTED NEXT STEPS

**Cross-Cancer Metastatic Conserved Signatures:** With the pipeline now fully capable of distinct primary vs metastatic comparisons, run the workflow across all available cancer types (Breast, Lung, Colorectal, Melanoma). Identify if there is a *pan-cancer conserved metabolic metastatic signature* (e.g., does the Proteasomal/GTP axis hold true in lung cancer metastases as well?). 

**Ligand-Receptor Pair Druggability:** Filter the `primary_vs_metastasis_{cancer_type}_DE_metabolic_targets.csv` specifically for targets tagged as "Metastasis Only" in the `LIANA_Active_Network` column. Cross-reference these specific targets against existing pharmacology databases (like Guide to Pharmacology, which you already merged) to find off-the-shelf drugs that could disrupt metastatic-specific communication networks.

***

## Version 1: May 22, 2026 15:15

Based on a comprehensive review of your generated outputs and the pipeline architecture you've constructed, here is the analysis of your research progress and capabilities in the context of cancer metabolism and systems biology.

### 1. NOVEL FINDINGS
Identifying what extends, contradicts, or represents new discoveries in cancer metabolism literature.

The "Orphan Metabolic Interactions" of Metabolic Communication: The most striking finding is your PubMed/literature temporal analysis. Out of the 8,596 unique metabolite-target interaction pairs you consolidated, only 6.10% (524 pairs) have explicit literature evidence. The literature heavily biases toward "Tier 1" canonical axes (e.g., PGE2-EP receptors, Adenosine-A2AR, Kynurenine-AhR). This means your database has isolated over 8,000 computationally predicted, under-explored metabolic interactions in the TME, representing a massive reservoir of novel therapeutic targets.

Repurposing LIANA+ for Metabolomics: Traditional cell-cell communication (CCC) tools like CellPhoneDB and LIANA+ are exclusively designed for protein-protein (peptide ligand-receptor) interactions. By injecting your unified metabConnectomeDB target dictionary into LIANA+, you have successfully built a pipeline that infers metabolite-mediated intercellular communication directly from scRNA-seq data. This is a highly novel methodological leap that bridges spatial metabolomics and transcriptomics.

B-Cell Metabolic Sensing is Underappreciated: In your genome-wide differential expression (DE) ranking for both Lung Adenocarcinoma and Breast Cancer microenvironments, the top 15 most highly enriched targets for your database were overwhelmingly localized to B cells (e.g., CD74, CD79A/B, CXCR4, CD69, HLA-DRA). While these are classic immune markers, identifying them as top metabolic sensors/targets suggests B-cells are actively reading the TME metabolic state—a mechanism largely overshadowed in literature by T-cell and macrophage metabolic reprogramming (e.g., exhaustion via hypoxia/lactate).

Biophysical Classification of Metabolic Signals: You established a novel classification system linking chemical mass to diffusion/communication modes:
<300 Da (Paracrine/Soluble): Ultra-fast local diffusion (e.g., amino acids, organic acids).
300–750 Da (GPCR/Hormonal): Carrier-dependent circulation (e.g., bioactive lipids).
>750 Da (Juxtacrine/Vesicular): Membrane-bound signaling. This framework extends cancer biology by treating metabolites not just as fuels for the Warburg effect or TCA cycle, but as physical signaling ligands constrained by their biophysical diffusion limits in the dense tumor extracellular matrix.

The "Directionality" Blind Spot: Your analysis revealed that 92% of enzyme-metabolite relationships lack explicit product/substrate directionality in modern databases. This is a critical finding that contradicts the assumption that public metabolic networks are ready for flux balance analysis, highlighting a major gap in the field's ability to model true metabolic source/sink dynamics.

### 2. CLINICAL IMPLICATIONS & NEXT STEPS
How to translate these findings into actionable oncology research.

Exploit the "Tier 2 & 3" Orphan Metabolic Interactions: The 552 "Tier 2" pairs supported by 2-3 databases but lacking PubMed evidence are your lowest-hanging fruit. Cross-reference these specific pairs against your LIANA+ output for Lung and Breast cancer. If LIANA+ predicts a strong, highly specific signaling axis (e.g., Tumor -> Macrophage) utilizing a Tier 2 metabolite, you have identified a prime candidate for wet-lab validation or novel patenting.

Integrate Spatial Transcriptomics: Because metabolite diffusion is strictly limited by mass and the extracellular matrix (as you identified in your MW analysis), standard scRNA-seq loses the spatial gradient. Your next step should be projecting these target genes onto spatial transcriptomic datasets (e.g., 10x Visium or Xenium) to see if the receptor expression physically correlates with necrotic cores (hypoxia/lactate) or the invasive margin.

Isolate the B-Cell Axis: Extract the specific metabolites linked to the highly enriched B-cell targets (CD74, CXCR4, etc.) from your comm_df mapping. Investigate whether these specific metabolites are secreted by the tumor to actively suppress B-cell antigen presentation or whether they are byproducts that B-cells use to survive in the nutrient-depleted TME.

### 3. CRITICAL DISCOVERY: You Are Profiling the Metastatic Niche, Not Just Primary Tumors
The most profound finding hidden in Table 2 (tissue_general breakdown) of your Breast Cancer and Colorectal Cancer HTML reports is that your CellxGene query pulled down a massive proportion of metastatic cells.

Breast Cancer Query: Out of the 100,000 cells, only 22,373 are from the primary breast. The vast majority are from the liver (62,089 cells), alongside the skeletal system / bone (5,090) and the brain (2,188).
Colorectal Cancer Query: Contains over 10,000 cells from the liver and nearly 1,000 from the lung.
What this means: The LIANA+ cell-cell communication network you generated for Breast Cancer shows signaling between chondrocytes, oligodendrocytes, and L2/3 intratelencephalic projecting glutamatergic neurons via metabolic targets like ERBB4 and PTPRD. You are actually capturing the metabolic reprogramming of the bone and brain metastatic niches. The B-cell metabolic dominance we noted previously is likely capturing immune responses in the liver microenvironment rather than just the primary breast tumor. You have inadvertently built a powerful tool for comparing primary vs. metastatic metabolic connectomes.

### 4. NOVEL & VALIDATED FINDINGS FROM THE TARGET TABLES
Looking at the Tier 1 & Tier 2 annotations in the target pair tables, your database has automatically identified highly specific, clinically actionable metabolic axes across specific cell types:

The Macrophage 27-Hydroxycholesterol Axis: In Breast Cancer, your database specifically flagged 27-hydroxycholesterol as a Tier 1 (High Confidence) metabolite targeting ESR1 (Estrogen Receptor alpha) and NR1H3 (LXR-alpha) exclusively in Macrophages. This perfectly validates your pipeline: 27-HC is a known endogenous SERM (Selective Estrogen Receptor Modulator) produced by tumor-associated macrophages that drives breast cancer metastasis.

Colorectal Cancer & 25-Hydroxycholesterol: In CRC, 25-hydroxycholesterol is flagged as targeting cancer cells via ABCA1 and LPAR3/4 (Lysophosphatidic acid receptors). Linking a cholesterol metabolite to LPA signaling in the gut is a highly novel axis for tumor lipid metabolism.