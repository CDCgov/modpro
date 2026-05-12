# MODPRO

MOD-eling PRO-teins Pipeline

The MODPRO pipeline was developed by the CDC's Influenza Division to analyze AlphaFold-predicted structures for antigenic drift of influenza hemagglutinin proteins. Viral strain fitness is predicted in terms of antigenic drift using a fitness landscape model (see repo referenced below) relative to a user-selected vaccine reference. It uses open-source programs along with custom CDC code. With minor tweaking, this program is usable by external partners but can also be run on the CDC cluster. The generated datasets offer insights into the evolutionary fitness of viral strains without the need for in-house sequencing or propagation efforts as long as high-quality sequencing is publicly available. The pipeline simplifies and visualizes complex sequence surveillance data, aiding leadership in making informed decisions regarding public health. It can be used to generate new training datasets and perform _ad hoc_ comparisons of biophysical characteristics between individual isolates or reference strains.

Structural bioinformatics pipeline orchestrated by a weekly (Mondays) NiFi cron job that queries ref_sys database top_model_priority table to identify high-priority HA variants not yet predicted by Openfold3 but targeted for experimentation by the laboratory. Pipeline then generates a run manifest after checking HDFS for duplicates and logs metadata before publishing manifest to GWA, triggering Snakemake structural prediction pipeline. The planned update includes database cleanup, automated variant processing, enhanced RBS epitope analysis with sialic acid binding prediction feature, and a strategic shift to HA1-focused relaxation and analytics. 

For each variant, amino acid sequences are retrieved from CDP, converted to fasta files, submitted for MSA generation, OpenFold3 structure prediction, and phenotypic analyses [glycosylation distance calculation, RBS epitopes, structural distance calculations etc] all the while capturing QC metrics such as full-length pLDDT, logging execution and CDP publication status. Full-length structures are archived, followed by HA1-focused processing that includes chain- and domain- trimming, Rosetta scoring, RMS/MaxSub evaluation, biophysical relaxation, with decision points enforcing QC thresholds and selecting improved structures based on logged improvement of energy score and RMSD. All intermediate and final metrics recorded in protein_modeling “logging” and “qc” tables. NiFi handles observability and quality control by ingesting checkpointed outputs from Snakemake, publishing validated files to HDFS, and writing SQL metadata and run results. The pipeline distinguishes clearly between compute and publication states, where “compute_status” tracks scientific execution outcomes (compute_status: PENDING, RUNNING, SUCCESS, QC_FAILED, FAILED) and “publication_status” tracks data delivery (publication_status:  NOT_READY, READY, POSTING, HDFS_POSTED, SQL_POSTED, COMPLETE, PUBLISH_FAILED), enabling robust retry, auditability, and reliable end-to-end processing.

Benchmarks:
Monomer/Trimer dataset of ~10 recent crystal structures for H1, H3, Hvic
Sialic Acid Binding Dataset <~40 known H5 sialic acid binding data, glycan array data>
Relaxed pocket scores— add in pipeline, comment out until passes QC
Minimal database for running openfold3 

Weekly Variant Structural Prediction Pipeline Flow
Updated Process Flow:
 -	Automated SQL queries against ref_sys database to identify high-priority variants (deduplication checks against existing protein_modeling entries)  
 -	Automated submission to OpenFold3 for structure prediction
 -	HA1 domain extraction, optimization by rosetta relaxation protocols for energy minimization (energy optimization enforcement for relaxed structures)
 -	Comprehensive logging and audit trails, error handling with retry mechanisms

Automated Variant Processing Pipeline: REFSYS ensures that limited laboratory resources (60-100 HI tests and 10-15 HINT tests per week per subtype) are allocated to specimens that will provide the maximum surveillance value.
The six-tier classification system is particularly elegant - it balances the need for comprehensive coverage with statistical rigor, ensuring that well-characterized variants receive inferred results while directing precious testing resources toward variants that need additional characterization or represent emerging threats.
The business case is compelling because it transforms what was likely a time-intensive, subjective process into an automated, statistically-driven system that improves both efficiency and quality while supporting CDC's broader public health mission and international surveillance partnerships. The updated pipeline should pull variants for weekly identification and processing of priority variants from REFSYS by automated SQL pull. 

Enhanced RBS Epitope Analysis: This epitope prediction tool represents integration of structural biology and immunology principles. The 15 Angstrom distance criterion is large enough to capture the full span of conformational epitopes while remaining within the physical reach of antibody binding domains. Inclusion of both surface accessibility and glycosylation analysis filter out sites without antibody accessibility. The centroid-based approach is especially valuable for RBS analysis since this region is the most functionally and immunologically important site for antigenic drift. 

Centroid-prediction tool pulls amino acids within 15 Angstroms with surface accessibility and glycosylation interaction. 
Benchmarks for RBS prediction tool: 

Sialic Acid Binding Site Analysis: Optimized Boltzmann scoring, benchmarks

Database Schema Optimization: Remove deprecated tables and implement new analytical structures and optimize join strategies 
HA1 Domain Focus: Transition from full HA to HA1-specific modeling
Targeted Tables for Removal:
unrelaxed_pdb_data - Contains unprocessed structural data
<convert to HA1 only, pLDDT average for whole and HA1>
secondary_structure - Legacy secondary structure annotations
test - Development/testing artifacts

Update Tables:  
 - 	calculated_rbs_epitope_residues - Centralized epitope analysis with geometric calculations. Update table from existing measurements to include a list of the predicted residues based on the criteria explained above. 
 - 	pipeline metadata - logging and weekly cron pull of variant runs
 - 	Sialic acid binding predictions [off mode until benchmarks passed]

HA1-Specific Modeling: Transitioning from full HA to HA1-focused analysis represents a strategic shift toward more targeted and efficient modeling because of increased functional relevance: HA1 contains the receptor binding domain and primary antigenic sites; better evolutionary focus as most antigenic drift occurs in the HA1 globular head domain. Full PDB will be stored in the database for ad hoc analysis. Critical receptor binding and neutralization sites are HA1-localized and enhanced accuracy better folding predictions for smaller, more stable domains.

Computational Efficiency: Smaller domain size enables higher resolution analysis as well as containing the region with the highest pLDDT values (highest confidence) for more reliable analysis.

This pipeline generates biophysical information on each structure and includes the following, in addition to more granular features :

  Qualitative Metrics
- Rosetta Energy Scores, total and per-residue (a measure of relative stability of protein in question)
- Root mean square difference - RMS and MaxSub (a scoring of the deviation in alignment between the protein and a reference structure)
- Relative surface accessibility per residue

  Biophysical Properties Calculators (CDC-generated)
- Glycosylation distance calculator
- RBS epitopes prediction calculator
- Structural distance calculator for generating fitness landscapes <https://git.biotech.cdc.gov/hxa6/fit_landscape>

## Full Workflow Diagram
```mermaid

flowchart TB
    %% =====================================================
    %% STYLES
    %% =====================================================
    classDef nifi fill:#EAF2FF,stroke:#2F6FED,stroke-width:1.5px,color:#1F2937;
    classDef model fill:#ECFDF3,stroke:#1F9D63,stroke-width:1.5px,color:#1F2937;
    classDef analytic fill:#F4FDF7,stroke:#159A5B,stroke-width:1.5px,color:#1F2937;
    classDef store fill:#FFF4E8,stroke:#D97706,stroke-width:1.5px,color:#1F2937;
    classDef decision fill:#FFFBEA,stroke:#CA8A04,stroke-width:1.5px,color:#1F2937;
    classDef report fill:#F3F0FF,stroke:#7C3AED,stroke-width:1.5px,color:#1F2937;

    %% =====================================================
    %% CONTROL PLANE
    %% =====================================================
    subgraph A["Structural Bioinformatics Workflow"]
        direction LR
        A1["Weekly Nifi Cron"]
        A2["ExecuteSQL<br/>Pull top ref_sys variants"]
        A3["Weekly REFSYS Pull<br/>Normalize and deduplicate"]
        A4["Post to GWA"]
        A5["Post to HDFS"]
        A6["Trigger Snakemake workflow from new GWA file"]
        A1 --> A2 --> A3 --> A4 --> A5 --> A6
    end

    %% =====================================================
    %% MODELING PLANE
    %% =====================================================
    subgraph B["Snakemake Workflow Part1"]
        direction LR

        subgraph B1["Input Configuration"]
            direction TB
            B1a["Ingest sequence"]
            B1b["Create seq.fasta"]
            B1c["Create msa_config.json"]
            B1d["Create seq_name.json"]
            B1a --> B1b --> B1c --> B1d
        end

        subgraph B2["Prediction"]
            direction TB
            B2a["Run MSAs"]
            B2b["Run OpenFold3"]
            B2c["Report HA1 pLDDT"]
            B2d{"QC threshold met?"}
            B2e["Re-run prediction"]
            B2a --> B2b --> B2c --> B2d
            B2d -- No --> B2e --> B2a
        end

        subgraph B3["HA1 Processing and Selection"]
            direction TB
            B3a["Remove B and C chains and HA2"]
            B3b["Archive full-length PDB to GWA"]
            B3c["Rosetta score for HA1"]
            B3d["RMS and MaxSub for HA1"]
            B3e["Rosetta relax"]
            B3f["Retain higher-quality HA1<br/>based on improved RMSD"]
            B3a --> B3b --> B3c --> B3d --> B3e --> B3f
        end

        B1d --> B2a
        B2d -- Yes --> B3a
    end

    %% =====================================================
    %% ANALYTICS PLANE
    %% =====================================================
    subgraph C["Snakemake Workflow Part 2"]
        direction LR
        C1["Glycosylation distance"]
        C2["GetContacts"]
        C3["FITscape"]
        C4["Centroid prediction and RBS residues"]
        C1 --> C2 --> C3 --> C4
    end

    %% =====================================================
    %% OBSERVABILITY / DELIVERY
    %% =====================================================
    subgraph D["Observability and Delivery"]
        direction LR
        D1["Logging"]
        D2["Run reporting"]
        D3["NiFi PostHDFS"]
        D4["Write SQL metadata and results"]
        D5["Summary outputs"]
        D1 --> D2 --> D3 --> D4 --> D5
    end

    %% =====================================================
    %% STORAGE
    %% =====================================================
    E1["GWA"]
    E2["HDFS"]
    E3["SQL Metadata and Results"]

    %% =====================================================
    %% FLOW LINKS
    %% =====================================================
    A4 --> E1
    A5 --> E2
    A6 --> B1a
    B3b --> E1
    B3f --> C1
    C4 --> D1
    D3 --> E2
    D4 --> E3
    D5 --> E2

    %% =====================================================
    %% CLASS ASSIGNMENTS
    %% =====================================================
    class A1,A2,A3,A4,A5,A6,D3,D4 nifi;
    class B1a,B1b,B1c,B1d,B2a,B2b,B2c,B2e,B3a,B3b,B3c,B3d,B3e,B3f model;
    class C1,C2,C3,C4 analytic;
    class E1,E2,E3 store;
    class B2d decision;
    class D1,D2,D5 report;
    
 ```
Pipeline logging and reporting is orchestrated by Nifi, which pulls from the Snakemake-generated JSON logs and posts to HDFS in near real-time. These can be accessed by SQL query of the CDP database: protein_modeling.modpro.logging.

DAGs for the Monday/Thursday Snakemake pipelines:

<svg width="800" height="1000" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M2 1L8 5L2 9" fill="none" stroke="#333" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
  </defs>
  
  <!-- Title -->

  <!-- Execution flow indicator -->
  <path d="M720 130 L720 400" stroke="#ccc" stroke-width="2" stroke-dasharray="4,4" fill="none"/>
  <text x="725" y="150" font-family="Arial, sans-serif" font-size="12" fill="#666" writing-mode="tb">Execution order</text>
</svg>
<img width="800" height="1000" alt="snakemake_dag_svg (1)" src="https://github.com/user-attachments/assets/7c7ea9b2-1c44-41b7-a5f7-eed8f960fd20" />

<img width="900" height="1200" alt="rosetta_ha1_pipeline_dag (1)" src="https://github.com/user-attachments/assets/af61d4ab-d675-4587-918b-5b0e8a3a2f69" />
![Uploading rosett<svg width="900" height="1200" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M2 1L8 5L2 9" fill="none" stroke="#333" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="decision-arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M2 1L8 5L2 9" fill="none" stroke="#d32f2f" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
  </defs>
  
  <!-- Title -->
  <text x="450" y="25" text-anchor="middle" font-family="Arial, sans-serif" font-size="18" font-weight="bold" fill="#333">Rosetta HA1 Processing Pipeline (Thursday Trigger)</text>
  
  <!-- Input trigger -->
  <rect x="390" y="50" width="120" height="44" rx="8" fill="#E1F5EE" stroke="#0F6E56" stroke-width="2"/>
  <text x="450" y="68" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#0F6E56">PDB detected</text>
  <text x="450" y="84" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#0F6E56">Thursday trigger</text>

  <!-- Chain processing (moved first) -->
  <rect x="370" y="130" width="160" height="56" rx="8" fill="#E6F1FB" stroke="#185FA5" stroke-width="2"/>
  <text x="450" y="150" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#185FA5">remove_chains</text>
  <text x="450" y="166" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#185FA5">Remove B,C,HA2</text>
  <text x="450" y="180" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#185FA5">Keep HA1 only</text>

  <!-- Baseline assessment (after chain removal) -->
  <rect x="370" y="220" width="160" height="56" rx="8" fill="#EEEDFE" stroke="#534AB7" stroke-width="2"/>
  <text x="450" y="240" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#534AB7">rosetta_baseline</text>
  <text x="450" y="256" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#534AB7">Pre-relax metrics</text>
  <text x="450" y="270" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#534AB7">RMS, MaxSub, pLDDT</text>

  <!-- Rosetta relax -->
  <rect x="370" y="310" width="160" height="56" rx="8" fill="#FAEEDA" stroke="#854F0B" stroke-width="2"/>
  <text x="450" y="330" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#854F0B">rosetta_relax</text>
  <text x="450" y="346" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#854F0B">Energy minimization</text>
  <text x="450" y="360" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#854F0B">Structure optimization</text>

  <!-- Post-relax scoring -->
  <rect x="370" y="400" width="160" height="56" rx="8" fill="#EEEDFE" stroke="#534AB7" stroke-width="2"/>
  <text x="450" y="420" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#534AB7">post_relax_score</text>
  <text x="450" y="436" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#534AB7">Post-relax metrics</text>
  <text x="450" y="450" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#534AB7">Energy assessment</text>

  <!-- Metrics comparison and selection -->
  <rect x="370" y="490" width="160" height="56" rx="8" fill="#FBEAF0" stroke="#993556" stroke-width="2"/>
  <text x="450" y="510" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#993556">select_best</text>
  <text x="450" y="526" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#993556">Compare energy scores</text>
  <text x="450" y="540" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#993556">Keep lower energy</text>

  <!-- Optional pocket relax -->
  <rect x="160" y="400" width="160" height="56" rx="8" fill="#FAECE7" stroke="#993C1D" stroke-width="2"/>
  <text x="240" y="420" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#993C1D">pocket_relax</text>
  <text x="240" y="436" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#993C1D">Optional targeted</text>
  <text x="240" y="450" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#993C1D">relaxation</text>

  <!-- Pocket relax metrics -->
  <rect x="580" y="490" width="160" height="56" rx="8" fill="#EEEDFE" stroke="#534AB7" stroke-width="2"/>
  <text x="660" y="510" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#534AB7">pocket_metrics</text>
  <text x="660" y="526" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#534AB7">Post-pocket relax</text>
  <text x="660" y="540" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#534AB7">assessment</text>

  <!-- Parallel structural analyses -->
  <rect x="50" y="580" width="160" height="56" rx="8" fill="#EAF3DE" stroke="#3B6D11" stroke-width="2"/>
  <text x="130" y="600" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#3B6D11">glyco_analysis</text>
  <text x="130" y="616" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#3B6D11">Glycosylation</text>
  <text x="130" y="630" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#3B6D11">distance</text>

  <rect x="240" y="580" width="160" height="56" rx="8" fill="#E1F5EE" stroke="#0F6E56" stroke-width="2"/>
  <text x="320" y="600" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#0F6E56">fitscape</text>
  <text x="320" y="616" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#0F6E56">Centroid prediction</text>
  <text x="320" y="630" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#0F6E56">RBS residues</text>

  <rect x="430" y="580" width="160" height="56" rx="8" fill="#FBEAF0" stroke="#993556" stroke-width="2"/>
  <text x="510" y="600" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#993556">epitope_analysis</text>
  <text x="510" y="616" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#993556">GetContacts</text>
  <text x="510" y="630" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#993556">analysis</text>

  <rect x="620" y="580" width="160" height="56" rx="8" fill="#fff3e0" stroke="#f57c00" stroke-width="2"/>
  <text x="700" y="600" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#f57c00">logging</text>
  <text x="700" y="616" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#f57c00">Step-by-step</text>
  <text x="700" y="630" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#f57c00">status tracking</text>

  <!-- Report generation -->
  <rect x="270" y="670" width="160" height="56" rx="8" fill="#EEEDFE" stroke="#534AB7" stroke-width="2"/>
  <text x="350" y="690" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#534AB7">generate_report</text>
  <text x="350" y="706" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#534AB7">Collect all metrics</text>
  <text x="350" y="720" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#534AB7">Write SQL metadata</text>

  <!-- Archive and upload -->
  <rect x="120" y="760" width="160" height="56" rx="8" fill="#F1EFE8" stroke="#5F5E5A" stroke-width="2"/>
  <text x="200" y="780" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#5F5E5A">archive_results</text>
  <text x="200" y="796" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#5F5E5A">Store optimized</text>
  <text x="200" y="810" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#5F5E5A">structures</text>

  <!-- HDFS upload -->
  <rect x="420" y="760" width="160" height="56" rx="8" fill="#E6F3FF" stroke="#0066cc" stroke-width="2"/>
  <text x="500" y="780" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#0066cc">hdfs_upload</text>
  <text x="500" y="796" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#0066cc">Post to CDP</text>
  <text x="500" y="810" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#0066cc">Update status</text>

  <!-- Final status -->
  <rect x="270" y="850" width="160" height="56" rx="8" fill="#e8f5e8" stroke="#4a9" stroke-width="2"/>
  <text x="350" y="870" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#4a9">pipeline_complete</text>
  <text x="350" y="886" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#4a9">All outputs</text>
  <text x="350" y="900" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#4a9">confirmed</text>

  <!-- Dependency arrows -->
  <!-- Input to chain removal -->
  <line x1="450" y1="94" x2="450" y2="128" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Chain removal to baseline -->
  <line x1="450" y1="186" x2="450" y2="218" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Baseline to relax -->
  <line x1="450" y1="276" x2="450" y2="308" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Relax to post-relax scoring -->
  <line x1="450" y1="366" x2="450" y2="398" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Post-relax scoring to selection -->
  <line x1="450" y1="456" x2="450" y2="488" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Optional pocket relax branch -->
  <path d="M410 428 L280 428" fill="none" stroke="#999" stroke-width="1.5" stroke-dasharray="3,3" marker-end="url(#arrow)"/>
  <text x="345" y="420" font-family="Arial, sans-serif" font-size="10" fill="#666">optional</text>

  <!-- Pocket relax to pocket metrics -->
  <path d="M400 428 L580 518" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Selection and optional pocket metrics to parallel analyses -->
  <path d="M410 546 L170 578" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M420 546 L320 578" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M480 546 L510 578" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M490 546 L700 578" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Optional pocket metrics to analyses -->
  <path d="M620 546 L170 578" fill="none" stroke="#999" stroke-width="1.5" stroke-dasharray="3,3" marker-end="url(#arrow)"/>
  <path d="M640 546 L320 578" fill="none" stroke="#999" stroke-width="1.5" stroke-dasharray="3,3" marker-end="url(#arrow)"/>
  <path d="M660 546 L510 578" fill="none" stroke="#999" stroke-width="1.5" stroke-dasharray="3,3" marker-end="url(#arrow)"/>

  <!-- All parallel analyses to report -->
  <path d="M170 636 L310 668" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M320 636 L340 668" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M510 636 L390 668" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M700 636 L430 668" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Report to final outputs -->
  <path d="M310 726 L240 758" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M390 726 L460 758" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Final convergence -->
  <path d="M240 816 L310 848" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M460 816 L390 848" fill="none" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>

  <!-- Data outputs -->
  <rect x="50" y="880" width="140" height="28" rx="4" fill="#f0f0f0" stroke="#999" stroke-width="1"/>
  <text x="120" y="897" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#333">pre_relax_energy_score</text>

  <rect x="200" y="880" width="130" height="28" rx="4" fill="#f0f0f0" stroke="#999" stroke-width="1"/>
  <text x="265" y="897" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#333">post_relax_energy_score</text>

  <rect x="340" y="880" width="140" height="28" rx="4" fill="#f0f0f0" stroke="#999" stroke-width="1"/>
  <text x="410" y="897" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#333">pocket_relax_metrics</text>

  <rect x="490" y="880" width="130" height="28" rx="4" fill="#f0f0f0" stroke="#999" stroke-width="1"/>
  <text x="555" y="897" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#333">predicted_RBS_residues</text>

  <rect x="630" y="880" width="120" height="28" rx="4" fill="#f0f0f0" stroke="#999" stroke-width="1"/>
  <text x="690" y="897" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#333">centroid_xyz</text>

  <rect x="760" y="880" width="100" height="28" rx="4" fill="#f0f0f0" stroke="#999" stroke-width="1"/>
  <text x="810" y="897" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#333">glyco_distances</text>

  <!-- Final outputs -->
  <rect x="150" y="930" width="120" height="28" rx="4" fill="#e8f5e8" stroke="#4a9" stroke-width="1"/>
  <text x="210" y="947" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#4a9">optimized_HA1.pdb</text>

  <rect x="290" y="930" width="120" height="28" rx="4" fill="#e8f5e8" stroke="#4a9" stroke-width="1"/>
  <text x="350" y="947" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#4a9">metrics_report.json</text>

  <rect x="430" y="930" width="120" height="28" rx="4" fill="#e8f5e8" stroke="#4a9" stroke-width="1"/>
  <text x="490" y="947" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#4a9">epitope_analysis.csv</text>

  <rect x="570" y="930" width="120" height="28" rx="4" fill="#e8f5e8" stroke="#4a9" stroke-width="1"/>
  <text x="630" y="947" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#4a9">hdfs_upload_status</text>

  <!-- Schedule annotation -->
  <rect x="50" y="50" width="200" height="36" rx="6" fill="#fff3e0" stroke="#f57c00" stroke-width="1"/>
  <text x="150" y="65" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" font-weight="bold" fill="#f57c00">Schedule: Thursday trigger</text>
  <text x="150" y="80" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#f57c00">Cron: 0 9 * * 4</text>

  <!-- Resource requirements -->
  <rect x="650" y="50" width="200" height="36" rx="6" fill="#e8f4f8" stroke="#1976d2" stroke-width="1"/>
  <text x="750" y="65" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" font-weight="bold" fill="#1976d2">Rosetta requirements</text>
  <text x="750" y="80" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#1976d2">CPU: 8 cores, Mem: 16GB</text>

  <!-- Legend -->
  <text x="50" y="990" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#333">Rule types</text>
  
  <rect x="50" y="1010" width="16" height="12" rx="2" fill="#E1F5EE" stroke="#0F6E56"/>
  <text x="74" y="1020" font-family="Arial, sans-serif" font-size="12" fill="#333">Input/Analysis</text>
  
  <rect x="170" y="1010" width="16" height="12" rx="2" fill="#E6F1FB" stroke="#185FA5"/>
  <text x="194" y="1020" font-family="Arial, sans-serif" font-size="12" fill="#333">Processing</text>
  
  <rect x="270" y="1010" width="16" height="12" rx="2" fill="#EEEDFE" stroke="#534AB7"/>
  <text x="294" y="1020" font-family="Arial, sans-serif" font-size="12" fill="#333">Assessment</text>
  
  <rect x="370" y="1010" width="16" height="12" rx="2" fill="#FAEEDA" stroke="#854F0B"/>
  <text x="394" y="1020" font-family="Arial, sans-serif" font-size="12" fill="#333">Rosetta compute</text>
  
  <rect x="500" y="1010" width="16" height="12" rx="2" fill="#FBEAF0" stroke="#993556"/>
  <text x="524" y="1020" font-family="Arial, sans-serif" font-size="12" fill="#333">QC/Metrics</text>
  
  <rect x="50" y="1040" width="16" height="12" rx="2" fill="#FAECE7" stroke="#993C1D"/>
  <text x="74" y="1050" font-family="Arial, sans-serif" font-size="12" fill="#333">Specialized relax</text>
  
  <rect x="190" y="1040" width="16" height="12" rx="2" fill="#EAF3DE" stroke="#3B6D11"/>
  <text x="214" y="1050" font-family="Arial, sans-serif" font-size="12" fill="#333">Structural analysis</text>
  
  <rect x="330" y="1040" width="16" height="12" rx="2" fill="#E6F3FF" stroke="#0066cc"/>
  <text x="354" y="1050" font-family="Arial, sans-serif" font-size="12" fill="#333">Data upload</text>
  
  <rect x="430" y="1040" width="16" height="12" rx="2" fill="#fff3e0" stroke="#f57c00"/>
  <text x="454" y="1050" font-family="Arial, sans-serif" font-size="12" fill="#333">Logging</text>
  
  <line x1="520" y1="1025" x2="540" y2="1025" stroke="#333" stroke-width="2" marker-end="url(#arrow)"/>
  <text x="550" y="1030" font-family="Arial, sans-serif" font-size="12" fill="#333">Required flow</text>
  
  <line x1="600" y1="1025" x2="620" y2="1025" stroke="#999" stroke-width="1.5" stroke-dasharray="3,3" marker-end="url(#arrow)"/>
  <text x="630" y="1030" font-family="Arial, sans-serif" font-size="12" fill="#333">Optional step</text>

  <!-- Database tables annotation -->
  <text x="50" y="1080" font-family="Arial, sans-serif" font-size="14" font-weight="bold" fill="#333">Database tables updated</text>
  <text x="50" y="1100" font-family="Arial, sans-serif" font-size="12" fill="#333">• modpro.QC (quality metrics, energy scores)</text>
  <text x="50" y="1115" font-family="Arial, sans-serif" font-size="12" fill="#333">• modpro.epitope (centroid, RBS residues, glycosylation)</text>
  <text x="50" y="1130" font-family="Arial, sans-serif" font-size="12" fill="#333">• modpro.logging (step-by-step pipeline status)</text>
</svg>a_ha1_pipeline_dag (1).svg…]()

URL for Nifi pipelines:
https://cdp-util-02.biotech.cdc.gov:8443/nifi/?processGroupId=416d35ee-136d-164d-8b00-8c9645776218&componentIds=

## Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the[Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md) and [Code of Conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md). For more information about CDC's privacy policy, please visit <http://www.cdc.gov/other/privacy.html>.

## Records Management Standard Notice

This repository is not a source of government records, it is to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov/).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).

## Attributions, License, and Limitations

## Limitations

This code is intended for external users and relies on code developed specifically for this purpose or from open source or publicly available code (getcontacts, Rosetta). This program is intended to create data that can be used for creating machine learning datasets, and is not intended to predict antigenic drift on its own or to any specific degree of confidence. Any downstream use of this data is the responsibility of the user and CDC and its affiliates take no responsibility for purposes beyond what is contained herein.

## License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later. This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version. This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details.

## Contributing Standard Notice

Anyone is encouraged to contribute to the repository by forking and submitting a pull request. (If you are new to GitHub, you might start with a basic tutorial.) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the Apache Software License v2 or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at <http://www.cdc.gov/other/privacy.html>.

## Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). All contributions to this repository will be released under the CC0 dedication. By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

## Points of contact

Maintainer Nicole Paterson (CDC) @patersonnicole <qxa4@cdc.gov>
Developer Nicholas Kovacs (CDC) @nicholasakovacs <nattila.kovacs@gmail.com>
Developer Brian Mann (CDC) @mannbrian <ytw4@cdc.gov>
