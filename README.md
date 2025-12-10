
# PBMC‑Derived Regulons Reveal Stress‑Defined Niches in Human Lymph Node 
### explore via 🔗[shiny app](https://yossohani.shinyapps.io/TF_LymphMap/)

<img width="1906" height="707" alt="image" src="https://github.com/user-attachments/assets/380f1476-c934-441c-b3ca-a8984fcb9e30" />


----

> **description:**  
> Integration of PBMC single‑cell RNA‑seq, pySCENIC regulon activity, and 10x Visium spatial transcriptomics to map activation and metabolic‑stress programs across immune niches in a human lymph node. plus a Shiny app for interactive exploration.

---

## 1. Project overview

This repository contains code  for the analysis and visualizations.


**Goal.**  
To determine whether a minimal set of transcription‑factor regulons, inferred from PBMC scRNA‑seq, can:

1. Reconstruct major immune compartments (B cells, memory CD4 T cells, dendritic cells) in a human lymph node, and  
2. Resolve them into distinct **regulatory niches** characterized by activation stress (AP‑1), oxidative/metabolic stress (Nrf2) and specialized B‑cell programs.

---

## 2. Data

### 2.1 PBMC single‑cell RNA‑seq

- Source: *[10X GENOMICS]*  
- Modalities: gene expression (UMI counts)  
- Preprocessing: QC, normalization, clustering, and cell‑type annotation performed with **Seurat**.


### 2.2 Human lymph node Visium dataset

- Platform: 10x Genomics Visium  
- Sample: human lymph node section  
- Source: *[10x genomics]*


---

## 3. Methods summary

1. **Regulon inference (PBMC)**
   - Tool: **pySCENIC**  
   - Input: PBMC scRNA‑seq  
   - Output: regulon activity matrix (cells × TFs),
     from which five immunologically relevant regulons were selected:
     **ATF3, FOSB, MAFF, NFE2L2 (Nrf2), IRF4**.

2. **PBMC → lymph node label transfer**
   - Tool: **Seurat v5** anchor‑based transfer (`FindTransferAnchors`, `TransferData`).  
   - Result: per‑spot probabilities for PBMC cell types; dominant label stored in `predicted.id`.

3. **Mapping regulon activity to Visium spots**
   - Regulon scores projected to Visium spots (e.g. via shared gene expression space ) and added as assay `"TF"`.

4. **TF‑based clustering of spots**
   - Features: regulon scores for ATF3, FOSB, MAFF, NFE2L2, IRF4.  
   - Workflow: `ScaleData` → `RunPCA` → `FindNeighbors` → `FindClusters`.  
   - Result: 7 TF clusters (`TF_cluster` ∈ 0–7) that segregate into DC‑, Memory CD4 T‑ and B‑cell–dominated niches.

5. **Differential expression and pathway analysis**
   - Differential genes between TF clusters with `FindMarkers`.  
   - Hallmark pathway enrichment with **fgsea** using MSigDB Hallmark sets.  
   - Key finding: B‑cell cluster 5 is enriched for **G2M, glycolysis, IL2–STAT5, IFN‑α**, while B‑cell cluster 4 is enriched for **CD74, HLA‑DRA, MS4A1, CD37, CXCL13**.

6. **Correlation analysis**
   - Spearman correlations of regulon activity across spots show:
     - Strong AP‑1 module: ATF3–FOSB (ρ ≈ 0.9).  
     - Strong antagonism: IRF4 vs NFE2L2 (ρ ≈ –0.75).  
     - Distinct MAFF‑biased stress program.

---

## 4. Repository structure

```text
.
├── scripts/
│ ├── Spatial_scripts ├── spatial_base.R
│                     ├── cor.R 
|                     ├── clust.R
   ├── Single_scripts ├── pyscenic
                      ├── base.R
                      ├── regulons.R
│                     ├── transpose.py
│ 
├── shiny_app/
│ ├── main_2.R.R # Shiny app ( ui.R + server.R)
│ └── pics # static images, css
└── README.md
