# scRNA-seq Analysis: PBMC 3k

End-to-end single-cell RNA sequencing analysis of 3,000 human peripheral 
blood mononuclear cells (PBMCs) using Python and Scanpy, with an interactive 
Streamlit dashboard for exploration.

## What this covers
- Quality control and filtering (dead cells, doublets, empty droplets)
- Normalization and log-transformation
- Feature selection (1,838 highly variable genes out of 32,738)
- PCA dimensionality reduction (50 components)
- Neighborhood graph construction (k=15 nearest neighbors)
- UMAP visualization
- Leiden clustering (6 distinct cell populations)
- Marker gene identification (Wilcoxon rank-sum test)
- Cell type annotation using published marker gene databases

## Cell populations identified
| Cluster | Cell Type | Key Marker Gene |
|---------|-----------|-----------------|
| 0 | CD4+ T cells | LDHB, CD3D |
| 1 | CD14+ Monocytes | FTL, CD14 |
| 2 | NK cells | NKG7, GZMA |
| 3 | B cells | CD74, CD79A |
| 4 | Dendritic cells | HLA-DPA1 |
| 5 | Platelets | PF4 |

## Dataset
PBMC 3k from 10x Genomics — the standard benchmark dataset for 
scRNA-seq methods. 2,700 peripheral blood mononuclear cells from 
a healthy donor.

## Stack
Python · Scanpy · AnnData · Matplotlib · Seaborn · Plotly · Streamlit

## Setup

```bash
# Clone the repo
git clone https://github.com/varghesebibin/scrna-pbmc3k.git
cd scrna-pbmc3k

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

## Run the analysis
Open Jupyter and run the notebook:
```bash
jupyter notebook notebooks/01_pbmc3k_analysis.ipynb
```

## Run the dashboard
```bash
streamlit run app.py
```

Then open your browser at `http://localhost:8501`

## Dashboard features
- Interactive UMAP colored by cell type
- Gene expression overlay — select any gene to see where it is expressed
- Cell type proportion bar chart
- Summary statistics

## ML techniques used
- Dimensionality reduction (PCA, UMAP)
- Feature selection (highly variable genes)
- Unsupervised clustering (Leiden algorithm)
- Statistical testing (Wilcoxon rank-sum)
- Graph-based methods (k-nearest neighbor graph)

## Limitations
- Cell type annotation based on top marker genes — automated tools 
  like CellTypist would provide more robust annotation
- Single timepoint from one healthy donor
- Resolution parameter (0.5) chosen as starting point — different 
  resolutions may reveal sub-populations

## Next steps
- Automated cell type annotation using CellTypist
- Sub-clustering to identify CD8+ T cells separately from CD4+ T cells
- Extension to spatial transcriptomics to map cell locations in tissue
