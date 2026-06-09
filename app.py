import streamlit as st
import scanpy as sc
import plotly.express as px
import pandas as pd
import scipy.sparse as sp

# ── Page config ──
st.set_page_config(page_title="PBMC 3k Analysis", layout="wide")

# ── Load data ──
@st.cache_data
def load_data():
    adata = sc.read_h5ad("data/pbmc3k_processed.h5ad")
    return adata

adata = load_data()

# ── Title ──
st.title("scRNA-seq Analysis: PBMC 3k")
st.markdown("Interactive dashboard for single-cell RNA sequencing analysis of 2,638 human blood cells.")

# ── Sidebar ──
st.sidebar.header("Controls")

# Gene selector
all_genes = adata.var_names.tolist()
default_genes = ['CD3D', 'CD14', 'NKG7', 'MS4A1', 'HLA-DPA1', 'PF4']
selected_gene = st.sidebar.selectbox("Select a gene to visualize", default_genes + sorted(all_genes))

# ── Build UMAP dataframe ──
umap_df = pd.DataFrame(
    adata.obsm['X_umap'],
    columns=['UMAP1', 'UMAP2']
)
umap_df['Cell Type'] = adata.obs['cell_type'].values
umap_df['Leiden Cluster'] = adata.obs['leiden'].values

# Add selected gene expression
if selected_gene in adata.var_names:
    gene_idx = list(adata.var_names).index(selected_gene)
    if sp.issparse(adata.X):
        umap_df[selected_gene] = adata.X[:, gene_idx].toarray().flatten()
    else:
        umap_df[selected_gene] = adata.X[:, gene_idx].flatten()

# ── Layout: two columns ──
col1, col2 = st.columns(2)

# ── Left: UMAP by cell type ──
with col1:
    st.subheader("UMAP — Cell Types")
    fig1 = px.scatter(
        umap_df, x='UMAP1', y='UMAP2',
        color='Cell Type',
        hover_data=['Cell Type', 'Leiden Cluster'],
        width=600, height=500
    )
    fig1.update_traces(marker=dict(size=3))
    fig1.update_layout(legend=dict(orientation="h", yanchor="bottom", y=-0.3))
    st.plotly_chart(fig1, use_container_width=True)

# ── Right: UMAP by gene expression ──
with col2:
    st.subheader(f"UMAP — {selected_gene} Expression")
    if selected_gene in umap_df.columns:
        fig2 = px.scatter(
            umap_df, x='UMAP1', y='UMAP2',
            color=selected_gene,
            color_continuous_scale='Viridis',
            hover_data=['Cell Type'],
            width=600, height=500
        )
        fig2.update_traces(marker=dict(size=3))
        st.plotly_chart(fig2, use_container_width=True)

# ── Cell type proportions ──
st.subheader("Cell Type Proportions")
cell_counts = adata.obs['cell_type'].value_counts()
fig3 = px.bar(
    x=cell_counts.index,
    y=cell_counts.values,
    labels={'x': 'Cell Type', 'y': 'Number of Cells'},
    color=cell_counts.index
)
fig3.update_layout(showlegend=False, height=400)
st.plotly_chart(fig3, use_container_width=True)

# ── Summary stats ──
st.subheader("Dataset Summary")
col3, col4, col5, col6 = st.columns(4)
col3.metric("Total Cells", f"{adata.n_obs:,}")
col4.metric("Total Genes", f"{adata.n_vars:,}")
col5.metric("Cell Types", f"{adata.obs['cell_type'].nunique()}")
col6.metric("Clusters", f"{adata.obs['leiden'].nunique()}")