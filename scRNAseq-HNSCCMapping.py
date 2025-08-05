# Step 1: Load what you need
import pandas as pd
import scanpy as sc

# Step 2: Open the full dataset and separate it into metadata and gene expression
# First 5 rows are metadata (like cell type info), everything else is gene expression
file_path = "Desktop/new/HN.csv"
df_full = pd.read_csv(file_path, index_col=0)

metadata = df_full.iloc[:5].T
expression = df_full.iloc[5:].T

# Step 3: Fix the cancer cell label column to be numeric for analysis
metadata['classified  as cancer cell'] = pd.to_numeric(metadata['classified  as cancer cell'], errors='coerce')

# Step 4: Create a Scanpy object with expression and attach metadata
adata = sc.AnnData(expression)
adata.obs = metadata

# Step 5: Convert cancer cell classification to readable labels
# We'll make a new column that says "cancer" or "non-cancer" based on values in the original column
def map_cancer_status(value):
    if str(value).startswith('1.00000'):
        return 'cancer'
    elif str(value).startswith('0'):
        return 'non-cancer'
    else:
        return 'unknown'

adata.obs['cancer_status'] = adata.obs['classified  as cancer cell'].apply(map_cancer_status).astype('category')

# Step 6: Normalize and log-transform the data (basic cleanup for single-cell analysis)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Step 7: Pick the most informative genes to focus on
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']]

# Step 8: Scale the data and reduce dimensions with PCA
sc.pp.scale(adata)
sc.tl.pca(adata)

# Step 9: Build a neighborhood graph (used by UMAP to learn structure)
sc.pp.neighbors(adata)

# Step 10: Run UMAP to visualize the cells in 2D space
sc.tl.umap(adata)
sc.pl.umap(adata, color='cancer_status')

# Step 11: Run differential expression to find marker genes for cancer vs non-cancer
sc.tl.rank_genes_groups(adata, groupby='cancer_status', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
