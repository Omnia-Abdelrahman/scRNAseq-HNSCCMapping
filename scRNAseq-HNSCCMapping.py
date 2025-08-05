# Step 1: Load what you need
import pandas as pd
import scanpy as sc

# Step 1: Load the full dataset
df_full = pd.read_csv("data/GSE103322_HNSCC_all_data.txt", sep="\t", index_col=0, low_memory=False)

# Step 2: Split metadata rows (first 5 rows) from gene expression
metadata = df_full.iloc[:5].T       # Transpose so each cell gets a row of metadata
expression = df_full.iloc[5:].T     # Transpose gene expression to cells × genes

# Step 3: Build AnnData with obs = metadata
adata = sc.AnnData(expression)
adata.obs = metadata

# Optional: check what metadata is now available
print(adata.obs.columns)

#just exploring data
adata.obs[['classified  as cancer cell']].head(10)
pd.set_option("display.max_columns", None)
df_full.iloc[:10, :10]

# Step 4: Fix the cancer cell label column to be numeric for analysis
metadata = df_full.iloc[:5].T
expression = df_full.iloc[5:].T
metadata['classified  as cancer cell'] = pd.to_numeric(metadata['classified  as cancer cell'], errors='coerce')

# Step 5: Create a Scanpy object with expression and attach metadata
adata = sc.AnnData(expression)
adata.obs = metadata

# Step 6: Convert cancer cell classification to readable labels
# We'll make a new column that says "cancer" or "non-cancer" based on values in the original column
def map_cancer_status(value):
    if str(value).startswith('1'):
        return 'cancer'
    elif str(value).startswith('0'):
        return 'non-cancer'
    else:
        return 'unknown'

adata.obs['cancer_status'] = adata.obs['classified  as cancer cell'].apply(map_cancer_status).astype('category')

#finding what we got
adata.obs['cancer_status'].value_counts()

#shorter 
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_pcs=50)
sc.tl.umap(adata)
sc.pl.umap(adata, color='cancer_status')

expression = df_full.iloc[5:].T.apply(pd.to_numeric, errors='coerce')
expression = expression.fillna(0).astype(float)

adata = sc.AnnData(expression)
adata.obs = metadata

def map_cancer_status(value):
    if str(value).startswith('1.0'):
        return 'cancer'
    elif str(value).startswith('0'):
        return 'non-cancer'
    else:
        return 'unknown'

adata.obs['cancer_status'] = adata.obs['classified  as cancer cell'].apply(map_cancer_status)
adata.obs['cancer_status'] = adata.obs['cancer_status'].astype('category')

sc.tl.rank_genes_groups(adata, groupby='cancer_status', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
########################

# OR  Normalize and log-transform the data (basic cleanup for single-cell analysis)
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
