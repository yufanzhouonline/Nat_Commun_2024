import pandas as pd
from cluster import *
from plot import *
from scintegrator import *
import sys, os
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
### read configuration file
mudi = MUDI('config.inh.json')

cell = Integrator(mudi.config)
### instance of class Cluster
cluster = Cluster(mudi.config)
cluster.read_contact()
cluster.call_pca()
cluster.CELL_LIST = cluster.CELL['cell_id'].values.tolist()

pca_df = pd.DataFrame({
    'PC1':cluster.RAW_PCA[1][:,0],
    'PC2':cluster.RAW_PCA[1][:,1],
    'Clust_Label':cluster.CELL['sample_type']
})

# visualize
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Clust_Label', linewidth=0, s=2)
outdir = cluster.config['output_dir']
plt.savefig(f'{outdir}/scHiC_MajorType.inh.png',dpi=300)
plt.close()

cell = Integrator(mudi.config)
cell.cell_label()
# replace cell.cell_cluster()
cell.PDCLUSTER = pd.DataFrame({
    'sample':cluster.CELL['cell_id'],
    'cell':cluster.CELL['sample_type'],
    'PC1':cluster.RAW_PCA[1][:,0],
    'PC2':cluster.RAW_PCA[1][:,1],
    'PC3':cluster.RAW_PCA[1][:,2],
    'path':cluster.CELL['path']
})
cell2cluster = dict(zip(cell.PDCLUSTER['cell'].unique(),np.arange(len(cell.PDCLUSTER['cell'].unique()))))
cell.PDCLUSTER['cluster'] = cell.PDCLUSTER['cell'].map(cell2cluster)
cell.PDCLUSTER['label'] = ['C' + str(i + 1) for i in cell.PDCLUSTER.cluster]
cell.CLUSTER_LIST = []
cell.CLUSTER_NUM = len(set(cell.PDCLUSTER.cluster))
for i in range(cell.CLUSTER_NUM):
    cell.CLUSTER_LIST.append(list(cell.PDCLUSTER[cell.PDCLUSTER.cluster == i]['sample']))

cell.common_binnum(sample_name = cell.CELL, cluster_list = cell.CLUSTER_LIST[0], binvalue_cutoff = 10)
cell.common_binset(sample_name = cell.CELL, cluster_list = cell.CLUSTER_LIST[0], binvalue_cutoff = 0)
cell.cads_score(sample_name = cell.CELL, cluster_list = cell.CLUSTER_LIST[0], binvalue_cutoff = 0)
cell.BIN_LIST
cell.BIN_SCORE

cell.cluster_gene_list(cutoff_list = [0] * cell.CLUSTER_NUM)

with open('cell.inh.pkl', 'rb') as f:
    cell = pickle.load(f)
cell2cluster = dict(zip(cell.PDCLUSTER['cell'].unique(),np.arange(len(cell.PDCLUSTER['cell'].unique()))))
allpd_list = cell.integration_score()
allpd_tmp = []
for i in range(len(cell.CLUSTER_LIST)):
    for j in range(len(cell.DEG.CLUSTER_NUMBER)):
        cur_df = allpd_list[i][j]
        cur_df['scHiC_clust'] = i
        cur_df['scRNA_clust'] = j
        allpd_tmp.append(cur_df)
allpd_df = pd.concat(allpd_tmp)

# assigned G group
group_df = pd.DataFrame(index=allpd_df['scHiC_clust'].unique(),columns=allpd_df['scRNA_clust'].unique())
group_mapping = {}
group_id = 1
for i in allpd_df['scHiC_clust'].unique():
    for j in allpd_df['scRNA_clust'].unique():
        group_mapping[(i,j)] = f'G{group_id}'
        group_df.loc[i,j] = f'G{group_id}'
        group_id += 1
allpd_df['Group'] = allpd_df.apply(lambda row: group_mapping[(row['scHiC_clust'], row['scRNA_clust'])], axis=1)
scRNAID2Type = { 0:'ID2', 1:'PV', 2:'SST', 3:'VIP'}
scHiCID2Type = {cell2cluster[type]:type for type in cell2cluster}
allpd_df['scRNA_clust_annot'] = allpd_df['scRNA_clust'].map(scRNAID2Type)
allpd_df['scHiC_clust_annot'] = allpd_df['scHiC_clust'].map(scHiCID2Type)
with pd.ExcelWriter("allpd_df.inh.xlsx") as writer:
    for cur_type,cur_df in allpd_df.groupby('Group'):
        cur_df.to_excel(writer, sheet_name=cur_type, index=False)



with open('cell.inh.pkl', 'wb') as outp:
    pickle.dump(cell, outp, pickle.HIGHEST_PROTOCOL)

# number
twocluster_table = pd.crosstab(allpd_df['scHiC_clust_annot'],allpd_df['scRNA_clust_annot'])
plt.figure(figsize=(10, 7))
sns.heatmap(twocluster_table, annot=True, fmt='d', cmap='Greens',linewidths=1)
plt.xlabel('scRNA_clust')
plt.ylabel('scHiC_clust')
plt.xticks(rotation=60)
plt.tight_layout()
plt.savefig('/data/kfang/Yufan/inter_genes_table.inh.png',dpi=300)
plt.close()
# score
score_table_3col = allpd_df[['scHiC_clust_annot','scRNA_clust_annot','integration']].groupby(['scHiC_clust_annot','scRNA_clust_annot']).agg({'integration':'mean'}).reset_index()
score_table = score_table_3col.pivot(index='scHiC_clust_annot',columns= 'scRNA_clust_annot',values='integration')
group_df = pd.DataFrame(group_df.values,index=twocluster_table.index.values,columns= twocluster_table.columns)
annot_df = group_df.applymap(str) + '(' + twocluster_table.applymap(str) + ')'
plt.figure(figsize=(10, 7))
sns.heatmap(score_table, annot=annot_df, fmt='s',cmap='Greens',linewidths=1)
plt.xlabel('scRNA_clust')
plt.ylabel('scHiC_clust')
plt.xticks(rotation=60)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig('/data/kfang/Yufan/inter_score_table.inh.png',dpi=300)
plt.close()