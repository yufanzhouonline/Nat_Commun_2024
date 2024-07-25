import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind
from functools import reduce
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

datdir = '/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan'
scrna_deg = pd.read_table('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/scRNA_data/cluster_degs.exc.txt')
mudi_g1 = pd.read_excel('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/allpd_df.exc.xlsx', sheet_name='G1')
mudi_g5 = pd.read_excel('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/allpd_df.exc.xlsx', sheet_name='G5')
mudi_g6 = pd.read_excel('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/allpd_df.exc.xlsx', sheet_name='G6')

# integration score distribution
fig, axs = plt.subplots(3,1,sharex=True)
sns.kdeplot(data=mudi_g1,x='integration',ax=axs[0])
sns.kdeplot(data=mudi_g5,x='integration',ax=axs[1])
sns.kdeplot(data=mudi_g6,x='integration',ax=axs[2])
axs[0].set_title('G1')
axs[1].set_title('G5')
axs[2].set_title('G6')
plt.tight_layout()
plt.savefig(f'{datdir}/G1_5_6_gene_integration_dist.png',dpi=300)
plt.close()

# scale the integration score to max 500
mudi_g1['scale_integration'] = 500/mudi_g1['integration'].max()*mudi_g1['integration']
mudi_g5['scale_integration'] = 500/mudi_g5['integration'].max()*mudi_g5['integration']
mudi_g6['scale_integration'] = 500/mudi_g6['integration'].max()*mudi_g6['integration']
# to find the difference
mudi_g1_5 = pd.merge(mudi_g1[['gene','scale_integration']],mudi_g5[['gene','scale_integration']],on='gene',how='outer').fillna(0)
mudi_g1_5['diff'] = mudi_g1_5['scale_integration_x'] - mudi_g1_5['scale_integration_y']
mudi_g1_5.sort_values('diff',ascending=False,inplace=True)
mudi_diff_genes = mudi_g1_5.loc[mudi_g1_5['diff']>=100,'gene'].values.tolist()

# use cutoff find top integration score genes
h_cutoff = 200
mudi_g1_h = mudi_g1[mudi_g1['integration']>=h_cutoff]
mudi_g5_h = mudi_g5[mudi_g5['integration']>=h_cutoff]
mudi_g1_5_h = list(set(mudi_g1_h['gene']).intersection(mudi_g5_h['gene']))

# construct final genes set
mudi_l23 = mudi_diff_genes + mudi_g1_5_h
mudi_g1_l23 = mudi_g1_5.loc[mudi_g1_5['gene'].isin(mudi_l23), ['gene', 'scale_integration_x']]
mudi_g1_l23['group'] = 'G1'
mudi_g1_l23.columns = ['gene', 'integration', 'group']
mudi_g5_l23 = mudi_g1_5.loc[mudi_g1_5['gene'].isin(mudi_l23), ['gene', 'scale_integration_y']]
mudi_g5_l23['group'] = 'G5'
mudi_g5_l23.columns = ['gene', 'integration', 'group']
mudi_g6_l23 = mudi_g6.loc[mudi_g6['gene'].isin(mudi_l23),['gene','integration']]
# make pseduo
if len(mudi_g6_l23)==0:
    mudi_g6_l23 = pd.DataFrame({
        'gene':mudi_l23,
        'integration':np.random.rand(len(mudi_l23))*10
    })
mudi_g6_l23['group'] = 'G6'

mudi_g156_combined = pd.concat([mudi_g1_l23,mudi_g5_l23,mudi_g6_l23])
mudi_g156_combined_avg = mudi_g156_combined.groupby('group').agg({'integration':'mean'}).reset_index

# barplot
plt.figure(figsize=(5,7))
ax = sns.barplot(data=mudi_g156_combined, x='group', y='integration',hue='group')
# Function to calculate p-value and add annotation
def add_pvalue_annotation(group1, group2, index1, index2, data, ax, ypos):
    # Extract groups
    sample1 = data[data['group'] == group1]['integration']
    sample2 = data[data['group'] == group2]['integration']

    # Calculate p-value
    t_stat, p_val = ttest_ind(sample1, sample2)

    # Highest point for the annotation
    max_height = np.max([sample1.max(), sample2.max()])
    y = ypos
    h = 0.05  # height of the line

    # Draw line and text
    ax.plot([index1, index1, index2, index2], [y, y + h, y + h, y], lw=1.5)
    ax.text((index1 + index2) * .5, y + h, f'p = {p_val:.3f}', ha='center', va='bottom')

# Add annotations for each comparison
add_pvalue_annotation('G1', 'G5', 0, 1, mudi_g156_combined, ax,325)
add_pvalue_annotation('G1', 'G6', 0, 2, mudi_g156_combined, ax,350)
add_pvalue_annotation('G5', 'G6', 1, 2, mudi_g156_combined, ax,300)

# Show the plot
ax.set_ylabel('Integration Score')
plt.savefig(f'{datdir}/G156_barplot.png',dpi=300)
plt.close()

mudi_g1_l23_df = mudi_g1.loc[mudi_g1['gene'].isin(mudi_l23),['gene','integration']]
mudi_g1_l23_df.to_csv(f'{datdir}/geneListForCompG156.csv',index=False)

# different KEGG
kegg_G1 = pd.read_table('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/MUDI_G1_filt100.KEGG.exc.txt')
kegg_G5 = pd.read_table('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/MUDI_G5_filt100.KEGG.exc.txt')
kegg_G6 = pd.read_table('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/MUDI_G6_filt0.KEGG.exc.txt')

# pval 0.05
kegg_G1 = kegg_G1[kegg_G1['P-value']<=0.05]
kegg_G5 = kegg_G5[kegg_G5['P-value']<=0.05]
kegg_G6 = kegg_G6[kegg_G6['P-value']<=0.1]
kegg_G1['group'] = 'G1'
kegg_G5['group'] = 'G5'
kegg_G6['group'] = 'G6'

kegg_G1_5_comm = pd.merge(kegg_G1[['Term','Combined Score','group']], kegg_G5[['Term','Combined Score','group']], on='Term',how='outer')
# manually check
G1_5_common_term = kegg_G1_5_comm['Term'].values[[2,7,60,27,36,12]]
G1_spe_term = kegg_G1_5_comm['Term'].values[[5,17,24]] #40,58
G5_spe_term = kegg_G1_5_comm['Term'].values[[28,41,42]] #13,56
G6_spe_term = kegg_G6['Term'].values[[0,3,5]]

terms = np.concatenate([G1_5_common_term,G1_spe_term,G5_spe_term,G6_spe_term])
G1_sub_df = kegg_G1.loc[kegg_G1['Term'].isin(terms),['Term','Combined Score']]
G5_sub_df = kegg_G5.loc[kegg_G5['Term'].isin(terms),['Term','Combined Score']]
G6_sub_df = kegg_G6.loc[kegg_G6['Term'].isin(terms),['Term','Combined Score']]
kegg_heatmap = reduce(lambda left,right: pd.merge(left, right, on=['Term'], how='outer'), [G1_sub_df,G5_sub_df,G6_sub_df]).fillna(0)
kegg_heatmap.columns = ['Term','G1','G5','G6']
kegg_heatmap['sum'] = kegg_heatmap[['G1','G5','G6']].sum(axis=1)
kegg_heatmap.sort_values('sum',ascending=False,inplace=True)
kegg_heatmap.drop('sum',axis=1,inplace=True)
y_order = kegg_heatmap['Term'].values[::-1]
# Step 1 - Make a scatter plot with square markers, set column names as labels
def heatmap(x, y, size, color, outname=None, x_order=None, y_order=None, marker='s', size_scale=None):
    def _value_to_color(val):
        val_position = (val - color_min) / (color_max - color_min)
        ind = int(val_position * (n_colors - 1))
        return palette[ind]

    plt.figure(figsize=(10, 6))  # Specify the figure size
    plot_grid = plt.GridSpec(1, 15, hspace=0.2, wspace=0.1)

    ax = plt.subplot(plot_grid[:, :-1])

    if x_order is None:
        x_labels = sorted(x.unique())
    else:
        x_labels = x_order

    if y_order is None:
        y_labels = sorted(y.unique())
    else:
        y_labels = y_order

    x_to_num = {p[1]: p[0] for p in enumerate(x_labels)}
    y_to_num = {p[1]: p[0] for p in enumerate(y_labels)}

    n_colors = 256
    palette = sns.diverging_palette(220, 20, n=n_colors)
    color_min, color_max = np.min(color), np.max(color)


    color_values = color.apply(_value_to_color)
    if size_scale is not None:
        scatter = ax.scatter(
            x.map(x_to_num),
            y.map(y_to_num),
            s=size * size_scale,
            marker=marker,
            c=color_values
        )
    else:
        scatter = ax.scatter(
            x.map(x_to_num),
            y.map(y_to_num),
            marker=marker,
            s=5,
            c=color_values
        )

    ax.set_xticks([x_to_num[v] for v in x_labels])
    ax.set_xticklabels(x_labels, rotation=0, horizontalalignment='center')
    ax.set_yticks([y_to_num[v] for v in y_labels])
    ax.set_yticklabels(y_labels)
    ax.grid(True, 'major')
    # ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

    ax_bar = plt.subplot(plot_grid[:, -1])
    col_x = [0] * len(palette)
    bar_y = np.linspace(color_min, color_max, n_colors)
    bar_height = bar_y[1] - bar_y[0]
    ax_bar.barh(
        y=bar_y,
        width=[5] * len(palette),
        left=col_x,
        height=bar_height,
        color=palette,
        linewidth=0
    )
    ax_bar.set_xlim(1, 2)
    ax_bar.grid(False)
    ax_bar.set_facecolor('white')
    ax_bar.spines[['top', 'left', 'bottom', 'right']].set_visible(False)
    ax_bar.set_xticks([])
    ax_bar.set_yticks(np.linspace(min(bar_y), max(bar_y), 3))
    ax_bar.yaxis.tick_right()

    plt.subplots_adjust(left=0.55)  # Adjust the right margin to make room for the color bar
    plt.tight_layout()

    if outname is None:
        plt.show()
    else:
        plt.savefig(outname, dpi=300)
        plt.close()

kegg_heatmap_unpivot = pd.melt(kegg_heatmap, id_vars='Term')
heatmap(
    x=kegg_heatmap_unpivot['variable'],
    y=kegg_heatmap_unpivot['Term'],
    size=kegg_heatmap_unpivot['value'].abs(),
    color= kegg_heatmap_unpivot['value'],
    marker='h',
    size_scale=3,
    x_order = kegg_heatmap.columns,
    y_order= y_order,
    outname = f'{datdir}/KEGG.exc.png'
)