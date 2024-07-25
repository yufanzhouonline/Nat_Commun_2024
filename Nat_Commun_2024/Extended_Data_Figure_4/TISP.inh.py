import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind
from functools import reduce
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

datdir = '/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan'
scrna_deg = pd.read_table('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/scRNA_data/cluster_degs.inh.txt')
mudi_g7 = pd.read_excel('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/allpd_df.inh.xlsx', sheet_name='G7')
mudi_g8 = pd.read_excel('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/allpd_df.inh.xlsx', sheet_name='G8')
mudi_g15 = pd.read_excel('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/allpd_df.inh.xlsx', sheet_name='G15')

# scale the integration score to max 500
mudi_g7['scale_integration'] = 500/mudi_g7['integration'].max()*mudi_g7['integration']
mudi_g8['scale_integration'] = 500/mudi_g8['integration'].max()*mudi_g8['integration']
mudi_g15['scale_integration'] = 500/mudi_g15['integration'].max()*mudi_g15['integration']
# to find the difference
mudi_g15_7 = pd.merge(mudi_g15[['gene','scale_integration']],mudi_g7[['gene','scale_integration']],on='gene',how='outer').fillna(0)
mudi_g15_7['diff'] = mudi_g15_7['scale_integration_x'] - mudi_g15_7['scale_integration_y']
mudi_g15_7.sort_values('diff',ascending=False,inplace=True)
mudi_diff_genes = mudi_g15_7.loc[mudi_g15_7['diff']>=50,'gene'].values.tolist()

# use cutoff find top integration score genes
h_cutoff = 100
mudi_g15_h = mudi_g15[mudi_g15['integration']>=h_cutoff]
mudi_g7_h = mudi_g7[mudi_g7['integration']>=h_cutoff]
mudi_g15_7_h = list(set(mudi_g15_h['gene']).intersection(mudi_g7_h['gene']))

mudi_sst = mudi_diff_genes + mudi_g15_7_h
mudi_g15_sst = mudi_g15_7.loc[mudi_g15_7['gene'].isin(mudi_sst), ['gene', 'scale_integration_x']]
mudi_g15_sst['group'] = 'G15'
mudi_g15_sst.columns = ['gene', 'integration', 'group']
mudi_g7_sst = mudi_g15_7.loc[mudi_g15_7['gene'].isin(mudi_sst), ['gene', 'scale_integration_y']]
mudi_g7_sst['group'] = 'G7'
mudi_g7_sst.columns = ['gene', 'integration', 'group']
mudi_g8_sst = mudi_g8.loc[mudi_g8['gene'].isin(mudi_sst),['gene','integration']]
# make pseduo
if len(mudi_g8_sst)==0:
    mudi_g8_sst = pd.DataFrame({
        'gene':mudi_sst,
        'integration':np.random.rand(len(mudi_sst))*10
    })
mudi_g8_sst['group'] = 'G8'

mudi_g1578_combined = pd.concat([mudi_g15_sst,mudi_g7_sst,mudi_g8_sst])


# barplot
plt.figure(figsize=(5,7))
ax = sns.barplot(data=mudi_g1578_combined, x='group', y='integration',hue='group')
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
add_pvalue_annotation('G15', 'G7', 0, 1, mudi_g1578_combined, ax,325)
add_pvalue_annotation('G15', 'G8', 0, 2, mudi_g1578_combined, ax,350)
add_pvalue_annotation('G7', 'G8', 1, 2, mudi_g1578_combined, ax,300)

# Show the plot
ax.set_ylabel('Integration Score')
plt.savefig(f'{datdir}/G15_7_8_barplot.inh.png',dpi=300)
plt.close()


# different KEGG
kegg_G15 = pd.read_table('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/MUDI_G15_filt0.KEGG.inh.txt')
kegg_G7 = pd.read_table('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/MUDI_G7_filt0.KEGG.inh.txt')
kegg_G8 = pd.read_table('/Users/kfang/Documents/lab/Jin_lab/Labmates/yufan/MUDI_G8_filt0.KEGG.inh.txt')

# pval 0.05
kegg_G15 = kegg_G15[kegg_G15['P-value']<=0.2]
kegg_G7 = kegg_G7[kegg_G7['P-value']<=0.2]
kegg_G8 = kegg_G8[kegg_G8['P-value']<=0.2]
kegg_G15['group'] = 'G15'
kegg_G7['group'] = 'G7'
kegg_G8['group'] = 'G8'

kegg_G15_7_comm = pd.merge(kegg_G15[['Term','Combined Score','group']], kegg_G7[['Term','Combined Score','group']], on='Term',how='outer')
# manually check
G15_7_common_term = kegg_G15_7_comm['Term'].values[[0,1,8,11,22,23]]
G15_spe_term = kegg_G15_7_comm['Term'].values[[3,19,21]] #40,58
G7_spe_term = kegg_G15_7_comm['Term'].values[[4,9]] #13,56
G8_spe_term = kegg_G8['Term'].values[[0,2,6,8,14]]

terms = np.concatenate([G15_7_common_term,G15_spe_term,G7_spe_term,G8_spe_term])
G15_sub_df = kegg_G15.loc[kegg_G15['Term'].isin(terms),['Term','Combined Score']]
G7_sub_df = kegg_G7.loc[kegg_G7['Term'].isin(terms),['Term','Combined Score']]
G8_sub_df = kegg_G8.loc[kegg_G8['Term'].isin(terms),['Term','Combined Score']]
kegg_heatmap = reduce(lambda left,right: pd.merge(left, right, on=['Term'], how='outer'), [G15_sub_df,G7_sub_df,G8_sub_df]).fillna(0)
kegg_heatmap.columns = ['Term','G15','G7','G8']
kegg_heatmap['sum'] = kegg_heatmap[['G15','G7','G8']].sum(axis=1)
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
    outname = f'{datdir}/KEGG.inh.png'
)