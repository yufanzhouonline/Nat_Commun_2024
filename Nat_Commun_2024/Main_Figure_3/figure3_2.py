import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.chdir('/data/yufan/schic/insulation/unnormalized')
os.getcwd()

res = 1000000
allchrno = ["chr" + str(i+1) for i in range(22)] + ['chrX']

samplename = pd.read_csv('/data/yufan/schic/schicnames.txt', header=0, sep="\t")
###Change the batch 3 to batch 6
samplename['Round'] = [6 if i==3 else i for i in samplename.Round]

###countcutoff must be > 1, otherwise cpu.hicluster_cpu will be error
countcutoff = 2

#MCF7
groupname = samplename[samplename.Sample=='MCF7'].copy()
groupname.reset_index(drop = True, inplace = True)

pdcluster = pd.read_csv('/data/yufan/schic/cluster/clusterlist001.txt',sep='\t',header=0)

np.array(pdcluster.loc[:,['PC1','PC2','PC3']])

###make plot
pc1 = 0
pc2 = 1
pc3 = 2

cellnum1 = pdcluster[pdcluster.label=='MCF7'].shape[0]
cellnum2 = pdcluster[pdcluster.label=='MCF7-T1M'].shape[0]
cellnum3 = pdcluster[pdcluster.label=='MCF7-TamR'].shape[0]
#cellnum4 = len(group4)

#cellnum4 = len(mcf7tbatch2)
print(cellnum1)
print(cellnum2)
print(cellnum3)
#print(cellnum4)

rawarray = np.array(pdcluster.loc[:,['PC1', 'PC2', 'PC3']])

###Figure PC1 and PC2
plt.figure()
p1, = plt.plot(rawarray[0:cellnum1,pc1],rawarray[0:cellnum1,pc2],'g.')
p2, = plt.plot(rawarray[cellnum1:(cellnum1 + cellnum2),pc1],rawarray[cellnum1:(cellnum1 + cellnum2),pc2],'b.')
p3, = plt.plot(rawarray[(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3),pc1],rawarray[(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3),pc2],'r.')

plt.legend([p1, p2, p3], ['MCF7', 'MCF7M1', 'MCF7TR'], loc='lower right')

plt.xlabel('PC1', fontsize=18)
plt.ylabel('PC2', fontsize=18)
plt.show()

###Figure PC2 and PC3
plt.figure()
p1, = plt.plot(rawarray[0:cellnum1, pc2], rawarray[0:cellnum1, pc3],'g.')
p2, = plt.plot(rawarray[cellnum1:(cellnum1 + cellnum2), pc2], rawarray[cellnum1:(cellnum1 + cellnum2), pc3],'b.')
p3, = plt.plot(rawarray[(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc2], rawarray[(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc3],'r.')
plt.legend([p1, p2, p3], ['MCF7', 'MCF7M1', 'MCF7TR'], loc='upper left')
plt.xlabel('PC2', fontsize=18)
plt.ylabel('PC3', fontsize=18)
plt.show()

###3D figure with PC1, PC2 and PC3
plt.figure()
ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
####separate the data point to three parts, marked them with different color
p1 = ax.scatter(rawarray[0:cellnum1, pc1], rawarray[0:cellnum1, pc2], rawarray[0:cellnum1, pc3], c='g')  ###draw the data points
p2 = ax.scatter(rawarray[cellnum1:(cellnum1 + cellnum2), pc1], rawarray[cellnum1:(cellnum1 + cellnum2), pc2], rawarray[cellnum1:(cellnum1 + cellnum2), pc3], c='b')  ###draw the data points
p3 = ax.scatter(rawarray[(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc1], rawarray[(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc2], rawarray[(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc3], c='r')  ###draw the data points

ax.set_zlabel('PC3')  ###axis
ax.set_ylabel('PC2')
ax.set_xlabel('PC1')
#ax.set_xticks(fontsize=10)
#ax.set_yticks(fontsize=10)
#ax.set_zticks(fontsize=10)
#Sometimes the default viewing angle is not optimal, 
#in which case we can use the view_init method to set the elevation and azimuthal angles. 
#In the following example, we'll use an elevation of 60 degrees (that is, 60 degrees above the x-y plane) 
#and an azimuth of 35 degrees (that is, rotated 35 degrees counter-clockwise about the z-axis):
#ax.view_init(60, 35)
ax.view_init(18, 161)
plt.legend([p1, p2, p3], ['MCF7', 'MCF7M1', 'MCF7TR'], loc='upper right')
plt.show()

pddf = pd.DataFrame(pdcluster.groupby(['cluster', 'label']).size())

pddf['cluster'] = ['C'+ str(pddf.index[i][0]+1) for i in range(pddf.shape[0])]
pddf['label'] = [pddf.index[i][1] for i in range(pddf.shape[0])]
pddf.columns = ['number', 'cluster', 'label']
pddf.reset_index(drop = True, inplace = True)
plt.barh(pddf.cluster, pddf.number, color='b')
plt.barh(pddf.cluster, pddf.number, color='r', left=pddf.number)
plt.show()
pickup = [0, 1, 2, 3, 5, 7, 10, 11, 13]
plt.barh(pddf[pddf.index.isin(pickup)].cluster, pddf[pddf.index.isin(pickup)].number, color='b')

###Make the stacked bar plot
pdplot = pd.DataFrame({'cluster':[], 'MCF7':[], 'MCF7M1':[], 'MCF7TR':[]})
pdplot = pdplot.astype({'cluster':int, 'MCF7':int, 'MCF7M1':int, 'MCF7TR':int})
for i in range(9):
	tempcount = []
	for j in ['MCF7', 'MCF7-T1M', 'MCF7-TamR']:
		pdcount = pdcluster[(pdcluster.cluster==i) & (pdcluster.label==j)].shape[0]
		print(i, j, pdcount)
		tempcount.append(pdcount)
	temppd = pd.DataFrame({'cluster':[i], 'MCF7':[tempcount[0]], 'MCF7M1':[tempcount[1]], 'MCF7TR':[tempcount[2]]})
	pdplot = pd.concat([pdplot, temppd],axis=0)

pdplot['label'] = ['C' + str(i+1) for i in pdplot.cluster]
pdplot.sort_values(by='cluster', ascending=False, inplace=True)

###make bar plot
p1 = plt.barh(pdplot.label, pdplot.MCF7, color='#00FF0088')
p2 = plt.barh(pdplot.label, pdplot.MCF7M1, color='#0000FF88', left=pdplot.MCF7)
p3 = plt.barh(pdplot.label, pdplot.MCF7TR, color='#FF000088', left=pdplot.MCF7 + pdplot.MCF7M1)
plt.legend([p1, p2, p3], ['MCF7', 'MCF7M1', 'MCF7TR'], loc='lower right')
plt.xlabel('# of single cells', fontsize=15)
plt.ylabel('Clusters', fontsize=15)
plt.show()

####clusters common bins
clusterc0 = pdcluster[pdcluster.cluster==0]['sample']
clusterc1 = pdcluster[pdcluster.cluster==1]['sample']
clusterc2 = pdcluster[pdcluster.cluster==2]['sample']
clusterc3 = pdcluster[pdcluster.cluster==3]['sample']
clusterc4 = pdcluster[pdcluster.cluster==4]['sample']
clusterc5 = pdcluster[pdcluster.cluster==5]['sample']
clusterc6 = pdcluster[pdcluster.cluster==6]['sample']
clusterc7 = pdcluster[pdcluster.cluster==7]['sample']
clusterc8 = pdcluster[pdcluster.cluster==8]['sample']


#common bins
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc5))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7c5) + list(mcf7t1mc5) + list(mcf7tamrc5))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc2))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7t1mc2) + list(mcf7tamrc2))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc7))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7c7) + list(mcf7tamrc7))].copy()
groupname = samplename[samplename.Name.isin(list(clusterc4))].copy()
groupname.reset_index(drop = True, inplace = True)
setlist = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	csvmat['twobin'] = csvmat.bin1 * 10000+ csvmat.bin2
	setlist.append(set(csvmat.twobin))

commonset = setlist[0]			###get the commonset of all cells to commonset
for i in range(1,len(setlist)):
	commonset = commonset.intersection(setlist[i])

###calculate the common set
commonset = sorted(list(commonset))
###add the bins of chr22:2846-2897
#commonset.extend([i*10000+i for i in range(2846,2898)])
print(len(commonset))

loci1 = pd.Series(commonset)//10000
loci2 = pd.Series(commonset)%10000
for i in loci1.index:
	print(i)
	if loci1[i] == loci2[i]:
		print('True')
	else:
		print('False')

binbed = pd.read_csv('/data/yufan/schic/batch6/hicpro2/output/hic_results/matrix/TA01/raw/1000000/TA01_1000000_abs.bed', header=None, sep='\t')
binbed.columns = ['chr', 'bin1', 'bin2', 'pos']
commonbin = binbed[binbed.pos.isin(loci1)].copy()

###RT06 chr14
###FOXA1: chr14:38,058,757-38,064,325
foxa1 = commonbin[commonbin.chr=='chr14'].copy()

tads = pd.read_csv('/data/yufan/schic/insulation/unnormalized/batch2/RP06/RP06_100000_chr14.is500001.ids200001.insulation.boundaries.bed', header=None, sep="\t", skiprows=1)
tads.columns = ['chr', 'startb', 'endb', 'name', 'value']

tads[(tads.startb>=37000000) & (tads.endb<39000000)]

hg19dim = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560]
chrlist = ['chr'+ str(i) for i in range(1,23)] + ['chrX']
chrsize = {chrlist[i]:hg19dim[i] for i in range(23)}
#{'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560}

tadsprefix = '/data/yufan/schic/insulation/unnormalized/batch'
sample = 'RP06'
chrno = 'chr14'
res = 100000

if str(res) == '40000':
    ispost = '520001'
    idspost = '240001'
elif str(res) == '200000':
    ispost = '600001'
    idspost = '400001'
elif str(res) == '500000':
    ispost = '1000001'
    idspost = '1000001'
elif str(res) == '1000000':
    ispost = '2000001'
    idspost = '2000001'
else:
    ispost = '500001'
    idspost = '200001'

def tadscount(clusterno):
	groupname = samplename[samplename.Name.isin(list(clusterno))].copy()
	cellsumlist = []
	celllenlist = []
	for i in groupname.index:
		tadscountlist = []
		tadslengthlist = []
		for j in chrlist:
			tadsfile = tadsprefix + str(groupname.Round[i]) + '/' + groupname.Name[i] + '/' + groupname.Name[i] + '_' + str(res) + '_' + chrno + '.is' + ispost + '.ids' + idspost + '.insulation.boundaries.bed'
			#print('Loading TADs file: ' + tadsfile)
			tadsbed = pd.read_csv(tadsfile, header=None, sep="\t", skiprows=1)
			tadsbed.columns = ['chr', 'startb', 'endb', 'name', 'value']
			tadscount = tadsbed.shape[0] + 1
			tadsmid = list(tadsbed.startb)
			tadsstart = [0] + tadsmid
			tadsend = tadsmid + [chrsize[chrno]]
			tadslength = np.array(tadsend)-np.array(tadsstart)
			tadscountlist.append(tadscount)
			tadslengthlist = tadslengthlist + list(tadslength)
		celltadssum = np.sum(tadscountlist)
		celltadslen = np.mean(tadslengthlist)
		print(groupname.Name[i], celltadssum, celltadslen)
		cellsumlist.append(celltadssum)
		celllenlist.append(celltadslen)
	return cellsumlist, celllenlist

cellsumc0, celllenc0 = tadscount(clusterc0)
cellsumc1, celllenc1 = tadscount(clusterc1)
cellsumc2, celllenc2 = tadscount(clusterc2)
cellsumc3, celllenc3 = tadscount(clusterc3)
cellsumc4, celllenc4 = tadscount(clusterc4)
cellsumc5, celllenc5 = tadscount(clusterc5)
cellsumc6, celllenc6 = tadscount(clusterc6)
cellsumc7, celllenc7 = tadscount(clusterc7)
cellsumc8, celllenc8 = tadscount(clusterc8)

###Make figures
all_data = [cellsumc0, cellsumc1, cellsumc2, cellsumc3, cellsumc4, cellsumc5, cellsumc6, cellsumc7, cellsumc8]
#all_data = [list(np.array(celllenc0)/1000), list(np.array(celllenc1)/1000), list(np.array(celllenc2)/1000), list(np.array(celllenc3)/1000), list(np.array(celllenc4)/1000), list(np.array(celllenc5)/1000), list(np.array(celllenc6)/1000), list(np.array(celllenc7)/1000), list(np.array(celllenc8)/1000)]
labels = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
medianprops = dict(linestyle='-', linewidth=1, color='red')

bplot = plt.boxplot(all_data, patch_artist=True, labels=labels, medianprops=medianprops, showfliers=False)  
#plt.title('Rectangular box plot')

#colors = ['darkred', 'darkorange', 'green', 'blue', 'c', 'm', 'y', 'lightblue', 'pink']
colors = ['#FF000088', '#00FF0088', '#0000FF88', '#44FFFF88', '#FF44FF88', '#FFFF4488', '#00888888', '#88008888', '#88880088']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)  

#for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
#	plt.setp(bplot[element], color=edge_color)


plt.xlabel('Clusters', fontsize=15)
plt.ylabel('# of TADs', fontsize=15)
#plt.ylabel('Size of TADs (Kb)', fontsize=15)
plt.show()

#Compute the Wilcoxon rank-sum statistic for two samples.
from scipy import stats

#alllen = celllenc0 + celllenc1 + celllenc2 + celllenc3 + celllenc4 + celllenc5 + celllenc6 + celllenc7 + celllenc8
alllen = celllenc0 + celllenc1 + celllenc2 + celllenc3 + celllenc4 + celllenc5 + celllenc6 + celllenc7

statvalue, pvalue = stats.ranksums(all_data[8], alllen)
print(pvalue)
