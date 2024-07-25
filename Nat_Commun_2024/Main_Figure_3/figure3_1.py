
import os
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from cluster import cpu
#from cluster import cputsne
from mpl_toolkits.mplot3d import Axes3D
import random

os.chdir('/data/yufan/schic/cluster/analysis06')
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
group1 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group1.append(groupname.Name[i])

#MCF7-T1M
groupname = samplename[samplename.Sample=='MCF7-T1M'].copy()
groupname.reset_index(drop = True, inplace = True)
group2 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group2.append(groupname.Name[i])

#MCF7-TamR
groupname = samplename[samplename.Sample=='MCF7-TamR'].copy()
groupname.reset_index(drop = True, inplace = True)
group3 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group3.append(groupname.Name[i])

###K562
k562list = ['102_K562-B', '110_K562-B', '123_K562-B', '141_K562-B',
	'156_K562-A', '215_K562-A', '217_K562-B', '220_K562-B', '221_K562-A', 
	'222_K562-B', '223_K562-B', '224_K562-B', '225_K562-A', '226_K562-B',
	'227_K562-A', '228_K562-B', '229_K562-A', '230_K562-A', '231_K562-A',
	'232_K562-A', '233_K562-A', '234_K562-A', '235_K562-A', '236_K562-A',
	'237_K562-A', '238_K562-A', '239_K562-A', '240_K562-A', '241_K562-A',
	'242_K562-A', '54_K562-B', '58_K562-B', '72_K562-B', '79_K562-B']
###	'K562_bulkA', 'K562_bulkB'
group4 = []
for i in k562list:
	runpath = "/data/yufan/epigenetics/venus/yufan/singlecells/k562/matrix/"
	rundata = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = i + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group4.append(i)

###Bing Ren WTC11C6
WTC11C6list = ['cell0' + str(i) if i < 10 else 'cell' + str(i) for i in range(1,95)]

group5 = []
for i in WTC11C6list:
	runpath = "/data/yufan/schic/public/ren4dn/hicpro/WTC11C6/output/hic_results/matrix/"
	rundata = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = 'WTC11C6_'+ i + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group5.append(i)

group5 = ['WTC11C6_' + i for i in group5]


###Bing Ren WTC11C6
WTC11C28list = ['cell0' + str(i) if i < 10 else 'cell' + str(i) for i in range(1,95)]

group6 = []
for i in WTC11C28list:
	runpath = "/data/yufan/schic/public/ren4dn/hicpro/WTC11C28/output/hic_results/matrix/"
	rundata = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = 'WTC11C28_'+ i + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group6.append(i)

group6 = ['WTC11C28_' + i for i in group6]

#network = group1 + group2 + group3
network = group1 + group2 + group3 + group4 + group5 + group6

label = network
hg19dim = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560]
chrom = [str(i+1) for i in range(22)] + ['X']
chromsize = {chrom[i]:hg19dim[i] for i in range(len(chrom))}
nc = 2
ndim = 20

###clusterlabel: cluster labels identified by KMeans, rawpca: the matrix transformed by PCA
#clusterlabel, rawpca = cpu.hicluster_cpu(network, chromsize, nc=nc, res=res, ncpus=5)
rawpca = cpu.hicluster_cpu(network, chromsize, nc=nc, res=res, ncpus=5)

###make plot
pc1 = 0
pc2 = 1
pc3 = 2

cellnum1 = len(group1)
cellnum2 = len(group2)
cellnum3 = len(group3)
cellnum4 = len(group4)
cellnum5 = len(group5)
cellnum6 = len(group6)

#cellnum4 = len(mcf7tbatch2)
print(cellnum1)
print(cellnum2)
print(cellnum3)
print(cellnum4)
print(cellnum5)
print(cellnum6)

plt.figure()
p1, = plt.plot(rawpca[1][0:cellnum1,pc1],rawpca[1][0:cellnum1,pc2],'g.')
p2, = plt.plot(rawpca[1][cellnum1:(cellnum1 + cellnum2),pc1],rawpca[1][cellnum1:(cellnum1 + cellnum2),pc2],'b.')
p3, = plt.plot(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3),pc1],rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3),pc2],'r.')
p4, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4),pc1],rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4),pc2],'y.')
p5, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5),pc1],rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5),pc2],'c.')
p6, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6),pc1],rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6),pc2],'k.')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper left')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower left')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower right')
plt.legend([p1, p2, p3, p4, p5, p6], ['MCF7', 'MCF7M1', 'MCF7TR', 'K562', 'WTC11C6', 'WTC11C28'], loc='upper center')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper center')
#plt.legend([p1, p2], ['MCF7', 'MCF7-TamR'], loc='upper right')
plt.xlabel('PC1', fontsize=18)
plt.ylabel('PC2', fontsize=18)
#plt.axis('square')
plt.show()

plt.figure()
p1, = plt.plot(rawpca[1][0:cellnum1, pc2],rawpca[1][0:cellnum1, pc3],'g.')
p2, = plt.plot(rawpca[1][cellnum1:(cellnum1 + cellnum2), pc2],rawpca[1][cellnum1:(cellnum1 + cellnum2), pc3],'b.')
p3, = plt.plot(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc2],rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc3],'r.')
p4, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc2],rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc3],'y.')
p5, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5),pc2],rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5),pc3],'c.')
p6, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6),pc2],rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6),pc3],'k.')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper left')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper center')
#plt.legend([p1, p2], ['MCF7', 'MCF7-TamR'], loc='upper right')
plt.legend([p1, p2, p3, p4, p5, p6], ['MCF7', 'MCF7M1', 'MCF7TR', 'K562', 'WTC11C6', 'WTC11C28'], loc='upper right')
plt.xlabel('PC2', fontsize=18)
plt.ylabel('PC3', fontsize=18)
#plt.axis('square')
plt.show()

###3D figure
plt.figure()
ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
####separate the data point to three parts, marked them with different color
p1 = ax.scatter(rawpca[1][0:cellnum1, pc1], rawpca[1][0:cellnum1, pc2], rawpca[1][0:cellnum1, pc3], c='g')  ###draw the data points
p2 = ax.scatter(rawpca[1][cellnum1:(cellnum1 + cellnum2), pc1], rawpca[1][cellnum1:(cellnum1 + cellnum2), pc2], rawpca[1][cellnum1:(cellnum1 + cellnum2), pc3], c='b')  ###draw the data points
p3 = ax.scatter(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc1], rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc2], rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc3], c='r')  ###draw the data points
p4 = ax.scatter(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc1], rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc2], rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc3], c='y')  ###draw the data points
p5 = ax.scatter(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5), pc1], rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5), pc2], rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5), pc3], c='c')  ###draw the data points
p6 = ax.scatter(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6), pc1], rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6), pc2], rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6), pc3], c='k')  ###draw the data points

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
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper left')
plt.legend([p1, p2, p3, p4, p5, p6], ['MCF7', 'MCF7M1', 'MCF7TR', 'K562', 'WTC11C6', 'WTC11C28'], loc='upper left', fontsize=8)

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper right')
plt.show()

###############################################
###Make the sub-cluster

from numpy import unique
from numpy import where
from sklearn.datasets import make_classification
from sklearn.mixture import GaussianMixture
from matplotlib import pyplot

# input of model
X = rawpca[1][:,0:3]

#from sklearn.cluster import KMeans
#model = KMeans(n_clusters=7)

import warnings
warnings.filterwarnings("ignore")

#Silhouette Coefficient
from sklearn.metrics import silhouette_score
scorelist = []
for i in range(20,1,-1):
	randomscore = 0
	for j in range(1, 1001, 1):
		#print('i: ' + str(i) + ', j: ' + str(j))
		kmeans = GaussianMixture(n_components=i,random_state=j).fit(X)
		randomscore = randomscore + silhouette_score(X, kmeans.predict(X))
	scorelist.append(randomscore/1000)
	print('The calinski_harabaz score for %d clustering is: %f'%(i, randomscore/1000))

plt.plot(range(20, 1, -1), scorelist, color='black')
plt.xlim(20, 1)
plt.xticks(range(20, 1, -1), range(20, 1, -1), fontsize=13)
plt.xlabel('Cluster #', fontsize=15)
plt.ylabel('Silhouette Coefficient', fontsize=15)
plt.show()

scorelist = []
for i in range(1, 1001, 1):
	kmeans = GaussianMixture(n_components=9,random_state=i).fit(X)
	score = silhouette_score(X, kmeans.predict(X))
	scorelist.append(score)
	print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))

maxseed = np.argmax(pd.Series(scorelist))

###4, 5, 6, 9, 10, 11, 13, 14 best 10
maxseed = 10
#Gaussian Mixture
model = GaussianMixture(n_components=9, random_state=maxseed)

#model = GaussianMixture(n_components=7, random_state=316)
#model = GaussianMixture(n_components=4, random_state=123)

################################

model.fit(X)

yhat = model.predict(X)

clusters = unique(yhat)

###PC1, PC2
for cluster in clusters:
	row_ix = where(yhat == cluster)
	pyplot.scatter(X[row_ix, 0], X[row_ix, 1])

pyplot.xlabel('PC1', fontsize=18)
pyplot.ylabel('PC2', fontsize=18)
pyplot.show()

###PC2, PC3
for cluster in clusters:
	row_ix = where(yhat == cluster)
	pyplot.scatter(X[row_ix, 1], X[row_ix, 2])

pyplot.xlabel('PC2', fontsize=18)
pyplot.ylabel('PC3', fontsize=18)
pyplot.show()

###3D figure
plt.figure()
ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
for cluster in clusters:
	row_ix = where(yhat == cluster)
	ax.scatter(X[row_ix, 0], X[row_ix, 1], X[row_ix, 2])


ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
ax.view_init(18, 161)
plt.show()

###Fowlkes-Mallows Score
labeltrue = [1] * cellnum1 + [2] * cellnum2 + [3] * cellnum3
from sklearn import metrics
scorelist = []
for i in range(20,1,-1):
	kmeans=GaussianMixture(n_components=i,random_state=316).fit(X)
#	score=metrics.adjusted_rand_score(labeltrue, kmeans.predict(X))
	score = metrics.fowlkes_mallows_score(labeltrue, kmeans.predict(X))
	scorelist.append(score)
	print('The adjusted rand score for %d clustering is: %f'%(i,score))

scorelist = []
for i in range(1, 1000, 1):
	kmeans = GaussianMixture(n_components=9,random_state=i).fit(X)
#	score = metrics.adjusted_rand_score(labeltrue, kmeans.predict(X))
	score = metrics.fowlkes_mallows_score(labeltrue, kmeans.predict(X))
	scorelist.append(score)
	print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))

maxseed = np.argmax(pd.Series(scorelist))
###maxseed: 181

###
clusterlabel = ['MCF7'] * cellnum1 + ['MCF7-T1M'] * cellnum2 + ['MCF7-TamR'] * cellnum3
pdcluster = pd.DataFrame({'sample': network, 'label': clusterlabel, 'cluster': yhat, 'PC1': X[:, 0], 'PC2': X[:, 1], 'PC3': X[:, 2]})

pdcluster.groupby(['cluster', 'label']).size()


mcf7c5 = pdcluster[(pdcluster.label=='MCF7') & (pdcluster.cluster==5)]['sample']
mcf7t1mc5 = pdcluster[(pdcluster.label=='MCF7-T1M') & (pdcluster.cluster==5)]['sample']
mcf7tamrc5 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==5)]['sample']

mcf7t1mc2 = pdcluster[(pdcluster.label=='MCF7-T1M') & (pdcluster.cluster==2)]['sample']
mcf7tamrc2 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==2)]['sample']

mcf7c7 = pdcluster[(pdcluster.label=='MCF7') & (pdcluster.cluster==7)]['sample']
mcf7tamrc7 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==7)]['sample']

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
groupname = samplename[samplename.Name.isin(list(clusterc7) + list(clusterc8))].copy()
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


###cluster name and label number
###SC1: 3, SC2: 6, SC3: 1, SC4: 0, SC5: 4, SC6: 5, SC7: 2
pdcluster['sclabel'] = pdcluster.apply(lambda row : \
	'sc1' if row['cluster']==3 else \
	'sc2' if row['cluster']==6 else \
	'sc3' if row['cluster']==1 else \
	'sc4' if row['cluster']==0 else \
	'sc5' if row['cluster']==4 else \
	'sc6' if row['cluster']==5 else \
	'sc7' if row['cluster']==2 else \
	'ERROR', axis=1)

pdcluster.groupby(['sclabel', 'label']).size()


pdcluster[pdcluster.sclabel=='sc3']

