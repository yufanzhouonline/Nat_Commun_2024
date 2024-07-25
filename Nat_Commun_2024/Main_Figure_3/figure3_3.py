
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.chdir('/data/yufan/schic/insulation/unnormalized')
os.getcwd()

allchrno = ["chr" + str(i+1) for i in range(22)] + ['chrX']

samplename = pd.read_csv('/data/yufan/schic/schicnames.txt', header=0, sep="\t")
###Change the batch 3 to batch 6
samplename['Round'] = [6 if i==3 else i for i in samplename.Round]

samplecount = []
samplecount.append(samplename[samplename['Sample'] == 'MCF7'].shape[0])
samplecount.append(samplename[samplename['Sample'] == 'MCF7-T1M'].shape[0])
samplecount.append(samplename[samplename['Sample'] == 'MCF7-TamR'].shape[0])

###Plot the sample count of single cell Hi-C
rects1 = plt.bar(['MCF7', 'MCF7M1', 'MCF7TR'], samplecount, color = ['#88224488', '#22884488', '#22448888'])
def autolabel(rects):
	for rect in rects:
		height = rect.get_height()
		#plt.text(rect.get_x()+rect.get_width() / 2 - 0.1, 1.03*height, '%s' % float(height))
		plt.text(rect.get_x()+rect.get_width() / 2 - 0.1, 1.03*height, str(height))

autolabel(rects1)
plt.ylim(0, 125)
plt.ylabel('# of cells', fontsize=15)
plt.show()

#MCF7
groupname = samplename[samplename.Sample=='MCF7'].copy()
groupname.reset_index(drop = True, inplace = True)

pdcluster = pd.read_csv('/data/yufan/schic/cluster/clusterlist001.txt',sep='\t',header=0)
pdcluster['clusterlabel'] = ['C' + str(i+1) for i in pdcluster.cluster]

####clusters common bins
clusterlist = []
for i in range(9):
	clusterlist.append(list(pdcluster[pdcluster.cluster==i]['sample']))


#common bins
def commonbin(clusterlist, res=1000000):
	groupname = samplename[samplename.Name.isin(clusterlist)].copy()
	groupname.reset_index(drop = True, inplace = True)
	setlist = []
	#res = 2000000
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
	#print(len(commonset))
	###CL1: 295, CL2: 2, CL3: 669 , CL4: 1, CL5: 268 , CL6: 8, CL7: 3, CL8: 2, CL9: 321
	###
	return len(commonset)

cellnumber = [25, 45, 8, 16, 54, 41, 26, 8, 8]
percell = []
for i in range(9):
	percell.append(commonbin(clusterlist[i], 1000000) / cellnumber[i])

#common bin set
def commonbinset(clusterlist, res=1000000):
	groupname = samplename[samplename.Name.isin(clusterlist)].copy()
	groupname.reset_index(drop = True, inplace = True)
	setlist = []
	#res = 2000000
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
	#print(len(commonset))
	###CL1: 295, CL2: 2, CL3: 669 , CL4: 1, CL5: 268 , CL6: 8, CL7: 3, CL8: 2, CL9: 321
	###
	return commonset

def commonbin(clusterno, res=1000000):
	commonset = commonbinset(clusterlist[clusterno], res=res)
	###loci of CADs
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
	commonbinbed = binbed[binbed.pos.isin(loci1)].copy()
	return commonbinbed

#import pdb

def isids(res):
	if str(res) == '40000':
		isidsdict = {'ispost':'520001', 'idspost':'240001'}
	elif str(res) == '200000':
		isidsdict = {'ispost':'600001', 'idspost':'400001'}
	elif str(res) == '300000':
		isidsdict = {'ispost':'600001', 'idspost':'600001'}
	elif str(res) == '400000':
		isidsdict = {'ispost':'800001', 'idspost':'800001'}
	elif str(res) == '500000':
		isidsdict = {'ispost':'1000001', 'idspost':'1000001'}
	elif str(res) == '1000000':
		isidsdict = {'ispost':'2000001', 'idspost':'2000001'}
	else:
		isidsdict = {'ispost':'500001', 'idspost':'200001'}
	return isidsdict

###Get the shifted TADs within CADs
###clusterno: Number of cluster
###commonbin: Common bin of cluster (CADs)
###res: Resolution for TADs
###cads: True, CADs, False: NADs
def cellshift(clusterno, commonbin, res=100000, cads=True):
	loadsample = samplename[samplename.Name.isin(pdcluster[pdcluster.cluster==clusterno]['sample'])].copy()
	loadsample.reset_index(drop = True, inplace = True)
	cellmeanlist = []
	for loadno in loadsample.index:
		shiftblist = []
		for chrno in allchrno:
			#for chrno in ['chr20']:
			#chrno = 'chr20'
			filepost = '_' + str(res) + '_' + chrno + '.is' + isids(res)['ispost'] + '.ids' + isids(res)['idspost'] + '.insulation.boundaries.bed'
			if res==300000 or res==400000:
				readfilepath = '/data/yufan/schic/res300k400k/insulation/batch' + str(loadsample.loc[loadno, 'Round']) + '/' + loadsample.loc[loadno, 'Name'] + '/' + loadsample.loc[loadno, 'Name'] + filepost
			else:
				readfilepath = '/data/yufan/schic/insulation/batch' + str(loadsample.loc[loadno, 'Round']) + '/' + loadsample.loc[loadno, 'Name'] + '/' + loadsample.loc[loadno, 'Name'] + filepost
			try:
				scmcf7tads = pd.read_csv(readfilepath, header=None, sep="\t", skiprows=1)
			except Exception:
				print('Read ERROR of file ' + readfilepath)
				continue
			scmcf7tads.columns = ['chr', 'startb', 'endb', 'name', 'value']
			#pdb.set_trace()
			scmcf7tads['midb'] = (scmcf7tads.startb + scmcf7tads.endb) // 2
			chrccd = commonbin[commonbin['chr']==chrno].copy()
			chrccd.reset_index(drop = True, inplace = True)
			scmcf7tads['incad'] = 0
			for i in scmcf7tads.index:
				withinccd = chrccd[(scmcf7tads.loc[i, 'midb'] > chrccd.bin1) & (scmcf7tads.loc[i, 'midb'] < chrccd.bin2)].shape[0]
				print('Cluster:', 'C' + str(clusterno+1), loadsample.loc[loadno, 'Name'], chrno, i, withinccd)
				if withinccd:
					scmcf7tads.loc[i, 'incad'] = 1
			if cads:
				scintads = scmcf7tads[scmcf7tads.incad==1].copy()
			else:
				scintads = scmcf7tads[scmcf7tads.incad==0].copy()
			scintads.reset_index(drop = True, inplace = True)
			if loadsample.loc[loadno, 'Sample'] == 'MCF7':
				if res==300000 or res==400000:
					mcf7tads = pd.read_csv('/data/yufan/schic/res300k400k/insulation/bulkcells/mcf7/mcf7p/mcf7p' + filepost, header=None, sep="\t", skiprows=1)
				else:
					mcf7tads = pd.read_csv('/data/yufan/schic/insulation/bulkcells/mcf7/mcf7p/mcf7p' + filepost, header=None, sep="\t", skiprows=1)
			elif loadsample.loc[loadno, 'Sample'] == 'MCF7-T1M':
				if res==300000 or res==400000:
					mcf7tads = pd.read_csv('/data/yufan/schic/res300k400k/insulation/bulkcells/mcf7m1/MCF7-T1M/MCF7-T1M' + filepost, header=None, sep="\t", skiprows=1)
				else:
					mcf7tads = pd.read_csv('/data/yufan/schic/insulation/bulkcells/mcf7m1/MCF7-T1M/MCF7-T1M' + filepost, header=None, sep="\t", skiprows=1)
			elif loadsample.loc[loadno, 'Sample'] == 'MCF7-TamR':
				if res==300000 or res==400000:
					mcf7tads = pd.read_csv('/data/yufan/schic/res300k400k/insulation/bulkcells/mcf7tr/MCF7TR/MCF7TR' + filepost, header=None, sep="\t", skiprows=1)
				else:
					mcf7tads = pd.read_csv('/data/yufan/schic/insulation/bulkcells/mcf7tr/MCF7TR/MCF7TR' + filepost, header=None, sep="\t", skiprows=1)
			else:
				print('ERROR')
				exit
			mcf7tads.columns = ['chr', 'startb', 'endb', 'name', 'value']
			mcf7tads['midb'] = (mcf7tads.startb + mcf7tads.endb) // 2
			for j in scintads.index:
				shiftboundary = np.min(np.abs(scintads.loc[j, 'midb'] - mcf7tads.midb))
				shiftblist.append(shiftboundary)
		cellmeanlist.append(np.mean(shiftblist))
	return cellmeanlist

reslist = [50000, 100000, 200000, 300000, 400000, 500000]
#reslist = [40000, 50000, 100000, 200000, 300000, 400000, 500000, 1000000]
shiftcadslist = []
shiftnadslist = []
for res in reslist:
	shiftcads = []
	for i in range(9):
		shiftcads.append(cellshift(i, commonbin(i), res=res, cads=True))
	shiftcadslist.append(shiftcads)
	shiftnads = []
	for i in range(9):
		shiftnads.append(cellshift(i, commonbin(i), res=res, cads=False))
	shiftnadslist.append(shiftnads)

#Compute the Wilcoxon rank-sum statistic for two samples.
from scipy import stats

###Print the mean
for i in range(len(reslist)):
	for j in range(9):
		allcads = [int(k) for k in shiftcadslist[i][j] if str(k)!='nan']
		allnads = [int(k) for k in shiftnadslist[i][j] if str(k)!='nan']
		meancads = round(np.mean(allcads), 0)
		meannads = round(np.mean(allnads), 0)
		diff = meancads - meannads
		statvalue, pvalue = stats.ranksums(allcads, allnads)
		print('Resolution:', reslist[i], 'Cluster:', 'C'+str(j+1), meancads, meannads, diff, pvalue, sep='\t')



valuelist1 = []
valuelist2 = []
stdlist1 = []
stdlist2 = []
for k in range(len(reslist)):
	value1 = [j for i in shiftcadslist[k] for j in i if str(j) != 'nan']
	value2 = [j for i in shiftnadslist[k] for j in i if str(j) != 'nan']
	valuelist1.append(value1)
	valuelist2.append(value2)
	stdlist1.append(np.std(value1))
	stdlist2.append(np.std(value2))


###Make figures for Standard deviation
label_list = [str(int(i/1000)) for i in reslist]    
num_list1 = [int(i/1000) for i in stdlist1]   
num_list2 = [int(i/1000) for i in stdlist2]   
x = range(len(num_list1))

rects1 = plt.bar(x, height=num_list1, width=0.4, alpha=0.8, color='#88448888', label="CADs")
rects2 = plt.bar([i + 0.4 for i in x], height=num_list2, width=0.4, color='#00AAAA88', label="NADs")
#plt.ylim(0, 50)   
plt.ylabel("SD of shifted boundaries of TADs (Kb)", fontsize=12)

plt.xticks([index + 0.2 for index in x], label_list)
plt.xlabel("TAD bin size (Kb)", fontsize=12)
#plt.title("Company")
plt.legend()   

for rect in rects1:
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2, height+1, str(height), ha="center", va="bottom")

for rect in rects2:
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2, height+1, str(height), ha="center", va="bottom")

plt.show()

###Make figures for mean of boxplot
boxdf = pd.DataFrame({'Resolution':[], 'Domain':[], 'Value':[]})
boxdf = boxdf.astype({'Resolution':int, 'Domain':str, 'Value':float})
for i in range(len(reslist)):
	for j in valuelist1[i]:
		tempdf = pd.DataFrame({'Resolution':[int(reslist[i]/1000)], 'Domain':['CADs'], 'Value':[j/1000]})
		boxdf = pd.concat([boxdf, tempdf], axis=0)
	for j in valuelist2[i]:
		tempdf = pd.DataFrame({'Resolution':[int(reslist[i]/1000)], 'Domain':['NADs'], 'Value':[j/1000]})
		boxdf = pd.concat([boxdf, tempdf], axis=0)

import seaborn as sns
sns.boxplot(x='Resolution', y='Value',data=boxdf, hue='Domain', palette='Set3', showfliers = True, fliersize = 5)
plt.legend(loc='upper left') #upper right, upper left, lower left, lower right, right, center left, center right, lower center, upper center, center
plt.xlabel('TAD bin size (Kb)', fontsize = 15)
plt.ylabel('Shifted boundaries of TADs (Kb)', fontsize = 15)
plt.ylim(-100, 1800)
plt.show()

for i in range(len(reslist)):
	statvalue, pvalue = stats.ranksums(valuelist1[i], valuelist2[i])
	print(reslist[i], pvalue, sep='\t')


###Each cluster in various resolution
cvaluelist1 = []
for i in range(9):
	cvalue = []
	for j in range(len(reslist)):
		cvalue.append([k for k in shiftcadslist[j][i] if str(k) != 'nan'])
	cvaluelist1.append(cvalue)

cvaluelist2 = []
for i in range(9):
	cvalue = []
	for j in range(len(reslist)):
		cvalue.append([k for k in shiftnadslist[j][i] if str(k) != 'nan'])
	cvaluelist2.append(cvalue)


for i in range(len(reslist)):
	diff = np.mean(cvaluelist1[4][i]) - np.mean(cvaluelist2[4][i])
	statvalue, pvalue = stats.ranksums(cvaluelist1[4][i], cvaluelist2[4][i])
	print(reslist[i], np.mean(cvaluelist1[4][i]), np.mean(cvaluelist2[4][i]), diff, pvalue, sep='\t')


clusterno = 0
for i in range(len(reslist)):
	valuelist1.append(cvaluelist1[clusterno][i])
	valuelist2.append()
	diff = np.mean(cvaluelist1[clusterno][i]) - np.mean(cvaluelist2[clusterno][i])
	statvalue, pvalue = stats.ranksums(cvaluelist1[clusterno][i], cvaluelist2[clusterno][i])
	print(reslist[i], np.mean(cvaluelist1[clusterno][i]), np.mean(cvaluelist2[clusterno][i]), diff, pvalue, sep='\t')



clusterno = 8
###Make figures for mean of boxplot
boxdf = pd.DataFrame({'Resolution':[], 'Domain':[], 'Value':[]})
boxdf = boxdf.astype({'Resolution':int, 'Domain':str, 'Value':float})
for i in range(len(reslist)):
	for j in cvaluelist1[clusterno][i]:
		tempdf = pd.DataFrame({'Resolution':[int(reslist[i]/1000)], 'Domain':['CADs'], 'Value':[j/1000]})
		boxdf = pd.concat([boxdf, tempdf], axis=0)
	for j in cvaluelist2[clusterno][i]:
		tempdf = pd.DataFrame({'Resolution':[int(reslist[i]/1000)], 'Domain':['NADs'], 'Value':[j/1000]})
		boxdf = pd.concat([boxdf, tempdf], axis=0)

import seaborn as sns
sns.boxplot(x='Resolution', y='Value',data=boxdf, hue='Domain', palette='Set3', showfliers=False)
plt.xlabel('TAD bin size (Kb)', fontsize = 15)
plt.ylabel('Shifted boundaries of TADs (Kb)', fontsize = 15)
plt.legend(loc='upper left') #upper right, upper left, lower left, lower right, right, center left, center right, lower center, upper center, center
plt.ylim(170, 820)
plt.show()







for i in range(len(reslist)):
	statvalue, pvalue = stats.ranksums(valuelist1[i], valuelist2[i])
	print(reslist[i], pvalue, sep='\t')



value1 = [j for i in shiftcadslist for j in i if str(j) != 'nan']
value2 = [j for i in shiftnadslist for j in i if str(j) != 'nan']

statvalue, pvalue = stats.ranksums(value1, value2)
print(pvalue)

print(pvalue)
#8.935337849460439e-05
np.mean(value1)
#168591.83988174552
np.mean(value2)
#201627.56869813104

valuelist1 = [value1]
valuelist2 = [value2]

valuelist1.append(value1)
valuelist2.append(value2)


###Make figures
all_data = [pd.Series(valuelist1[2])//1000, pd.Series(valuelist2[2])//1000]
#all_data = [list(np.array(celllenc0)/1000), list(np.array(celllenc1)/1000), list(np.array(celllenc2)/1000), list(np.array(celllenc3)/1000), list(np.array(celllenc4)/1000), list(np.array(celllenc5)/1000), list(np.array(celllenc6)/1000), list(np.array(celllenc7)/1000), list(np.array(celllenc8)/1000)]
labels = ['CADs', 'NADs']
medianprops = dict(linestyle='-', linewidth=1, color='red')
#bplot = plt.boxplot(all_data, patch_artist=True, labels=labels, medianprops=medianprops) 
bplot = plt.boxplot(all_data, patch_artist=True, labels=labels, medianprops=medianprops, showfliers=False)  
#plt.title('Rectangular box plot')

#colors = ['darkred', 'darkorange', 'green', 'blue', 'c', 'm', 'y', 'lightblue', 'pink']
colors = ['#FF000088', '#00FF0088', '#0000FF88', '#44FFFF88', '#FF44FF88', '#FFFF4488', '#00888888', '#88008888', '#88880088']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)  

#for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
#	plt.setp(bplot[element], color=edge_color)


#plt.xlabel('Clusters', fontsize=15)
plt.ylabel('Shifted boundaries of TADs (Kb)', fontsize=15)
#plt.ylabel('Size of TADs (Kb)', fontsize=15)
#plt.ylim(180, 270)
plt.show()

np.std(valuelist1[2])

np.std(valuelist2[2])

statvalue, pvalue = stats.ranksums(valuelist1[2], valuelist2[2])
print(pvalue)
