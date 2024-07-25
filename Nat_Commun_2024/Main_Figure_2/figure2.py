###make plot of genomic distance and relative contact probability

import os
import warnings
warnings.filterwarnings("ignore")
from itertools import combinations

# import semi-core packages
import matplotlib.pyplot as plt
from matplotlib import colors
#%matplotlib inline
plt.style.use('seaborn-poster')
import numpy as np
import pandas as pd

# import open2c libraries
import bioframe

import cooler
import cooltools
import cooltools.expected

samplename = pd.read_csv('/data/yufan/schic2/schicnames2.txt', header=0, sep='\t', skiprows=0)

def getcvd(prefix = '/data/yufan/schic2/contactprobability/mcf7', cell = 'RP01'):
    clr = cooler.Cooler(prefix + '/' + cell + '.cool::/')
    hg19_chromsizes = bioframe.fetch_chromsizes('hg19')
    hg19_cens = bioframe.fetch_centromeres('hg19')
    hg19_arms = bioframe.core.construction.add_ucsc_name_column(bioframe.make_chromarms(hg19_chromsizes,  hg19_cens))
    hg19_arms = hg19_arms[hg19_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)
    # cvd == contacts-vs-distance
    cvd = cooltools.expected.diagsum(
        clr=clr,
        view_df=hg19_arms,
        transforms={'balanced': lambda p: p['count']*p['weight1']*p['weight2']}
    )
    return clr, cvd

################################################
###MCF7

###combined single cells
combined_clr, combined_cvd = getcvd(prefix='/data/yufan/schic2/contactprobability/combined', cell='mcf7')

###single cell
cwd = '/data/yufan/schic2/contactprobability/mcf7'
os.chdir(cwd)
os.getcwd()

cellname = samplename[samplename.Sample=='MCF7']['Name']
for cell in cellname:
    print('Treat single cell: ' + cell)
    sc_cvd = combined_cvd.copy()
    clr, cvd = getcvd(prefix = '/data/yufan/schic2/contactprobability/mcf7', cell = cell)
    sc_cvd['n_valid'] = cvd['n_valid'] * cellname.shape[0] + combined_cvd['n_valid']
    sc_cvd['count.sum'] = cvd['count.sum'] * cellname.shape[0] + combined_cvd['count.sum']
    sc_cvd['balanced.sum'] = cvd['balanced.sum'] * cellname.shape[0] + combined_cvd['balanced.sum']
    ############################
    ###single cell
    # Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
    sc_lb_cvd, sc_lb_slopes, sc_lb_distbins = cooltools.expected.logbin_expected(sc_cvd)
    # The resulting table contains P(s) curves for each individual region.
    # Aggregating these curves into a single genome-wide curve is involving too, 
    # so we created a separate function for this too.
    sc_lb_cvd_agg, sc_lb_slopes_agg = cooltools.expected.combine_binned_expected(
        sc_lb_cvd,
        binned_exp_slope=sc_lb_slopes
    )
    sc_lb_cvd_agg['s_bp'] = sc_lb_cvd_agg['diag.avg'] * clr.binsize 
    sc_lb_slopes_agg['s_bp'] = sc_lb_slopes_agg['diag.avg'] * clr.binsize
    #################
    plt.loglog(
        sc_lb_cvd_agg['s_bp'],
        sc_lb_cvd_agg['balanced.avg'] / sc_lb_cvd_agg.loc[0, 'balanced.avg'],
        '-',
        markersize=5,
        color='#88FF88FF',
        linewidth=0.5
    )

############################
###combined plot
# Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(combined_cvd)
# The resulting table contains P(s) curves for each individual region.
# Aggregating these curves into a single genome-wide curve is involving too, 
# so we created a separate function for this too.
lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(
    lb_cvd,
    binned_exp_slope=lb_slopes
)
lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * combined_clr.binsize 
lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * combined_clr.binsize

plt.loglog(
    lb_cvd_agg['s_bp'],
    #lb_cvd_agg['balanced.avg'],
    lb_cvd_agg['balanced.avg'] / lb_cvd_agg.loc[0, 'balanced.avg'],
    '-',
    markersize=5,
    color='black',
    linewidth=2
)


plt.xlim(10**(4),10**(9))
plt.ylim(10**(-4),10**(2/3))
#plt.xlabel('Genomic distance (bp)')
#plt.ylabel('Relative contact probability')
plt.show()


################################################
###MCF7M1

###combined single cells
combined_clr, combined_cvd = getcvd(prefix='/data/yufan/schic2/contactprobability/combined', cell='mcf7m1')

###single cell
cwd = '/data/yufan/schic2/contactprobability/mcf7m1'
os.chdir(cwd)
os.getcwd()

cellname = samplename[samplename.Sample=='MCF7-T1M']['Name']
for cell in cellname:
    print('Treat single cell: ' + cell)
    sc_cvd = combined_cvd.copy()
    clr, cvd = getcvd(prefix = '/data/yufan/schic2/contactprobability/mcf7m1', cell = cell)
    sc_cvd['n_valid'] = cvd['n_valid'] * cellname.shape[0] + combined_cvd['n_valid']
    sc_cvd['count.sum'] = cvd['count.sum'] * cellname.shape[0] + combined_cvd['count.sum']
    sc_cvd['balanced.sum'] = cvd['balanced.sum'] * cellname.shape[0] + combined_cvd['balanced.sum']
    ############################
    ###single cell
    # Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
    sc_lb_cvd, sc_lb_slopes, sc_lb_distbins = cooltools.expected.logbin_expected(sc_cvd)
    # The resulting table contains P(s) curves for each individual region.
    # Aggregating these curves into a single genome-wide curve is involving too, 
    # so we created a separate function for this too.
    sc_lb_cvd_agg, sc_lb_slopes_agg = cooltools.expected.combine_binned_expected(
        sc_lb_cvd,
        binned_exp_slope=sc_lb_slopes
    )
    sc_lb_cvd_agg['s_bp'] = sc_lb_cvd_agg['diag.avg'] * clr.binsize 
    sc_lb_slopes_agg['s_bp'] = sc_lb_slopes_agg['diag.avg'] * clr.binsize
    #################
    plt.loglog(
        sc_lb_cvd_agg['s_bp'],
        sc_lb_cvd_agg['balanced.avg'] / sc_lb_cvd_agg.loc[0, 'balanced.avg'],
        '-',
        markersize=5,
        color='#88FFFFFF',
        linewidth=0.5
    )

############################
###combined plot
# Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(combined_cvd)
# The resulting table contains P(s) curves for each individual region.
# Aggregating these curves into a single genome-wide curve is involving too, 
# so we created a separate function for this too.
lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(
    lb_cvd,
    binned_exp_slope=lb_slopes
)
lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * combined_clr.binsize 
lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * combined_clr.binsize

plt.loglog(
    lb_cvd_agg['s_bp'],
    #lb_cvd_agg['balanced.avg'],
    lb_cvd_agg['balanced.avg'] / lb_cvd_agg.loc[0, 'balanced.avg'],
    '-',
    markersize=5,
    color='black',
    linewidth=2
)


plt.xlim(10**(4),10**(9))
plt.ylim(10**(-4),10**(2/3))
#plt.xlabel('Genomic distance (bp)')
#plt.ylabel('Relative contact probability')
plt.show()

################################################
###MCF7TR

###combined single cells
combined_clr, combined_cvd = getcvd(prefix='/data/yufan/schic2/contactprobability/combined', cell='mcf7tr')

###single cell
cwd = '/data/yufan/schic2/contactprobability/mcf7tr'
os.chdir(cwd)
os.getcwd()

cellname = samplename[samplename.Sample=='MCF7-TamR']['Name']
for cell in cellname:
    print('Treat single cell: ' + cell)
    sc_cvd = combined_cvd.copy()
    clr, cvd = getcvd(prefix = '/data/yufan/schic2/contactprobability/mcf7tr', cell = cell)
    sc_cvd['n_valid'] = cvd['n_valid'] * cellname.shape[0] + combined_cvd['n_valid']
    sc_cvd['count.sum'] = cvd['count.sum'] * cellname.shape[0] + combined_cvd['count.sum']
    sc_cvd['balanced.sum'] = cvd['balanced.sum'] * cellname.shape[0] + combined_cvd['balanced.sum']
    ############################
    ###single cell
    # Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
    sc_lb_cvd, sc_lb_slopes, sc_lb_distbins = cooltools.expected.logbin_expected(sc_cvd)
    # The resulting table contains P(s) curves for each individual region.
    # Aggregating these curves into a single genome-wide curve is involving too, 
    # so we created a separate function for this too.
    sc_lb_cvd_agg, sc_lb_slopes_agg = cooltools.expected.combine_binned_expected(
        sc_lb_cvd,
        binned_exp_slope=sc_lb_slopes
    )
    sc_lb_cvd_agg['s_bp'] = sc_lb_cvd_agg['diag.avg'] * clr.binsize 
    sc_lb_slopes_agg['s_bp'] = sc_lb_slopes_agg['diag.avg'] * clr.binsize
    #################
    plt.loglog(
        sc_lb_cvd_agg['s_bp'],
        sc_lb_cvd_agg['balanced.avg'] / sc_lb_cvd_agg.loc[0, 'balanced.avg'],
        '-',
        markersize=5,
        color='#FF88FFFF',
        linewidth=0.5
    )

############################
###combined plot
# Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(combined_cvd)
# The resulting table contains P(s) curves for each individual region.
# Aggregating these curves into a single genome-wide curve is involving too, 
# so we created a separate function for this too.
lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(
    lb_cvd,
    binned_exp_slope=lb_slopes
)
lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * combined_clr.binsize 
lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * combined_clr.binsize

plt.loglog(
    lb_cvd_agg['s_bp'],
    #lb_cvd_agg['balanced.avg'],
    lb_cvd_agg['balanced.avg'] / lb_cvd_agg.loc[0, 'balanced.avg'],
    '-',
    markersize=5,
    color='black',
    linewidth=2
)


plt.xlim(10**(4),10**(9))
plt.ylim(10**(-4),10**(2/3))
#plt.xlabel('Genomic distance (bp)')
#plt.ylabel('Relative contact probability')
plt.show()

##########################################################################
################################################
###Three combined plot

###combined single cells MCF7
combined_clr, combined_cvd = getcvd(prefix='/data/yufan/schic2/contactprobability/combined', cell='mcf7')
############################
###combined plot
# Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(combined_cvd)
# The resulting table contains P(s) curves for each individual region.
# Aggregating these curves into a single genome-wide curve is involving too, 
# so we created a separate function for this too.
lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(
    lb_cvd,
    binned_exp_slope=lb_slopes
)
lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * combined_clr.binsize 
lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * combined_clr.binsize

plt.loglog(
    lb_cvd_agg['s_bp'],
    #lb_cvd_agg['balanced.avg'],
    lb_cvd_agg['balanced.avg'] / lb_cvd_agg.loc[0, 'balanced.avg'],
    '-',
    markersize=5,
    color='#88FF88',
    linewidth=2
)

###combined single cells MCF7M1
combined_clr, combined_cvd = getcvd(prefix='/data/yufan/schic2/contactprobability/combined', cell='mcf7m1')
############################
###combined plot
# Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(combined_cvd)
# The resulting table contains P(s) curves for each individual region.
# Aggregating these curves into a single genome-wide curve is involving too, 
# so we created a separate function for this too.
lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(
    lb_cvd,
    binned_exp_slope=lb_slopes
)
lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * combined_clr.binsize 
lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * combined_clr.binsize

plt.loglog(
    lb_cvd_agg['s_bp'],
    #lb_cvd_agg['balanced.avg'],
    lb_cvd_agg['balanced.avg'] / lb_cvd_agg.loc[0, 'balanced.avg'],
    '-',
    markersize=5,
    color='#88FFFF',
    linewidth=2
)

###combined single cells MCF7TR
combined_clr, combined_cvd = getcvd(prefix='/data/yufan/schic2/contactprobability/combined', cell='mcf7tr')
############################
###combined plot
# Logbin-expected aggregates P(s) curves per region over exponentially increasing distance bins.
lb_cvd, lb_slopes, lb_distbins = cooltools.expected.logbin_expected(combined_cvd)
# The resulting table contains P(s) curves for each individual region.
# Aggregating these curves into a single genome-wide curve is involving too, 
# so we created a separate function for this too.
lb_cvd_agg, lb_slopes_agg = cooltools.expected.combine_binned_expected(
    lb_cvd,
    binned_exp_slope=lb_slopes
)
lb_cvd_agg['s_bp'] = lb_cvd_agg['diag.avg'] * combined_clr.binsize 
lb_slopes_agg['s_bp'] = lb_slopes_agg['diag.avg'] * combined_clr.binsize

plt.loglog(
    lb_cvd_agg['s_bp'],
    #lb_cvd_agg['balanced.avg'],
    lb_cvd_agg['balanced.avg'] / lb_cvd_agg.loc[0, 'balanced.avg'],
    '-',
    markersize=5,
    color='#FF88FF',
    linewidth=2
)

plt.xlim(10**(4),10**(9))
plt.ylim(10**(-4),10**(2/3))
#plt.xlabel('Genomic distance (bp)')
#plt.ylabel('Relative contact probability')
plt.show()
