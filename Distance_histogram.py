import os 
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set_style("white")
from scipy.optimize import minimize
import os
from scipy.optimize import curve_fit
from matplotlib import rcParams
import matplotlib.font_manager as font_manager
import os
import pickle
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")
import numpy as np
######################################################################################################################
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'] ; i  = 0
######################### Plotters #######################################################################################
rcParams['lines.linewidth'] = 3.0
rcParams['xtick.labelsize'] = 20
rcParams['xtick.major.width'] = 0.8
rcParams['xtick.major.pad'] = 5
rcParams['xtick.major.size'] = 12.0
rcParams['xtick.major.width'] = 4.0
rcParams['ytick.labelsize'] = 20
rcParams['ytick.major.width'] = 0.8
rcParams['ytick.major.size'] = 12.0
rcParams['ytick.major.width'] = 4.0
rcParams['ytick.major.pad'] = 5
rcParams['axes.linewidth']  = 4.0
rcParams['axes.labelsize'] = 20
rcParams['ytick.right'] = True 
rcParams['xtick.top'] = True
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['figure.figsize'] = 12,10
######################################################################################################################
cwd = os.getcwd()
dir = os.listdir(cwd)
#print(dir)
fold = []
for a in dir:
    if a[0]=='1':
        fold.append(a)


f = open('kinetic_input.dat', 'r')
lines = f.readlines()
for line in lines:
    if 'cutoff' in line:
        D = float(line.split()[-1].replace('\n',''))
    if 'slope' in line:
        m = float(line.split()[-1].replace('\n',''))
print(m,D)
Dees = []
'''
for f in fold:
    os.chdir('{0}/{1}'.format(cwd,f))
    try:
        f = open('Distance_statistics.txt','r')
    except:
        continue
    lines = f.readlines()
    items = lines[0].split()
    for i in items:
        Dees.append(float(i))
'''
f = open('Static_distances.txt','r')
lines = f.readlines()
items = lines[0].split()
for i in items:
    Dees.append(float(i))
d = np.linspace(0,max(Dees),10000) ; dec = []
for dis in d:
    #dec.append(0.5*(math.erf(3.0699*(2*D-dis)+1)))
    #dec.append(0.5*(math.erf(m*(2*D-dis))+1))
    dec.append(1/(1 + np.exp(m*(dis-2*D))))    
    #dec.append(np.tanh(m*(dis-2*D)))

fig, ax1   = plt.subplots()
#print(Dees)
kwargs = dict(alpha=0.5, bins=100, density=True, stacked=True)
kwargs = dict(hist_kws={'alpha':0.5}, kde_kws={'linewidth':2})

#plt.hist(Dees, bins=20)
#plt.hist(Dees, **kwargs, color='g', label='Ideal')
sns.distplot(Dees, color= 'k', label="Compact", **kwargs)
plt.gca().set(ylabel='Normalized probability', xlabel = r'Cu(I) pair distance d/ $\AA$')
plt.xlim(0, 100)
plt.ylim(0, 0.1)
ax2 = ax1.twinx()
ax2.plot(np.array(d), np.array(dec), linewidth = 1.25, linestyle = 'dashed', color = 'r')
ax2.set_ylabel(r'Decay factor $\lambda$(d)', fontsize = 20, color = 'r')
ax2.set_ylim(0,1.05)
ax2.spines['right'].set_color('red')
ax2.yaxis.label.set_color('red')
ax2.tick_params(axis='y', colors='red')
os.chdir(cwd)
plt.savefig('Pair-wise static distribution.png', dpi=600, orientation='portrait',bbox_inches='tight', pad_inches=0.1)
plt.close()
#plt.legend();
#plt.show()

Dees = []
f = open('Distance_statistics.txt','r')
lines = f.readlines()
items = lines[0].split()
for i in items:
    Dees.append(float(i))
d = np.linspace(0,max(Dees),10000) ; dec = []
for dis in d:
    #dec.append(0.5*(math.erf(3.0699*(2*D-dis)+1)))
    #dec.append(0.5*(math.erf(m*(2*D-dis))+1))
    dec.append(1/(1 + np.exp(m*(dis-2*D))))    
    #dec.append(np.tanh(m*(dis-2*D)))

fig, ax1   = plt.subplots()
#print(Dees)
kwargs = dict(alpha=0.5, bins=100, density=True, stacked=True)
kwargs = dict(hist_kws={'alpha':0.5}, kde_kws={'linewidth':2})

#plt.hist(Dees, bins=20)
#plt.hist(Dees, **kwargs, color='g', label='Ideal')
sns.distplot(Dees, color= 'k', label="Compact", **kwargs)
plt.gca().set(ylabel='Normalized probability', xlabel = r'Cu(I) pair distance d/ $\AA$')
plt.xlim(0, 100)
plt.ylim(0, 0.1)
ax2 = ax1.twinx()
ax2.plot(np.array(d), np.array(dec), linewidth = 1.25, linestyle = 'dashed', color = 'r')
ax2.set_ylabel(r'Decay factor $\lambda$(d)', fontsize = 20, color = 'r')
ax2.set_ylim(0,1.05)
ax2.spines['right'].set_color('red')
ax2.yaxis.label.set_color('red')
ax2.tick_params(axis='y', colors='red')
os.chdir(cwd)
plt.savefig('Pair-wise kinetic distribution.png', dpi=600, orientation='portrait',bbox_inches='tight', pad_inches=0.1)
plt.close()
#plt.legend();
#plt.show()
