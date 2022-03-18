import os 
import subprocess 
import numpy as np 
from random import seed 
from random import random
from timeit import default_timer as timer
from .Zmc import Zmc 
from .monkeypatch import monkeypatch_class
import math
import matplotlib.pyplot as plt
import statistics

@monkeypatch_class(Zmc)
def EWMA(self,t,x):
    nsize = self.NCu
    EWMA = []
    lam = 0.10 #Decay function 
    for i in range(0,len(x)):
        if i == 0:
            temp = float(x[i])
        else:
            temp = lam*float(x[i-1])+(1-lam)*float(x[i])   
        EWMA.append(temp)
    nums = []
    for a in x:
        nums.append(float(a))
    #plt.plot(t,EWMA) #,t,EWMA)
    #plt.ylim(0,10)
    #plt.show()
    #print(len(EWMA)
    #return EWMA
    sig = np.std(nums)
    L = 1
    LCL = []
    UCL = []
    for i in range(0,len(x)):
        UCL.append(EWMA[-1] + L*sig*np.sqrt(lam*(1-np.power((1-lam),i))/(2-lam)))
        LCL.append(EWMA[-1] - L*sig*np.sqrt(lam*(1-np.power((1-lam),i))/(2-lam)))
    print('Plotting things here')
    #plt.plot(t,UCL,'-')
    #plt.plot(t,LCL,'-')
    #plt.plot(t,EWMA,'*')
    #plt.show()
    cf = [] 
    for i in range(0,len(x)):
        if(EWMA[i]<UCL[i]):
            first = i
            break
    #This tool does not take care that the decay signal lies between UCL and LCL so I am going to tech this with a checker
    pc = 0.10
    #print(first)
    self.cutoff = first
    supsum=0
    flag = 0
    covlist = []
    conv = []
    if(first/len(x)<pc):
        print("Simulation sufficiently converged. Calculating averages")
    else:
        print("Simulation not converged")
        return False
    for i in range(first,len(x)):
        supsum = supsum + nums[i]/nsize
        flag = flag + 1
        covlist.append(nums[i]/(nsize))
        conv.append(nums[i])
    conv = np.array(conv)
    sd = statistics.stdev(conv)
    fin_nums = []
    for i in range(0,len(x)):
        if EWMA[i] < UCL[i] and EWMA[i] > LCL[i]:
            fin_nums.append(nums[i])
    print('The modified fraction is:',statistics.mean(fin_nums)/nsize)
    print('The median fraction is:', statistics.median(nums[first :])/nsize)
    return(supsum/flag,flag,covlist,sd)

@monkeypatch_class(Zmc)
def get_fraction(self):
    ss,nc,cl,std = self.EWMA(self.t,self.tw)
    print('The unmodified fraction is:', ss)
   
@monkeypatch_class(Zmc)
def get_rate(self):
    model = np.polyfit(self.t[self.cutoff :],self.rev[self.cutoff :],1)
    print('The rate is:',model[0])
    