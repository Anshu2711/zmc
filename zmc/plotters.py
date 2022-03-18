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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

@monkeypatch_class(Zmc)
def history_plotter(self):
    os.system('rm -r Figures')
    os.system('mkdir Figures')
    L = self.xl
    VecStart_x = [0,0,L,L,0,0,L,L,0,0,0,0]
    VecStart_y = [0,L,0,L,0,0,0,0,0,0,L,L]
    VecStart_z = [0,0,0,0,0,L,0,L,0,L,0,L]
    VecEnd_x = [0,0,L,L,0,0,L,L,L,L,L,L]
    VecEnd_y = [0,L,0,L,L,L,L,L,0,0,L,L]
    VecEnd_z  =[L,L,L,L,0,L,0,L,0,L,0,L]
    xs = [] ; ys = [] ; zs = []
    for c in self.Cu_list:
        [x,y,z] = c.position
        xs.append(x) ; ys.append(y) ; zs.append(z)
    flag = 1
    ss = [48]*self.NCu
    for x,y,z in zip(self.snap_times,self.positions,self.snap_events):
        #print(x) ; print(y)
        cs = []
        for p in y:
            if p == 1:
                cs.append('b')
            else:
                cs.append('r')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax._axis3don = False
        for i in range(12):
            ax.plot([VecStart_x[i], VecEnd_x[i]], [VecStart_y[i],VecEnd_y[i]],zs=[VecStart_z[i],VecEnd_z[i]], color = 'k')
        ax.scatter(xs,ys,zs,color=cs,marker='o',s=ss)
        plt.title('No. of events = {:03d} Time = {:1.5f} seconds'.format(z,x))
        plt.savefig('Figures/Figure {0}.png'.format(flag))
        flag += 1
        ax.view_init(elev=10, azim=10)
        plt.close()
