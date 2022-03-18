import os 
import subprocess 
import numpy as np 
from random import seed 
from random import random
from timeit import default_timer as timer
import math
import matplotlib.pyplot as plt
import statistics 
from collections import Counter

f = open('Minimum_distance_frequencies.txt','r')
lines = f.readlines()
xs = [] ; ys = [] ; ays = []
for line in lines:
    items = line.split()
    xs.append(int(items[1]))
    ys.append(int(items[-1].replace(']','')))
    if ys[-1] < 0:
        ays.append(ys[-1])

dlow = 36*365
dhigh = 36*365 + 35 
cagenums = [] ; events = []
print(max(ays))
central_3x3_cage_nums = [294,375,456,293,374,455,292,373,454,285,366,447,284,365,446,283,364,445,276,357,438,275,356,437,274,355,436]
central_to_corner = {285:6, 366:15, 447:24, 284:5, 365:14, 446:23, 283:4, 364:13, 445:22, 274:1, 355:10, 436:19, 275:2, 356:11, 437:20, 276:3, 357:12, 438:21, 292:7, 373:16, 454:25, 293:8, 374:17, 455:26, 294:9, 375:18, 456:27}
corner_inds_x1 = [] ; corner_inds_x2 = [] ; corner_inds_x3 = [] 
for x,y in zip(xs,ys):
    '''
    if dlow <= x <= dhigh:
        nums += 1
    ''' 
    c = int(x/36)
    if c in central_3x3_cage_nums:
        cagenums.append(c)
        if y < 0:
            events.append(y)
        #print(x-(36*c))
    
        if central_to_corner[c] in [1,2,3,4,5,6,7,8,9]:
            fac = abs(y)/y
            corner_inds_x1.append(fac*(36*(central_to_corner[c]-1) + (x - 36*c)))
        if central_to_corner[c] in [10,11,12,13,14,15,16,17,18]:
            fac = abs(y)/y
            corner_inds_x2.append(fac*(36*(central_to_corner[c]-1) + (x - 36*c)))
        if central_to_corner[c] in [19,20,21,22,23,24,25,26,27]:
            fac = abs(y)/y
            corner_inds_x3.append(fac*(36*(central_to_corner[c]-1) + (x - 36*c)))

print(corner_inds_x1)
print(corner_inds_x2)
print(corner_inds_x3)
    

ucage = Counter(cagenums)

for uc in ucage.keys():
    print(uc, ucage[uc])
#print(len(xs),nums)

print(events)



'''
plt.figure()
#plt.plot(np.array(xs), np.array(ys), 'r', alpha = 0.5)
plt.stem(np.array(xs), np.array(ys), 'r', markerfmt = ' ')
plt.show()
'''
