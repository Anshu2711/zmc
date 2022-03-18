import random 
import numpy as np
import matplotlib.pyplot as plt
import math
NI = 100 
NII = 0 
Ntot = NI + NII
ro_Cu = 0.001
Ap = 1.0E+010
Eap = 0.45
Ar = 1.0E+010
Ear = 0.35 
kB = 8.61733034e-05 
T = 473 
pO2 = 2800
pNH3 = 1.0
pNO = 1.0
kp = Ap*np.exp(-1*Eap/(kB*T))*pO2 
kr = Ar*np.exp(-1*Ear/(kB*T))*pNH3*pNO  
print('kp is:', kp)
print('kr is:', kr)
ones = []
twos = []
t = 0 
counter = 0
random.seed(12311654)
time = []
Vtot = Ntot/ro_Cu

def combo(a,b):
    ret = math.factorial(a)/(math.factorial(b)*math.factorial(a-b))
    return ret
    
while t < 1e-3 and NI >= -1 and NII >= -1:
    rr = kr*(NII) 
    rp = kp*NI*(NI-1)/Vtot #NI   #*NI  # combo(NI,2)   
    Rtot = rr + rp
    frac1 = rr/Rtot 
    frac2 = rp/Rtot
    r1 = random.uniform(0,1)
    if r1 > frac1:
        NII = NII + 2 #  1  2
        NI = NI - 2   #  1 2
    else:
        NI = NI + 1
        NII = NII - 1 
    r2 = random.uniform(0,1)
    t = t - np.log(r2)/Rtot
    ones.append(NI)
    twos.append(NII)
    time.append(t)
    counter += 1 
    if counter == 5000:
        print('Done!')
        break
    #print(t,r1,NI*(NI-1),NI,NII)
#print(NI[-1])
time = np.array(time)
ones = np.array(ones)
plt.plot(time,ones)
plt.show()
