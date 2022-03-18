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

##### This file contains the core kMC algorithms used by ZMC.  ###############################################
'''
The first algorithm below is implemented in ZACROS. 
Basically one generates a queue that contains occurence times of all candidate events. The event with the least occurence time 
is executed first. The occurence time is updated. 
There is a possibility of implementing alternate algorithms that can be considered for future versions. 
'''
class Event:
    def __init__(self,name,sites,time):
        self.name = name
        self.sites = sites
        self.time = time
        
class Copper:
    def __init__(self,lab,state,pos):
        self.label = lab
        self.state = state 
        self.position = pos

def getKeysByValue(dictOfElements, valueToFind):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] == valueToFind:
            listOfKeys.append(item[0])
    return listOfKeys

def create_pairs(Culist):
    out = []
    for i in range(0,len(Culist)):
        for j in range(i+1,len(Culist)):
            out.append([Culist[i],Culist[j]])
    return out

def create_all_pairs(Culist):
    out = []
    for i in range(0,len(Culist)):
        for j in range(0,len(Culist)):
            if j == i:
                continue
            out.append([Culist[i],Culist[j]])
    return out

def distance(p1,p2):
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

def normalized_list(arr):
    #This function transforms the array into a normalized cumulative array which can be used for kMC 
    Norm_sum = sum(arr)
    size = len(arr)
    ret = []
    for i in range(0,size):
        temp = arr[0:i+1]
        ret.append(sum(temp)/Norm_sum)
    return ret

def find_ind(arr,r):
    Sumo = sum(arr)
    Sumi = 0
    size = len(arr)
    if r > 0.5:
        Sum = sum(arr)
        for i in range(size-1,-1,-1):
            if Sum <= r*Sumo:
                #print(Sum, '>', r*Sumo)
                return i + 1
            Sum -= arr[i]
    else:
        for i in range(0,size):
            Sumi += arr[i]
            if Sumi >= r*Sumo:
                #print(Sumi, '>', r*Sumo)
                return i
    return 0

def find_ind_null(arr,r, m):
    arr = [a/m for a in arr]
    Sum = sum(arr)
    Sumi = 0
    size = len(arr)
    i = -1 
    s = 0 
    while s < r:
        i+=1
        if i == size:
            return None
        s+=arr[i]
    return i
    #print('Random', r, 'Sum', Sum)
    '''
    for i in range(size-1,-1,-1):
        if Sum <= r:
            return i
        Sum -= arr[i]
        if Sum == 0:
            return 0
    if r > 0.5:
        Sum = sum(arr)
        for i in range(size-1,-1,-1):
            if Sum <= r:
                #print(Sum, '>', r*Sumo)
                return i + 1
            Sum -= arr[i]
    else:
        for i in range(0,size):
            Sumi += arr[i]
            if Sumi >= r:
                #print(Sumi, '>', r*Sumo)
                return i
    return 0
    '''
    
def delete(distlist,k):
#   Delete a species from the neighbor list
    i = len(distlist)
    while i > 0:
        if distlist[i-1][0] == k or distlist[i-1][1]== k:
            distlist.pop(i-1)
        i-=1
    return distlist

def residual(arr, N):
    pairs = 0
    while(arr):
        #print(distlist)
        ireacted = arr[0][0:2]
        distlist = delete(arr,ireacted[0])
        distlist = delete(arr,ireacted[1])
        pairs += 1   
    #print(distlist)
    #print('Pairs',pairs)
    #print('Pairs',pairs)
    r = (N-2.*pairs)/N
    return r
    

@monkeypatch_class(Zmc)
def read_outputs(self):
    k_B = 8.61733034E-5
    print('The outputs are:','\n')
    print('Pairing prefactor =',self.Ap)
    print('Pressure O2 =',self.P_O2)
    print('Activation factor =',self.Eap)
    print('Temperature =',self.dec_mode)
    print('Dimensions=',self.xl,self.yl,self.zl)
    print('CuI=',self.NCuI) 
    print('CuII=',self.NCuII)
    #print('Rate constant is',self.Ap*self.P_O2*np.exp(-1*self.Eap/(k_B*self.T))*self.expo_function(1.0))

@monkeypatch_class(Zmc)
def expo_function(self,d):
    return self.A_expo*np.exp(-1*self.B_expo*d)

@monkeypatch_class(Zmc)
def sigmo_function(self,d):
    import math 
    m = self.Slo_sigmo 
    #decay = 0.5*(math.erf(3.0699*(self.Cut_sigmo-d)+1))         #This line will play a critical role in determining the decay function. This line will be added to the main ZMC input file
    decay = 1/(1 + np.exp(m*(d-2*self.Cut_sigmo)))
    #print('d, sigmo = ',d,decay)
    if decay == 0:
        return 0  
    else:
        return decay

@monkeypatch_class(Zmc)
def step_function(self,d):
    import numpy
    #ret = np.heaviside(d-self.Cut_step,1E-15)
    #print('d, step = ',d,ret)
    if d > self.Cut_step:
        return 0
    else:
        return 1

@monkeypatch_class(Zmc)
def calc_oxi_escape_time(self,d, i = None):
    k_B = 8.61733034E-5
    if self.dec_mode=='exponential':
        dec = self.expo_function(d) 
    elif self.dec_mode == 'sigmoidal':
        dec = self.sigmo_function(d) 
    elif self.dec_mode == 'step':
        dec = self.step_function(d)
    r = self.Ap*self.P_O2*np.exp(-1*self.Eap/(k_B*self.T))*dec
    if i != None:
        #print('None type')
        if r < 1E-20:
            r = 0
        #print('Hell')
        #print('Oxidation',r)
        return r        #This will return the rate of oxidation instead of escape time. 
    #seed(int(self.rseed))
    else:
        #print('Nope')
        return -1*np.log(1-random())/r
    
@monkeypatch_class(Zmc)
def calc_oxi_rate(self,d):
    k_B = 8.61733034E-5
    if self.dec_mode=='exponential':
        dec = self.expo_function(d) 
    elif self.dec_mode == 'sigmoidal':
        dec = self.sigmo_function(d) 
    elif self.dec_mode == 'step':
        dec = self.step_function(d)
    r = self.Ap*self.P_O2*np.exp(-1*self.Eap/(k_B*self.T))*dec
    return r

@monkeypatch_class(Zmc)
def calc_red_escape_time(self, i = None):
    k_B = 8.61733034E-5
    r = self.Ar*self.P_NO*self.P_NH3*np.exp(-1*self.Ear/(k_B*self.T))
    if i != None:
        if r < 1E-20:
            r = 0 
        #print('Hell')
        #print('Reduction',r,'from',self.Ar,',',self.P_NO,',',self.P_NH3,',',self.Ear)
        return r        #This will return the rate of oxidation instead of escape time. 
    #seed(int(self.rseed))
    return -1*np.log(1-random())/r

@monkeypatch_class(Zmc)
def initialize_event_list(self):
    # Creating dictionaries of the label and oxidation states 
    e_list = []
    self.Cudict = {}
    for i in range(0,self.NCu):
        self.Cudict.update({self.Cu_list[i].label : self.Cu_list[i].state})
    CuI_labels = getKeysByValue(self.Cudict, 1) 
    CuII_labels = getKeysByValue(self.Cudict, 2)
    self.CuI_list = CuI_labels 
    #self.CuII_list = CuII_labels
    if self.kin_mode == 'steady_state':
        #CuI_labels = getKeysByValue(self.Cudict, 1)  
        if len(CuI_labels) > 1:
            pairs = create_pairs(CuI_labels)
        #print(pairs[10][1])
            for i in range(0,len(pairs)):
                p1 = pairs[i][0]
                p2 = pairs[i][1]
                d = distance(self.Cu_list[p1-1].position,self.Cu_list[p2-1].position)
                time = self.calc_oxi_escape_time(d)*len(self.Cu_list)
                e_list.append(Event('oxidation',[p1,p2],time))
        #CuII_labels = getKeysByValue(self.Cudict, 2)
        if len(CuII_labels) > 0:
            if self.red_mode == 'single':
                for i in range(0,len(CuII_labels)):
                    time = self.calc_red_escape_time()
                    e_list.append(Event('reduction',CuII_labels[i],time))
        return e_list
    elif self.kin_mode == 'transient_reduction':
        #CuII_labels = getKeysByValue(self.Cudict, 2)
        if len(CuII_labels) > 0:
            if self.red_mode == 'single':
                for i in range(0,len(CuII_labels)):
                    time = self.calc_red_escape_time()
                    e_list.append(Event('reduction',CuII_labels[i],time))
        return e_list
    elif self.kin_mode == 'transient_oxidation':
        #CuI_labels = getKeysByValue(self.Cudict, 1)  
        if len(CuI_labels) > 1:
            pairs = create_pairs(CuI_labels)
        #print(pairs[10][1])
            for i in range(0,len(pairs)):
                p1 = pairs[i][0]
                p2 = pairs[i][1]
                d = distance(self.Cu_list[p1-1].position,self.Cu_list[p2-1].position)
                time = self.calc_oxi_escape_time(d)*len(self.Cu_list)
                e_list.append(Event('oxidation',[p1,p2],time))
        return e_list 
        
@monkeypatch_class(Zmc)
def calculate(self):
    #self.events = self.initialize_event_list() 
    #self.init_snap()
    #self.init_spec()
    #self.init_proc()
    ###### Now diving into the kMC algorithm ############# 
    if self.kin_mode == 'steady_state' and self.red_mode == 'single':
        print('Shall we ?')
        #self.calc_steady_state_standard_periodic()
        self.calc_steady_state_standard_periodic()
    elif self.kin_mode == 'transient_oxidation':
        print('We shall')
        self.calc_trans_ox()
    else:
        print('Unfortunately, this version has not been implemented yet. Will implement in a latter version')
        exit()
    '''
    elif self.kin_mode == 'transient_reduction' and self.red_mode == 'single':
        self.calc_trans_red()
    elif self.kin_mode == 'transient_oxidation' and self.red_mode == 'single':
        self.calc_trans_oxi()
    '''
   
  
@monkeypatch_class(Zmc)
def calc_steady_state(self):
    t_start = timer()
    #print(self.max_events)  #This is the algorithm of least escape times 
    while True: #self.simevents < self.max_events or self.simutime < self.max_time:
        t_end = timer()
        if (t_end-t_start) > self.wtime:
            break
        '''
        for e in self.events:
            print(e.name, ",", e.sites, ",", e.time)
        '''
        #### First we collect escape times and find the mininum with index ###### 
        exec_time = min([e.time for e in self.events])
        min_ind = [e.time for e in self.events].index(exec_time)
        ##### Update the time and event counter ################################
        self.simutime += exec_time
        self.simevents += 1
        ##### Extract the name of the event and take actions accordingly #######
        if self.events[min_ind].name=='oxidation':
            #print('oxidation')
            [a,b] = self.events[min_ind].sites 
            self.Cu_list[a-1].state = 2 
            self.Cu_list[b-1].state = 2
            self.NCuI -= 2 
            self.NCuII += 2 
            self.n_O2 -= 1
            self.oxievents += 1 
            self.events.pop(min_ind) 
            elim_ind = []
            #print(a,b)
            for e in self.events:   
                if e.name == 'oxidation':
                    if a in e.sites or b in e.sites:
                        elim_ind.append(self.events.index(e))
            #print(elim_ind)
            self.events = [self.events[x] for x in range(0,len(self.events)) if x not in elim_ind]
            self.refresh_reduction(a)
            self.refresh_reduction(b)
        elif self.events[min_ind].name=='reduction':
            #print('reduction')
            c = self.events[min_ind].sites
            #print(c)
            #self.Cu_list[c-1].state = 1 
            self.NCuII -= 1 
            self.NCuI += 1 
            self.n_NO -= 2
            self.n_NH3 -= 2 
            self.redevents += 1 
            self.Cu_list[c-1].state = 1 
            self.events.pop(min_ind)
            self.refresh_oxidation(c)    
        if self.simevents == self.max_events:
            break
        #self.events = self.initialize_event_list()       #Generate a fresh events list with updated states as alternative 1. 
        print(self.simutime,"NCuI",self.NCuI)
        #print(self.simutime, " , ", self.n_NO)
        

@monkeypatch_class(Zmc)
def calc_steady_state_standard(self):
    t_start = timer()
    self.initialize_rate_list()
    #print(self.max_events)  #This is the standard kMC algorithm. 
    while True: #self.simevents < self.max_events or self.simutime < self.max_time:
        t_end = timer()
        if (t_end-t_start) > self.wtime:
            break       
        #### First we generate a random number and get the index of the executable event ###### 
        ind = find_ind(self.rate_nums_norm,random())
        #### Subsequently, we advance the simulation clock time by the rate of the event
        try:
            exec_time = random()/self.rate_nums[ind]
        except:
            print('No more oxidations!')
            break                      #This is a contingency for the case when all ions are Copper 1's and are unable to be oxidized. 
        ##### Update the time and event counter ####################################
        self.simutime += exec_time
        self.simevents += 1
        sites = self.site_list[ind]  #<----This refers to the sites where the event will occur
        #print(exec_time, sites)
        #### Extract the event name and take the actions accordingly ###########
        if sites[1]==0:
            #This is a reduction event  
            red_site = sites[0]
            self.Cu_list[red_site-1].state = 1 
            self.NCuII -= 1 
            self.NCuI += 1 
            self.n_NO -= 2
            self.n_NH3 -= 2 
            self.redevents += 1 
            #self.Cu_list[c-1].state = 1 
        else:
            #This is an oxidation event 
            [a,b] = sites 
            self.Cu_list[a-1].state = 2 
            self.Cu_list[b-1].state = 2
            self.NCuI -= 2 
            self.NCuII += 2 
            self.n_O2 -= 1
            self.oxievents += 1 
        #print('Refresh time:',self.NCuII)
        self.refresh_rate_list()
        print(self.simutime,self.NCuII)
        #print(self.simutime, self.NCuII)
        if self.simevents == self.max_events:
            print('Events done')
            break
        '''
        if self.simevents == 3:
            print(exec_time,",",self.site_list[ind])
            break
        '''

@monkeypatch_class(Zmc)
def calc_steady_state_standard_periodic(self):
    import random
    t_start = timer()
    self.minimum_image_convention()
    random.seed(int(self.rseed))
    #print('Net coppers:',len(self.Cu_list))
    #print(self.max_events)  #This is the standard kMC algorithm.
    print('Volume=',self.vol)
    if self.seed_mode == 'CHA':
        from ase.io import read 
        #atoms = read('/afs/crc.nd.edu/user/a/agoswami/zmc/zmc/POSCAR_CHA')
        atoms = read('/afs/crc.nd.edu/user/a/agoswami/zmc/zmc/POSCAR_CHA')
        del atoms[[atom.index for atom in atoms if atom.symbol=='O']]           # Delete all O atoms 
        atomsrep = atoms.repeat([self.al,self.bl,self.cl])
    time = [0]
    ones = [self.NCuI]
    twos = [self.NCuII]
    events = [0]
    oxidevents = [0]
    redevents = [0]
    first = []
    for c in self.Cu_list:
        first.append(c.state)
    self.positions = [first]
    self.snap_times = [0]
    self.snap_events = [0]
    n_int = 1
    n_int_snap = 1
    simtime = [] ; os = [] ; distances =  [] 
    print('The label is', self.spec_int.label)
    while True: #self.simevents < self.max_events or self.simutime < self.max_time:
        t_end = timer()
        if (t_end-t_start) > self.wtime:
            break     
        #### If there are no pairable coppers present ########################################
        if len(self.rate_dict) == 0:
            print('No more pairable coppers here')
            break
        #### First we generate a random number and get the index of the executable event ###### 
        r1 = random.uniform(0,1)
        #print(self.site_list)
        #print(self.rate_nums)
        ind = find_ind(self.rate_nums,r1)
        #t_imd = timer() ; print('Time taken is',t_imd-t_start)
        #print('Found index at:', timer()-t_start)
        #print(ind)
        #print(r1,ind)
        if ind == None:
            #print(self.site_list)
            #print(r1)
            print('No more pairable coppers present')
            break
        #### Subsequently, we advance the simulation clock time by the rate of the event
        if self.rate_nums[ind] != 0:
            r2 = random.uniform(0,1)
            exec_time = -1*np.log(r2)/sum(self.rate_nums)
        else:
            print('No more oxidations!')
            break    #This is a contingency for the case when all ions are Copper 1's and are unable to be oxidized.  
        ##### Update the time and event counter ####################################
        self.simutime += exec_time
        self.simevents += 1
        sites = list(self.rate_dict)[ind]
        self.sites = sites
        self.ind = ind
        #sites = self.site_list[ind]  #<----This refers to the sites where the event will occur
        #print(sites)
        #print(exec_time, sites)
        '''
        for e in self.Cu_list:
            print(e.label, ",", e.state)
        '''
        '''
        ########## Applying a random selection algorithm here ##################################
        sites = random.choice(self.site_list)
        r = random.uniform(0,1)
        exec_time = -1*np.log(r)/sum(self.rate_nums)
        self.simutime += exec_time
        self.simevents += 1
        '''
        #### Extract the event name and take the actions accordingly ##########################
        if sites[1]==0:
            #This is a reduction event  
            red_site = sites[0]
            self.Cu_list[red_site-1].state = 1 
            self.NCuII -= 1 
            self.NCuI += 1 
            self.n_NO -= 1
            self.n_NH3 -= 1 
            self.redevents += 1 
            self.frequency[red_site].append(-1*self.simevents)
            #self.refresh_rlist(red_site, 0, ind)
            #print(red_site)
            #self.Cu_list[c-1].state = 1 
        else:
            #This is an oxidation event 
            [a,b] = sites 
            self.Cu_list[a-1].state = 2 
            self.Cu_list[b-1].state = 2
            self.NCuI -= 2 
            self.NCuII += 2 
            self.n_O2 -= 1
            self.oxievents += 1 
            self.frequency[a].append(1*self.simevents)  ; self.frequency[b].append(1*self.simevents)
            if self.statistics == True:
                if self.seed_mode == 'CHA':
                    #distances.append(self.frac_mic_distance(self.Cu_list[a-1].position, self.Cu_list[b-1].position))
                    distances.append(atomsrep.get_distance(self.p2ind[tuple(self.Cu_list[a-1].position)],self.p2ind[tuple(self.Cu_list[b-1].position)],mic = True))
                else:
                    distances.append(self.mic_distance(self.Cu_list[a-1].position, self.Cu_list[b-1].position))
            #self.refresh_rlist(a, b, ind)
            #print(a,b)
        #print('Refresh time:',self.NCuII)
        '''
        for e in self.Cu_list:
            print(e.label, ",", e.state)
        '''
        #Going for something daring here. When 100 events are up, the time interval of recording events/species will be the next highest multiple of 10 of simulation time. 
        if self.spec_int.label == 'time':
            if self.simevents == 100:
                b = math.floor(math.log10(self.simutime)) + 1
                ones = [] ; time = [] ; twos = [] ; events = [] ; oxidevents = [] ; redevents = [] ; n_int = 0
                self.spec_int.interval = 10**b 
                self.max_time = 100*self.spec_int.interval
                print('The simulation time:',self.simutime,'for which species interval:',self.spec_int.interval)
        temp_pos = []
        if self.snap_int.label == 'event':
            if self.simevents >= n_int_snap*self.snap_int.interval:    
                while self.simevents >= n_int_snap*self.snap_int.interval:
                    for c in self.Cu_list:
                        temp_pos.append(c.state)
                    n_int_snap += 1 
                self.snap_events.append(self.simevents)
                self.snap_times.append(self.simutime)
                self.positions.append(temp_pos)
        else:
            if self.simutime >= n_int_snap*self.snap_int.interval:
                while self.simutime >= n_int_snap*self.snap_int.interval:
                    for c in self.Cu_list:
                        temp_pos.append(c.state)
                    n_int_snap += 1
                self.snap_events.append(self.simevents)
                self.snap_times.append(self.simutime)
                self.positions.append(temp_pos)
        if self.spec_int.label == 'event':
            if self.simevents >= n_int*self.spec_int.interval:
                ones.append(self.NCuI)
                time.append(self.simutime)
                twos.append(self.NCuII)
                events.append(self.simevents)
                oxidevents.append(self.oxievents)
                redevents.append(self.redevents)
                n_int += 1
                print(n_int, timer()-t_start)
        else:
            if self.simutime >= n_int*self.spec_int.interval:
                while self.simutime >= n_int*self.spec_int.interval:
                    ones.append(self.NCuI)
                    time.append(n_int*self.spec_int.interval)                 #time.append(self.simutime)
                    twos.append(self.NCuII)
                    events.append(self.simevents)
                    oxidevents.append(self.oxievents)
                    redevents.append(self.redevents)
                    n_int += 1       
                    print(n_int, timer()-t_start)
        #self.refresh_mic()
        self.refresh_rlist()
        #print(r2,self.simutime,self.NCuI)
        #print('This',self.simevents, self.NCuI, timer()-t_start)
        #print(self.simevents)
        simtime.append(self.simutime)
        os.append(self.NCuI/self.NCu)
        if self.simevents == self.max_events or self.simutime >= self.max_time:
            #print(self.rate_nums)
            print('Events done')
            break
        '''
        if self.simevents == 3:
            print(exec_time,",",self.site_list[ind])
            break
         '''
    print('The total time taken is:',self.simutime)
    wall_time = timer() - t_start 
    f = open('raw_outputs.txt','w')
    f.write('Events           Time                      CuI                   CuII\n')
    for a,b,c,d in zip(time,ones,twos,events):
        f.write(' ' +str('{:03d}'.format(d)))
        f.write('           ')
        f.write(str('{:0.12f}'.format(a)))
        f.write('                   ')
        f.write(str('{:03d}'.format(b)))
        f.write('                   ')
        f.write(str('{:03d}'.format(c)))
        f.write('\n')
    f.close()
    f = open('event_outputs.txt','w')
    f.write('Events           Time                      Oxi_events                   Red_events\n')
    for a,b,c,d in zip(time,oxidevents,redevents,events):
        f.write(' ' +str('{:03d}'.format(d)))
        f.write('           ')
        f.write(str('{:0.5f}'.format(a)))
        f.write('                      ')
        f.write(str('{:03d}'.format(b)))
        f.write('                          ')
        f.write(str('{:03d}'.format(c)))
        f.write('\n')
    f.close()
    if self.statistics == True:
        f = open('Distance_statistics.txt','w')
        for d in distances:
            f.write(str(d))
            f.write('     ')
        f.close()
        f = open('Minimum_distance_frequencies.txt','w')
        for k in self.min_distance.keys():
            f.write('{0}                     {1}                   {2}\n'.format(k, self.indices[k-1], self.frequency[k])) # self.min_distance[k]
        f.close()
    self.t = time 
    self.o = ones 
    self.tw = twos
    self.ev = events 
    self.oev = oxidevents
    self.rev = redevents
    f.close()
    f = open('simulation_outputs.txt','w')
    f.write('The simulation time taken by the code is {0} seconds\n'.format(self.simutime))
    f.write('The total calculation time taken by the code is {0} seconds\n'.format(wall_time))
    f.write('The total kMC events simulated by the code is {0}'.format(self.simevents))
    f.close()
    #print(self.snap_times)
    #print(self.positions)
    #plt.plot(simtime,os)
    #plt.show()

    
@monkeypatch_class(Zmc)
def initialize_rate_list(self):
    self.Cudict = {}
    #self.r_list = {}         # Initialize a dictionary that maps site of events ----> rate of reaction
    self.site_list = []
    self.rate_nums = []
    print('Start')
    for i in range(0,self.NCu):
        self.Cudict.update({self.Cu_list[i].label : self.Cu_list[i].state})
    CuI_labels = getKeysByValue(self.Cudict, 1)
    CuII_labels = getKeysByValue(self.Cudict, 2)
    pairs = create_all_pairs(CuI_labels)
    for p in pairs:
        p1 = p[0]
        p2 = p[1]
        d = distance(self.Cu_list[p1-1].position,self.Cu_list[p2-1].position)
        rate = self.calc_oxi_escape_time(d,1)
        self.rate_nums.append(rate)
        self.site_list.append([p1,p2])
        #self.r_list[(p1,p2)] = rate 
        #self.r_list.update({[p1,p2]:rate})    #       [(p1,p2),rate]
    for l in CuII_labels:
        rate = self.calc_red_escape_time(1)
        self.rate_nums.append(rate)
        self.site_list.append([l,0])
        #self.r_list[(l,0)] = rate
    print('Finish')
    self.rate_nums_norm = normalized_list(self.rate_nums)   #Returns a cumulative normalized array suitable for comparing with kMC
    print('Super finish')

@monkeypatch_class(Zmc)
def refresh_rate_list(self):
    self.Cudict = {}
    #self.r_list = {}         # Initialize a dictionary that maps site of events ----> rate of reaction
    self.site_list = []
    self.rate_nums = []
    for i in range(0,self.NCu):
        #print(self.Cu_list[i].label,self.Cu_list[i].state)
        self.Cudict.update({self.Cu_list[i].label : self.Cu_list[i].state})
    CuI_labels = getKeysByValue(self.Cudict, 1)
    #CuII_labels = getKeysByValue(self.Cudict, 2)
    for c in self.Cu_list:
        #print(c.label, c.state)
        if c.state == 1:        #This is copper 1 
            for j in CuI_labels:
                if j == c.label:
                    continue
                d = distance(self.Cu_list[j-1].position,c.position)
                self.rate_nums.append(self.calc_oxi_escape_time(d,1))
                self.site_list.append([c.label,j])
        else:                   #This is copper 2 
            self.rate_nums.append(self.calc_red_escape_time(1))
            self.site_list.append([c.label,0])
    #print(self.site_list)
    self.rate_nums_norm = normalized_list(self.rate_nums)   #Returns a cumulative normalized array suitable for comparing with kMC
    #print(self.rate_nums_norm)
    
@monkeypatch_class(Zmc)
def period_image_Cu(self):
    self.per_list = []
    for c in self.Cu_list:
        [x,y,z] = c.position 
        points = [[0,0,0],[self.xl,0,0],[0,self.yl,0],[0,0,self.zl],[self.xl,self.yl,0],[0,self.yl,self.zl],[self.xl,0,self.zl],[self.xl,self.yl,self.zl],[0,0,z],[self.xl,0,z],[0,self.yl,z],[self.xl,self.yl,z],[0,y,0],[self.xl,y,0],[0,y,self.zl],[self.xl,y,self.zl],[x,0,0],[x,0,self.zl],[x,self.yl,0],[x,self.yl,self.zl],[0,y,z],[self.xl,y,z],[x,0,z],[x,self.yl,z],[x,y,0],[x,y,self.zl]]
        displace = [[1,1,1],[-1,1,1],[1,-1,1],[1,1,-1],[-1,-1,1],[1,-1,-1],[-1,1,-1],[-1,-1,-1],[1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0],[1,0,1],[-1,0,1],[1,0,-1],[-1,0,-1],[0,1,1],[0,1,-1],[0,-1,1],[0,-1,-1],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
        index = 0 
        for p in points:
            d = distance(c.position,p)
            if d < self.cutoff:
                xnew = x + displace[index][0]*self.xl
                ynew = y + displace[index][1]*self.yl
                znew = z + displace[index][2]*self.zl
                self.per_list.append(Copper(lab = c.label,state=c.state,pos=[xnew,ynew,znew]))
            index += 1
    self.per_Cu_list = self.Cu_list + self.per_list
    
@monkeypatch_class(Zmc)
def initiate_pair_rate_list(self):
    if self.dec_mode != 'sigmoidal' and self.dec_mode != 'step':
        print('Pairlist function can only work for sigmoidal or step dependences')
        exit()
    self.site_list = []
    self.rate_nums = []
    if self.dec_mode == 'sigmoidal':
        self.cutoff = self.Cut_sigmo
    else:
        self.cutoff = self.Cut_step
    self.per_Cu_list = [] 
    self.pairs = {}
    self.period_image_Cu()
    for c in self.Cu_list:
        nearest_neigh = {}               #A list that will keep track of nearest neighbours for a given Cu 
        for u in self.per_Cu_list:
            d = distance(c.position,u.position)
            if d == 0:
                continue
            if d < self.cutoff:
                try:
                    if d < nearest_neigh[u.label]:
                        nearest_neigh[u.label] = d
                except:
                    nearest_neigh[u.label] = d 
                '''
                if nearest_neigh[u.label] == None:
                    nearest_neigh[u.label] = d  
                else:
                '''   
        nn_new = {}                                                   # Reminder: This method does not sort the distances themselves. 
        for n in sorted(nearest_neigh):
            nn_new.update({n:nearest_neigh[n]})              #Sorts all the pairs in ascending order of the site numbers
        #nearest_neigh = list(set(nearest_neigh))
        #nearest_neigh.sort()
        self.pairs[c.label] = nn_new     
    for i in range(0,len(self.Cu_list)):
        nn = self.pairs[i+1]
        if len(nn)==0:
            continue
        for n in nn:
            if i+1 not in self.pairs[n]:
                print(i+1,'is not self-consistent')
                break
    print(self.pairs)
    #print(self.pairs[1][3],'==',self.pairs[3][1])
    print('Full power consistent!')
    print(len(self.Cu_list))
    for i in range(0,len(self.pairs)):
        if self.Cu_list[i].state == 1:      #If the copper is copper 1. 
            if len(self.pairs[i+1])==0:
                pass
            else:
                for p in self.pairs[i+1]:
                    if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                        continue
                    else:
                        '''
                        if p < i+1:
                            continue
                        '''
                        d = self.pairs[i+1][p]
                        self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))
                        self.site_list.append([i+1,p])
        elif self.Cu_list[i].state == 2:
            self.rate_nums.append(self.calc_red_escape_time(1))
            self.site_list.append([i+1,0])
    #print(self.site_list)
  
@monkeypatch_class(Zmc)            
def refresh_pair_rate_list(self):
   ### Now to refresh the rate list again. 
    self.site_list = []
    self.rate_nums = []
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    for i in range(0,len(self.pairs)):
        if self.Cu_list[i].state == 1:      #If the copper is copper 1. 
            if len(self.pairs[i+1])==0:
                pass
            else:
                for p in self.pairs[i+1]:
                    #print(p)
                    if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                        continue
                    else:
                        '''
                        if p < i+1:
                            continue
                        '''
                        d = self.pairs[i+1][p]
                        #oxid_rates.append(self.calc_oxi_escape_time(d,1))
                        #oxid_sites.append([i+1,p])
                        self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))
                        self.site_list.append([i+1,p])
        elif self.Cu_list[i].state == 2:
            #redu_rates.append(self.calc_red_escape_time(1))
            #redu_sites.append([i+1,0])                                   #<----Uncomment lines 563,564,568,569,572,573 for diff algo
            self.rate_nums.append(self.calc_red_escape_time(1))
            self.site_list.append([i+1,0])
    #self.rate_nums =  redu_rates + oxid_rates 
    #self.site_list =  redu_sites + oxid_sites 
    #print(self.site_list)   
       
@monkeypatch_class(Zmc)
def refresh_reduction(self,Z):
    ### This function adds new potential reduction events to the event list. 
    self.events.append(Event('reduction',Z,self.calc_red_escape_time()))
    self.CuI_list.remove(Z)
    
@monkeypatch_class(Zmc)
def refresh_oxidation(self,X):
    ### This function adds new potential oxidation events to the event list. 
    for c in self.CuI_list:
        d = distance(self.Cu_list[X-1].position,self.Cu_list[c-1].position)
        if X > c:
            self.events.append(Event('oxidation',[c,X],self.calc_oxi_escape_time(d)*len(self.Cu_list)))
        else:
            self.events.append(Event('oxidation',[X,c],self.calc_oxi_escape_time(d)*len(self.Cu_list)))
    self.CuI_list.append(X)
   
@monkeypatch_class(Zmc)
def mic_distance(self,a,b):
    a = np.array(a)
    b = np.array(b)
    #x1 = a[0] 
    #y1 = a[1]
    #z1 = a[2]
    #x2 = b[0]
    #y2 = b[1]
    #z2 = b[2]
    #xnew = x1 - x2 ; xnew = xnew - self.xl*np.rint(xnew/self.xl)
    #ynew = y1 - y2 ; ynew = ynew - self.yl*np.rint(ynew/self.yl)
    #znew = z1 - z2 ; znew = znew - self.zl*np.rint(znew/self.zl)
    v = a - b; v = v - self.xl * np.rint(v/self.xl)
    '''
    xnew = abs(x1-x2) - math.floor(abs(x1-x2)/(0.5*self.xl))*self.xl
    ynew = abs(y1-y2) - math.floor(abs(y1-y2)/(0.5*self.yl))*self.yl
    znew = abs(z1-z2) - math.floor(abs(z1-z2)/(0.5*self.zl))*self.zl
    '''
    '''
    dold = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    dnew = np.sqrt((xnew)**2 + (ynew)**2 + (znew)**2)
    print('Old: ',dold,' vs New: ',dnew)
    '''
    return np.linalg.norm(v)
    
@monkeypatch_class(Zmc)
def frac_mic_distance(self,f1,f2):
    #print('#->',f1,'<-#')
    #print('#->',f2,'<-#')
    [ax,ay,az] = [13.6750000000000007 ,   0.0000000000000000   , 0.0000000000000000]
    [bx,by,bz] = [-6.8374999999999977 ,  11.8428973967521998   , 0.0000000000000000]
    [cx,cy,cz] = [0.0000000000000000  ,  0.0000000000000000  , 14.7669999999999995]
    adist = abs(f1[0] - f2[0]) - int(2*abs(f1[0] - f2[0])) 
    bdist = abs(f1[1] - f2[1]) - int(2*abs(f1[1] - f2[1]))
    cdist = abs(f1[2] - f2[2]) - int(2*abs(f1[2] - f2[2]))
    #print(xdist,ydist,zdist)
    delx = self.al*adist*ax + self.bl*bdist*bx + self.cl*cdist*cx 
    dely = self.al*adist*ay + self.bl*bdist*by + self.cl*cdist*cy 
    delz = self.al*adist*az + self.bl*bdist*bz + self.cl*cdist*cz 
    return np.sqrt((delx)**2 + (dely)**2 + (delz)**2)

@monkeypatch_class(Zmc)
def minimum_image_convention(self):    #Since the Coppers wont be moving, the only distances that need to be measured is from MIC. 
    #First we compute all possible pair distances 
    self.site_list = []
    self.rate_nums = []
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    self.pairs = {}
    if self.seed_mode == 'CHA':
        from ase.io import read 
        #atoms = read('/afs/crc.nd.edu/user/a/agoswami/zmc/zmc/POSCAR_CHA')
        atoms = read('/afs/crc.nd.edu/user/a/agoswami/zmc/zmc/POSCAR_CHA')
        del atoms[[atom.index for atom in atoms if atom.symbol=='O']]           # Delete all O atoms 
        atomsrep = atoms.repeat([self.al,self.bl,self.cl])
        self.vol = atomsrep.get_volume()
    else:
        self.vol = self.xl*self.yl*self.zl
    #self.period_image_Cu()
    #self.mic_parellel()
    max_distance = 0 
    min_distance = 1000
    d_list = [] 
    for c in self.Cu_list:
        #print(c.label)
        nearest_neigh = {} 
        dees = []
        for u in self.Cu_list:             #A list that will keep track of nearest neighbours for a given Cu 
            if self.seed_mode == 'CHA':
                #d = self.mic_cha_distance(c.position,u.position)
                d = atomsrep.get_distance(self.p2ind[tuple(c.position)],self.p2ind[tuple(u.position)],mic = True)
            else:
                d = self.mic_distance(c.position,u.position)
            #print(d)
            if u.label==c.label:
                continue
            if d < min_distance:
                min_distance = d 
            if d > max_distance:
                max_distance = d     
            nearest_neigh[u.label] = d
            dees.append(d)
            if u.label > c.label:
                d_list.append(d)
            #print(c.label, u.label, d)
        #Sorts all the pairs in ascending order of the site numbers
        #nearest_neigh = list(set(nearest_neigh))
        #nearest_neigh.sort()
        self.pairs[c.label] = nearest_neigh
        self.min_distance[c.label] = min(dees)
        self.frequency[c.label] = []                       # This initializes the dictionary of frequencies of every ion which can pair if in +1 state. 
    #print(nearest_neigh)
    #print(self.pairs)
    print('The minimum distance is: ',min_distance)
    print('The maximum distance is: ',max_distance) 
    if self.statistics == True:
        f = open('Static_distances.txt','w')
        for d in d_list:
            f.write(str(d))
            f.write('     ')
        f.close()
        '''
        f = open('Minimum_distances.txt','w')
        for k in self.min_distance.keys():
            f.write('{0}                     {1}\n'.format(k, self.avg_distance[k]))
        f.close()
        '''
    self.rate_dict = {}
    #self.oxi_dict = {}
    #self.red_dict = {}
    for i in range(0,len(self.pairs)):
        if self.Cu_list[i].state == 1:                #If the copper is copper 1. 
            for p in self.pairs[i+1]:
                if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                    continue
                else:
                    '''
                    if p < i+1:
                        continue
                    '''
                    d = self.pairs[i+1][p]
                    if self.calc_oxi_rate(d) != 0:
                        #self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))   #<-- This factor after / is to be added.
                        self.site_list.append([i+1,p])
                        self.rate_dict.update({tuple([i+1,p]):self.calc_oxi_escape_time(d,1)/self.vol})
                        #self.rate_dict[self.calc_oxi_escape_time(d,1)/self.vol] = [i+1,p]
                        '''
                        oxid_rates.append(self.calc_oxi_escape_time(d,1)/self.vol)
                        oxid_sites.append([i+1,p])
                        '''
                        #print(d)
        elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            self.site_list.append([i+1,0])
            self.rate_dict.update({tuple([i+1,0]):self.calc_red_escape_time(1)})
            #self.rate_dict[self.calc_red_escape_time(1)] = [i+1,0]
            '''
            redu_rates.append(self.calc_red_escape_time(1))
            redu_sites.append([i+1,0])
            '''
    #print('The rate_dict is:',self.rate_dict)
    #self.rate_dict.update(self.oxi_dict)
    #self.rate_dict.update(self.red_dict)
    self.rate_nums = list(self.rate_dict.values())
    '''
    self.rate_nums =  redu_rates + oxid_rates 
    self.site_list =   redu_sites + oxid_sites
    self.redlength = len(redu_rates)
    self.oxilength = len(oxid_rates)
    '''
    #print(self.oxilength)
    #print(self.site_list)
    #print(self.rate_nums)
    #print(self.oxilength)
    #print('Redu_size check',len(redu_rates),len(redu_sites))
    #print('Oxid_size check',len(oxid_rates),len(oxid_sites))
    #print(self.site_list)
          
@monkeypatch_class(Zmc)
def refresh_mic(self): 
    print(self.rate_dict)
    self.rate_dict = {}
    for i in range(0,len(self.pairs)):
        if self.Cu_list[i].state == 1:                #If the copper is copper 1. 
            for p in self.pairs[i+1]:
                if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                    continue
                else:
                    '''
                    if p < i+1:
                        continue
                    '''
                    d = self.pairs[i+1][p]
                    if self.calc_oxi_rate(d) != 0:
                        #self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))   #<-- This factor after / is to be added.
                        #self.site_list.append([i+1,p])
                        self.rate_dict.update({tuple([i+1,p]):self.calc_oxi_escape_time(d,1)/self.vol})
                        self.rate_dict.update({tuple([p,i+1]):self.calc_oxi_escape_time(d,1)/self.vol})
                        #self.rate_dict[self.calc_oxi_escape_time(d,1)/self.vol] = [i+1,p]
                        '''
                        oxid_rates.append(self.calc_oxi_escape_time(d,1)/self.vol)
                        oxid_sites.append([i+1,p])
                        '''
                        #print(d)
        elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            #self.site_list.append([i+1,0])
            self.rate_dict.update({tuple([i+1,0]):self.calc_red_escape_time(1)})
            #self.rate_dict[self.calc_red_escape_time(1)] = [i+1,0]
            '''
            redu_rates.append(self.calc_red_escape_time(1))
            redu_sites.append([i+1,0])
            '''
    #print('The rate_dict is:',self.rate_dict)
    self.rate_nums = list(self.rate_dict.values())
    '''
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    self.site_list = []
    self.rate_nums = []
    self.vol = self.xl*self.yl*self.zl
    for i in range(0,len(self.pairs)):
        if self.Cu_list[i].state == 1:                #If the copper is copper 1. 
            for p in self.pairs[i+1]:
                if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                    continue
                else:
                    
                    if p < i+1:
                        continue
                   
                    d = self.pairs[i+1][p]
                    if self.calc_oxi_escape_time(d,1) != 0:
                        #self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))
                        #self.site_list.append([i+1,p])
                        oxid_rates.append(self.calc_oxi_escape_time(d,1)/self.vol)
                        oxid_sites.append([i+1,p])
        elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            #self.site_list.append([i+1,0])
            redu_rates.append(self.calc_red_escape_time(1))
            redu_sites.append([i+1,0])
    self.rate_nums =   redu_rates + oxid_rates 
    self.site_list =   redu_sites + oxid_sites
    #print('Redu_size check',len(redu_rates),len(redu_sites))
    #print('Oxid_size check',len(oxid_rates),len(oxid_sites))
    #print(self.site_list)
    '''

@monkeypatch_class(Zmc)
def refresh_rlist(self):
    #print(self.rate_dict)
    #print(self.sites)
    if self.sites[1]!=0:
        #Heads up: This is oxidation land 
        [a,b] = self.sites 
        for x in list(self.rate_dict):
            if a in x or b in x:
                del self.rate_dict[x]
        self.rate_dict.update({tuple([a,0]):self.calc_red_escape_time(1)})   #Updating the dictionary with reduction steps
        self.rate_dict.update({tuple([b,0]):self.calc_red_escape_time(1)})
    else:
        #Heads up: This is reduction land 
        a = self.sites[0]
        del self.rate_dict[(a,0)]
        for p in self.pairs[a]:
            if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                continue
            else:
                '''
                if p < i+1:
                    continue
                '''
                d = self.pairs[a][p]
                if self.calc_oxi_escape_time(d,1) != 0:               #Updating the dictionary with oxidation steps 
                    self.rate_dict.update({tuple([a,p]):self.calc_oxi_escape_time(d,1)/self.vol})
                    self.rate_dict.update({tuple([p,a]):self.calc_oxi_escape_time(d,1)/self.vol})
    self.rate_nums = list(self.rate_dict.values())
    #print(self.rate_dict)
'''       
Notes: In the matlab Zeolite_modelling.m file, the microkinetic results are obtained. The final value of the coverage will depend on whether I treat the ODE as dx/dt = k1*(1-x) - k2*x*x or dx/dt = k1*(1-x) - 2*k2*x*x. The first case will solve the microkinetic equation for all unique Copper 1 pairs in the box. As in (1,2),(1,3)...(2,3),(2,4)... It wont consider (2,1) in the pair list. The latter case will solve the microkinetic equation for all 2 BODY COMBINATIONS. Basically (2,1) will also be included.
For the first case, we will include the lines of code if p < i+1: continue in minimum_image_convention and refresh_mic functions above. For the second case, we will OMIT those lines of code. 
EDIT: 16th June. The rate of pairing is to be normalized relative to VOLUME OF THE SUPERCELL, NOT THE TOTAL COPPERS. 
NOTES: 27/06/2020: I recognize that the code is taking so much time to run especially in constructing the pair lists in between events. That can be optimized easily by changing how to maintain the pair lists. Basically, I construct the list as l = oxid + red. For an oxidation event, I take the two I sites, pop the events in the list and add the 2 new sites into the reduction list. 
For a reduction event, I take the one site, pop the reduction event, add all potential oxidation pairs. Let's see if this can be done. 
'''

@monkeypatch_class(Zmc)
def calc_trans_ox(self):
    import random
    t_start = timer()
    print('Here')
    self.minimum_image_convention()
    random.seed(int(self.rseed))
    #print('Net coppers:',len(self.Cu_list))
    #print(self.max_events)  #This is the standard kMC algorithm.
    time = [0]
    ones = [self.NCuI/self.NCu]
    twos = [self.NCuII/self.NCu]
    events = [0]
    oxidevents = [0]
    redevents = [0]
    n_int = 1
    #print(self.site_list)
    while True: #self.simevents < self.max_events or self.simutime < self.max_time:
        t_end = timer()
        if (t_end-t_start) > self.wtime:
            break       
        #### First we generate a random number and get the index of the executable event ###### 
        r1 = random.uniform(0,1)
        #print(self.site_list)
        #print(self.rate_nums)
        ind = find_ind(self.rate_nums,r1)
        #print('Found index at:', timer()-t_start)
        #print(ind)
        #print(r1,ind)
        if len(self.rate_nums)==0: # or len(self.rate_nums)==1:
            print('All coppers oxidized!')
            break
        '''
        if ind == None:
            #print(self.site_list)
            #print(r1)
            print('No more pairable coppers present')
            break
        '''
        #### Subsequently, we advance the simulation clock time by the rate of the event
        if self.rate_nums[ind] != 0:
            r2 = random.uniform(0,1)
            exec_time = -1*np.log(r2)/sum(self.rate_nums)
        else:
            print('No more oxidations!')
            break    #This is a contingency for the case when all ions are Copper 1's and are unable to be oxidized.  
        ##### Update the time and event counter ####################################
        self.simutime += exec_time
        self.simevents += 1
        sites = self.site_list[ind]  #<----This refers to the sites where the event will occur
        #print(sites)
        #print(exec_time, sites)
        '''
        for e in self.Cu_list:
            print(e.label, ",", e.state)
        '''
        '''
        ########## Applying a random selection algorithm here ##################################
        sites = random.choice(self.site_list)
        r = random.uniform(0,1)
        exec_time = -1*np.log(r)/sum(self.rate_nums)
        self.simutime += exec_time
        self.simevents += 1
        '''
        #### Extract the event name and take the actions accordingly ##########################
        if sites[1]==0:
            #This is a reduction event  
            red_site = sites[0]
            self.Cu_list[red_site-1].state = 1 
            self.NCuII -= 1 
            self.NCuI += 1 
            self.n_NO -= 1
            self.n_NH3 -= 1 
            self.redevents += 1 
            #self.refresh_rlist(red_site, 0, ind)
            #print(red_site)
            #self.Cu_list[c-1].state = 1 
        else:
            #This is an oxidation event 
            [a,b] = sites 
            self.Cu_list[a-1].state = 2 
            self.Cu_list[b-1].state = 2
            self.NCuI -= 2 
            self.NCuII += 2 
            self.n_O2 -= 1
            self.oxievents += 1 
            #self.refresh_rlist(a, b, ind)
            #print(a,b)
        #print('Refresh time:',self.NCuII)
        '''
        for e in self.Cu_list:
            print(e.label, ",", e.state)
        '''
        #print(r2,self.simutime,self.NCuI)
        if self.spec_int.label == 'event':
            if self.simevents >= 1: #n_int*self.spec_int.interval:
                ones.append(self.NCuI)
                time.append(self.simutime)
                twos.append(self.NCuII)
                events.append(self.simevents)
                oxidevents.append(self.oxievents)
                redevents.append(self.redevents)
                n_int += 1
        else:
            if self.simutime >= n_int*self.spec_int.interval:
                ones.append(self.NCuI)
                time.append(n_int*self.spec_int.interval)                 #time.append(self.simutime)
                twos.append(self.NCuII)
                events.append(self.simevents)
                oxidevents.append(self.oxievents)
                redevents.append(self.redevents)
                n_int += 1        
        self.refresh_mic_TOX()
        print('This',self.simutime, self.NCuI, timer()-t_start)
        if self.simevents == self.max_events:
            #print(self.rate_nums)
            print('Events done')
            break
        '''
        if self.simevents == 3:
            print(exec_time,",",self.site_list[ind])
            break
         '''
    f = open('raw_outputs_TOX.txt','w')
    f.write('Events           Time                      CuI                   CuII\n')
    for a,b,c,d in zip(time,ones,twos,events):
        f.write(' ' +str('{:03d}'.format(d)))
        f.write('           ')
        f.write(str('{:0.3e}'.format(a)))
        f.write('                   ')
        f.write(str('{:0.3f}'.format(b)))
        f.write('                   ')
        f.write(str('{:0.3f}'.format(c)))
        f.write('\n')
    f.close()
    f = open('event_outputs_TOX.txt','w')
    f.write('Events           Time                      Oxi_events                   Red_events\n')
    for a,b,c,d in zip(time,oxidevents,redevents,events):
        f.write(' ' +str('{:03d}'.format(d)))
        f.write('           ')
        f.write(str('{:0.3e}'.format(a)))
        f.write('                      ')
        f.write(str('{:03d}'.format(b)))
        f.write('                          ')
        f.write(str('{:03d}'.format(c)))
        f.write('\n')
    self.t = time 
    self.o = ones 
    self.tw = twos
    self.ev = events 
    self.oev = oxidevents
    self.rev = redevents
    f.close()
    #plt.plot(time,ones)
    #plt.show()
    
@monkeypatch_class(Zmc)
def refresh_mic_TOX(self):    
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    self.site_list = []
    self.rate_nums = []
    if self.seed_mode == 'CHA':
        from ase.io import read 
        #atoms = read('/afs/crc.nd.edu/user/a/agoswami/zmc/zmc/POSCAR_CHA')
        atoms = read('/afs/crc.nd.edu/user/a/agoswami/zmc/zmc/POSCAR_CHA')
        del atoms[[atom.index for atom in atoms if atom.symbol=='O']]           # Delete all O atoms 
        atomsrep = atoms.repeat([self.al,self.bl,self.cl])
        self.vol = atomsrep.get_volume()
    else:
        self.vol = self.xl*self.yl*self.zl
    for i in range(0,len(self.pairs)):
        if self.Cu_list[i].state == 1:                #If the copper is copper 1. 
            for p in self.pairs[i+1]:
                if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                    continue
                else:
                    d = self.pairs[i+1][p]
                    if self.calc_oxi_rate(d) != 0:
                        '''
                        if p < i+1:
                            continue
                        '''
                        #self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))
                        #self.site_list.append([i+1,p])
                        oxid_rates.append(self.calc_oxi_rate(d)/self.vol)
                        oxid_sites.append([i+1,p])
        #elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            #self.site_list.append([i+1,0])
            #redu_rates.append(self.calc_red_escape_time(1))
            #redu_sites.append([i+1,0])
    self.rate_nums =   redu_rates + oxid_rates 
    self.site_list =   redu_sites + oxid_sites
    #print(self.site_list)
    #print('Redu_size check',len(redu_rates),len(redu_sites))
    #print('Oxid_size check',len(oxid_rates),len(oxid_sites))
    #print(self.site_list)    
    
@monkeypatch_class(Zmc)
def calc_resid_fraction(self):
    if self.kin_mode!='transient_oxidation':
        print('This cannot work for transient oxidation case. Change the mode in kinetic_input.dat') 
    self.site_list = []
    self.rate_nums = []
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    self.pairs = {}
    self.vol = self.xl*self.yl*self.zl
    #self.period_image_Cu()
    for c in self.Cu_list:
        #print(c.label)
        nearest_neigh = {} 
        for u in self.Cu_list:             #A list that will keep track of nearest neighbours for a given Cu 
            d = self.mic_distance(c.position,u.position)
            #print(d)
            if u.label==c.label:
                continue
            nearest_neigh[u.label] = d
            #print(c.label, u.label, d)
        #Sorts all the pairs in ascending order of the site numbers
        #nearest_neigh = list(set(nearest_neigh))
        #nearest_neigh.sort()
        self.pairs[c.label] = nearest_neigh   
    #print(nearest_neigh)
    #print(self.pairs)
    for i in range(0,len(self.pairs)):
        if self.Cu_list[i].state == 1:                #If the copper is copper 1. 
            for p in self.pairs[i+1]:
                if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                    continue
                else:
                    if p < i+1:
                        continue
                    d = self.pairs[i+1][p]
                    if self.calc_oxi_rate(d) != 0:
                        '''
                        if p < i+1:
                            continue
                        '''
                        #self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))   #<-- This factor after / is to be added.
                        #self.site_list.append([i+1,p])
                        oxid_rates.append(self.calc_oxi_rate(d)/self.vol)
                        oxid_sites.append([i+1,p])
                        #print(d)
        elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            #self.site_list.append([i+1,0])
            redu_rates.append(self.calc_red_escape_time(1))
            redu_sites.append([i+1,0])
    self.rate_nums =   redu_rates + oxid_rates 
    self.site_list =   redu_sites + oxid_sites
    self.redlength = len(redu_rates)
    self.oxilength = len(oxid_rates)
    #print(self.oxilength)
    #print(self.site_list)
    #print(self.rate_nums)
    #print(self.oxilength)
    #print('Redu_size check',len(redu_rates),len(redu_sites))
    #print('Oxid_size check',len(oxid_rates),len(oxid_sites))
    #print(self.site_list)
    origin1 = self.site_list
    print('The number of pairs is:',len(self.site_list))
    #print(self.pairs)
    import random
    random.seed(int(self.rseed))
    print(self.NCu)
    flag = 0 ; resid = []
    while flag < 100:
        origin = random.sample(origin1,len(self.site_list)) ; 
        resid.append(residual(origin, self.NCu))
        flag += 1 
        print(flag)
    import statistics 
    res_frac = statistics.mean(resid)
    print(res_frac)
    f = open('Resid_frac.txt','w')
    for r in resid:
        f.write(str(r))
        f.write('\n')
    f.write('The residual fraction is {0}'.format(res_frac))

@monkeypatch_class(Zmc)
def minimum_image_convention_null(self):    #Since the Coppers wont be moving, the only distances that need to be measured is from MIC. 
    '''
    self.Cu_list_mic = []
    for c1 in self.Cu_list:
        [x,y,z] = c1.position 
        xnew = x - math.floor(x/(0.5*self.xl))*self.xl
        ynew = y - math.floor(y/(0.5*self.yl))*self.xl
        znew = z - math.floor(z/(0.5*self.xl))*self.xl
        self.Cu_list_mic.append(Copper(c1.label,c1.state,pos=[xnew,ynew,znew]))
        #print(xnew,ynew,znew)
    for c1 in self.Cu_list_mic:
        for c2 in self.Cu_list_mic:
            [x1,y1,z1] = c1.position
            [x2,y2,z2] = c2.position 
            if abs(x1-x2)<0.5*self.xl or abs(y1-y2)<0.5*self.yl or abs(z1-z2)<0.5*self.zl:
                print(c1.label,' is not following mic!')
    '''
    #First we compute all possible pair distances 
    self.site_list = []
    self.rate_nums = []
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    self.pairs = {}
    sum_rates = []
    self.vol = self.xl*self.yl*self.zl
    #self.period_image_Cu()
    for c in self.Cu_list:
        #print(c.label)
        nearest_neigh = {} 
        for u in self.Cu_list:             #A list that will keep track of nearest neighbours for a given Cu 
            d = self.mic_distance(c.position,u.position)
            #print(d)
            if u.label==c.label:
                continue
            nearest_neigh[u.label] = d
            #print(c.label, u.label, d)
        #Sorts all the pairs in ascending order of the site numbers
        #nearest_neigh = list(set(nearest_neigh))
        #nearest_neigh.sort()
        self.pairs[c.label] = nearest_neigh   
    #print(nearest_neigh)
    #print(self.pairs)
    for i in range(0,len(self.pairs)):
        if self.Cu_list[i].state == 1:                #If the copper is copper 1. 
            for p in self.pairs[i+1]:
                temp = 0
                if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                    continue
                else:
                    '''
                    if p < i+1:
                        continue
                    '''
                    d = self.pairs[i+1][p]
                    if self.calc_oxi_rate(d) != 0:
                        '''
                        if p < i+1:
                            continue
                        '''
                        #self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))   #<-- This factor after / is to be added.
                        #self.site_list.append([i+1,p])
                        oxid_rates.append(self.calc_oxi_escape_time(d,1)/self.vol)
                        oxid_sites.append([i+1,p])
                        temp += self.calc_oxi_escape_time(d,1)/self.vol
                        #print(d)
                sum_rates.append(temp)
        elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            #self.site_list.append([i+1,0])
            redu_rates.append(self.calc_red_escape_time(1))
            redu_sites.append([i+1,0])
    self.rate_nums =   redu_rates + oxid_rates 
    self.site_list =   redu_sites + oxid_sites
    self.redlength = len(redu_rates)
    self.oxilength = len(oxid_rates)
    #print(sum_rates)
    print(redu_rates)
    self.tau_max = max(sum_rates)
    print('Tau max:',self.tau_max)
    #print(self.oxilength)
    #print(self.site_list)
    #print(self.rate_nums)
    #print(self.oxilength)
    #print('Redu_size check',len(redu_rates),len(redu_sites))
    #print('Oxid_size check',len(oxid_rates),len(oxid_sites))
    #print(self.site_list)
    
@monkeypatch_class(Zmc)
def calc_steady_state_standard_null(self):
    import random
    t_start = timer()
    self.minimum_image_convention_null()
    random.seed(int(self.rseed))
    #print('Net coppers:',len(self.Cu_list))
    #print(self.max_events)  #This is the standard kMC algorithm.
    time = [0]
    ones = [self.NCuI]
    twos = [self.NCuII]
    events = [0]
    oxidevents = [0]
    redevents = [0]
    n_int = 1
    while True: #self.simevents < self.max_events or self.simutime < self.max_time:
        t_end = timer()
        if (t_end-t_start) > self.wtime:
            break       
        #### First we generate a random number and get the index of the executable event ###### 
        r1 = random.uniform(0,1)
        site= math.floor(r1*self.NCu) + 1 
        #print(site)
        if self.Cu_list[site-1].state == 1:
            event = 'Oxidation'
            #print('It is a one')
            oxid_sites = []
            oxid_rates = []
            for p in self.pairs[site]:
                if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                    continue
                else:
                    '''
                    if p < i+1:
                        continue
                    '''
                    d = self.pairs[site][p]
                    if self.calc_oxi_rate(d) != 0:
                        oxid_rates.append(self.calc_oxi_rate(d)/(self.vol))
                        oxid_sites.append([site,p]) 
            r2 = random.uniform(0,1)
            ind = find_ind_null(oxid_rates,r2,sum(oxid_rates)) 
            #ind = find_ind_null(oxid_rates,r2,self.tau_max) 
            if ind is None:
                print('NULL')
                continue
            else:
                #print(ind)
                r3 = random.uniform(0,1)
                rate = oxid_rates[ind]
                exec_time = -1*np.log(r3)/(rate)
                self.simutime += exec_time
                self.simevents += 1
                [a,b] = oxid_sites[ind]  #<----This refers to the sites where the event will occur
                self.Cu_list[a-1].state = 2 
                self.Cu_list[b-1].state = 2
                self.NCuI -= 2 
                self.NCuII += 2 
                self.n_O2 -= 1
                self.oxievents += 1 
        else:
            event = 'Reduction'
            #print('It is a two')
            rate = self.calc_red_escape_time(1)
            r3 = random.uniform(0,1)
            exec_time = -1*np.log(r3)/rate
            self.simutime += exec_time
            self.simevents += 1
            self.Cu_list[site-1].state = 1 
            self.NCuII -= 1 
            self.NCuI += 1 
            self.n_NO -= 1
            self.n_NH3 -= 1 
            self.redevents += 1 
        #print(self.tau_max)
        print('This',self.simutime, self.NCuI, timer()-t_start)
        #print(oxid_rates)
        #print(oxid_sites)
        if self.simevents == 100:
            break
        continue
        #print(self.site_list)
        #print(self.rate_nums)
        ind = find_ind(self.rate_nums,r1)
        #print('Found index at:', timer()-t_start)
        #print(ind)
        #print(r1,ind)
        if ind == None:
            #print(self.site_list)
            #print(r1)
            print('No more pairable coppers present')
            break
        #### Subsequently, we advance the simulation clock time by the rate of the event
        if self.rate_nums[ind] != 0:
            r2 = random.uniform(0,1)
            exec_time = -1*np.log(r2)/sum(self.rate_nums)
        else:
            print('No more oxidations!')
            break    #This is a contingency for the case when all ions are Copper 1's and are unable to be oxidized.  
        ##### Update the time and event counter ####################################
        self.simutime += exec_time
        self.simevents += 1
        sites = self.site_list[ind]  #<----This refers to the sites where the event will occur
        #print(sites)
        #print(exec_time, sites)
        '''
        for e in self.Cu_list:
            print(e.label, ",", e.state)
        '''
        '''
        ########## Applying a random selection algorithm here ##################################
        sites = random.choice(self.site_list)
        r = random.uniform(0,1)
        exec_time = -1*np.log(r)/sum(self.rate_nums)
        self.simutime += exec_time
        self.simevents += 1
        '''
        #### Extract the event name and take the actions accordingly ##########################
        if sites[1]==0:
            #This is a reduction event  
            red_site = sites[0]
            self.Cu_list[red_site-1].state = 1 
            self.NCuII -= 1 
            self.NCuI += 1 
            self.n_NO -= 1
            self.n_NH3 -= 1 
            self.redevents += 1 
            #self.refresh_rlist(red_site, 0, ind)
            #print(red_site)
            #self.Cu_list[c-1].state = 1 
        else:
            #This is an oxidation event 
            [a,b] = sites 
            self.Cu_list[a-1].state = 2 
            self.Cu_list[b-1].state = 2
            self.NCuI -= 2 
            self.NCuII += 2 
            self.n_O2 -= 1
            self.oxievents += 1 
            #self.refresh_rlist(a, b, ind)
            #print(a,b)
        #print('Refresh time:',self.NCuII)
        '''
        for e in self.Cu_list:
            print(e.label, ",", e.state)
        '''
        if self.spec_int.label == 'event':
            if self.simevents >= n_int*self.spec_int.interval:
                ones.append(self.NCuI)
                time.append(self.simutime)
                twos.append(self.NCuII)
                events.append(self.simevents)
                oxidevents.append(self.oxievents)
                redevents.append(self.redevents)
                n_int += 1
        else:
            if self.simutime >= n_int*self.spec_int.interval:
                ones.append(self.NCuI)
                time.append(n_int*self.spec_int.interval)                 #time.append(self.simutime)
                twos.append(self.NCuII)
                events.append(self.simevents)
                oxidevents.append(self.oxievents)
                redevents.append(self.redevents)
                n_int += 1        
        self.refresh_mic()
        #print(r2,self.simutime,self.NCuI)
        
        print('This',self.simutime, self.NCuI, timer()-t_start)
        if self.simevents == self.max_events:
            #print(self.rate_nums)
            print('Events done')
            break
        '''
        if self.simevents == 3:
            print(exec_time,",",self.site_list[ind])
            break
         '''
    f = open('raw_outputs.txt','w')
    f.write('Events           Time                      CuI                   CuII\n')
    for a,b,c,d in zip(time,ones,twos,events):
        f.write(' ' +str('{:03d}'.format(d)))
        f.write('           ')
        f.write(str('{:0.3e}'.format(a)))
        f.write('                   ')
        f.write(str('{:03d}'.format(b)))
        f.write('                   ')
        f.write(str('{:03d}'.format(c)))
        f.write('\n')
    f.close()
    f = open('event_outputs.txt','w')
    f.write('Events           Time                      Oxi_events                   Red_events\n')
    for a,b,c,d in zip(time,oxidevents,redevents,events):
        f.write(' ' +str('{:03d}'.format(d)))
        f.write('           ')
        f.write(str('{:0.3e}'.format(a)))
        f.write('                      ')
        f.write(str('{:03d}'.format(b)))
        f.write('                          ')
        f.write(str('{:03d}'.format(c)))
        f.write('\n')
    self.t = time 
    self.o = ones 
    self.tw = twos
    self.ev = events 
    self.oev = oxidevents
    self.rev = redevents
    f.close()
    #plt.plot(time,ones)
    #plt.show()
