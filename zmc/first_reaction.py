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

@monkeypatch_class(Zmc)
def minimum_image_convention_FR(self):    #Since the Coppers wont be moving, the only distances that need to be measured is from MIC. 
    #First we compute all possible pair distances 
    import random
    self.site_list = []
    self.rate_nums = []
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    self.pairs = {}
    self.vol = self.xl*self.yl*self.zl
    #self.period_image_Cu()
    #self.mic_parellel()
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
    random.seed(int(self.rseed))
    self.rate_dict = {}
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
                        #self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))   #<-- This factor after / is to be added.
                        #self.site_list.append([i+1,p])
                        #temp = -2*math.log(random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)
                        self.rate_dict.update({tuple([i+1,p]):-1*math.log(random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)})
                        self.rate_dict.update({tuple([p,i+1]):-1*math.log(random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)})
                        #self.rate_dict[self.calc_oxi_escape_time(d,1)/self.vol] = [i+1,p]
                        '''
                        oxid_rates.append(self.calc_oxi_escape_time(d,1)/self.vol)
                        oxid_sites.append([i+1,p])
                        '''
                        #print(d)
        elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            #self.site_list.append([i+1,0])
            temp = -1*math.log(1-random.uniform(0,1))/self.calc_red_escape_time(1)
            self.rate_dict.update({tuple([i+1,0]):temp})
            #self.rate_dict.update({tuple([0,i+1]):-1*math.log(random.uniform(0,1))/self.calc_red_escape_time(1)})
            #self.rate_dict[self.calc_red_escape_time(1)] = [i+1,0]
            '''
            redu_rates.append(self.calc_red_escape_time(1))
            redu_sites.append([i+1,0])
            '''
    #print('The rate_dict is:',self.rate_dict)
    #self.rate_dict.update(self.oxi_dict)
    #self.rate_dict.update(self.red_dict)
    self.rate_nums = list(self.rate_dict.values())
    self.state = random.getstate()
                          
@monkeypatch_class(Zmc)
def calc_steady_state_standard_periodic_FR(self):
    import random
    t_start = timer()
    self.minimum_image_convention_FR()
    #random.seed(int(self.rseed))
    #print('Net coppers:',len(self.Cu_list))
    #print(self.max_events)  #This is the standard kMC algorithm.
    time = [0]
    ones = [self.NCuI]
    twos = [self.NCuII]
    events = [0]
    oxidevents = [0]
    redevents = [0]
    n_int = 1
    simtime = [] ; os = []
    while True: #self.simevents < self.max_events or self.simutime < self.max_time:
        t_end = timer()
        if (t_end-t_start) > self.wtime:
            break       
        #### First we generate a random number and get the index of the executable event ###### 
        if len(self.rate_dict) == 0:
            print('No more pairable coppers here')
            break
        min_time =  min(self.rate_nums)
        ind = self.rate_nums.index(min_time)
        #print(min_time)
        #### Subsequently, we advance the simulation clock time by the rate of the event
        '''
        if self.rate_nums[ind] != 0:
            r2 = random.uniform(0,1)
            exec_time = -1*np.log(r2)/sum(self.rate_nums)
        else:
            print('No more oxidations!')
            break    #This is a contingency for the case when all ions are Copper 1's and are unable to be oxidized. 
        '''
        ##### Update the time and event counter ####################################
        self.simutime += min_time
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
        if sites[1]==0 or sites[0]==0:
            #This is a reduction event  
            red_site = sites[0] + sites[1]
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
        '''
        if self.simevents == 1:
            power = math.floor(math.log(self.simutime)) + 1 
            self.spec_int.interval = 10**power        #Automatically setting the time interval 
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
        #self.refresh_mic()
        #self.refresh_rlist_FR()
        self.refresh_rlist_scratch()
        #print(r2,self.simutime,self.NCuI)
        print('This',self.simutime, self.NCuI, timer()-t_start)
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
    plt.plot(simtime,os)
    plt.show()

@monkeypatch_class(Zmc)
def refresh_rlist_scratch(self):
    import random
    import numpy as np
    self.rate_dict = {}
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
                        #self.rate_nums.append(self.calc_oxi_escape_time(d,1)/len(self.Cu_list))   #<-- This factor after / is to be added.
                        #self.site_list.append([i+1,p])
                        #temp = -2*math.log(random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)
                        self.rate_dict.update({tuple([i+1,p]):-1*math.log(1-random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)})
                        self.rate_dict.update({tuple([p,i+1]):-1*math.log(1-random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)})
                        #self.rate_dict[self.calc_oxi_escape_time(d,1)/self.vol] = [i+1,p]
                        '''
                        oxid_rates.append(self.calc_oxi_escape_time(d,1)/self.vol)
                        oxid_sites.append([i+1,p])
                        '''
                        #print(d)
        elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            #self.site_list.append([i+1,0])
            temp = -1*math.log(1-random.uniform(0,1))/self.calc_red_escape_time(1)
            self.rate_dict.update({tuple([i+1,0]):temp})
            #self.rate_dict.update({tuple([0,i+1]):-1*math.log(random.uniform(0,1))/self.calc_red_escape_time(1)})
            #self.rate_dict[self.calc_red_escape_time(1)] = [i+1,0]
            '''
            redu_rates.append(self.calc_red_escape_time(1))
            redu_sites.append([i+1,0])
            '''
    #print('The rate_dict is:',self.rate_dict)
    #self.rate_dict.update(self.oxi_dict)
    #self.rate_dict.update(self.red_dict)
    self.rate_nums = list(self.rate_dict.values())

@monkeypatch_class(Zmc)
def refresh_rlist_FR(self):
    import random
    import numpy as np
    random.setstate(self.state)
    #print(self.rate_dict)
    #self.seed = self.rseed + 1
    #random.seed(int(self.seed))
    #print(self.sites)
    #if self.sites[1]!=0:
    if 0 not in self.sites:
        #Heads up: This is oxidation land 
        [a,b] = self.sites 
        for x in list(self.rate_dict):
            if a in x or b in x:
                del self.rate_dict[x]
        self.rate_dict.update({tuple([a,0]):-1*np.log(random.uniform(0,1))/self.calc_red_escape_time(1)})   #Updating the dictionary with reduction steps
        self.rate_dict.update({tuple([b,0]):-1*np.log(random.uniform(0,1))/self.calc_red_escape_time(1)})
        #self.rate_dict.update({tuple([0,a]):-1*math.log(random.uniform(0,1))/self.calc_red_escape_time(1)})
        #self.rate_dict.update({tuple([0,b]):-1*math.log(random.uniform(0,1))/self.calc_red_escape_time(1)})
    else:
        #Heads up: This is reduction land 
        a = self.sites[0] + self.sites[1]
        del self.rate_dict[(a,0)]
        #del self.rate_dict[(0,a)]
        for p in self.pairs[a]:
            if self.Cu_list[p-1].state == 2:            #Ignore the pair if one of the coppers is a 2. 
                continue
            else:
                if p < a:
                    continue
                d = self.pairs[a][p]
                if self.calc_oxi_escape_time(d,1) != 0:               #Updating the dictionary with oxidation steps 
                    #temp = math.log(random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)
                    #temp = -1*math.log(random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)
                    self.rate_dict.update({tuple([a,p]):-1*np.log(random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)})
                    self.rate_dict.update({tuple([p,a]):-1*np.log(random.uniform(0,1))*self.vol/self.calc_oxi_escape_time(d,1)})
    self.rate_nums = list(self.rate_dict.values())
    self.state = random.getstate()
    #print(self.rate_dict)                         