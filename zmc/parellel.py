#This file is to take advantage of python's multiprocessing modules to see if I can speed up calculations 

import multiprocessing as mp
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
from multiprocessing import Process, Manager, Value, Queue

pair_list = {}
counter = 0


@monkeypatch_class(Zmc)
def create_pairs(self,c, uni_dict):
    #print(c.label)
    nearest_neigh = {} 
    for u in self.Cu_list:             #A list that will keep track of nearest neighbours for a given Cu 
        d = self.mic_distance(c.position,u.position)
        #print(d)
        if u.label==c.label:
            continue
        nearest_neigh[u.label] = d
        #uni_rate.append(self.calc_oxi_escape_time(d,1)/self.vol)
        #uni_list.append(temp)
        #print(c.label, u.label, d)
    #Sorts all the pairs in ascending order of the site numbers
    #nearest_neigh = list(set(nearest_neigh))
    #nearest_neigh.sort()
    uni_dict[c.label] = nearest_neigh 
    #print('New',pair_list)
    
@monkeypatch_class(Zmc)
def create_rate_lists(self,c):
    #print(c.label)
    import random
    nearest_neigh = {} 
    if c.state == 2:
        string = '{0}-0'.format(c.label)
        #uni_dict[random.uniform(0,1)/self.calc_red_escape_time(1)] = '{0}-0'.format(c.label)
        nearest_neigh[random.uniform(0,1)/self.calc_red_escape_time(1)] =  string
        #uni_dict += nearest_neigh
        #uni_dict[c.label] = nearest_neigh
        #nearest_neigh['sum'] = self.calc_red_escape_time(1)
        #s.value += self.calc_red_escape_time(1)  
        #uni_dict[tuple([c.label,0])] = nearest_neigh
        #summer.put(self.calc_red_escape_time(1))
        #uni_dict['sum'] += self.calc_red_escape_time(1)
        #uni_list.append([c.label,0])
        return nearest_neigh
    s = 0 
    for u in self.Cu_list:             #A list that will keep track of nearest neighbours for a given Cu 
        d = self.mic_distance(c.position,u.position)
        #print(d)
        if u.label==c.label:
            continue
        if u.state == 1 and c.state == 1 and self.calc_oxi_escape_time(d,1)!=0:
            string = '{0}-{1}'.format(c.label, u.label)
            #uni_list.append([c.label,u.label])
            #uni_dict[random.uniform(0,1)/(self.calc_oxi_escape_time(d,1)/self.vol)] = '{0}-{1}'.format(c.label,u.label) #d u.label
            #uni_dict[tuple([c.label,u.label])] = (self.calc_oxi_escape_time(d,1)/self.vol)
            nearest_neigh[random.uniform(0,1)/(self.calc_oxi_escape_time(d,1)/self.vol)] =  string
            #uni_dict['sum'] += self.calc_oxi_escape_time(d,1)/self.vol 
            #uni_dict['sum'] += self.calc_red_escape_time(1)
            #summer.put(self.calc_oxi_escape_time(d,1)/self.vol)
            #s.value += self.calc_oxi_escape_time(d,1)/self.vol
        #uni_rate.append(self.calc_oxi_escape_time(d,1)/self.vol)
        #uni_list.append(temp)
        #print(c.label, u.label, d)
    #Sorts all the pairs in ascending order of the site numbers
    #nearest_neigh = list(set(nearest_neigh))
    #nearest_neigh.sort()
    #nearest_neigh['sum'] = s
    #uni_dict += nearest_neigh
    return nearest_neigh
    uni_dict[c.label] = nearest_neigh     
    #print('New',pair_list)

@monkeypatch_class(Zmc)
def create_sums(self,P):
    self.rate_list['Net_sum'] += self.rate_list[P+1]['sum']
        
@monkeypatch_class(Zmc)
def create_rate_dict(self,P):
    self.rate_list.update(P)
        
@monkeypatch_class(Zmc)
def create_lists_nex(self,P,uni_rate,uni_site):   
    for p in self.pairs[P+1]:
        if self.Cu_list[p-1].state == 1:            #Ignore the pair if one of the coppers is a 2.            
            if p < P+1:
                continue      
            d = self.pairs[P+1][p]
            if self.calc_oxi_rate(d) != 0:
                uni_rate.append(self.calc_oxi_escape_time(d,1)/self.vol)
                uni_site.append([P+1,p])

                
#This is just a reference function to compare non-parellel and parellel capabilities
@monkeypatch_class(Zmc)
def mic_normal(self):
    self.site_list = []
    self.rate_nums = []
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    self.pairs = {}
    self.vol = self.xl*self.yl*self.zl
    t_start = timer()
    print('Starting now')
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
    t_end = timer()
    print('The normal time taken for the pair lists is:',t_end-t_start)
    
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
                        oxid_rates.append(self.calc_oxi_escape_time(d,1)/self.vol)
                        oxid_sites.append([i+1,p])
                        #print(d)
        elif self.Cu_list[i].state == 2:
            #self.rate_nums.append(self.calc_red_escape_time(1))
            #self.site_list.append([i+1,0])
            redu_rates.append(self.calc_red_escape_time(1))
            redu_sites.append([i+1,0])
    t_end2 = timer()
    print('The normal time taken for the pair lists is:',t_end2-t_end)
    
@monkeypatch_class(Zmc)
def mic_parellel(self, num = 0):
    import random
    import concurrent.futures
    self.site_list = []
    self.rate_nums = []
    oxid_sites = []
    oxid_rates = []
    redu_sites = []
    redu_rates = []
    manager = Manager()
    d = manager.dict()
    r = manager.dict()
    self.vol = self.xl*self.yl*self.zl
    t_start = timer()
    print('Starting now')
    #self.period_image_Cu()
    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(self.create_pairs, [(c,d) for c in self.Cu_list])
    #results = pool.starmap(self.create_pairs_nex, [(c,d,r,s) for c in self.Cu_list])
    #pool.close()
    self.pairs = d
    pool.close()
    pool.join()
    t_end = timer()
    print('The parellel time taken for the pair lists is:',t_end-t_start)
    '''
    #processes =[mp.Process(target=self.create_pairs, args=(c,d)) for c in self.Cu_list]
    #for p in processes: p.start()
    #for p in processes: p.join()
    t_end = timer()
    #self.pairs = results
    #s = Queue()
    #r['Net_sum'] = 0
    print('The parellel time taken for the pair lists is:',t_end-t_start) #It seems only the dictionary gets updated
    '''
    '''
    if num == 0:
        pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(12)
    #pool = mp.Pool(1)
    #results = [pool.apply(self.create_rate_lists, args=(c,r)) for c in self.Cu_list]
    results = pool.starmap_async(self.create_rate_lists, [(c,r) for c in self.Cu_list]).get()
    '''
    '''
    results = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for c in self.Cu_list:
            future = executor.submit(self.create_rate_lists, c)
            results.append(future.result())
    #self.rate_list = r 
    t_end2 = timer()
    #random.seed(int(self.rseed))
    print('The parellel time taken for the rate pair lists is:',t_end2-t_end)
    #pool.close()
    #pool.join()
    #print(results)
    self.rate_list = {}
    for res in results:
        self.rate_list.update(res)
    t_end3 = timer()
    print('The combination time for rate pair lists is:',t_end3-t_end2)
    self.min_time = min(self.rate_list.keys())
    self.site = self.rate_list[self.min_time]
    #print(min_time)
    print(self.rate_list[self.min_time])
    t_end4 = timer()
    print('The minimum finding time is:',t_end4-t_end3)
    #print(self.rate_list)
    '''

@monkeypatch_class(Zmc)
def calc_steadyscr_parellel(self):
    import random
    t_start = timer()
    self.mic_parellel()
    random.seed(int(self.rseed))
    print('Shall we begin?')
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
        self.simutime += self.min_time
        self.simevents += 1
        sites = [int(self.site.split('-')[0]),int(self.site.split('-')[1])]
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
        print('This',self.simutime, self.NCuI, timer()-t_start)
        self.mic_parellel(1)
        #print(r2,self.simutime,self.NCuI)
        if self.simevents == self.max_events:
            #print(self.rate_nums)
            print('Events done')
            break
            
@monkeypatch_class(Zmc)
def mic_parellel_refresh(self):
    pool = mp.Pool(12)
    print('Here')
    results = pool.map(self.create_rate_lists, [c for c in self.Cu_list])
    print('and here')
    pool.close()
    pool.join()
    pool.terminate()
    self.rate_list = {}
    for res in results:
        self.rate_list.update(res)
    self.min_time = min(self.rate_list.keys())
    self.site = self.rate_list[self.min_time]
