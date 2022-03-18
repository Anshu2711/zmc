import os 
import subprocess 
import numpy as np 
from random import seed 
from random import random
from timeit import default_timer as timer

########### Defining a Zmc class that will initialize by reading all the input files #############################
class interval:
    def __init__(self):
        self.label = 'Event' 
        self.interval = 10
        
class Copper:
    def __init__(self,lab,state,pos):
        self.label = lab
        self.state = state 
        self.position = pos

def error_write(nline):
    f = open('error_output.txt','w')
    f.write(nline)
    f.close()

def detect_num(line,char):
        try:
            return float(line.split()[-1].replace('\n',''))
        except:
            error_write('The {0} isnt a numeric value. Please re-enter'.format(char))
            return False 
        
class Zmc:
    '''
    Contains attributes corresponding to the various input parameters
    Reaction conditions: r_seed, T, P_O2, P_NH3, P_NO 
    Box description: lx, ly, lz, N, x_Cu1 
    Kinetic Parameters: A_pair, E_pair, A_red, E_red, f(d) 
    Simulation parameters: N_max, t_final, t_wall 
    Recording parameter: event_intval, time_intval 
    '''
 
    def __init__(self):
        '''
        Initialize by reading the input files
        '''
        self.rseed = None # 10000
        self.T = None #473
        self.P_O2 = None
        self.P_NH3 = None 
        self.P_NO = None 
        #self.P_tot = None #1 
        #self.x_O2 = None #0.1
        #self.x_NH3 = None #0.4
        #self.x_NO = None #0.3 
        self.max_events = None #100
        self.max_time = None #360 
        self.wtime = None #50
        self.xl = None #10.0
        self.yl = None #10.0
        self.zl = None #10.0
        self.NCu = None #20 
        self.xCuI = None #0.5
        self.Ar = None #5.0E+08 
        self.Ear = None #0.35 
        self.Ap = None #4.0E+09
        self.Eap = None #0.97
        self.seed_mode = None #'auto'
        self.autostat = 0 
        self.manstat = 0
        self.kin_mode = None
        self.red_mode = None #'single'
        self.dec_mode = None #'exponential'
        self.dec_counters = [0,0,0]  #Lists to validate the expo/sigmo/step dependence. 
        self.expo_count = 0 
        self.sig_count = 0 
        self.step_count = 0 
        self.steady_count = 0 
        self.tred_count = 0
        self.tox_count = 0
        self.simutime = 0
        self.simevents = 0
        self.oxievents = 0
        self.redevents = 0
        self.Cut_sigmo = None 
        self.Cut_step = None
        self.n_O2 = 0 
        self.n_NO = 0
        self.n_NH3 = 0 
        self.Cu_list = []
        self.snap_int = interval()
        self.proc_int = interval()
        self.spec_int = interval()
        if 'error_output.txt' in os.listdir():
            os.system('rm error_output.txt')
        if 'raw_outputs.txt' in os.listdir():
            os.system('rm raw_outputs.txt')
        if 'event_outputs.txt' in os.listdir():
            os.system('rm event_outputs.txt')
        stat1 = self.read_simu_input()
        stat2 = self.read_box_input()
        stat3 = self.read_kin_input()
        if stat1 == False or stat2 == False or stat3 == False:
            print('kMC calculation exitting') 
            exit()
        self.read_outputs()
                  
       
    def read_simu_input(self):
        if 'simulation_input.dat' not in os.listdir():
            error_write('simulation_input.dat not detected. Please create the required file in the directory')
            return False
        lines = open('simulation_input.dat','r').readlines()
        for line in lines:
            if line.startswith('#'):
                continue               #<-----Comment 
            elif line == '\n':
                continue               #<-----Empty line 
            elif '#' in line:
                line = line.split('#')[0]   #<------ Hash tag in the line which is ignored while parsing. 
                if len(line) == 0 or line == '':
                    continue
            if 'random_seed' in line:
                self.rseed = detect_num(line,'rseed')
                if self.rseed == False: return False
            elif 'temperature' in line:
                self.T = detect_num(line,'temperature')
                if self.T == False: return False
            elif 'O2_pressure' in line:
                self.P_O2 = detect_num(line,'O2_pressure')
                if self.P_O2 == False: return False 
            elif 'NO_pressure' in line:
                self.P_NO = detect_num(line,'NO_pressure')
                if self.P_NO == False: return False
            elif 'NH3_pressure' in line:
                self.P_NH3 = detect_num(line,'NH3_pressure')
                if self.P_NH3 == False: return False
            elif 'snapshots' in line: 
                if 'event' in line:
                    self.snap_int.label = 'event' 
                    self.snap_int.interval = detect_num(line,'snapshot')
                    if self.snap_int.interval == False: return False
                    self.snap_int.interval = int(self.snap_int.interval)
                elif 'time' in line:
                    self.snap_int.label = 'time'
                    self.snap_int.interval = detect_num(line,'snapshot')
                    if self.snap_int.interval == False: return False
                else:
                    error_write('Snapshot interval specified which is neither event or time. Please rectify in simulation_input.dat')
            elif 'process_statistics' in line: 
                if 'event' in line:
                    self.proc_int.label = 'event' 
                    self.proc_int.interval = detect_num(line,'process_statistics')
                    if self.proc_int.interval == False: return False
                    self.proc_int.interval = int(self.proc_int.interval)
                elif 'time' in line:
                    self.proc_int.label = 'time'
                    self.proc_int.interval = detect_num(line,'process_statistics')
                    if self.proc_int.interval == False: return False
                else:
                    error_write('Process statistics interval specified which is neither event or time. Please rectify in simulation_input.dat')
            elif 'species_numbers' in line: 
                if 'event' in line:
                    self.spec_int.label = 'event' 
                    self.spec_int.interval = detect_num(line,'species_numbers')
                    if self.spec_int.interval == False: return False
                    self.spec_int.interval = int(self.spec_int.interval)
                elif 'time' in line:
                    self.spec_int.label = 'time'
                    self.spec_int.interval = detect_num(line,'species_numbers')
                    if self.spec_int.interval == False: return False
                else:
                    error_write('Specie number interval specified which is neither event or time. Please rectify in simulation_input.dat')
            elif 'max_events' in line: 
                self.max_events = detect_num(line,'max_events')
                if self.max_events == False: return False
            elif 'max_time' in line: 
                self.max_time = detect_num(line,'max_time')
                if self.max_time == False: return False  
            elif 'wall_time' in line: 
                self.wtime = detect_num(line,'wall_time')
                if self.wtime == False: return False
            elif len(line.strip())==0:
                continue
            else:
                error_write('Unknown line {0} detected. Please re-enter in simulation_input.dat'.format(line))
                return False 
        msg = ['Random seed','Temperature','O2_pressure','NO_pressure','NH3_pressure','Snapshot interval','Process statistics interval','Specie number interval','Max events','Max time','Wall time']
        var =[self.rseed,self.T,self.P_O2,self.P_NO,self.P_NH3,self.snap_int.interval,self.proc_int.interval,self.spec_int.interval,self.max_events, self.max_time,self.wtime]
        ind = 0
        for v in var:
            if v == None:
                error_write('{0} has not been specified in simulation_input.dat. Please recheck and re-enter'.format(msg[ind]))
                return False 
            ind += 1 
    
    def read_box_input(self):
        if 'box_input.dat' not in os.listdir():
            error_write('box_input.dat not detected. Please create the required file in the directory')
            return False
        #lines = open('box_input.dat','r').readlines()
        with open('box_input.dat','r') as fi:
            lines = fi.readlines()
        for line in lines:
            if line.startswith('#'):
                continue               #<-----Comment 
            elif line == '\n':
                continue               #<-----Empty line 
            elif '#' in line:
                line = line.split('#')[0]   #<------ Hash tag in the line which is ignored while parsing. 
                if len(line) == 0 or line == '':
                    continue
            if 'x_length' in line:
                self.xl = detect_num(line,'x_length')
                if self.xl == False: return False
            elif 'y_length' in line:
                self.yl = detect_num(line,'y_length')
                if self.yl == False: return False
            elif 'z_length' in line:
                self.zl = detect_num(line,'z_length')
                if self.zl == False: return False
            elif 'N_Copper' in line:
                self.NCu = detect_num(line,'N_Copper')
                self.autostat += 1
                if self.NCu == False: return False
            elif 'Cu1_frac' in line:
                self.xCuI = detect_num(line,'Cu1_frac')
                self.autostat += 1
                if self.xCuI is False: return False
            elif 'Seed_mode' in line: 
                self.seed_mode = line.split()[-1].replace('\n','')
                if self.seed_mode != 'auto' and self.seed_mode != 'manual':
                    error_write('Unknown seed_mode {0} detected. Please re-edit the box_input.dat file'.format(self.seed_mode))
                    return False
            elif len(line.split())==4:               #<-- Checks if 4 numbers are entered as a first pass. Advanced checking in seed_box_manual
                for li in line.split(): 
                    try:
                        temp = float(li)
                    except:
                        error_write('Unknown line {0} detected. Please re-edit the box_input.dat file'.format(line))
                        return False
                self.manstat += 1
            elif len(line.strip())==0:
                continue
            else:
                error_write('Unknown line {0} detected. Please re-edit the box_input.dat file'.format(line))
                return False
        ############# Now beginning the seeding of the 3D box ########################################################
        if self.autostat == 2 and self.manstat !=0:         #<--Checks if entries for BOTH auto and manual have been entered. 
            if self.seed_mode == 'auto':
                error_write('Manual seeding entries detected. Only auto seeding entries N_Copper and Cu1_frac needed!')
                return False
            elif self.seed_mode == 'manual':
                error_write('Auto seeding entries detected. Only manual seeding entries of Cu oxidation state and positions needed')
                return False
        if self.seed_mode == 'auto':
            stat = self.seed_box_auto() 
            if stat == False: return False
        elif self.seed_mode == 'manual':
            stat = self.seed_box_manual()
            if stat == False: return False
        else:
            error_write('Seed mode not specified. Please re-enter in box_input.dat')
            return False
    
    def read_kin_input(self):
        if 'kinetic_input.dat' not in os.listdir():
            error_write('kinetic_input.dat not detected. Please create the required file in the directory')
            return False
        lines = open('kinetic_input.dat','r').readlines()
        for line in lines:
            if line.startswith('#'):
                continue               #<-----Comment 
            elif line == '\n':
                continue               #<-----Empty line 
            elif '#' in line:
                line = line.split('#')[0]   #<------ Hash tag in the line which is ignored while parsing. 
                if len(line) == 0 or line == '':
                    continue
            if '_expo' in line:
                self.dec_counters[0] += 1
                continue
            elif 'Sigmoid' in line:
                self.dec_counters[1] += 1 
                continue
            elif 'Step' in line:
                self.dec_counters[2] += 1
                continue
            elif 'Kinetic_mode' in line:
                self.kin_mode = line.split()[-1].replace('\n','')
                if self.kin_mode != 'steady_state' and self.kin_mode != 'transient_reduction' and self.kin_mode != 'transient_oxidation':
                    error_write('Unknown kinetic_mode {0} detected. Please re-edit the kinetic_input.dat file'.format(self.kin_mode))
                    return False
            elif 'A_reduction' in line:
                self.Ar = detect_num(line,'A_reduction')
                if self.Ar == False: return False 
                self.steady_count += 1
                self.tred_count += 1 
            elif 'Ea_reduction' in line:
                self.Ear = detect_num(line,'A_reduction')
                if self.Ear is False: return False 
                self.steady_count += 1
                self.tred_count += 1
            elif 'Reduce_mode' in line:
                self.red_mode = line.split()[-1].replace('\n','')
                if self.red_mode != 'single' and self.red_mode != 'dual':
                    error_write('Unknown reduce_mode {0} detected. Please re-edit the kinetic_input.dat file'.format(self.red_mode))
                    return False
                self.steady_count += 1
                self.tred_count += 1
            elif 'A_pairing' in line:
                self.Ap = detect_num(line,'A_pairing')
                if self.Ap == False: return False
                self.steady_count += 1
                self.tox_count += 1
            elif 'Ea_pairing' in line:
                self.Eap = detect_num(line,'A_pairing')
                if self.Eap is False: return False 
                self.steady_count += 1
                self.tox_count += 1
            elif 'Decay' in line:
                self.dec_mode = line.split()[-1].replace('\n','')
                if self.dec_mode != 'exponential' and self.dec_mode != 'sigmoidal' and self.dec_mode != 'step':
                    error_write('Unknown reduce_mode {0} detected. Please re-edit the kinetic_input.dat file'.format(self.dec_mode))
                    return False
                self.steady_count += 1
                self.tox_count += 1
            elif len(line.strip())==0:
                continue
            else:
                error_write('Unknown line {0} detected. Please re-edit the kinetic_input.dat file'.format(line))
                return False 
        if self.kin_mode == 'steady_state':
            stat = self.set_steady_state()
            if stat == False: return False
        elif self.kin_mode == 'transient_oxidation':
            stat = self.set_trans_oxidation()
            if stat == False: return False
        elif self.kin_mode == 'transient_reduction':
            stat = self.set_trans_reduction()
            if stat == False: return False
        if self.dec_counters!=[2,0,0] and self.dec_counters!=[0,1,0] and self.dec_counters!=[0,0,1]:    
            if self.dec_mode == 'exponential':
                error_write('Entries relevant to sigmoidal/step decay detected. Remove all those entries from kinetic_input.dat.')
                return False
            elif self.dec_mode == 'sigmoidal':
                error_write('Entries relevant to exponential/step decay detected. Remove all those entries from kinetic_input.dat.')
                return False
            elif self.dec_mode == 'step':
                error_write('Entries relevant to exponential/sigmoidal decay detected. Remove all those entries from kinetic_input.dat.')
                return False
    
    def seed_box_auto(self):
        if self.autostat != 2:
            error_write('Auto-seeding components N_Copper and Cu1_frac not entered. Please check and re-enter in box_input.dat')
            return False
        if self.xCuI == None or self.NCu == None or self.rseed == None:
            error_write('CuI fraction/Total Copper/Random seed not entered. Please check and re-enter in box_input.dat')
            return False
        if self.xCuI > 1.0 or self.xCuI < 0:
            error_write('The CuI fraction has to be between 0 and 1!. Please check and re-enter in box_input.dat')
            return False
        if self.xCuI == 0:
            print('Shoot!')
        self.NCuI = int(self.xCuI*self.NCu) 
        self.NCuII = int(self.NCu - self.NCuI)
        self.NCu = self.NCuI + self.NCuII
        if self.xl == None or self.yl == None or self.zl == None:
            error_write('x_length/y_length/z_length not entered. Please check and re-enter in box_input.dat')
            return False
        maxi_posi = np.array([self.xl,self.yl,self.zl])
        seed(int(self.rseed))
        for i in range(1,self.NCuI+1):
            rand_arr = np.array([random(), random(), random()])
            posi = np.multiply(rand_arr,maxi_posi)
            self.Cu_list.append(Copper(lab=i,state=1,pos=posi))
        for i in range(self.NCuI+1,int(self.NCu)+1):
            rand_arr = np.array([random(), random(), random()])
            posi = np.multiply(rand_arr,maxi_posi)
            self.Cu_list.append(Copper(lab=i,state=2,pos=posi))
    
    def seed_box_manual(self):
        if self.manstat==0:
            error_write('Please recheck the manual entries in box_input.dat. If not entered, please do so.')
            return False
        if self.xl == None or self.yl == None or self.zl == None:
            error_write('x_length/y_length/z_length not entered. Please check and re-enter in box_input.dat')
            return False
        lines = open('box_input.dat','r').readlines()
        self.NCuI = 0
        self.NCuII = 0
        Cu_list1 = []
        Cu_list2 = []
        flag1 = 0
        num1s = 0 
        for line in lines:
            if len(line.split())==4:
                ox_state = (float((line.split()[0])))
                if ox_state == 1:
                    num1s += 1 
        flag2 = num1s                 
        for line in lines:
            if len(line.split())==4:
                ox_state = (float((line.split()[0])))
                if ox_state - int(ox_state) < 1E-15:
                    ox_state = int(ox_state)
                else:
                    error_write('The oxidation state of Copper cannnot be a decimal here! Please recheck in box_input.dat')
                    return False 
                X_pos = float(line.split()[1])
                Y_pos = float(line.split()[2])
                Z_pos = float(line.split()[3])
                posi = np.array([X_pos, Y_pos, Z_pos])
                if X_pos > self.xl or Y_pos > self.yl or Z_pos > self.zl:
                    error_write('The position specified are outside the box. Please check the entries again in box_input.dat')
                    return False 
                if ox_state != 1 and ox_state !=2:
                    error_write('Invalid oxidation state {0} specified in box_input.dat. Please recheck'.format(ox_state))
                    return False 
                elif ox_state == 1:
                    flag1 += 1 
                    Cu_list1.append(Copper(lab=flag1,state=ox_state,pos=posi))
                    self.NCuI += 1 
                elif ox_state == 2:
                    flag2 += 1 
                    self.NCuII += 1
                    Cu_list2.append(Copper(lab=flag2,state=ox_state,pos=posi))
        self.Cu_list = Cu_list1 + Cu_list2
        self.NCu = self.NCuI + self.NCuII    #This automatically fills the Cu list first with O.S. 1 and then O.S. 2. 
        
    def set_steady_state(self):
        if self.steady_count!=6:
            error_write('All the inputs required for steady state calculations not entered. Please recheck and re-enter')
            return False
        if self.dec_mode == 'exponential':
            stat = self.set_expo_dependence()
            if stat == False: return False
        elif self.dec_mode == 'sigmoidal':
            stat = self.set_sig_dependence()
            if stat == False: return False
        elif self.dec_mode == 'step':
            stat = self.set_step_dependence()
            if stat == False: return False
        
    def set_trans_reduction(self):
        if self.Ap != None or self.Eap != None or self.dec_mode != None:
            error_write('For transient reduction, oxidation parameters A_pairing, Ea_pairing, decay mode dont need to be specified. Please recheck and re-enter in kinetic_input.dat')
            return False
        if self.NCuI != 0:
            error_write('For transient reduction, all Coppers have to be in +2 oxidation state initially. Please recheck and re-enter in box_input.dat')
            return False
        
    def set_trans_oxidation(self):
        if self.Ar != None or self.Ear != None or self.red_mode != None:
            error_write('For transient oxidation, reduction parameters A_reduction, Ea_reduction, reduce_mode dont need to be specified. Please recheck and re-enter in kinetic_input.dat')
            return False
        if self.NCuII != 0:
            error_write('For transient oxidation, all Coppers have to be in +1 oxidation state initially. Please recheck and re-enter in box_input.dat')
            return False
        if self.dec_mode == 'exponential':
            stat = self.set_expo_dependence()
            if stat == False: return False
        elif self.dec_mode == 'sigmoidal':
            print('Yes!')
            stat = self.set_sig_dependence()
            if stat == False: return False
        elif self.dec_mode == 'step':
            stat = self.set_step_dependence()
            if stat == False: return False
                                
    def set_expo_dependence(self):
        if self.dec_counters[0]!=2:
            error_write('Either A_expo or B_expo have not been specified. Please check and re-enter in kinetic_input.dat')
            return False
        lines = open('kinetic_input.dat','r').readlines()
        for line in lines:
            if 'A_expo' in line:
                self.A_expo = detect_num(line,'A_expo')
                if self.A_expo is False: return False 
                self.expo_count +=1
            elif 'B_expo' in line:
                self.B_expo = detect_num(line,'A_reduction')
                if self.B_expo is False: return False
                self.expo_count += 1 
        if self.expo_count != 2:
            error_write('Exponential dependence components A_expo, B_expo not recognized separately. Please check and re-enter in kinetic_input.dat')
            return False
        
    def set_sig_dependence(self):
        if self.dec_counters[1]!=1:
            error_write('Sigmoidal cutoff has not been specified. Please check and re-enter in kinetic_input.dat')
            return False
        lines = open('kinetic_input.dat','r').readlines()
        for line in lines:
            if 'Sigmoid_cutoff' in line:
                self.Cut_sigmo = detect_num(line,'Sigmoid_cutoff')
                if self.Cut_sigmo == False: return False 
        if self.Cut_sigmo >= self.xl or self.Cut_sigmo >= self.yl or self.Cut_sigmo >= self.zl:
            error_write('Sigmoid cutoff specified is greater than/equal to box dimensions. Please recheck the input at kinetic_input.dat.')
            return False

    def set_step_dependence(self):
        if self.dec_counters[2]!=1:
            error_write('Step cutoff has not been specified. Please check and re-enter in kinetic_input.dat')
            return False
        lines = open('kinetic_input.dat','r').readlines()
        for line in lines:
            if 'Step_cutoff' in line:
                self.Cut_step = detect_num(line,'Step_cutoff')
                if self.Cut_step == False: return False 
        if self.Cut_step > self.xl or self.Cut_step > self.yl or self.Cut_step > self.zl:
            error_write('Sigmoid cutoff specified is greater than box dimensions. Please recheck the input at kinetic_input.dat.')
            return False
    
    def init_snap(self):
        f = open('Box_state.txt','w')
        f.write('##### 3D Coordinates of the seeded coppers written below#######\n')
        f.write('S.No.     X      Y     Z\n')
        i = 0
        for cop in self.Cu_list:
            i += 1
            posit = cop.position
            f.write('{0}     {1}      {2}      {3}\n'.format(i,posit[0],posit[1],posit[2]))
        f.close()
        
    def init_spec(self):
        f = open('Specie_numbers.txt','w')
        f.write('####### Specie numbers of all the participating species########\n')
        f.write('S.No     Time         Events      CuI       CuII       O2       NO       NH3\n')
        f.close()
    
    def init_proc(self):
        f = open('Process_statistics.txt','w')
        f.write('###### Process counters of all the occuring processes###########\n')
        f.write('S.No     Time        Events        Pairing          Reduction\n')
        f.close()
        
        
        