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
                if self.kin_mode != 'steady' and self.kin_mode != 'transient_reduction' and self.kin_mode != 'transient_oxidation':
                    error_write('Unknown kinetic_mode {0} detected. Please re-edit the kinetic_input.dat file'.format(self.kin_mode))
                    return False
            elif 'A_reduction' in line:
                self.Ar = detect_num(line,'A_reduction')
                if self.Ar == False: return False 
            elif 'Ea_reduction' in line:
                self.Ear = detect_num(line,'A_reduction')
                if self.Ear == False: return False 
            elif 'Reduce_mode' in line:
                self.red_mode = line.split()[-1].replace('\n','')
                if self.red_mode != 'single' and self.red_mode != 'dual':
                    error_write('Unknown reduce_mode {0} detected. Please re-edit the kinetic_input.dat file'.format(self.red_mode))
                    return False
            elif 'A_pairing' in line:
                self.Ap = detect_num(line,'A_pairing')
                if self.Ap == False: return False
            elif 'Ea_pairing' in line:
                self.Eap = detect_num(line,'A_pairing')
                if self.Eap == False: return False 
            elif 'Decay' in line:
                self.dec_mode = line.split()[-1].replace('\n','')
                if self.dec_mode != 'exponential' and self.dec_mode != 'sigmoidal' and self.dec_mode != 'step':
                    error_write('Unknown reduce_mode {0} detected. Please re-edit the kinetic_input.dat file'.format(self.dec_mode))
                    return False
            elif len(line.strip())==0:
                continue
            else:
                error_write('Unknown line {0} detected. Please re-edit the kinetic_input.dat file'.format(line))
                return False 
        ############# Now beginning the functional dependence description #################################################
        if self.dec_mode == 'exponential':
            stat = self.set_expo_dependence()
            if stat == False: return False
        elif self.dec_mode == 'sigmoidal':
            stat = self.set_sig_dependence()
            if stat == False: return False
        elif self.dec_mode == 'step':
            stat = self.set_step_dependence()
            if stat == False: return False
        #print(self.dec_counters)
        #print(self.dec_counters!=[2,0,0])
        #print(self.dec_counters!=[0,1,0])
        #print(self.dec_counters!=[0,0,1])
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