# -*- coding: utf-8 -*-
"""
Created on Fri Apr 05 19:43:42 2013

@author: Dr Dave
"""

#WARNING!!!!
#THIS IS CUSTOMIZED FOR ANTHONY'S ESR EXPERIMENT
#DO NOT USE FOR A QUBIT EXPERIMENT

from numpy import *
from slab.instruments import InstrumentManager
from slab.instruments import Alazar, AlazarConfig
from pulse_sequences import *
from slab import dsfit
from slab.dsfit import gaussfunc
from slab.datamanagement import load_array
import time
from scipy import optimize
#These classes do all the necessary data taking to determine 
#the static qubit parameters  (EC, EJ, fcavity, g)

#class which characterizes the qubit-cavity system
#Eventually integrate with a method for fitting to get parameters


def _qubit_freq_v_flux(p,x):
        
    #p[0]: max qubit frequency
    #p[1]: flux voltage offset
    #p[2]: flux to voltage tuning   
    #p[3]: "EC:
    
    return p[0]*((abs(cos((array(x)-p[1])*p[2])))+(4*p[3]/p[0])**2)**0.5/(1+(4*p[3]/p[0])**2)**0.5
    
def _cavity_freq_v_flux(p,x,qubit_fitvals):
        
    #p[0]: resonator frequency
    #p[1]: g
    #p[2]: phi to flux
    #f_qubit: qubit frequencies (same length as f_qubit)
    
    qubit_fitvals[2] = p[2] 
    f_qubit = _qubit_freq_v_flux(qubit_fitvals,x)
    
    return ((p[0]+f_qubit)/2.0-sqrt((f_qubit-p[0])**2/4.0+p[1]**2))*(f_qubit>=p[0]) + ((p[0]+f_qubit)/2.0+sqrt((f_qubit-p[0])**2/4.0+p[1]**2))*(f_qubit<p[0])



class cQED:
    
    def __init__(self):
        
        self.config = dict()    
        
        #qubit and cavity parameters
        self.config['EC'] = 1
        self.config['EJ'] = 1
        self.config['flux_volt_offset'] = 1.55
        self.config['flux_to_phi'] = 1
        self.config['max_qubit_freq'] = 7.15e9
        self.config['f0'] = 5759900000 #resonator bare frequency
        self.config['g']=67e6
        
        #from pulse probe...this is at the maximum flux point
        self.config['max_omega02_2'] = 0
        
        #from vacuum rabi (2 x n array with flux points)
        self.fread = []
        
        #from pulse probe (near the maximum flux point, 2 x n array with flux points)
        self.omega01 = []
        
        
        #from ramsey and t1 (4 x n array: 1: flux, 2: freq, 3: t1, 4: t2)
        self.qubit_props = []
   
       
    def load_qubit_char_data(self, vac_rabi_obj, pp_data_obj, omega02_2_max):
    
        #load fread data from vac_rabi_obj
        self.fread = zeros((2,len(vac_rabi_obj.flux_pts)))
        self.fread[0] = vac_rabi_obj.flux_pts
        self.fread[1] = vac_rabi_obj.vacuum_rabi_fitdata[0]
                
        #load omega data from pp_data_obj
        self.omega01 = zeros((2,len(pp_data_obj.flux_pts)))
        self.omega01[0] = pp_data_obj.flux_pts
        self.omega01[1] = pp_data_obj.pulse_probe_fitdata
        
        self.config['max_omega02_2'] = omega02_2_max
     
    #get qubit frequency from fit parameters
    def get_qubit_freq(self,flux_pts,n_levels=1):
        
        #construct matrix in charge space
        n_trunc = 10

        qubit_freqs = zeros((n_levels,len(flux_pts))) 
        
        EC = self.config['EC']/1.0e9
        EJ = self.config['EJ']/1.0e9
        
        for ii in range(len(flux_pts)):
            qubit_ham = zeros((2*n_trunc+1,2*n_trunc+1))
            
            for i,n in enumerate(range(-n_trunc,(n_trunc+1))):
                qubit_ham[i,i] = 4*EC*(n**2)
                   
                if i<(2*n_trunc):
                    qubit_ham[(i+1),i] = EJ*abs(cos((flux_pts[ii]-self.config['flux_volt_offset'])*self.config['flux_to_phi']))/2
                    qubit_ham[i,(i+1)] = EJ*abs(cos((flux_pts[ii]-self.config['flux_volt_offset'])*self.config['flux_to_phi']))/2
                    
            w = linalg.eigvalsh(qubit_ham)
            
            #sort w
            w = sort(w)
            
            for i in range(n_levels):            
                qubit_freqs[i,ii] = (w[i+1]-w[i])*1.0e9
        
        return qubit_freqs
        
        #return self.max_freq/sqrt(1+0.4**2)*sqrt(abs(cos(2*pi/4*(flux_pt-self.flux_offset)))+0.4**2) 
        
    #get the read frequency from fit parameters (1 photon shift, dispersive only and not including counter rotating terms)
    def get_read_freq(self,flux_pts):
        
        
        #get qubit frequencies
        q_freqs = self.get_qubit_freq(flux_pts,2)
        
        read_freqs = zeros((4,len(flux_pts))) 
        
        use_CR_terms = 0;
        
        for ii in range(len(flux_pts)):
            qubit_ham = zeros((9,9))
            
            qubit_ham[0,0] = 0 #|0,0>
            qubit_ham[0,4] = self.config['g']*use_CR_terms #CR term
            qubit_ham[4,0] = self.config['g']*use_CR_terms #CR term
            
            qubit_ham[1,1] = self.config['f0'] #|0,1>
            qubit_ham[1,6] = sqrt(2)*self.config['g']*use_CR_terms #CR term
            qubit_ham[6,1] = sqrt(2)*self.config['g']*use_CR_terms #CR term
            qubit_ham[1,2] = self.config['g']
            qubit_ham[2,1] = self.config['g']
            qubit_ham[2,2] = q_freqs[0,ii] #|1,0>
            qubit_ham[2,7] = 2*self.config['g']*use_CR_terms #CR term
            qubit_ham[7,2] = 2*self.config['g']*use_CR_terms #CR term
            
            qubit_ham[3,3] = 2*self.config['f0'] #|0,2>
            qubit_ham[3,4] = sqrt(2)*self.config['g']
            qubit_ham[4,3] = sqrt(2)*self.config['g']
            qubit_ham[4,4] = q_freqs[0,ii]+self.config['f0']  #|1,1>
            qubit_ham[4,8] = sqrt(6)*self.config['g']*use_CR_terms #CR term
            qubit_ham[8,4] = sqrt(6)*self.config['g']*use_CR_terms #CR term
            qubit_ham[4,5] = sqrt(2)*self.config['g']
            qubit_ham[5,4] = sqrt(2)*self.config['g']
            qubit_ham[5,5] = q_freqs[0,ii]+q_freqs[1,ii]  #|2,0>
            
            qubit_ham[6,6] = q_freqs[0,ii]+2*self.config['f0']  #|1,2>
            
            qubit_ham[7,7] = q_freqs[0,ii]+q_freqs[1,ii]+self.config['f0']  #|2,1>
            qubit_ham[8,8] = q_freqs[0,ii]+q_freqs[1,ii]+2*self.config['f0']  #|2,2>
            
            w,v = linalg.eigh(qubit_ham)
            
            #need to get eigenvalues corresponding to |g,0>,|g,1>,|e,0>,|e,1>
            
            v = transpose(v)            
            
            for jj in range(len(w)):
                
                if max(abs(v[jj]))==abs(v[jj,0]):
                    eg0 = w[jj]
                    
                if max(abs(v[jj]))==abs(v[jj,1]):
                    eg1 = w[jj]
                    
                if max(abs(v[jj]))==abs(v[jj,2]):
                    ee0 = w[jj]
                    
                if max(abs(v[jj]))==abs(v[jj,4]):
                    ee1 = w[jj]
            
            #sort w
            w = sort(w)
            #print (w[1]-w[0]) - 
            #print w

            #shift of cavity with qubit in ground state          
            read_freqs[0,ii] = (eg1-eg0)
            #shift of cavity with qubit in excited state
            read_freqs[1,ii] = (ee1-ee0)
            
            #shift of qubit with 0 photons in cavity
            read_freqs[2,ii] = (eg0-ee0)
            
            #shift of qubit with 1 photon in cavity
            read_freqs[3,ii] = (eg1-ee1)
            
                
        return read_freqs
    
    #get read frequency from vacuum rabi data
    def get_read_freq2(self,flux_pts):
        
        return interp(flux_pts,self.fread[0],self.fread[1])
    
        
    #compute from already loaded fread vs flux and omega01, omega02
    def compute_qubit_parameters(self, plotter, output_traces=False):
        

        if plotter is not None:
            plotter.init_plot("Frequency Fit", rank=1, accum=False)
            plotter.init_plot("Read Cavity Fit", rank=1, accum=False)
            plotter.init_plot("Qubit Freqs", rank=1, accum=False)

        #find max qubit frequency
        max_freq = self.omega01[1,0]
        max_flux = self.omega01[0,0]

        for i in range(len(self.omega01[0])):
            if self.omega01[1,i]>max_freq:
                max_freq = self.omega01[1,i]
                max_flux = self.omega01[0,i]
                
        print max_freq, max_flux
        
        fitvals = fitgeneral(self.omega01[0],self.omega01[1],_qubit_freq_v_flux,[max_freq,max_flux,1.0,1.0e8])
        
        #fitvals[1] = 0.42
        
        self.config['max_qubit_freq'] = fitvals[0] 
        self.config['flux_volt_offset'] = fitvals[1] 
        self.config['flux_to_phi'] = fitvals[2]*1.01
        
        fitvals[3] = (fitvals[0]-self.config['max_omega02_2'])*2
        
        if plotter is not None:
            print fitvals
            plotter.plot((concatenate((self.omega01[0],self.omega01[0])),concatenate((_qubit_freq_v_flux(fitvals,self.omega01[0]),self.omega01[1]))),"Frequency Fit")
        
        if output_traces:
            qubit_v_flux_fit = _qubit_freq_v_flux(fitvals,self.omega01[0])
            
        #fit read cavity data
        fitvals2 = fitgeneral(self.fread[0],self.fread[1],lambda p,x: _cavity_freq_v_flux(p,x,fitvals),[mean(self.fread[1])*1.001,200.0e6,self.config['flux_to_phi']])
                
        
        if output_traces:
            cavity_v_flux_fit = _cavity_freq_v_flux(fitvals2,self.fread[0],fitvals)        
        
        self.config['EC'] = (fitvals[0]-self.config['max_omega02_2'])*2    
        self.config['EJ'] = (self.config['max_qubit_freq']+self.config['EC'])**2/8/self.config['EC']
        self.config['f0'] = fitvals2[0]
        self.config['g'] = fitvals2[1]
        self.config['flux_to_phi'] = fitvals2[2]    
        
        #print self.config['EJ'],self.config['EC'],self.get_qubit_freq([-0.447])[0], self.get_read_freq([-0.45])[0]
        
        if plotter is not None:
            print fitvals2
            fitvals2[1] = 75.0e6
            plotter.plot((concatenate((self.fread[0], self.fread[0])),concatenate((_cavity_freq_v_flux(fitvals2,self.fread[0],fitvals),self.fread[1]))),"Read Cavity Fit")        
            x2 = linspace(-0.5,1.0,500) 
            #print shape(_qubit_freq_v_flux(fitvals,x2))
            #print shape(self.get_qubit_freq(x2))
            plotter.plot((concatenate((x2,x2)),concatenate((_qubit_freq_v_flux(fitvals,x2),self.get_qubit_freq(x2)[0]))),"Qubit Freqs") 
            
        if output_traces:
            return qubit_v_flux_fit,cavity_v_flux_fit
    
    
    def save_data(self,exp_path,file_prefix,overwrite=False):
        
        #print expt_path+fname
        try:
            if overwrite:
                f=SlabFile(exp_path+file_prefix+'.h5','w')
            else:
                fname=get_next_filename(exp_path,prefix=file_prefix,suffix='.h5')
                print fname
                f=SlabFile(exp_path+fname,'w')
        except NameError:
            print "Error opening file for saving cqed data"
            return
            
        #save
        f.create_dataset('fread',data=self.fread)
        f.create_dataset('omega01', data=self.omega01)
        f.create_dataset('qubit_props', data= self.qubit_props)
        f.attrs['qubit_data_type'] = "cqed"
        
        f.save_settings(self.config)
        
        f.close()        
        
    def load_data(self,exp_path,filename,settings_only=True):
        
        try:
            f=SlabFile(exp_path+filename,'r')
        except:
            print "Error opening file, will use default settings"
            return
            
        #check that this is a vacuum rabi file
        if f.attrs['qubit_data_type'] != 'cqed':
            print "Wrong data file! Not cqed data"
            f.close()
            return
        
        #load
        if not settings_only:
            
            self.fread = load_array(f,'fread') 
            self.omega01 = load_array(f,'omega01') 
            self.qubit_props = load_array(f,'qubit_props')
                
        
        self.config = f.load_settings()
        
        f.close()
        
#For fitting filters and two qubit "J"
class multipole:
    
    def __init__(self):
        
        self.config = dict()    
        
        #qubit and cavity parameters
        self.config['Qubit1'] = cQED()
        self.config['Qubit2'] = cQED()
        self.config['use_both_qubits'] = True 
        self.config['qubit1_crossflux'] = 0.0
        self.config['qubit2_crossflux'] = 0.0
        self.config['num_filters'] = 3
        self.config['filter_freq'] = 7e9
        self.config['g_filter_filter'] = 150e6
        self.config['g_filter_q1'] = 150e6
        self.config['g_filter_q2'] = 150e6
        self.config['J_12'] = 0 #direct J coupling
        
    def fit_multipole(self,flux1,filter_lines):
        
        #will vary filter_freq,g_filter_filter,g_filter_g1,and flux conversions for qubit 1
        #assume only 1 qubit
        
        (self.config['Qubit1']).config['flux_to_phi']
        (self.config['Qubit1']).config['flux_volt_offset']

        start_params = zeros(6)
        
        start_params[0] = self.config['filter_freq']
        start_params[1] = self.config['g_filter_filter']
        start_params[2] = self.config['g_filter_q1']
        start_params[3] = self.config['g_filter_q2']
        start_params[4] = (self.config['Qubit1']).config['flux_to_phi']
        start_params[5] = (self.config['Qubit1']).config['flux_volt_offset']
        
        bestfitparams, success = optimize.leastsq(self.__fit_err_func__, start_params, args=(flux1,filter_lines))
        
        self.config['filter_freq'] = bestfitparams[0]
        self.config['g_filter_filter'] = bestfitparams[1]
        self.config['g_filter_q1'] = bestfitparams[2]
        self.config['g_filter_q2'] = bestfitparams[3]
        (self.config['Qubit1']).config['flux_to_phi'] = bestfitparams[4]
        (self.config['Qubit1']).config['flux_volt_offset'] = bestfitparams[5]        
        
        return bestfitparams        
        
    def __fit_err_func__(self,fit_params,flux1,filter_lines):
        
        self.config['filter_freq'] = fit_params[0]
        self.config['g_filter_filter'] = fit_params[1]
        self.config['g_filter_q1'] = fit_params[2]
        self.config['g_filter_q2'] = fit_params[3]
        (self.config['Qubit1']).config['flux_to_phi'] = fit_params[4]
        (self.config['Qubit1']).config['flux_volt_offset'] = fit_params[5]

        c = self.multipole_spect(flux1,[0.0])
        c = c.T
        
        for i in range(4):
            c[i] = c[i] - filter_lines[i]
            
        return c.flatten()
    

    def multipole_spect(self,flux1,flux2):
        
        #flux1 and flux2 are 1D arrays of length n (or 1 if the flux was fixed)
        if len(flux1)==1:
            flux1 = ones(size(flux2))*flux1
            
        if len(flux2)==1:
            flux2 = ones(size(flux1))*flux2
            
        #print self.config['num_filters']
        num_qubits = 1
        if self.config['use_both_qubits']:
            num_qubits += 1
           
            
        num_levels = num_qubits+self.config['num_filters']
        spect_levels = zeros((len(flux1),num_levels))
        
        qubit1_energy = (self.config['Qubit1']).get_qubit_freq(flux1+array(flux2)*(self.config['qubit1_crossflux']))
        qubit1_energy = qubit1_energy[0]
            
        if self.config['use_both_qubits']:
            qubit2_energy = (self.config['Qubit2']).get_qubit_freq(flux2+array(flux1)*(self.config['qubit2_crossflux']))
            qubit2_energy = qubit2_energy[0]
        
        #print qubit1_energy,qubit2_energy
        
        #print shape(qubit1_energy),shape(qubit2_energy)
        
        for i in range(len(flux1)):
            
            h_matrix = zeros((num_levels,num_levels))
            
            h_matrix[0,0] = qubit1_energy[i]
            h_matrix[1,0] = self.config['g_filter_q1']
            h_matrix[0,1] = self.config['g_filter_q1']
            
            if self.config['use_both_qubits']: 
                h_matrix[num_levels-1,num_levels-1] = qubit2_energy[i]
                h_matrix[num_levels-2,num_levels-1] = self.config['g_filter_q2']
                h_matrix[num_levels-1,num_levels-2] = self.config['g_filter_q2']
            
                h_matrix[0,num_levels-1] = self.config['J_12']
                h_matrix[num_levels-1,0] = self.config['J_12']
            
            for j in range(self.config['num_filters']):
                h_matrix[j+1,j+1] = self.config['filter_freq']
                if j!=0:
                    h_matrix[j,j+1] = self.config['g_filter_filter']
                    h_matrix[j+1,j] = self.config['g_filter_filter']
                    
            w = linalg.eigvalsh(h_matrix)
            
            spect_levels[i] = w
            
                        
        return spect_levels
    
    def save_data(self,exp_path,file_prefix,overwrite=False):
        
        #print expt_path+fname
        try:
            if overwrite:
                f=SlabFile(exp_path+file_prefix+'.h5','w')
            else:
                fname=get_next_filename(exp_path,prefix=file_prefix,suffix='.h5')
                print fname
                f=SlabFile(exp_path+fname,'w')
        except NameError:
            print "Error opening file for saving cqed data"
            return
            
        #save
        f.attrs['qubit_data_type'] = "multipole"
        
        f.save_settings(self.config)
        
        f.close()        
        
    def load_data(self,exp_path,filename,settings_only=True):
        
        try:
            f=SlabFile(exp_path+filename,'r')
        except:
            print "Error opening file, will use default settings"
            return
            
        #check that this is a vacuum rabi file
        if f.attrs['qubit_data_type'] != 'multipole':
            print "Wrong data file! Not multipole data"
            f.close()
            return
        
        #load
        if not settings_only:
            
            self.fread = load_array(f,'fread') 
            self.omega01 = load_array(f,'omega01') 
            self.qubit_props = load_array(f,'qubit_props')
                
        
        self.config = f.load_settings()
        
        f.close()    
        
#generic class for taking data where the parameter 
#that is swept is a frequency (e.g. vacuum rabi, spectroscopy)

class frequency_sweep_exp():
    
    def __init__(self):
        
        self.config = dict(
            monitor_pulses=False,
            read_delay=50,
            card_delay=256, # 256*2+100
            read_length=2048,
            acq_length=1024,
            meas_range = [1.0, 1.0],
            total_length=8192,
            front_buffer=500,
            generator_enable_delay=200, #150
            
            #Qubit drive properties
            drive_freq=[5.5e9,7.1405e9],
            drive_sideband = [False,False],
            drive_sideband_freq = [50.0e6,50.0e6],
            drive_pwr=[13.0,13.0], #note that this doesn't do anything when driving with the BNC
            drive=['LB3','LB1'],
            shaped_pulse=True,
            pulsed_drive=True,
            
            #Qubit read properties
            read_freq=[5.6e9,5e9],
            lo_power=[16.0,16.0],
            read_pwr=[16.0,16.0],
            RF=['RF2','RF1'],
            LO=['RF2','RF1'],
            PHASE=['LBPHASE1','LBPHASE2'],
            read_phase = [70.0,70.0],
            IFreq=[10.0e6,10.0e6], #for heterodyne

            #flux properties
            flux_instr=[['YOKO2',1],['YOKO1',3]],
            flux=[0.0,0.0],
            
            qubit_enable = [True,True], #1 or 2 qubit experiment
    
            #awg properties
            #first awg is the main "trigger" AWG
            seq_file = "S:\\_Data\\sequence_files\temp.awg",
            AWG=['TEK','AGILENT811'],
            awg_amp = [[2.0,2.0,2.0,2.0],[1.0,1.0]],
            awg_offset = [[0.0,0.0,0.0,0.0],[0.0,0.0]],
            awg_marker_amps = [[1.0,1.0,1.0,1.0,2.0,1.0,1.0,1.0],[1.0,1.0,1.0,1.0]],
            awg_clocks = [1.2,4.2], #GSa/s (Note: these need to be set manually)
            awg_delay = [0,0.0,0.0,0.0,0,0], #These are the delays between the various channels (note: there is an arbitrary overall delay and the delays are defined w.r.t. the LOCAL clock)
            load_agilent = True, #if this is false then the agilent will not load under any condition!
            start_agilent = True, #if this is false then the agilent outputs will be set to off
            awg_pulse_offset = [0.0,0.0,0.0,0.0], #this is the offset when the drive pulse is on

            exp_id = '', #filled out by the specified experiment (rabi, ramsey, etc.)
            
            #acquisition properties
            num_avgs = -1, #-1 is run continuously (not an option here),
            data_window = [[0,1024],[0,1024]], #part of the acqusition to use for determining the qubit state,
            save_fulldata = False, #this saves the full trace
            
            cw_meas = False, #cw measurement does not pulse the generators
            use_awg = True,  #use the awg for pulsing (if not cw_meas) and triggering...if set to false then card is triggered using an alternate source
            
        )
        
    def __init_card_and_gens__(self):
        
        im = InstrumentManager()
        
        if self.config['num_avgs'] <= 0:
            raise NameError("Number of Averages must be a Positive Number")
        
        #if self.config['LO'][0]==self.config['RF'][0]:
        homodyne_meas=True
        print "Homodyne measurement"
        coup = 'DC'
        #else:
        #    homodyne_meas=False
        #    coup = 'AC'
            
        if self.config['num_avgs']==1:
            rec_per_buff = 1
        else:
            rec_per_buff = 512
        
       
        #load card
        card_config={'clock_edge': 'rising', 'clock_source': 'reference',  
            'trigger_coupling': 'DC',  'trigger_operation': 'or', 
            'trigger_source2': 'disabled','trigger_level2': 1.0, 'trigger_edge2': 'rising', 
            'trigger_level1': 0.6, 'trigger_edge1': 'rising', 'trigger_source1': 'external', 
            'ch1_enabled': True,  'ch1_coupling':coup,'ch1_range':self.config['meas_range'][0], 'ch1_filter': False, 
            'ch2_enabled': True, 'ch2_coupling': coup,'ch2_range':self.config['meas_range'][1], 'ch2_filter': False,            
            'bufferCount': 10,'recordsPerBuffer': rec_per_buff, 'trigger_delay': 0, 'timeout': 1000,
            'samplesPerRecord':self.config['acq_length'],'recordsPerAcquisition': self.config['num_avgs'], 'sample_rate': 1000000
            }
        
        scope_settings= AlazarConfig(card_config)
    
        card=Alazar(scope_settings)
        card.configure(scope_settings)
        
        RF = [None,None]
        LO = [None,None]
        drive = [None,None]        
        
        for j in range(2):
            
            #setup measurement and drive instruments for both qubits            
            if self.config['qubit_enable'][j]:
                
#                RF[j]=im[self.config['RF'][j]]
#                RF[j].set_output(True)
#                if not self.config['cw_meas']:
#                    RF[j].set_ext_pulse()
#                    RF[j].set_mod(True)
#                else:
#                    RF[j].set_mod(False)
#                RF[j].set_power(self.config['read_pwr'][j])
#                
#                if not homodyne_meas:
#                    LO[j]=im[self.config['LO'][j]]
#                    LO[j].set_output(True)
#                    LO[j].set_internal_pulse(10e-6)
#                    LO[j].set_mod(True)        
#                    LO[j].set_power(self.config['lo_power'][j])
#                
#                #set phase
#                if self.config['PHASE'][j] != '':
#                    phase_instr = im[self.config['PHASE'][j]]
#                    phase_instr.set_working_frequency(self.config['read_freq'][j])
#                    phase_instr.set_phase(self.config['read_phase'][j])                
                
                #setup drive
                drive[j] = im[self.config['drive'][j]]
                
                if self.config['drive'][j][0:2] == "RF": 
                    #labbrick
                    if not self.config['cw_meas']:
                        drive[j].set_pulse_ext(mod=True)
                        #Internal pulse must be set to false!
                    else:
                        drive[j].set_pulse_ext(mod=False)
                    drive[j].set_mod(False)
                else:
                    #Agilent
                    if not self.config['cw_meas']:
                        drive[j].set_ext_pulse()
                    else:
                        drive[j].set_ext_pulse()
                        
                drive[j].set_power(self.config['drive_pwr'][j])
                drive[j].set_frequency(self.config['drive_freq'][j])
                drive[j].set_output(True)
                
        return card, drive, RF, LO, homodyne_meas

#take a vacuum rabi trace        
class vacuum_rabi(frequency_sweep_exp):
    
    
    def __init__(self):
        
        frequency_sweep_exp.__init__(self)    
        
        self.config['exp_id'] = 'VAC_RABI'
        
        self.vacuum_rabi_data = []
        self.vacuum_rabi_fitdata = []
        self.flux_pts = []
        self.freq_pts = []
        
        #read sweep settings
        self.config['start_freq'] = [4.19e9,4.65e9]
        self.config['end_freq'] = [4.75e9,4.75e9]
        self.config['num_fpts'] = 200
        self.config['num_avgs'] = 20000*2
        
        #flux settings
        self.config['flux_start'] = [-0.5,-0.5]
        self.config['flux_end'] = [1.0,1.0]
        self.config['flux_pts'] = 50
        
        #do a pulse before the vacuum rabi?
        self.config['do_pulse'] = [False,False]
        self.config['pulse_length'] = [15.0,15.0]
        self.config['pulse_height'] = [0.1,0.1]
        
        #this is the read power that goes into the fridge
        #this doesn't set any property and is just for record keeping
        self.config['read_pwr2'] = [-50.0,-50.0]
        
        
    
     #take a 2D rabi data
    def take_read_data(self,plotter=None):
        
        #initialize the card and generator
        card, drive, RF, LO, homodyne_meas = self.__init_card_and_gens__()
            
        #define instrument manager    
        im=InstrumentManager()        
        
        #turn off the drive if we are not doing a pulse
        for j in range(2):
            if (drive[j] is not None) and (not self.config['do_pulse']):
                drive[j].set_output(False)
        
        if self.config['use_awg']:

            #stop the awg
            for i in range(len(self.config['AWG'])):
                awg=im[self.config['AWG'][i]]
                awg.stop()
                     
            #setup awg
                  
            #pulse sequence
            awg_pulse = pulse_sequence(1, self.config['total_length']*self.config['awg_clocks'][0],TEK_clk=self.config['awg_clocks'][0],Agilent_clk=self.config['awg_clocks'][1],delays=self.config['awg_delay'])
        
            awg_pulse.marker[awg_pulse.channel_index['card_trig']][0].square_pulse((self.config['front_buffer']+self.config['card_delay']),self.config['acq_length'])        
        
            if not self.config['cw_meas']:
                for j in range(2):
                    if self.config['qubit_enable']:
                        
                        awg_pulse.marker[awg_pulse.channel_index['read_trig'][j]][0].square_pulse((self.config['front_buffer']+self.config['read_delay']),self.config['read_length'])

                        #do a pulse before the vacuum rabi measurement (eg. to measure in the excited state)
                        if self.config['do_pulse'][j]:
                            pulse_loc = self.config['front_buffer'] - self.config['pulse_length'][j]
                            awg_pulse.marker[awg_pulse.channel_index['drive_trig'][j]][0].square_pulse2(pulse_loc,2*self.config['generator_enable_delay']+2*self.config['pulse_length'][j])
                            awg_pulse.analogwf[awg_pulse.channel_index['drive_I'][j]][0].gauss_pulse(pulse_loc,self.config['pulse_length'][j]/2,self.config['pulse_height'][j])

                        
                    
            awg_pulse.load_into_awg(self.config['seq_file'],self.config['AWG'][0])
            
            #set to sequence mode and set and the channel amplitudes and offsets
            awg=im[self.config['AWG'][0]]
            awg.prep_experiment()
            awg.set_amps_offsets(self.config['awg_amp'][0],self.config['awg_offset'][0],self.config['awg_marker_amps'][0])
            awg.run()
        
        #copied
        if plotter is not None:
            for j in range(2):
                if self.config['qubit_enable'][j]:
                    plotter.init_plot("Scope" + str(j),rank=1,accum=False)
                    plotter.init_plot("Readout"+ str(j),rank=1,accum=False)
                    plotter.init_plot("2DData"+ str(j),rank=2,accum=False)               
        
        self.flux_pts = [[],[]]
        self.freq_pts = [[],[]]
        
        for j in range(2):        
        
            self.flux_pts[j] = linspace(self.config['flux_start'][j],self.config['flux_end'][j],self.config['flux_pts'])
            self.freq_pts[j] = linspace(self.config['start_freq'][j],self.config['end_freq'][j],self.config['num_fpts']) 
        
        if self.config['save_fulldata']:
            self.vacuum_rabi_data = zeros((2,self.config['flux_pts'],self.config['num_fpts'],self.config['acq_length']))
        else:
            self.vacuum_rabi_data = zeros((2,self.config['flux_pts'],self.config['num_fpts']))
            
        #flux
        flux1 = im[self.config['flux_instr'][0][0]]
        flux2 = im[self.config['flux_instr'][1][0]]
        ch_pts = zeros((2,self.config['acq_length']))
        flux1.set_mode('voltage')
        flux2.set_mode('voltage')
      
        for i in range(self.config['flux_pts']):
                        
            flux1.set_volt(self.flux_pts[0][i],self.config['flux_instr'][0][1])
            flux2.set_volt(self.flux_pts[1][i],self.config['flux_instr'][1][1])
            
            phase_instr = im[self.config['PHASE'][0]]
            phase_instr.set_working_frequency(self.config['read_freq'][0])
            phase_instr.set_phase(20.0*i)
                    
            for j in range(self.config['num_fpts']):
        
                
                #set read frequency
                for k in range(2):
                    if self.config['qubit_enable'][k]:
                        RF[k].set_frequency(self.freq_pts[k][j])
                        if not homodyne_meas:
                            LO[k].set_frequency(self.freq_pts[k][j]+self.config['IFreq'][k])
                        
                        #let the generators settle to the new frequencies
                        while (not RF[k].get_settled()) or (not homodyne_meas and not LO[k].get_settled()):
                            pass
                            
                               
                tpts,ch_pts[0],ch_pts[1]=card.acquire_avg_data()
                
                for k in range(2):
                    if self.config['qubit_enable'][k]:
                        if self.config['save_fulldata']:
                            self.vacuum_rabi_data[k][i][j] = ch_pts[k]
                        else:
                            
                            if not homodyne_meas:
                                amp,_,_,_ = heterodyne(tpts,ch_pts[k],None,self.config['IFreq'][k],anti_alias=True)
                            else:
                                amp = mean(ch_pts[k][self.config['data_window'][k][0]:self.config['data_window'][k][1]])
                           
                            self.vacuum_rabi_data[k][i][j] = amp
                           
                            if plotter is not None:
                                plotter.plot(ch_pts[k],'Scope'+str(k))
                                plotter.plot((self.freq_pts[k][0:j],self.vacuum_rabi_data[k][i][0:j]),"Readout"+str(k))
                            
                            
                if plotter is not None:
                    for k in range(2):
                        if self.config['qubit_enable'][k] and not self.config['save_fulldata']:
                            plotter.plot((self.vacuum_rabi_data[k]),"2DData"+str(k))
       
            print im[self.config['PHASE'][0]].get_phase()
        
    def save_data(self,exp_path,file_prefix,overwrite=False):
        
        #print expt_path+fname
        try:
            if overwrite:
                f=SlabFile(exp_path+file_prefix+'.h5','w')
            else:
                fname=get_next_filename(exp_path,prefix=file_prefix,suffix='.h5')
                print fname
                f=SlabFile(exp_path+fname,'w')
        except NameError:
            print "Error opening file for saving vacuum rabi data"
            return
            
        #save
        f.create_dataset('vacuum_rabi_data',data=self.vacuum_rabi_data)
        f.create_dataset('vacuum_rabi_fitdata', data=self.vacuum_rabi_fitdata)
        f.create_dataset('flux_pts', data= self.flux_pts)
        f.create_dataset('freq_pts', data= self.freq_pts)
        f.attrs['qubit_data_type'] = "vacuum_rabi"
        
        #convert config dictionary into saveable format (get rid of nested lists)
        b = self.config.copy()
        for i in b.keys():
            if asarray(b[i]).dtype=="object":
                for j in range(len(b[i])):
                    b[i+str(j)] = b[i][j]
                    
                b[i] = []        
        
        f.save_settings(b)
        
        f.close()        
        
    def load_data(self,exp_path,filename,settings_only=True):
        
        try:
            f=SlabFile(exp_path+filename,'r')
        except:
            print "Error opening file, will use default settings"
            return
            
        #check that this is a vacuum rabi file
        if f.attrs['qubit_data_type'] != 'vacuum_rabi':
            print "Wrong data file! Not vacuum rabi data"
            f.close()
            return
        
        #load
        if not settings_only:
            self.vacuum_rabi_data = load_array(f,'vacuum_rabi_data')
            self.vacuum_rabi_fitdata = load_array(f,'vacuum_rabi_fitdata')
            self.flux_pts = load_array(f,'flux_pts')
            self.freq_pts = load_array(f,'freq_pts')
        
            
        self.config = f.load_settings()
        
        f.close()
        
        
    def fit_data(self,plotter=None,gauss_fit=True):
        
        #fit flux data to get f0 and Q
        if not plotter is None:
            plotter.init_plot("FitData",rank=1,accum=False)
            plotter.init_plot("PeakFreq",rank=1,accum=False)
        
        #Get the maximum
        freq_pts_fit = linspace(self.config['start_freq'],self.config['end_freq'],1000)
        
        self.vacuum_rabi_fitdata = zeros((2,len(self.flux_pts))) 
        self.vacuum_rabi_fitdata[1,:] = 1.0e6
        max_values = zeros(len(self.flux_pts))
        
        for i in range(len(self.flux_pts)):
            trans_data = interp(freq_pts_fit,self.freq_pts,self.vacuum_rabi_data[i])
            trans_max = trans_data[0]
            freq_max = freq_pts_fit[0]            
            for j in range(len(trans_data)):
                if trans_data[j]>trans_max:
                    trans_max = trans_data[j]
                    freq_max = freq_pts_fit[j]
            
            max_values[i] = trans_max
            self.vacuum_rabi_fitdata[0,i] = freq_max
        
        if not plotter is None:
            plotter.plot((self.flux_pts,self.vacuum_rabi_fitdata),"PeakFreq")
            
        #Do a more advanced gaussian fit
        if gauss_fit:
            for i in range(len(self.flux_pts)):
                fitvals = fitgeneral(self.freq_pts,self.vacuum_rabi_data[i],gaussfunc,[self.vacuum_rabi_data[i,0],max_values[i],self.vacuum_rabi_fitdata[0,i],self.vacuum_rabi_fitdata[1,i]])
                
                if not plotter is None:
                    plotter.plot((self.freq_pts,gaussfunc(fitvals,self.freq_pts)),"FitData")
                self.vacuum_rabi_fitdata[1,i] = fitvals[3]
                self.vacuum_rabi_fitdata[0,i] = fitvals[2]
            
            if not plotter is None:
                plotter.plot((self.flux_pts,self.vacuum_rabi_fitdata),"PeakFreq")
    
    def display_data(self,plotter,run_slice):
        
        #move towards matplotlib eventually        
        
        plotter.init_plot("ReadFreq",rank=1,accum=False)
        plotter.init_plot("ReadWidth",rank=1,accum=False)
        plotter.init_plot("Readout",rank=1,accum=False)
        plotter.init_plot("2DData",rank=2,xpts=self.freq_pts,ypts=self.flux_pts,accum=False)        
        
        if len(self.vacuum_rabi_data) < run_slice:
            print "Cannot display this slice, out of bounds"
        else:
            plotter.plot((self.freq_pts,self.vacuum_rabi_data[run_slice]),"Readout")
            
        plotter.plot(self.vacuum_rabi_data,"2DData")
        plotter.plot((self.flux_pts,self.vacuum_rabi_fitdata[0,:]),"ReadFreq")
        plotter.plot((self.flux_pts,self.vacuum_rabi_fitdata[1,:]),"ReadWidth")
    
    def cavity_peak(self,flux_pts):
        pass
        #use data to get the cavity peak (eg for a pulse probe spectroscopy)
        #plus pts can be an array
        cavity_peaks = interp(flux_pts,self.flux_pts,self.vacuum_rabi_fitdata[0,:])
        return cavity_peaks
    
    
#do pulse probe spectroscopy to get qubit frequency
class pulse_probe_spectroscopy(frequency_sweep_exp):
    
    
    def __init__(self):
        
        frequency_sweep_exp.__init__(self)
        
        self.config['exp_id'] = 'PP'
        
        self.pulse_probe_data = []
        #this a fit to omega01 versus flux
        self.pulse_probe_fitdata = []
       
        self.flux_pts = []
        
        #this has the same dimensions as pulse_probe_data beause the frequency array is not
        #constant
        self.freq_pts = []
        
        self.read_freq_pts = []
        
        #for fast flux pulse probe
        self.config['flux_pulse_length'] = [0.0,0.0]
        self.config['flux_pulse_delay'] = [0.0,0.0]
        self.config['do_fast_flux'] = False
        self.config['dc_flux'] = [0.0,0.0]
        self.config['add_flux_noise'] = False #adds flux noise from the Agilent AWG    
        
        self.config['drive_end_pwr'] = [0.0,0.0]
        self.config['num_drive_pulses'] = 1
        self.config['wait_time'] = 100.0
       
     #take a 2D rabi data
    def take_spec_data(self,plotter,read_freqs=[None,None],drive_freqs=[None,None]):
        
        
        use_fixed_range = [True,True]     
        
        #initialize the card and generator
        card, drive, RF, LO, homodyne_meas = self.__init_card_and_gens__()
       
        for j in range(2):
            
            if self.config['qubit_enable'][j]:            
            
                #if read_freqs is none then we just use the specified range
                if read_freqs[j] is None:
                    read_freqs[j] = linspace(self.config['read_freq_start'][j],self.config['read_freq_end'][j],self.config['flux_pts'])
                else:
                    if len(read_freqs[j]) != self.config['flux_pts']:
                        raise NameError("Read frequencies must have the same number of entries as flux points")
                
                if drive_freqs[j] is None:
                    use_fixed_range[j] =True
                else:
                    use_fixed_range[j] = False
                    if len(drive_freqs[j]) != self.config['flux_pts']:
                        raise NameError("Drive frequencies must have the same number of entries as flux points")
         
    
        
        self.flux_pts = zeros((2,self.config['flux_pts']))
        for j in range(2):
            self.flux_pts[j] = linspace(self.config['flux_start'][j],self.config['flux_end'][j],self.config['flux_pts'])
        self.read_freq_pts = read_freqs
        self.pulse_probe_data = zeros((2,self.config['flux_pts'],self.config['num_fpts']))
        self.freq_pts = zeros((2,self.config['flux_pts'],self.config['num_fpts']))
          
            
        #define instrument manager    
        im=InstrumentManager()

        for i in range(len(self.config['AWG'])):
            awg=im[self.config['AWG'][i]]
            awg.stop()
       
        if not self.config['cw_meas']:
            
            #Pulsed experiment!
            awg_clk_ratio = self.config['awg_clocks'][1]/self.config['awg_clocks'][0]
        
            #Load the waveforms for the card, measurement and drive pulses
            #pulse sequence
            awg_pulse = pulse_sequence(1,self.config['total_length'],ceil(self.config['total_length']*awg_clk_ratio))
            
            pulse_wait = -self.config['drive_pulse_length'][0]*0 - self.config['read_length']*0
                       
            max_wait = self.config['read_pulse_wait']+self.config['front_buffer']+self.config['generator_enable_delay']+self.config['wait_time']*(self.config['num_drive_pulses']-1.0)+self.config['num_drive_pulses']*max(self.config['drive_pulse_length'])+max(self.config['flux_pulse_length'])

            #trigger the card         
            awg_pulse.marker[awg_pulse.channel_index['card_trig']][0].square_pulse((max_wait+pulse_wait+self.config['card_delay']),self.config['acq_length'])
            
            #trigger the flux AWG
            if self.config['do_fast_flux']:
                awg_pulse.marker[awg_pulse.channel_index['flux_trig']][0].square_pulse(1,100)   
        
            for j in range(2):
                
                if self.config['qubit_enable'][j]:
            
                    if self.config['shaped_pulse']:
                        
                            for k in range(self.config['num_drive_pulses']):
                                #awg_pulse.analogwf[awg_pulse.channel_index['drive1_I']][0].gauss_pulse((max_wait-self.config['drive_pulse_length']),self.config['drive_pulse_length']/2.,self.config['drive_pulse_height'],False)
                                awg_pulse.analogwf[awg_pulse.channel_index['drive_I'][j]][0].smooth_square((max_wait-self.config['drive_pulse_length'][j]/2.-k*(self.config['drive_pulse_length'][j]+self.config['wait_time'])),20,self.config['drive_pulse_length'][j],self.config['drive_pulse_height'][j])
                                #awg_pulse.analogwf[awg_pulse.channel_index['drive_I'][j]][0].smooth_square_w_chirp((max_wait-self.config['drive_pulse_length'][j]/2.),20,self.config['drive_pulse_length'][j],self.config['drive_pulse_height'][j],-0.001,0.001)
                            awg_pulse.marker[awg_pulse.channel_index['drive_trig'][j]][0].square_pulse2((max_wait-self.config['drive_pulse_length'][j]/2.-(self.config['drive_pulse_length'][j]+self.config['wait_time'])*(self.config['num_drive_pulses']-1.0)/2),(2*self.config['generator_enable_delay']+1*self.config['drive_pulse_length'][j]+(self.config['drive_pulse_length'][j]+self.config['wait_time'])*(self.config['num_drive_pulses']-1.0)))
                            
                    else:
                        awg_pulse.marker[awg_pulse.channel_index['drive_trig'][j]][0].square_pulse2((max_wait-self.config['drive_pulse_length'][j]/2),self.config['drive_pulse_length'][j])
                    
                    awg_pulse.marker[awg_pulse.channel_index['read_trig'][j]][0].square_pulse((max_wait+pulse_wait+self.config['read_delay']),self.config['read_length'])
                
            
            
            awg_pulse.load_into_awg(self.config['seq_file'],self.config['AWG'][0])
            
            #if doing the fast flux, load those waveforms into the awg (note: don't sequence)
            if self.config['do_fast_flux']:
                
                for i in range(self.config['flux_pts']):
                    
                    for j in range(2):
                        if self.config['qubit_enable'][j]:
                            ind1 = awg_pulse.add_analog_pulse(awg_pulse.channel_index['flux'][j])
                            awg_pulse.analog_seq_table[awg_pulse.channel_index['flux'][j]][i] = ind1
                            
                            awg_pulse.analogwf[awg_pulse.channel_index['flux'][j]][ind1].square_pulse((max_wait-self.config['flux_pulse_length'][j]-self.config['flux_pulse_delay'][j])*awg_clk_ratio-self.config['awg_delay'][j],self.config['flux_pulse_length'][j]*awg_clk_ratio,self.flux_pts[j][i])
                            
                            #do a negative pulse at the end
                            awg_pulse.analogwf[awg_pulse.channel_index['flux'][j]][ind1].square_pulse((self.config['total_length']-2*self.config['flux_pulse_length'][j])*awg_clk_ratio,self.config['flux_pulse_length'][j]*awg_clk_ratio,-self.flux_pts[j][i])
                            
                      
                #load into agilent, but do not sequence!
                awg_pulse.load_full_into_agilent(self.config['AWG'][1],False)
                
            if self.config['add_flux_noise'] and self.config['do_fast_flux']:
                raise NameError("Can't do fast flux and add flux noise simultaneously")
                
            if self.config['add_flux_noise']:
                
                awg=im[self.config['AWG'][1]]
                for j in range(2):
                    awg.select_channel(j+1)
                    awg.write(':SOUR:FUNC:MODE FIX') #standard waveshape mode
                    awg.write(':SOUR:FUNC:SHAP NOIS') #noise
                    awg.write(':SOUR:FREQ:CW 10e3') #highest frequency
                    #awg.write(':SOUR:VOLT:LEV:AMPL 0.1') #set noise amplitude
                    awg.write(':ENAB 1') #CW mode
                    awg.write(':INIT:CONT 1') #CW mode
                    
                    
                
            #set to sequence mode and set and the channel amplitudes and offsets
            for i in range(len(self.config['AWG'])):
                awg=im[self.config['AWG'][i]]
                awg.set_amps_offsets(self.config['awg_amp'][i],self.config['awg_offset'][i],self.config['awg_marker_amps'][i])
                
                if i==0:  
                    #don't switch the agilent to sequence mode!
                    awg.prep_experiment()
                
                if not (i==1 and self.config['do_fast_flux']==False):
                    awg.run()
                    
                if i==1 and self.config['add_flux_noise']:
                    for j in range(2):
                        awg.select_channel(j+1)
                        awg.set_output(True)
                    
                
        
        else:
            
            #if we are doing a CW measurement without the awg we don't need to do anything! The card is being triggered by a BNC generator            
            
            if self.config['use_awg']:
                
                #pulse sequence
                awg_pulse = pulse_sequence(1,self.config['total_length'])
                awg_pulse.marker[awg_pulse.channel_index['card_trig']][0].square_pulse((1+self.config['card_delay']),self.config['acq_length'])
                awg_pulse.load_into_awg(self.config['seq_file'],self.config['AWG'][0])  
                
                #set to sequence mode and set and the channel amplitudes and offsets
                for i in range(len(self.config['AWG'])):
                    awg=im[self.config['AWG'][i]]
                    awg.set_amps_offsets(self.config['awg_amp'][i],self.config['awg_offset'][i],self.config['awg_marker_amps'][i])
                    
                    if i==0:  
                        #don't switch the agilent to sequence mode!
                        awg.prep_experiment()
                    
                    awg.run()                
                
                       
        #copied
        if not plotter is None:
            for j in range(2):
                if self.config['qubit_enable'][j]:
                    plotter.init_plot("Scope"+str(j),rank=1,accum=False)
                    plotter.init_plot("Readout"+str(j),rank=1,accum=False)
                    plotter.init_plot("2DData"+str(j),rank=2,accum=False)
            
               
        
               
        #flux
        flux1 = im[self.config['flux_instr'][0][0]]
        flux2 = im[self.config['flux_instr'][1][0]]
  
        if self.config['do_fast_flux']:
            #set the DC flux voltage
            flux1.set_volt(self.config['dc_flux'][0],self.config['flux_instr'][0][1])
            flux2.set_volt(self.config['dc_flux'][1],self.config['flux_instr'][1][1])
        
        drive_freq_pts = zeros((2,self.config['num_fpts']))
        
        if use_fixed_range:
            for j in range(2):
                drive_freq_pts[j] = linspace(self.config['freq_start'][j],self.config['freq_end'][j],self.config['num_fpts'])
        
        try:        
        
            for i in range(self.config['flux_pts']):
                            
    
                if not self.config['do_fast_flux']:
                    flux1.set_volt(self.flux_pts[0][i],self.config['flux_instr'][0][1])
                    flux2.set_volt(self.flux_pts[1][i],self.config['flux_instr'][1][1])
#                    phase_instr = im[self.config['PHASE'][0]]
#                    phase_instr.set_working_frequency(self.config['read_freq'][0])
#                    phase_instr.set_phase(i*180.0/self.config['flux_pts'])
                    time.sleep(1)
                
                #change read frequency
                for j in range(2):
                    if self.config['qubit_enable'][j]:
                        
                        RF[j].set_frequency(self.read_freq_pts[j][i])
                        if not homodyne_meas:
                            LO[j].set_frequency(self.read_freq_pts[j][i]+self.config['IFreq'][j])
                        
                        #drive frequency array
                        if not use_fixed_range:
                            drive_freq_pts[j] = [linspace(drive_freqs[j][i]-self.config['freq_range'][j]/2,drive_freqs[j][i]+self.config['freq_range'][j]/2,self.config['num_fpts'])]
                                  
                
                        #add to array
                        self.freq_pts[j,i] = drive_freq_pts[j]
                        
                        #change drive power
                        drive[j].set_power(self.config['drive_pwr'][j]+(self.config['drive_end_pwr'][j]-self.config['drive_pwr'][j])*i/self.config['flux_pts'])
                
                #if doing fast flux change agilent AWG sequence!
                if self.config['do_fast_flux']:
                    flux_awg = im[self.config['AWG'][1]]
                    flux_awg.set_to_trace(i+2,i+2)
                
                for j in range(self.config['num_fpts']):
                   
                    #change drive freq
                    for k in range(2):
                        if self.config['qubit_enable'][k]:
                            drive[k].set_frequency(drive_freq_pts[k][j])
                            
                            if (self.config['drive'][k] == "RF1" or self.config['drive'][k] == "RF2"):
                                 while not drive[k].get_settled():
                                    pass
                            else:
                                #time to let the generator settled
                                #was set to 200ms!
                                time.sleep(.05)
                            
                    tpts,ch1_pts,ch2_pts=card.acquire_avg_data()
                    
                    if not homodyne_meas:
                        amp1,phase1,amp2,phase2 = heterodyne(tpts,ch1_pts,ch2_pts,self.config['IFreq'])
                    else:
                        amp1 = mean(ch1_pts)
                        amp2 = mean(ch2_pts)
                    
                    ch_pts = [ch1_pts,ch2_pts]
                    self.pulse_probe_data[0,i,j] = amp1
                    self.pulse_probe_data[1,i,j] = amp2
                    
                    if not plotter is None:
                        for k in range(2):
                            if self.config['qubit_enable'][k]:
                                plotter.plot(ch_pts[k],'Scope'+str(k))
                                plotter.plot((drive_freq_pts[k,0:j],self.pulse_probe_data[k,i,0:j]),"Readout"+str(k))
                 
                if not plotter is None:
                    for k in range(2):
                        if self.config['qubit_enable'][k]:
                                plotter.plot((self.pulse_probe_data[k]),"2DData"+str(k))
                                
#                print im[self.config['PHASE'][0]].get_phase()
                
        except (KeyboardInterrupt, SystemExit):
            print "Exiting"
              
        
    def save_data(self,exp_path,file_prefix,overwrite=False):
        
        #print expt_path+fname
        try:
            if overwrite:
                f=SlabFile(exp_path+file_prefix+'.h5','w')
            else:
                fname=get_next_filename(exp_path,prefix=file_prefix,suffix='.h5')
                print fname
                f=SlabFile(exp_path+fname,'w')
        except NameError:
            print "Error opening file for saving pulse probe data"
            return
            
        #save
        f.create_dataset('pulse_probe_data',data=self.pulse_probe_data)
        f.create_dataset('pulse_probe_fitdata', data=self.pulse_probe_fitdata)
        f.create_dataset('flux_pts', data= self.flux_pts)
        f.create_dataset('freq_pts', data= self.freq_pts)
        f.create_dataset('read_freq_pts', data= self.read_freq_pts)
        
        
        f.attrs['qubit_data_type'] = "pulse_probe"
        
        #convert config dictionary into saveable format (get rid of nested lists)
        b = self.config
        for i in b.keys():
            if asarray(b[i]).dtype=="object":
                for j in range(len(b[i])):
                    b[i+str(j)] = b[i][j]
                    
                b[i] = []
        
        f.save_settings(self.config)
        
        f.close()        
        
    def load_data(self,exp_path,filename,settings_only=True):
        
        try:
            f=SlabFile(exp_path+filename,'r')
        except:
            print "Error opening file, will use default settings"
            return
            
        #check that this is a vacuum rabi file
        if f.attrs['qubit_data_type'] != 'pulse_probe':
            print "Wrong data file! Not pulse probe spectroscopy data"
            f.close()
            return
        
        #load
        if not settings_only:
            self.pulse_probe_data = load_array(f,'pulse_probe_data')
            self.pulse_probe_fitdata = load_array(f,'pulse_probe_fitdata')
            self.flux_pts = load_array(f,'flux_pts')
            self.freq_pts = load_array(f,'freq_pts')
            self.read_freq_pts = load_array(f,'read_freq_pts')
            
            
        self.config = f.load_settings()
        
        f.close()
        
        
    def fit_data(self, plotter, gauss_fit=True):
        
        #fit flux data to get omega01
        #NOTE: Looking for a dip in the transmission!
        if not plotter is None:
            plotter.init_plot("FitData",rank=1,accum=False)
            plotter.init_plot("PeakFreq",rank=1,accum=False)
        
        self.pulse_probe_fitdata = zeros(len(self.flux_pts)) 
        
        min_values = zeros(len(self.flux_pts))
        
        for i in range(len(self.flux_pts)):
            
            trans_min = self.pulse_probe_data[i,0]
            freq_min = self.freq_pts[i,0]            
            for j in range(len(self.pulse_probe_data[i])):
                if self.pulse_probe_data[i,j]<trans_min:
                    trans_min = self.pulse_probe_data[i,j]
                    freq_min = self.freq_pts[i,j]
            
            min_values[i] = trans_min
            self.pulse_probe_fitdata[i] = freq_min
        
        if not plotter is None:
            plotter.plot((self.flux_pts,self.pulse_probe_fitdata),"PeakFreq")
            
        #Do a more advanced gaussian fit
        if gauss_fit:
            for i in range(len(self.flux_pts)):
                fitvals = fitgeneral(self.freq_pts[i],self.pulse_probe_data[i],gaussfunc,[mean(self.pulse_probe_data[i]), min_values[i]-mean(self.pulse_probe_data[i]),self.pulse_probe_fitdata[i],1.0e6])
                
                if not plotter is None:
                    plotter.plot((self.freq_pts,gaussfunc(fitvals,self.freq_pts[i])),"FitData")
                self.pulse_probe_fitdata[i] = fitvals[2]
                
            if not plotter is None:
                plotter.plot((self.flux_pts,self.pulse_probe_fitdata),"PeakFreq")
    
    def display_data(self,plotter,run_slice):
        
        #put in code to pad        
        
        plotter.init_plot("Omega01",rank=1,accum=False)
        plotter.init_plot("Readout",rank=1,accum=False)
        plotter.init_plot("2DData",rank=2,accum=False)        
        
        if len(self.pulse_probe_data) < run_slice:
            print "Cannot display this slice, out of bounds"
        else:
            plotter.plot((self.freq_pts[run_slice-1],self.pulse_probe_data[run_slice-1]),"Readout")
            
        plotter.plot((self.pulse_probe_data),"2DData")
        if len(self.pulse_probe_fitdata) != 0:
            plotter.plot((self.flux_pts,self.pulse_probe_fitdata),"Omega01")
            
    def cavity_peak(self,flux_pt):
        
        #use data to get the cavity peak (eg for a pulse probe spectroscopy)
        #plus pts can be an array
        cavity_peaks = interp(flux_pts,self.flux_pts,self.vacuum_rabi_fitdata[0,:])
        return cavity_peaks
        
#do number splitting to get photon temperature
class number_splitting:
    
    
    #these are both 1D data
    splitting_data = []
    freq_pts = []
            
    config = dict()
    
    config['IFreq'] = 10e6
    
    #AWG settings
    config['read_delay']=50
    config['acq_time']=512
    config['card_delay']=256
    config['read_length']=1024
    config['total_length']=8192
    config['generator_enable_delay'] = 50
    
    
    #used freq start and freq end for a fixed drive range
    #these are the drive parameters
    config['freq_start'] = 6.5e9
    config['freq_end'] = 7.5e9
    config['freq_range'] = 0.2e9
    config['num_fpts'] = 200
    config['num_avgs'] = 20000
    config['drive_pwr'] = 10
    config['drive_pulse_height'] = 0.2
    config['drive_pulse_length'] = 500
    config['shaped_drive'] = True
    
    #specifications of the *loading* read pulse
    config['load_freq'] = 6.5e9
    config['load_pulse_length'] = 150
    config['load_pwr'] = -20
    config['load_wait'] = 50 #time between the end of the load pulse and the start of the drive pulse
    config['no_load'] = False    
    
    #specifications of the read pulse
    config['read_pulse_freq'] = 5.760e9
    config['lo_power'] = 16
    config['rf_power'] = -5 #this is the read pulse power!
    
    #flux is constat for this measurement
    config['flux_value'] = -0.5
    
    config['drive'] = 'LB1'
    config['RF'] = 'RF2' #read RF
    config['LO'] = 'RF1'
    config['flux'] = 'SRS'
    config['LOAD_RF'] = 'RF2' #RF to load photons (can be the same as read)
        
    #temperature of the fridge (must be inputted manually)
    config['temperature'] = 20e-3
    
     #take a 2D rabi data
    def take_splitting_data(self,plotter):
        
            
       
        #load card
        card_config={'clock_edge': 'rising', 'clock_source': 'reference',  
            'trigger_coupling': 'DC',  'trigger_operation': 'or', 
            'trigger_source2': 'disabled','trigger_level2': 1.0, 'trigger_edge2': 'rising', 
            'trigger_level1': 0.6, 'trigger_edge1': 'rising', 'trigger_source1': 'external', 
            'ch1_enabled': True,  'ch1_coupling': 'AC','ch1_range':0.1, 'ch1_filter': False, 
            'ch2_enabled': True, 'ch2_coupling': 'AC','ch2_range': 4, 'ch2_filter': False,            
            'bufferCount': 10,'recordsPerBuffer': 100, 'trigger_delay': 0, 'timeout': 1000,
            'samplesPerRecord':self.config['acq_time'],'recordsPerAcquisition': self.config['num_avgs'], 'sample_rate': 1000000
            }
        
        scope_settings= AlazarConfig(card_config)
    
        card=Alazar(scope_settings)
        card.configure(scope_settings)
            
        #define instrument manager    
        im=InstrumentManager()        
        
        #Setup transmission measurement equipement (RF and LO)
        LO=im[self.config['LO']]
        RF=im[self.config['RF']]
        
        if self.config['RF']==self.config['LOAD_RF']:
            print "Load and read sources are the same. Only a single frequency and power will be used."
            single_RF = True
        else:
            RF_load = im[self.config['LOAD_RF']]
            single_RF = False
           
        LO.set_output(True)
        RF.set_output(True)
        LO.set_mod(True)
        RF.set_mod(True)
        LO.set_power(self.config['lo_power'])
        RF.set_power(self.config['rf_power'])
        
        #set frequency
        RF.set_frequency(self.config['read_pulse_freq'])
        LO.set_frequency(self.config['read_pulse_freq']+self.config['IFreq'])
            
        
        if not single_RF:
            #assume loading RF is a lab brick
            RF_load.set_pulse_ext(mod=True)
            #Internal pulse must be set to false!
            RF_load.set_mod(False) 
            RF_load.set_power(self.config['load_pwr'])
            RF_load.set_frequency(self.config['load_freq'])
            RF_load.set_output(True)
        
        #setup drive
        drive = im[self.config['drive']]
        drive.set_pulse_ext(mod=True)
        #Internal pulse must be set to false!
        drive.set_mod(False)    
        drive.set_power(self.config['drive_pwr'])
        drive.set_output(True)
        
            
        #setup awg
        awg=im['AWG']
        
        awg.reset()
        
        #pulse sequence
        awg_pulse = pulse_sequence(self.config['total_length'])
        
        
        load_pulse_center = 1+self.config['generator_enable_delay']+self.config['load_pulse_length']/2
        
        drive_pulse_center = 1+self.config['generator_enable_delay']+self.config['load_pulse_length']+self.config['drive_pulse_length']+self.config['load_wait']
        
        if not self.config['no_load']:
            if single_RF:
                awg_pulse.read_trigger.square_pulse2(load_pulse_center,self.config['load_pulse_length'])
            else:
                awg_pulse.drive2_trigger.square_pulse2(load_pulse_center,self.config['load_pulse_length'])
            
        if self.config['shaped_drive']:
            awg_pulse.drive.gauss_pulse(drive_pulse_center,self.config['drive_pulse_length']/2.,self.config['drive_pulse_height'],False)
            awg_pulse.drive_trigger.square_pulse2(drive_pulse_center,2*self.config['generator_enable_delay']+2*self.config['drive_pulse_length'])
            drive_pulse_end = drive_pulse_center+self.config['drive_pulse_length']
        else:
            awg_pulse.drive_trigger.square_pulse2(drive_pulse_center,self.config['drive_pulse_length'])
            drive_pulse_end = drive_pulse_center+self.config['drive_pulse_length']/2
            
        awg_pulse.card_trigger.square_pulse(drive_pulse_end+self.config['card_delay'],self.config['read_length'])
        awg_pulse.read_trigger.square_pulse(drive_pulse_end+self.config['read_delay'],self.config['read_length'])
        awg_pulse.load_into_awg(awg,1,not single_RF)
        
        #awg_pulse.plot_pulses(plotter)
        
        #start
        awg.select_channel(1)
        awg.set_output(True)
        awg.set_mode("USER")    
        awg.set_amplitude(2.0)
        awg.set_offset(0.0)
        awg.define_sequence_step(1,1)
        
        awg.select_channel(2)
        awg.set_output(True)
        awg.set_mode("USER")    
        awg.set_amplitude(2.0)
        awg.set_offset(0.0)
        awg.set_trigger(src='EXT')
        
        #copied
        if not plotter is None:
            plotter.init_plot("Scope",rank=1,accum=False)
            plotter.init_plot("Readout",rank=1,accum=False)
            
        
        self.splitting_data = zeros(self.config['num_fpts'])
        self.freq_pts = linspace(self.config['freq_start'],self.config['freq_end'],self.config['num_fpts'])
        
        #flux
        flux = im[self.config['flux']]
        flux.set_volt(self.config['flux_value'])
        
       
        for i in range(len(self.freq_pts)):
                        
           
            #change drive freq
            drive.set_frequency(self.freq_pts[i])
            time.sleep(.2)
            tpts,ch1_pts,ch2_pts=card.acquire_avg_data()
            amp,_,_,_ = heterodyne(tpts,ch1_pts,ch2_pts,self.config['IFreq'])
            self.splitting_data[i] = amp
            if not plotter is None:
                plotter.plot(ch1_pts,'Scope')
                plotter.plot((self.freq_pts,self.splitting_data),"Readout")
                          
          
        
    def save_data(self,exp_path,file_prefix,overwrite=False):
        
        #print expt_path+fname
        try:
            if overwrite:
                f=SlabFile(exp_path+file_prefix+'.h5','w')
            else:
                fname=get_next_filename(exp_path,prefix=file_prefix,suffix='.h5')
                print fname
                f=SlabFile(exp_path+fname,'w')
        except NameError:
            print "Error opening file for saving splitting data"
            return
            
        #save
        f.create_dataset('splitting_data',data=self.splitting_data)
        f.create_dataset('freq_pts', data= self.freq_pts)
        
        f.attrs['qubit_data_type'] = "splitting"
        
        
        f.save_settings(self.config)
        
        f.close()        
        
    def load_data(self,exp_path,filename,settings_only=True):
        
        try:
            f=SlabFile(exp_path+filename,'r')
        except:
            print "Error opening file, will use default settings"
            return
            
        #check that this is a vacuum rabi file
        if f.attrs['qubit_data_type'] != 'splitting':
            print "Wrong data file! Not splitting data"
            f.close()
            return
        
        #load
        if not settings_only:
            self.splitting_data = load_array(f,'splitting_data')
            self.freq_pts = load_array(f,'freq_pts')
            
        self.config = f.load_settings()
        
        f.close()
        
        
    def fit_data(self, plotter, gauss_fit=True):
        
        #fit flux data to get omega01
        #NOTE: Looking for a dip in the transmission!
        if not plotter is None:
            plotter.init_plot("FitData",rank=1,accum=False)
            plotter.init_plot("PeakFreq",rank=1,accum=False)
        
        self.pulse_probe_fitdata = zeros(len(self.flux_pts)) 
        
        min_values = zeros(len(self.flux_pts))
        
        for i in range(len(self.flux_pts)):
            
            trans_min = self.pulse_probe_data[i,0]
            freq_min = self.freq_pts[i,0]            
            for j in range(len(self.pulse_probe_data[i])):
                if self.pulse_probe_data[i,j]<trans_min:
                    trans_min = self.pulse_probe_data[i,j]
                    freq_min = self.freq_pts[i,j]
            
            min_values[i] = trans_min
            self.pulse_probe_fitdata[i] = freq_min
        
        if not plotter is None:
            plotter.plot((self.flux_pts,self.pulse_probe_fitdata),"PeakFreq")
            
        #Do a more advanced gaussian fit
        if gauss_fit:
            for i in range(len(self.flux_pts)):
                fitvals = fitgeneral(self.freq_pts[i],self.pulse_probe_data[i],gaussfunc,[mean(self.pulse_probe_data[i]), min_values[i]-mean(self.pulse_probe_data[i]),self.pulse_probe_fitdata[i],1.0e6])
                
                if not plotter is None:
                    plotter.plot((self.freq_pts,gaussfunc(fitvals,self.freq_pts[i])),"FitData")
                self.pulse_probe_fitdata[i] = fitvals[2]
                
            if not plotter is None:
                plotter.plot((self.flux_pts,self.pulse_probe_fitdata),"PeakFreq")
    
    def display_data(self,plotter):
        
        #put in code to pad        
        
        plotter.init_plot("Readout",rank=1,accum=False)
        
        plotter.plot((self.freq_pts,self.splitting_data),"Readout")
         
#do flux sideband spectroscopy 
class sideband_spectroscopy:
    
    
    sideband_data = []    
   
    sideband_fitdata = []
    
    #this has the same dimensions as pulse_probe_data beause the frequency array is not
    #constant
    freq_pts = []
    
     
    config = dict()
    
    config['IFreq'] = 10e6
    
    #AWG settings
    config['read_delay']=50
    config['card_delay']=256
    config['read_length']=1024
    config['total_length']=4096
    config['generator_enable_delay'] = 50
    
    #these are the first drive pulse
    config['drive_pwr'] = 10
    config['drive_pulse_height'] = 0.2
    config['drive_pulse_length'] = 500
    config['drive_freq'] = 5.429e9
    config['shaped_drive'] = True
    
    config['lo_power'] = 16
    config['rf_power'] = -30
    config['read_freq'] = 5.76e9
    
    #two "modes" -> if this is false try and see photon coming out of resonator from sideband
    #if true then do a read pulse and look at state of qubit
    config['do_read'] = False
    
    #these are the sideband drive properties
    config['flux_drive_height'] = 0.2
    config['flux_drive_length'] = 1000
    config['freq_start'] = 250e6
    config['freq_end'] = 350e6
    config['num_fpts'] = 200
    config['num_avgs'] = 20000
    config['flux_drive_offset'] = 0
    
    config['flux_val'] = 0.18
    
    config['drive'] = 'LB1'
    config['RF'] = 'RF2'
    config['LO'] = 'RF1'
    config['flux'] = 'SRS'
        
    
     #take a 2D rabi data
    def take_sideband_data(self,plotter,):
        
            
        #load card
        card_config={'clock_edge': 'rising', 'clock_source': 'reference',  
            'trigger_coupling': 'DC',  'trigger_operation': 'or', 
            'trigger_source2': 'disabled','trigger_level2': 1.0, 'trigger_edge2': 'rising', 
            'trigger_level1': 0.6, 'trigger_edge1': 'rising', 'trigger_source1': 'external', 
            'ch1_enabled': True,  'ch1_coupling': 'AC','ch1_range':4, 'ch1_filter': False, 
            'ch2_enabled': True, 'ch2_coupling': 'AC','ch2_range': 4, 'ch2_filter': False,            
            'bufferCount': 10,'recordsPerBuffer': 100, 'trigger_delay': 0, 'timeout': 1000,
            'samplesPerRecord':self.config['read_length'],'recordsPerAcquisition': self.config['num_avgs'], 'sample_rate': 1000000
            }
        
        scope_settings= AlazarConfig(card_config)
    
        card=Alazar(scope_settings)
        card.configure(scope_settings)
            
        #define instrument manager    
        im=InstrumentManager()        
        
        #Setup transmission measurement equipement (RF and LO)
        LO=im[self.config['LO']]
        RF=im[self.config['RF']]
           
        LO.set_output(True)
        RF.set_output(True)
        LO.set_mod(True)
        RF.set_mod(True)
        LO.set_power(self.config['lo_power'])
        RF.set_power(self.config['rf_power'])
        
        #setup drive
        drive = im[self.config['drive']]
        drive.set_pulse_ext(mod=True)
        #Internal pulse must be set to false!
        drive.set_mod(False)    
        drive.set_power(self.config['drive_pwr'])
        drive.set_output(True)
        drive.set_frequency(self.config['drive_freq'])
        
            
        #setup awg
        awg=im['AWG']
        
        
        
        #copied
        if not plotter is None:
            plotter.init_plot("Scope",rank=1,accum=False)
            plotter.init_plot("Readout",rank=1,accum=False)
            
           
        RF.set_frequency(self.config['read_freq'])
        LO.set_frequency(self.config['read_freq']+self.config['IFreq'])
        
        self.sideband_data = zeros(self.config['num_fpts'])
        self.freq_pts = linspace(self.config['freq_start'], self.config['freq_end'],self.config['num_fpts'])
        
        
        #flux
        flux = im[self.config['flux']]
        flux.set_volt(self.config['flux_val'])
            
        
        for i in range(self.config['num_fpts']):
                        
            
            #change lo freq 
            if not self.config['do_read']:     
                LO.set_frequency(self.freq_pts[i]+self.config['drive_freq']+self.config['IFreq'])
            
            awg.reset()
        
            #pulse sequence
            awg_pulse = pulse_sequence(self.config['total_length'])
            
            pulse_wait = 0
            max_wait = 1+self.config['generator_enable_delay']+2*self.config['drive_pulse_length']+2*self.config['flux_drive_length']
            drive_pulse_center = max_wait-self.config['drive_pulse_length']-2*self.config['flux_drive_length']
            flux_drive_center = max_wait-self.config['flux_drive_length']-self.config['flux_drive_offset']
    
            if self.config['shaped_drive']:
                awg_pulse.drive.gauss_pulse(drive_pulse_center,self.config['drive_pulse_length']/2.,self.config['drive_pulse_height'],False)
                awg_pulse.drive_trigger.square_pulse2(drive_pulse_center,2*self.config['generator_enable_delay']+2*self.config['drive_pulse_length'])
                
            else:
                awg_pulse.drive_trigger.square_pulse2(max_wait-self.config['drive_pulse_length']/2,self.config['drive_pulse_length'])
                
            awg_pulse.card_trigger.square_pulse(flux_drive_center+self.config['flux_drive_length']+self.config['card_delay'],self.config['read_length'])
            if self.config['do_read']:        
                awg_pulse.read_trigger.square_pulse(flux_drive_center+self.config['flux_drive_length']+pulse_wait+self.config['read_delay'],self.config['read_length'])
        
            #sideband pulse
            awg_pulse.read.gauss_pulse_with_freq(flux_drive_center, self.config['flux_drive_length']/2, self.config['flux_drive_height'], self.freq_pts[i]/1.0e9, True)       
        
            print self.freq_pts[i]/1.0e9        
        
            #load into AWG
            awg_pulse.load_into_awg(awg,1,True)
                
        
            #awg_pulse.plot_pulses(plotter)
            
            #start
            awg.select_channel(1)
            awg.set_output(True)
            awg.set_mode("USER")    
            awg.set_amplitude(2.0)
            awg.set_offset(0.0)
            awg.define_sequence_step(1,1)
            
            awg.select_channel(2)
            awg.set_output(True)
            awg.set_mode("USER")    
            awg.set_amplitude(2.0)
            awg.set_offset(0.0)
            awg.set_trigger(src='EXT')            
                
            
           
            #change drive freq
            tpts,ch1_pts,ch2_pts=card.acquire_avg_data()
            amp,_,_,_ = heterodyne(tpts,ch1_pts,ch2_pts,self.config['IFreq'])
            self.sideband_data[i] = amp
            if not plotter is None:
                plotter.plot(ch1_pts,'Scope')
                plotter.plot((self.freq_pts,self.sideband_data),"Readout")
             
        
    def save_data(self,exp_path,file_prefix,overwrite=False):
        
        #print expt_path+fname
        try:
            if overwrite:
                f=SlabFile(exp_path+file_prefix+'.h5','w')
            else:
                fname=get_next_filename(exp_path,prefix=file_prefix,suffix='.h5')
                print fname
                f=SlabFile(exp_path+fname,'w')
        except NameError:
            print "Error opening file for saving sideband data"
            return
            
        #save
        f.create_dataset('sideband_data',data=self.sideband_data)
        f.create_dataset('sideband_fitdata', data=self.sideband_fitdata)
        f.create_dataset('freq_pts', data= self.freq_pts)
        
        f.attrs['qubit_data_type'] = "sideband"
        
        
        f.save_settings(self.config)
        
        f.close()        
        
    def load_data(self,exp_path,filename,settings_only=True):
        
        try:
            f=SlabFile(exp_path+filename,'r')
        except:
            print "Error opening file, will use default settings"
            return
            
        #check that this is a vacuum rabi file
        if f.attrs['qubit_data_type'] != 'sideband':
            print "Wrong data file! Not sideband spectroscopy data"
            f.close()
            return
        
        #load
        if not settings_only:
            self.sideband_data = load_array(f,'sideband_data')
            self.sideband_fitdata = load_array(f,'sideband_fitdata')
            self.freq_pts = load_array(f,'freq_pts')
            
        self.config = f.load_settings()
        
        f.close()
        
        
    def fit_data(self, plotter, gauss_fit=True):
        
        #fit flux data to get omega01
        #NOTE: Looking for a dip in the transmission!
        pass
    
    def display_data(self,plotter,run_slice):
        
        #put in code to pad        
        
        plotter.init_plot("Omega01",rank=1,accum=False)
        plotter.init_plot("Readout",rank=1,accum=False)
        plotter.init_plot("2DData",rank=2,accum=False)        
        
        if len(self.pulse_probe_data) < run_slice:
            print "Cannot display this slice, out of bounds"
        else:
            plotter.plot((self.freq_pts[run_slice-1],self.pulse_probe_data[run_slice-1]),"Readout")
            
        plotter.plot((self.pulse_probe_data),"2DData")
        if len(self.pulse_probe_fitdata) != 0:
            plotter.plot((self.flux_pts,self.pulse_probe_fitdata),"Omega01")
            
    def cavity_peak(self,flux_pt):
        pass
        #use data to get the cavity peak (eg for a pulse probe spectroscopy)
        

#do esr experiment
class esr_exp(frequency_sweep_exp):
    
    
    def __init__(self):
        
        frequency_sweep_exp.__init__(self)
        
        self.config['exp_id'] = 'ESR'
        
        self.esr_data = []
        
       
        self.tau_pts = []
        
                
        self.config['use_awg'] = True        
        
        #custom esr options
        self.config['b_field'] = 0.0
        self.config['tau'] = 10.0 #in us
        self.config['led_delay'] = 2000 #in us
        self.config['led_pulse_length'] = 2.0 #in ms
        self.config['pulse_height'] = [1.0,1.0]
        self.config['pulse_length'] = [100,100] #in ns
        self.config['pulse_phase'] = [0,pi/2]
        self.config['master_trigger'] = self.config['led_pulse_length'] + self.config['led_delay']*1e-3 + 3*self.config['tau'] #ms
        self.config['switch_buffer'] = 2000 #ns
        
        self.config['tau_start'] = 1.0
        self.config['tau_end'] = 2.0
        self.config['tau_pts'] = 100
        
        self.config['trigger_instr'] = 'BNC_1'
        self.config['awg_trigger'] = 'BNC_2'
        self.config['led_trigger'] = 'BNC_3'
        
        #there is only 1 "qubit"
        self.config['qubit_enable'][1] = False
        
        #save single shot
        self.config['save_single_shot'] = False
       
       
    def take_esr_data(self,plotter):
        
        #initialize the card and generator
        card, drive, RF, LO, homodyne_meas = self.__init_card_and_gens__()
            
        #define instrument manager    
        im=InstrumentManager()  
        
        #setup triggers
        
        #master trigger
        mast_trigger = im[self.config['trigger_instr']]
        mast_trigger.set_function('PULSE')
        mast_trigger.set_pulse_width(100e-9)
        mast_trigger.set_period(self.config['master_trigger']*1.0e-3)
        mast_trigger.set_output(False)
        
        #awg trigger
        awg_trigger = im[self.config['awg_trigger']]        
        awg_trigger.set_function('PULSE')
        awg_trigger.set_pulse_width(100e-9)
        awg_trigger.set_period(self.config['tau']*1.0e-6)
        awg_trigger.set_burst_cycles(3)
        awg_trigger.set_burst_mode('triggered')
        awg_trigger.set_burst_state('on')
        awg_trigger.set_trigger_source('EXT')
        
        #led trigger
        led_trigger = im[self.config['led_trigger']]
        led_trigger.set_function('PULSE')
        led_trigger.set_pulse_width(self.config['led_pulse_length']*1.0e-3)
        led_trigger.set_burst_cycles(1)
        led_trigger.set_burst_mode('triggered')
        led_trigger.set_amplitude(5.0)
        led_trigger.set_offset(0.0)
        led_trigger.set_burst_state('on')
        
       
        #Load the waveforms for the card, measurement and drive pulses
        #pulse sequence
        awg_seq = pulse_sequence(3,int(ceil(self.config['total_length']*self.config['awg_clocks'][0])),int(ceil(self.config['total_length']*self.config['awg_clocks'][1])),self.config['awg_clocks'][0],self.config['awg_clocks'][1],self.config['awg_delay'])
            
        pulse_center = self.config['total_length']/2.0
        
        #setup channels
        awg_seq.channel_index['drive_I'][0] = 4
        awg_seq.channel_index['drive_Q'][0] = 5
        
        awg_seq.channel_index['card_trig'] = 8
        awg_seq.channel_index['led_trigger'] = 9
        awg_seq.channel_index['hittite_switch'] = 10
        
        #three pulses
        for j in range(3):
            
            #add new pulses
            ind1a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][0])
            awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][0]][j] = ind1a
        
            ind1b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][0])
            awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][0]][j] = ind1b
        
            ind2 = awg_seq.add_marker_pulse(awg_seq.channel_index['card_trig'])
            awg_seq.marker_seq_table[awg_seq.channel_index['card_trig']][j] = ind2
            ind3 = awg_seq.add_marker_pulse(awg_seq.channel_index['led_trigger'])
            awg_seq.marker_seq_table[awg_seq.channel_index['led_trigger']][j] = ind3
            ind4 = awg_seq.add_marker_pulse(awg_seq.channel_index['hittite_switch'])
            awg_seq.marker_seq_table[awg_seq.channel_index['hittite_switch']][j] = ind4
                     
            if j<2:
                
                #microwave pulses
                awg_seq.analogwf[awg_seq.channel_index['drive_I'][0]][ind1a].square_pulse2(pulse_center,self.config['pulse_length'][j],self.config['pulse_height'][j]*cos(self.config['pulse_phase'][j]))
                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][0]][ind1b].square_pulse2(pulse_center,self.config['pulse_length'][j],self.config['pulse_height'][j]*sin(self.config['pulse_phase'][j]))
                
                #hittite switch
                awg_seq.marker[awg_seq.channel_index['hittite_switch']][ind4].square_pulse2(pulse_center,self.config['pulse_length'][j]+self.config['switch_buffer'])
                
                
            if j==2:
                
                #trigger the card         
                awg_seq.marker[awg_seq.channel_index['card_trig']][ind2].square_pulse2(pulse_center,self.config['acq_length'])
            
                #led switch
                awg_seq.marker[awg_seq.channel_index['led_trigger']][ind3].square_pulse2(self.config['total_length']-100.0,100.0)
                
            
           
                
        #load the agilent awg
        for i in [1]:
            #load into agilent, but do not sequence!
            awg_seq.load_full_into_agilent(self.config['AWG'][1],True)
            awg=im[self.config['AWG'][i]]
            awg.prep_experiment()
            awg.set_amps_offsets(self.config['awg_amp'][i],self.config['awg_offset'][i],self.config['awg_marker_amps'][i])
            awg.run()
            
            #pause until the awg is loaded
            while not awg.query("*OPC?"):
                pass
            
        #turn on the master trigger
            mast_trigger.set_output(True)
        ###TO DO
        
        #copied
        if plotter is not None:
            for j in range(2):
                
                plotter.init_plot("Scope" + str(j),rank=1,accum=False)
                plotter.init_plot("Readout"+ str(j),rank=1,accum=False)
                plotter.init_plot("2DData"+ str(j),rank=2,accum=False)               
        
        
        self.tau_pts = linspace(self.config['tau_start'],self.config['tau_end'],self.config['tau_pts'])
        
        if self.config['save_single_shot']:
            self.esr_data = zeros((2,self.config['tau_pts'],self.config['acq_length'],self.config['num_avgs']))
        elif self.config['save_fulldata']:
            self.esr_data = zeros((2,self.config['tau_pts'],self.config['acq_length']))
        else:
            self.esr_data = zeros((2,self.config['tau_pts']))
            
        
        ch_pts = zeros((2,self.config['acq_length']))
        
        for i in range(self.config['tau_pts']):
            
            #turn off master trigger
            ##TO DO
            awg_trigger.set_period(self.tau_pts[i]*1.0e-6)
            
            #turn on master trigger 
            ###TO DO
            
            tpts,ch_pts[0],ch_pts[1]=card.acquire_avg_data()
                
            for k in range(2):
                
                if self.config['save_fulldata']:
                    self.esr_data[k][i] = ch_pts[k]
                else:
                    
                    
                    amp = mean(ch_pts[k][self.config['data_window'][k][0]:self.config['data_window'][k][1]])
                    
                    self.esr_data[k][i] = amp
                   
                    if plotter is not None:
                        plotter.plot(ch_pts,'Scope'+str(k))
                        plotter.plot((self.tau_pts[0:j],self.esr_data[k][0:j]),"Readout"+str(k))
                            
                            
            if plotter is not None:
                for k in range(2):
                    if self.config['qubit_enable'][k] and not self.config['save_fulldata']:
                        plotter.plot((self.esr_data[k]),"2DData"+str(k))
                        
    def save_data(self,exp_path,file_prefix,overwrite=False):
        
        #print expt_path+fname
        try:
            if overwrite:
                f=SlabFile(exp_path+file_prefix+'.h5','w')
            else:
                fname=get_next_filename(exp_path,prefix=file_prefix,suffix='.h5')
                print fname
                f=SlabFile(exp_path+fname,'w')
        except NameError:
            print "Error opening file for saving esr data"
            return
            
        #save
        f.create_dataset('esr_data',data=self.esr_data)
        f.create_dataset('tau_pts', data= self.tau_pts)
        
        f.attrs['qubit_data_type'] = "esr"
        
        
        #convert config dictionary into saveable format (get rid of nested lists)
        b = self.config.copy()
        for i in b.keys():
            if asarray(b[i]).dtype == "object":
                for j in range(len(b[i])):
                    b[i + str(j)] = b[i][j]

                b[i] = []


        #load keys (allow for new keys to be added later)
        f.save_settings(b)
        
        f.close()        
       
           
if __name__=="__main__":
    
    expt_path="S:\\_Data\\130405 - Updated Qubit Characterization\\static qubit data\\"    
    prefix="vac_rabi"
    
    #take_vacuum_rabi_data_script(expt_path,prefix)
        