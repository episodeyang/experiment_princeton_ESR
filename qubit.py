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

import util as util
from liveplot import LivePlotClient
from data_cache import dataCacheProxy

class mckay_master():

    def __init__(self):

        self._config = dict(
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

        # im = InstrumentManager()

        # if self._config['num_avgs'] <= 0:
        #     raise NameError("Number of Averages must be a Positive Number")
        #
        # #if self._config['LO'][0]==self._config['RF'][0]:
        homodyne_meas = True
        # print "Homodyne measurement"
        # coup = 'DC'
        # #else:
        # #    homodyne_meas=False
        # #    coup = 'AC'

        # if self._config['num_avgs']==1:
        #     rec_per_buff = 1
        # else:
        #     rec_per_buff = 512


        rec_per_buff = 512
        coup = 'DC'
        # load card
        card__config={'clock_edge': 'rising', 'clock_source': 'reference',
            'trigger_coupling': 'DC',  'trigger_operation': 'or',
            'trigger_source2': 'disabled','trigger_level2': 1.0, 'trigger_edge2': 'rising',
            'trigger_level1': 0.6, 'trigger_edge1': 'rising', 'trigger_source1': 'external',
            'ch1_enabled': True,  'ch1_coupling':coup,'ch1_range':self._config['meas_range'][0], 'ch1_filter': False,
            'ch2_enabled': True, 'ch2_coupling': coup,'ch2_range':self._config['meas_range'][1], 'ch2_filter': False,
            'bufferCount': 10,'recordsPerBuffer': rec_per_buff, 'trigger_delay': 0, 'timeout': 1000,
            'samplesPerRecord':self._config['acq_length'],'recordsPerAcquisition': self._config['num_avgs'], 'sample_rate': 1000000
            }

        scope_settings = AlazarConfig(card__config)

        card = Alazar(scope_settings)
        card.configure(scope_settings)

        RF = [None, None]
        LO = [None, None]
        drive = [None, None]

        return card, drive, RF, LO, homodyne_meas

class esrExperiment(mckay_master):

    def attach_instruments(self):
        self.im = InstrumentManager()
        self.na = self.im['NWA']
        self.fridge = self.im['FRIDGE']
        self.rf = self.im['RF_1']
        self.masterTrig = self.im['BNC_1']
        self.awgTrig = self.im['BNC_2']
        self.led = self.im['BNC_3']
        self.awg = self.im['AWG']
        self.srs = self.im['SRS']
        self.zMain = self.im['ZMain']
        self.xShim = self.im['XShim']
        self.yShim = self.im['YShim']
        self.zShim = self.im['ZShim']
        self.alazar = Alazar()

    def __init__(self, expt_path=None, prefix=None, alazarConfig=None, newDataFile=False):
        
        mckay_master.__init__(self)
        
        self._config['exp_id'] = 'ESR'
        
        self.esr_data = []

        self.tau_pts = []

        self._config['use_awg'] = True
        
        #custom esr options
        self._config['b_field'] = 0.0
        self._config['tau'] = 10.0 #in us
        self._config['led_delay'] = 2000 #in us
        self._config['led_pulse_length'] = 2.0 #in ms
        self._config['pulse_height'] = [1.0,1.0]
        self._config['pulse_length'] = [100,100] #in ns
        self._config['pulse_phase'] = [0,pi/2]
        self._config['master_trigger'] = self._config['led_pulse_length'] + self._config['led_delay']*1e-3 + 3*self._config['tau'] #ms
        self._config['switch_buffer'] = 2000 #ns
        
        self._config['tau_start'] = 1.0
        self._config['tau_end'] = 2.0
        self._config['tau_pts'] = 100
        
        self._config['trigger_instr'] = 'BNC_1'
        self._config['awg_trigger'] = 'BNC_2'
        self._config['led_trigger'] = 'BNC_3'
        
        #there is only 1 "qubit"
        self._config['qubit_enable'][1] = False
        
        #save single shot
        self._config['save_single_shot'] = False

        # Ge's code
        if expt_path != None and prefix != None:
            self.expt_path = expt_path
            self.prefix = prefix

        self.note_maxLength = 79
        self.config = lambda: None

        if alazarConfig != None:
            self.plotter = LivePlotClient()
            self.attach_instruments()

            # self.na.set_default_state()
            # self.na.params = naParams
            # self.na.update = self.updateNWA
            self.na.set_trigger_source('bus')

            self.xShim.Remote()
            self.xShim.set_voltage(10.0)
            self.xShim.set_output()
            self.yShim.Remote()
            self.yShim.set_voltage(10.0)
            self.yShim.set_output()
            self.zShim.Remote()
            self.zShim.set_voltage(10.0)
            self.zShim.set_output()

            # # self.lb.set_output(False)
            # self.lb.set_pulse_ext(False)
            # self.lb.set_mod(False)
            # self.lb.set_power(0)

            self.alazar.configure(AlazarConfig(alazarConfig))
            # self.alazarConfig = alazarConfig
            # self.alazar.config = AlazarConfig(alazarConfig)
            # self.alazar.configure()

            self.configNWA()
            self.na.take_one = self.na_take_one;

            #this is the dataCache attached to the experiment.
            self.dataCache = dataCacheProxy(self, newFile=newDataFile)
            self.filename = self.dataCache.filename
        else:
            print "using esrExperiment as a method proxy."
        self.count = -1
        self.t0 = time.time()

    def take_esr_data(self,plotter):
        
        #initialize the card and generator
        card, drive, RF, LO, homodyne_meas = self.__init_card_and_gens__()
            
        #define instrument manager    
        im = InstrumentManager()
        
        #setup triggers
        
        #master trigger
        mast_trigger = im[self._config['trigger_instr']]
        mast_trigger.set_function('PULSE')
        mast_trigger.set_pulse_width(100e-9)
        mast_trigger.set_period(self._config['master_trigger']*1.0e-3)
        mast_trigger.set_output(False)
        
        #awg trigger
        awg_trigger = im[self._config['awg_trigger']]
        awg_trigger.set_function('PULSE')
        awg_trigger.set_pulse_width(100e-9)
        awg_trigger.set_period(self._config['tau']*1.0e-6)
        awg_trigger.set_burst_cycles(3)
        awg_trigger.set_burst_mode('triggered')
        awg_trigger.set_burst_state('on')
        awg_trigger.set_trigger_source('EXT')
        
        #led trigger
        led_trigger = im[self._config['led_trigger']]
        led_trigger.set_function('PULSE')
        led_trigger.set_pulse_width(self._config['led_pulse_length']*1.0e-3)
        led_trigger.set_burst_cycles(1)
        led_trigger.set_burst_mode('triggered')
        led_trigger.set_amplitude(5.0)
        led_trigger.set_offset(0.0)
        led_trigger.set_burst_state('on')
        
       
        #Load the waveforms for the card, measurement and drive pulses
        #pulse sequence
        awg_seq = pulse_sequence(3,int(ceil(self._config['total_length']*self._config['awg_clocks'][0])),int(ceil(self._config['total_length']*self._config['awg_clocks'][1])),self._config['awg_clocks'][0],self._config['awg_clocks'][1],self._config['awg_delay'])
            
        pulse_center = self._config['total_length']/2.0
        
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
                awg_seq.analogwf[awg_seq.channel_index['drive_I'][0]][ind1a].square_pulse2(pulse_center,self._config['pulse_length'][j],self._config['pulse_height'][j]*cos(self._config['pulse_phase'][j]))
                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][0]][ind1b].square_pulse2(pulse_center,self._config['pulse_length'][j],self._config['pulse_height'][j]*sin(self._config['pulse_phase'][j]))
                
                #hittite switch
                awg_seq.marker[awg_seq.channel_index['hittite_switch']][ind4].square_pulse2(pulse_center,self._config['pulse_length'][j]+self._config['switch_buffer'])
                
                
            if j==2:
                
                #trigger the card         
                awg_seq.marker[awg_seq.channel_index['card_trig']][ind2].square_pulse2(pulse_center,self._config['acq_length'])
            
                #led switch
                awg_seq.marker[awg_seq.channel_index['led_trigger']][ind3].square_pulse2(self._config['total_length']-100.0,100.0)

        #load the agilent awg
        for i in [1]:
            #load into agilent, but do not sequence!
            awg_seq.load_full_into_agilent(self._config['AWG'][1],True)
            awg=im[self._config['AWG'][i]]
            awg.prep_experiment()
            awg.set_amps_offsets(self._config['awg_amp'][i],self._config['awg_offset'][i],self._config['awg_marker_amps'][i])
            awg.run()
            
            #pause until the awg is loaded
            while not awg.query("*OPC?"):
                pass
            
        #turn on the master trigger
            mast_trigger.set_output(True)
        ###TO DO
        
        #copied
        # if plotter is not None:
        #     for j in range(2):
        #
        #         plotter.init_plot("Scope" + str(j),rank=1,accum=False)
        #         plotter.init_plot("Readout"+ str(j),rank=1,accum=False)
        #         plotter.init_plot("2DData"+ str(j),rank=2,accum=False)
        
        self.tau_pts = linspace(self._config['tau_start'],self._config['tau_end'],self._config['tau_pts'])
        
        if self._config['save_single_shot']:
            self.esr_data = zeros((2,self._config['tau_pts'],self._config['acq_length'],self._config['num_avgs']))
        elif self._config['save_fulldata']:
            self.esr_data = zeros((2,self._config['tau_pts'],self._config['acq_length']))
        else:
            self.esr_data = zeros((2,self._config['tau_pts']))
            
        
        ch_pts = zeros((2,self._config['acq_length']))
        
        for i in range(self._config['tau_pts']):
            
            #TODO: turn off master trigger
            awg_trigger.set_period(self.tau_pts[i]*1.0e-6)
            
            #TODO: turn on master trigger
            
            tpts,ch_pts[0],ch_pts[1]=card.acquire_avg_data()
                
            for k in range(2):
                
                if self._config['save_fulldata']:
                    self.esr_data[k][i] = ch_pts[k]
                else:
                    
                    
                    amp = mean(ch_pts[k][self._config['data_window'][k][0]:self._config['data_window'][k][1]])
                    
                    self.esr_data[k][i] = amp
                   
                    if plotter is not None:
                        plotter.plot(ch_pts,'Scope'+str(k))
                        plotter.plot((self.tau_pts[0:j],self.esr_data[k][0:j]),"Readout"+str(k))

            # if plotter is not None:
            #     for k in range(2):
            #         if self._config['qubit_enable'][k] and not self._config['save_fulldata']:
            #             plotter.plot((self.esr_data[k]),"2DData"+str(k))
                        
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
        
        
        #convert _config dictionary into saveable format (get rid of nested lists)
        b = self._config.copy()
        for i in b.keys():
            if asarray(b[i]).dtype == "object":
                for j in range(len(b[i])):
                    b[i + str(j)] = b[i][j]

                b[i] = []


        #load keys (allow for new keys to be added later)
        f.save_settings(b)
        f.close()

    def configNWA(self, params=None):
        if params != None:
            print "now load parameters for network analyzer"
            self.na.set_averages(params['avg'])
            self.na.set_power(params['power'])
            self.na.set_center_frequency(params['center'])
            self.na.set_span(params['span'])
            self.na.set_ifbw(params['ifbw'])
            self.na.set_sweep_points(params['sweep_pts'])
        self.na.set_trigger_source('BUS')
        self.na.set_trigger_average_mode(True)
        self.na.set_timeout(10000)
        self.na.set_format('MLOG')

    def na_take_one(self, plotName='na spectrum'):
        """Setup Network Analyzer to take a single averaged trace and grab data,
        returning fpts, mags, phases"""
        self.na.clear_averages()
        self.na.trigger_single()
        self.na.averaging_complete()
        ans = self.na.read_data()
        if plotName == None:
            plotName = 'na spectrum';
        self.plotter.append_z(plotName, ans[1])
        return ans

if __name__=="__main__":
    
    expt_path="S:\\_Data\\130405 - Updated Qubit Characterization\\static qubit data\\"    
    prefix="vac_rabi"
    
    #take_vacuum_rabi_data_script(expt_path,prefix)
