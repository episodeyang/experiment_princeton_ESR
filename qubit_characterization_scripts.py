# -*- coding: utf-8 -*-
"""
Created on Sat Apr 06 23:34:06 2013

@author: Dr Dave
"""
from qubit import *
from pulse_sequences import *

# flux1_instr = "YOKO"
#flux2_instr = "SRS"

def qubit_pulses_man(exp_path):
    pulse_type = 18

    #1 rabi, 2 T1, 3 Ramsey, 4 ramsey1rabi2, 5 qubitflux, 6 state tomography, 
    #7 pulse test, 8 rabi2. 9 fast flux spectroscopy, 10 tomography basis states
    #11 filter tests, #12 mixer test, #13 process tomography, #14 sideband rabi, #15 blue sideband, #16 tek7k test
    #17 rabi ef, #18 rabi with enhanced readout

    qubit_enable = [True, False]
    #flux = [0.073,0.549] --> qubit 1 = 6.404GHz, qubit2 = 6.202GHz
    #flux = [0.07222, 0.53777] #[0.094, 0.525] #with bias tee #0.598
    #flux = [-0.05, 0.524]
    #flux = [0.04, 0.5024]
    #'''Flux set by current in Amps!'''
    #flux = [0.8685e-3,4.4520e-3] #with divider #Optimized location (6.2322e9, 6.462e9)
    #flux = [1.875e-4, 10.02e-4] #w/o divider'
    #flux = [1.88e-2*2.498,10.02e-2*2.476] #volts
    '''Flux set by voltage in Volts!'''
    flux = [0.0692, 0.355]  #for sideband, qubit_freq = [5.528e9, 5.528e9+4.1918e9]
    #flux = [0.0505,0.6625]
    seq_file = exp_path + "sequence_files\\"

    #if this is -1 then use the predicted qubit freq
    qubit_freq = [6.02e9, 5.0e9 + 4.1904e9]  #for sideband, flux = [0.139,0.6865]
    #qubit_freq = [6.802e9,6.8995e9]
    #qubit_freq = [7.68e9, 6.4071e9]
    delta_freq = [30e6, 40e6]
    delta_freq = []
    read_pwr_type = 0  #0 - med power, 1 - low power

    read_freq = 0
    #med power read (put about 40dB on attenuator)\
    #phase 45, 100

    if read_pwr_type == 0:
        #read_freq = [4.1939e9,4.6546e9]
        read_freq = [4.193e9, 4.65365e9]  #for sideband, flux = [0.139,0.6865]
        #read_freq = [4.191e9,4.6534e9]#[4.19192e9,4.6534e9] 
        #read_freq[0] = 4.194e9+0.5e6

    #low power read (put about 50dB on attenuator 1, 60dB on 2N)\
    #phase 0, 100
    if read_pwr_type == 1:
        read_freq = [4.193e9, 4.653e9 + 0.0e6]

    #high power
    #phase 0, 90
    #correlations?
    #20db on attenuator
    #read_freq = [4.1973e9,4.6584e9] #high power

    do_homodyne = True
    do_sequential = False

    #file name with the cavity and qubit characterizations
    if (read_freq == -1 or qubit_freq == -1):
        fname_qubit_char = "0000_cqed_data.h5"
        qubit1 = cQED()
        qubit1.load_data(exp_path, fname_qubit_char, False)

    if pulse_type == 1:
        pulse_exp = rabi()
    elif pulse_type == 2:
        pulse_exp = t1()
    elif pulse_type == 3:
        pulse_exp = ramsey()
    elif pulse_type == 4:
        pulse_exp = ramsey1rabi2()
    elif pulse_type == 5:
        pulse_exp = qubit_fastflux()
    elif pulse_type == 6 or pulse_type == 10:
        pulse_exp = gate_w_statetomography()
    elif pulse_type == 7:
        pulse_exp = pulse_tests()
    elif pulse_type == 8:
        pulse_exp = rabi2()
    elif pulse_type == 9:
        pulse_exp = qubit_fastflux_spectroscopy()
    elif pulse_type == 11:
        pulse_exp = filter_tests()
    elif pulse_type == 12:
        pulse_exp = mixer_test()
    elif pulse_type == 13:
        pulse_exp = gate_w_process_tomography()
    elif pulse_type == 14:
        pulse_exp = sideband_rabi()
    elif pulse_type == 15:
        pulse_exp = blue_sideband()
    elif pulse_type == 16:
        pulse_exp = awg_test()
    elif pulse_type == 17:
        pulse_exp = rabi_ef()
    elif pulse_type == 18:
        pulse_exp = rabi_enhanced_readout()
    else:
        raise NameError("Invalid Experiment Type")


    #setup standard options    
    #-------------------------
    pulse_exp.config['pulsed_drive'] = [True, True]

    if read_freq == -1:
        pulse_exp.config['read_freq'] = qubit1.get_read_freq([flux])[0] + 3.0e6
    else:
        pulse_exp.config['read_freq'] = read_freq

    if qubit_freq == -1:
        pulse_exp.config['drive_freq'] = qubit1.get_qubit_freq([flux])[0, 0] + 96.4e6
    else:
        pulse_exp.config['drive_freq'] = qubit_freq

    pulse_exp.config['flux'] = flux
    pulse_exp.config['flux_type'] = 'voltage'
    pulse_exp.config['qubit_enable'] = qubit_enable
    pulse_exp.config['meas_range'] = [1.0, 1.0]
    pulse_exp.config['drive_pwr'] = [10.0, 10.0]

    if do_homodyne:
        pulse_exp.config['LO'] = pulse_exp.config['RF']
        pulse_exp.config['read_pwr'] = [16.0, 16.0]

    pulse_exp.config['save_each_run'] = False
    pulse_exp.config['save_each_run_fulldata'] = False

    if pulse_exp.config['save_each_run'] or pulse_exp.config['save_each_run_fulldata']:
        pulse_exp.config['num_avgs'] = 5000

    #for now with Lab brick, need to set marker channel 3, marker 1 to 2.0V
    pulse_exp.config['awg_marker_amps'][0][4] = 2.0
    pulse_exp.config['awg_marker_amps'][0][2] = 2.0
    pulse_exp.config['drive'][1] = 'LB2'  #'LB1'
    pulse_exp.config['drive'][0] = 'LB1'
    #pulse_exp.config['drive'][0] = 'RF1'
    #pulse_exp.config['RF'][1] = 'LB2'
    if read_pwr_type == 0:
        pulse_exp.config['data_window'] = [[200, 850], [200, 950]]
        #pulse_exp.config['data_window'] = [[600,700],[550,650]]
        #pulse_exp.config['data_window'] = [[400,500],[400,500]]
        #pulse_exp.config['data_window'] = [[200,1550],[200,1500]]
    else:
        pulse_exp.config['data_window'] = [[200, 1000], [200, 1000]]

    print pulse_exp.config['data_window']

    pulse_exp.config['read_delay'] = -100 / 1.2 * 0  #100
    pulse_exp.config['card_delay'] = 300 / 1.2 + pulse_exp.config['read_delay']  #300+
    pulse_exp.config['read_length'] = 2048
    pulse_exp.config['acq_length'] = 1024

    pulse_exp.config['awg_amp'][0] = [0.7, 0.7, 0.8, 0.8]

    #old settings (qubit 1- 6.2, qubit 2-6.4)
    #BNC-845 offsets    
    #pulse_exp.config['awg_offset'][0][0] = 0.008 + 0.000#-0.081 #0.008 #6GHz
    #pulse_exp.config['awg_offset'][0][1] = 0.015-0.000 #0.098 #0.016 #6GHz

    #settings for qubits at 6.2,6.4    
    pulse_exp.config['awg_offset'][0][0] = 0.004 + 0.000  #-0.081 #0.008 #6GHz
    pulse_exp.config['awg_offset'][0][1] = 0.017 - 0.000  #0.098 #0.016 #6GHz
    pulse_exp.config['awg_offset'][0][2] = -0.002 + 0.000  #-0.002 #6GHz
    pulse_exp.config['awg_offset'][0][3] = 0.021 - 0.000  #0.021 #6GHz

    #settings for blue sideband, LB3(qubit 1) at 5.0GHz, LB1(sideband) at 9.1904GHz
    pulse_exp.config['awg_offset'][0][0] = -0.003
    pulse_exp.config['awg_offset'][0][1] = -0.029
    pulse_exp.config['awg_offset'][0][2] = 0.015
    pulse_exp.config['awg_offset'][0][3] = 0.011

    #    pulse_exp.config['awg_offset'][0][0] = 0.004#-0.081 #0.008 #6GHz
    #    pulse_exp.config['awg_offset'][0][1] = -0.015 #0.098 #0.016 #6GHz
    #    pulse_exp.config['awg_offset'][0][2] = 0.000 #-0.002 #6GHz
    #    pulse_exp.config['awg_offset'][0][3] = 0.000 #0.021 #6GHz

    #    pulse_exp.config['awg_offset'][0][0] = 0.008 + 0.000#-0.081 #0.008 #6GHz
    #    pulse_exp.config['awg_offset'][0][1] = 0.015-0.000 #0.098 #0.016 #6GHz
    #    pulse_exp.config['awg_offset'][0][2] = 0.002+0.000 #-0.002 #6GHz
    #    pulse_exp.config['awg_offset'][0][3] = 0.015 #0.021 #6GHz

    #pulse_exp.config['awg_pulse_offset'] = [-0.081,0.098,0,0]
    #fine adjust offset
    pulse_exp.config['awg_pulse_offset'] = [0 * -0.001, 0 * -0.001, 0 * -0.001, 0 * -0.001]

    pulse_exp.config['drive_sideband'] = [False, False]
    pulse_exp.config['drive_sideband_freq'] = [50.0e6, -50.0e6]

    if pulse_exp.config['drive_sideband'][1]:
        #+50MHz (sideband) changes the offsets for qubit 2 mixer
        pulse_exp.config['awg_offset'][0][2] = 0.001  #6GHz
        pulse_exp.config['awg_offset'][0][3] = 0.017  #6GHz

    #pulse_exp.config['awg_delay'][4:5] = [300,300]

    #--------------------


    #setup pulse specific options  
    #-------------------------
    if pulse_type == 1:

        seq_file = seq_file + "rabi.awg"

        pulse_exp.config['PHASE'] = ['LBPHASE1', '']
        pulse_exp.config['read_phase'] = [40.0, 0.0]
        pulse_exp.config['read_freq'] = [read_freq[0], read_freq[1]]

        pulse_exp.config['rabi_vary'] = "width"

        #pulse_exp.config['drive_sideband'] = [True,True]
        #pulse_exp.config['drive_sideband_freq'] = [50.0e6,-50.0e6]

        pulse_exp.config['rabi_height'] = [0.2, 0.5]
        pulse_exp.config['rabi_length'] = [20.0, 20.0]
        pulse_exp.config['drag_pulse'] = [False, False]
        pulse_exp.config['height_scaling'] = [1.0, 1.0]

        if pulse_exp.config['rabi_vary'] == "height":
            pulse_exp.xdata = linspace(0, 0.6, 61)
        elif pulse_exp.config['rabi_vary'] == "width":
            pulse_exp.xdata = linspace(5, 300, 60)
        else:
            raise NameError("Invalid rabi vary")

        #uncomment to test generator null
        #pulse_exp.config['generator_enable_delay'] = 3000

        pulse_exp.config['monitor_pulses'] = False
        pulse_exp.config['shaped_pulse'] = True

        pulse_exp.config['drive_angle'] = pi / 2 * 0

        pulse_exp.config['start_agilent'] = False

        pulse_exp.config['num_avgs'] = 100000



    elif pulse_type == 2:

        seq_file = seq_file + "t1.awg"

        pulse_exp.config['PHASE'] = ['LBPHASE1', '']
        pulse_exp.config['read_phase'] = [90.0, 0.0]

        pulse_exp.config['t1_pulse_height'] = [0.5, 0.28]
        pulse_exp.config['t1_pulse_length'] = [20, 20]

        pulse_exp.xdata = linspace(100, 10000, 100)
        pulse_exp.config['total_length'] = 8192 * 3

        pulse_exp.config['start_agilent'] = False


    elif pulse_type == 3:

        ramsey_freq = +50.0e6

        seq_file = seq_file + "ramsey.awg"

        #pulse_exp.config['read_delay'] = -100+200
        #pulse_exp.config['card_delay'] = 200+300+pulse_exp.config['read_delay']
        #pulse_exp.config['read_length'] = 2048

        #pulse_exp.config['drive_sideband'] = [True,True]
        #pulse_exp.config['drive_sideband_freq'] = [150.0e6,150.0e6]

        pulse_exp.config['PHASE'] = ['LBPHASE1', '']
        pulse_exp.config['read_phase'] = [0.0, 0.0]

        #optimized!
        pulse_exp.config['ramsey_pulse_height'] = [0.5, 0.5]
        pulse_exp.config['ramsey_pulse_height2'] = pulse_exp.config['ramsey_pulse_height']
        pulse_exp.config['ramsey_pulse_length'] = [10.0, 10.0]
        pulse_exp.config['phase_advance_rate'] = [0.1 * 0, 0.1 * 0]
        pulse_exp.config['ramsey_time'] = [30.0, 30.0]
        #pulse_exp.config['ramsey_vary'] = 'time2'

        pulse_exp.config['drive_freq'] = array(pulse_exp.config['drive_freq']) + ramsey_freq

        #ram_obj.config['ramsey_pulse_height'] = 0.256
        #ram_obj.config['ramsey_pulse_length']= 50

        pulse_exp.config['spin_echo'] = False
        pulse_exp.config['echo_pulse_length'] = array(pulse_exp.config['ramsey_pulse_length']) * 2
        pulse_exp.config['echo_pulse_height'] = array(pulse_exp.config['ramsey_pulse_height'])

        start_late = 0.0
        if pulse_exp.config['spin_echo']:
            pulse_exp.config['phase_advance_rate'] = [0.02, 0.02]
            start_late = 1.0

        #ram_obj.ramsey_time = linspace(10,1000,100)
        #note: make sure the starting time leaves enough space for the pulses
        pulse_exp.config['ramsey_vary'] = 'time'
        pulse_exp.xdata = linspace(max(pulse_exp.config['ramsey_time']) + start_late * 200, 100, 100)

        if pulse_exp.config['ramsey_vary'] == 'time2':
            pulse_exp.config['ramsey_time'][1] = 1000.

        pulse_exp.config['monitor_pulses'] = False
        pulse_exp.config['total_length'] = 8192 * 2
        pulse_exp.config['awg_amp'][0] = [0.8 * 1.05, 0.72 * 1.05, 0.8, 0.95]

        pulse_exp.config['do_fast_flux'] = False
        pulse_exp.config['ramsey_fast_flux'] = [-0.5 * 0, -0.5 * 0]
        pulse_exp.config['flux_wait_time_percent'] = [0.25,
                                                      0.25]  #fraction of the wait time that that the flux pulse is applied
        pulse_exp.config['flux_up_time'] = [25,
                                            25]  #time to the fast flux pulse height (for smooth square...ignored if doing trap pulses)
        pulse_exp.config['trap_flux_pulse'] = [True, True]  #Do trapezoidal pulse
        pulse_exp.config['trap_slope'] = [-0.02, -0.02]  #Slope of the trap pulse
        pulse_exp.config['awg_amp'][1] = [0.15, 0.15]

        #Remember to turn this on if we are doing the fast flux!
        pulse_exp.config['start_agilent'] = False


    elif pulse_type == 4:

        seq_file = seq_file + "ramsey1rabi2.awg"

        pulse_exp.config['read_delay'] = -100 + 200
        pulse_exp.config['card_delay'] = 200 + 100 + pulse_exp.config['read_delay']
        pulse_exp.config['read_length'] = 2048

        #optimized!
        pulse_exp.config['ramsey_pulse_height'] = 0.3
        pulse_exp.config['ramsey_pulse_length'] = 24.0
        pulse_exp.config['phase_advance_rate'] = 0.002
        pulse_exp.config['ramsey_time'] = 200.0

        #rabi height for qubit 2
        pulse_exp.xdata = linspace(0, 0.9, 50)

        pulse_exp.config['monitor_pulses'] = False
        pulse_exp.config['total_length'] = 8192
        pulse_exp.config['awg_amp'][0] = [0.8, 0.9, 1.4, 1.8]
        pulse_exp.config['awg_offset'][0][0] = 0.004  #6GHz
        pulse_exp.config['awg_offset'][0][1] = 0.015  #6GHz
        pulse_exp.config['awg_offset'][0][2] = 0.005  #6GHz
        pulse_exp.config['awg_offset'][0][3] = 0.014  #6GHz
        #ram_obj.config['awg_offset'][0] = 0.011 #8GHz
        #ram_obj.config['awg_offset'][1] = 0.003 #8GHz
        #ram_obj.config['awg_offset'][2] = 0.11

        pulse_exp.config['do_fast_flux'] = True
        pulse_exp.config['flux_height'] = [-0.9, -0.9]
        pulse_exp.config['flux_time'] = [40.0, 20.0]  #total time for the flux pulse
        pulse_exp.config['flux_up_time'] = [30,
                                            30]  #time to the fast flux pulse height (for smooth square...ignored if doing trap pulses)
        pulse_exp.config['trap_flux_pulse'] = [False, False]  #Do trapezoidal pulse
        pulse_exp.config['trap_slope'] = [-0.015, -0.015]  #Slope of the trap pulse

        pulse_exp.config['awg_amp'][1] = [1.5, 1.5]

    elif pulse_type == 5 or pulse_type == 6 or pulse_type == 9 or pulse_type == 10 or pulse_type == 13:

        #common options
        #optimized!
        pulse_exp.config['first_pulse_height'] = [0.169, 0.229]
        pulse_exp.config['first_pulse_length'] = [10.0, 10.0]
        pulse_exp.config['second_pulse_height'] = [0.16, 0.229]
        pulse_exp.config['second_pulse_length'] = [10.0, 10.0]
        pulse_exp.config['second_pulse_phase'] = [0, 0]
        #pulse_exp.config['pulse_wait'] = 180.0

        #tomography settings
        pulse_exp.config['tomo_pulse_length'] = 10.0

        #first list is the X(I),Y(Q) pulse for qubit 1...second list for qubit 2
        pulse_exp.config['tomo_pulse_height'] = [[0.16, 0.2131], [0.2389, 0.2706]]
        pulse_exp.config['tomo_pulse_height_neg'] = [[0.16, 0.219], [0.2127, 0.2386]]

        #first list is the X(I),Y(Q) pulse for qubit 1...second list for qubit 2
        #no offsets
        #pulse_exp.config['tomo_pulse_height'] = [[0.165,0.2166],[0.224,0.253]]
        #pulse_exp.config['tomo_pulse_height_neg'] = pulse_exp.config['tomo_pulse_height']

        #pulse used to calibrate the read out (I channel!)
        pulse_exp.config['tomo_pi_pulse'] = [0.334, 0.464]

        pulse_exp.config['tomo_wait'] = 5.0

        pulse_exp.config['monitor_pulses'] = False
        pulse_exp.config['total_length'] = 8192

        #pulse_exp.config['drag_pulse'] = [False,False]

        #uncomment if using TEK for flux
        #pulse_exp.config['awg_amp'][0][1] = 2.0
        #pulse_exp.config['awg_offset'][0][1] = 0.0

        do_second_pulse_rabi = False

        #this is for characterizing the second pulse fidelity with state tomography        
        if do_second_pulse_rabi:
            pulse_exp.config['first_pulse_height'] = [0.0, 0.0]
            #pulse_exp.config['second_pulse_height'][1] = 0.469
            #pulse_exp.config['second_pulse_height'][0] = 0.0
            pulse_exp.config['second_pulse_height'][0] = 0.1695 * 0  #I - pi/2
            #pulse_exp.config['second_pulse_height'][0] = 0.1672 #Q - pi/2
            pulse_exp.config['second_pulse_height'][0] = 0.334  #I - pi
            pulse_exp.config['second_pulse_phase'][0] = pi / 2 * 0



        #flux pulse properties
        pulse_exp.config['do_fast_flux'] = True
        pulse_exp.config['flux_height'] = [-0.99, -0.99]
        pulse_exp.config['pedestal_height'] = [-0.083, -0.]  #[-0.097-0.1,-0.112*0]
        pulse_exp.config['flux_pulse_type'] = [3, 3]  #2: trap pulse, 3: tan pulse
        pulse_exp.config['trap_time'] = [62.5 + 10.0, 5.8 + 60.0]  #total time for the flux pulse
        #pulse_exp.config['trap_time'] = [120.0,57.0] #total time for the flux pulse

        #extra long pulse to look at cphase vs time
        pulse_exp.config['trap_time'] = [300., 70.]

        #trapezoid flux options
        pulse_exp.config['trap_slope_up'] = (array(pulse_exp.config['flux_height']) - array(
            pulse_exp.config['pedestal_height'])) / 20  #Slope of the trap pulse
        pulse_exp.config['trap_slope_down'] = pulse_exp.config['trap_slope_up']

        #tan pulse options 
        #pulse_exp.config['tan_pulse_up_duration'] = [[25.0,5.0],[8.0,8.0]]
        pulse_exp.config['tan_pulse_up_duration'] = [[25.0, 5.0], [25.0, 5.0]]
        pulse_exp.config['tan_pulse_slopes'] = [[0.2, 0.2], [0.2, 0.2]]
        pulse_exp.config['tan_pulse_midpt'] = [-0.6, -0.6]

        #60ns pulse comp
        pulse_exp.config['comp_pulse_height'] = [[-0.0095, -0.013, -0.016], [-0.004, -0.006, -0.008]]
        pulse_exp.config['comp_pulse_length'] = [[200.0, 1500.0], [300.0, 1500.0]]

        pulse_exp.config['comp_pulse_height'] = [[-0.002, -0.014, -0.015], [-0.000, -0.004, -0.005]]
        pulse_exp.config['comp_pulse_length'] = [[90.0, 1500.0], [50.0, 1500.0]]

        pulse_exp.config['comp_pulse_height'] = [[-0.001, -0.007, -0.012, -0.013, -0.014], [0.002, -0.004, -0.005]]
        pulse_exp.config['comp_pulse_length'] = [[29.0, 35.0, 150.0, 1400.0], [50.0, 1500.0]]

        pulse_exp.config['comp_pulse_height'] = [[-0.000, -0.004, -0.011, -0.014, -0.015, -0.016],
                                                 [0.002, -0.002, -0.008, -0.012]]
        pulse_exp.config['comp_pulse_length'] = [[13.0, 23.0, 40.0, 150.0, 1400.0], [50.0, 150.0, 1200.0]]


        #pulse_exp.config['comp_pulse_height'] = [[-0.000,-0.007,-0.0015,-0.015],[0.002,-0.004,-0.005]]
        #pulse_exp.config['comp_pulse_length'] = [[35.0,40.0,1400.0],[50.0,1500.0]]

        #pulse_exp.config['comp_pulse_height'] = [[-0.009,-0.012,-0.015,-0.016],[0.002,-0.004,-0.005]]
        #pulse_exp.config['comp_pulse_length'] = [[50.0,150.0,1400.0],[50.0,1500.0]]

        #72ns pulse compensation
        #pulse_exp.config['comp_pulse_height'] = [[-0.005,-0.037,-0.042],[-0.003,-0.017,-0.02]]

        #300ns pulse compensation
        pulse_exp.config['comp_pulse_height'] = [[-0.12, -0.135, -0.15], [-0.01, -0.03, -0.04]]
        pulse_exp.config['comp_pulse_length'] = [[200., 1500.], [300., 1500.]]

        pulse_exp.config['comp_pulse_height'] = array(pulse_exp.config['comp_pulse_height'])

        #120ns pulse compensation
        #pulse_exp.config['comp_pulse_height'] = [[-0.036,-0.043,-0.048],[-0.01,-0.03,-0.04]]


        #pulse_exp.config['comp_pulse_height'] = [[-0.03,-0.045,-0.05],[-0.01,-0.015,-0.02]]
        #pulse_exp.config['comp_pulse_height'] = [[-0.000,-0.0,-0.0],[0.0,0.0,0.0]]        
        #pulse_exp.config['comp_pulse_length'] = [[100.0,1400.0],[300.0,1500.0]]

        #pulse_exp.config['comp_pulse_height'][1] = array(pulse_exp.config['comp_pulse_height'][1])*0        

        pulse_exp.config['pedestal_time'] = [5.0, 0.0]
        #pulse_exp.config['pedestal_time'] = [10.0,5.0]

        pulse_exp.config[
            'pulse_wait'] = 350.  # pulse_exp.config['trap_time'][0] + 2*pulse_exp.config['pedestal_time'][0]+ max(pulse_exp.config['first_pulse_length'])+ max(pulse_exp.config['second_pulse_length'])+5.0
        #pulse_exp.config['pulse_wait'] += 20.0        
        print pulse_exp.config['pulse_wait']

        pulse_exp.config['awg_amp'][1] = [2.0, 2.0]

        pulse_exp.config['flux_height'][1] /= 2.0  #3.9
        pulse_exp.config['tan_pulse_midpt'][1] /= 2.0  #3.9

        pulse_exp.config['flux_height'][1] *= -1.0  #3.9
        pulse_exp.config['tan_pulse_midpt'][1] *= -1.0  #3.9

        #pulse_exp.config['first_pulse_height'][1]=0
        #pulse_exp.config['second_pulse_height'] = [0.0,0.0]


        if pulse_type == 5:

            seq_file = seq_file + "qubitfastflux.awg"

            pulse_exp.config['load_agilent'] = True
            pulse_exp.config['start_agilent'] = True
            pulse_exp.config['num_avgs'] = 100000

            flux_step = 2  #1,2,3,4

            #step 1 is to vary the pedestal height for 
            #qubit 1 to zero out the single qubit phase
            if flux_step == 1:
                #pulse_exp.config['pedestal_height'][1] = 0
                #pulse_exp.config['flux_height'][1] = 0
                pulse_exp.config['first_pulse_height'][1] = 0.45
                pulse_exp.config['second_pulse_height'][1] = 0.0
                pulse_exp.config['fastflux_vary'] = 'pedestal1'
                #pulse_exp.xdata = linspace(1.0,30.0,50)
                pulse_exp.xdata = linspace(pulse_exp.config['pedestal_height'][0] + 0.02,
                                           pulse_exp.config['pedestal_height'][0] - 0.02, 50)

            #step 2 is to vary the length of the qubit 2 pulse to change the conditional phase of qubit 2
            #after doing a pi pulse
            elif flux_step == 2:
                pulse_exp.config['qubit2_ramsey'] = False
                if not pulse_exp.config['qubit2_ramsey']:
                    pulse_exp.config['first_pulse_height'][1] = 0.0 + 0.0 * pulse_exp.config['first_pulse_height'][1]
                    pulse_exp.config['second_pulse_height'][1] = 0.0
                pulse_exp.config['fastflux_vary'] = 'length2'
                pulse_exp.xdata = linspace(pulse_exp.config['trap_time'][1] - 10,
                                           pulse_exp.config['trap_time'][1] + 150, 100)
                #pulse_exp.config['save_each_run'] = True

            #step 3 is to vary the height of the first pulse for qubit 2 (i.e. a Rabi on qubit 2)
            elif flux_step == 3:
                pulse_exp.config['second_pulse_height'][1] = 0.0
                pulse_exp.config['fastflux_vary'] = 'pulse_height2'
                pulse_exp.xdata = linspace(0.0, 0.5, 40)

            #step 4 to take single shot data
            elif flux_step == 4:
                pulse_exp.config['second_pulse_height'][1] = 0.0
                #pulse_exp.config['first_pulse_height'][1] = 0.0
                #pulse_exp.config['second_pulse_height'][0] = 0.0
                pulse_exp.config['fastflux_vary'] = 'pulse_height1'
                pulse_exp.xdata = linspace(0.4, 0.5, 5)
                pulse_exp.config['save_each_run'] = True
                pulse_exp.config['num_calib_seqs'] = 0


        elif pulse_type == 6 or pulse_type == 10 or pulse_type == 13:

            seq_file = seq_file + "statetomo.awg"



            #set to true for no flux pulsing!
            null_flux = False
            if null_flux:
                pulse_exp.config['flux_height'] = [0, 0]
                pulse_exp.config['pedestal_height'] = [0, 0]
                pulse_exp.config['comp_pulse_height'] = [0, 0]

            pulse_exp.config['load_agilent'] = True
            pulse_exp.config['start_agilent'] = True

            pulse_exp.config['num_avgs'] = 10000

            pulse_exp.config['expanded_tomo'] = True

            do_calib = False

            if pulse_type == 10:

                do_calib = False
                pulse_exp.config['tomo_pi_pulse'] = [0.245, 0.296]
                pulse_exp.config['tomo_vary'] = 'basis'
                pulse_exp.config['num_avgs'] = 1000000
                pulse_exp.config['start_agilent'] = False
                pulse_exp.config['load_agilent'] = False
                pulse_exp.config['save_each_run_fulldata'] = False
                pulse_exp.config['save_each_run'] = True

            elif pulse_type == 6:
                if do_calib:
                    pulse_exp.config['load_agilent'] = True
                    pulse_exp.config['start_agilent'] = True
                    #calibrate the single qubit pulses
                    #calib1: calibrate the I channel tomography pulses
                    #calib2: calibrate the Q channel tomography pulses
                    #calib3: calibrate the second gate pulse
                    #calib3b: calibrate the second gate pulse (Q channel - for process tomography)
                    #calib3c: second gate pulse 0,+,- versus time after the flux pulse                    
                    #calib4: calibrate the first gate pulse
                    #calib4b: calibrate the first gate pulse with no flux
                    #calib4c: claibrate the first gate pulse (Q channel - for process tomography)

                    pulse_exp.config['tomo_vary'] = 'calib4'
                    #pulse_exp.xdata = linspace(0,0.8,20)
                    pulse_exp.xdata = linspace(-0.6, 0.6, 100)

                    if pulse_exp.config['tomo_vary'] == 'calib3c':
                        pulse_exp.xdata = linspace(0, 100, 10)
                        pulse_exp.config['tomo_wait'] = pulse_exp.xdata[-1]
                    #pulse_exp.xdata = list(linspace(0,0.45,24))
                    #pulse_exp.xdata.append(0.263) #qubit 1 pi
                    #pulse_exp.xdata.append(0.339) #qubit 2 pi    

                    #pulse_exp.config['flux_height'][0] = 0.0

                    #pulse_exp.xdata = (list(linspace(0.235,0.25,24)))
                    #pulse_exp.xdata.insert(0,0.0)
                    #pulse_exp.xdata = [0,0.128,0.263,0.165,0.339,0.0]            
                    #pulse_exp.config['save_each_run'] = True
                    pulse_exp.config['num_avgs'] = 20000

                else:

                    #uncomment to do a trivial pi/2 on both
                    #pulse_exp.config['first_pulse_height'] = [0.0,0.0]
                    #pulse_exp.config['first_pulse_height'][0] = 0.0
                    #pulse_exp.config['first_pulse_height'][1] = 0.43
                    #pulse_exp.config['second_pulse_height'] = [0.172,0.223]
                    #pulse_exp.config['second_pulse_height'] = [0.0,0.0]
                    #pulse_exp.config['second_pulse_phase'] = [0,0]

                    pulse_exp.config['tomo_vary'] = 'tomo'
                    pulse_exp.config['save_each_run'] = True
                    pulse_exp.config['num_avgs'] = 50000

            elif pulse_type == 13:

                pulse_exp.config['num_avgs'] = 200000
                pulse_exp.config['first_pulse_pi'] = [0.3497, 0.4251]  #pi pulse for the first pulse (I channel)
                pulse_exp.config['first_pulse_Q'] = [0.1946, 0.2491]
                pulse_exp.config['exp_path'] = exp_path


        elif pulse_type == 9:

            after_pulse = True
            before_pulse = False
            seq_file = seq_file + "qubitfastfluxspect.awg"
            pulse_exp.config['num_avgs'] = 10000

            #pulse_exp.config['flux_height'][1] = -0.2

            #uncomment for filter adiabaticity test settings
            #            pulse_exp.config['flux_height'] = [-0.99,-0.0]
            #            pulse_exp.config['pedestal_height'] = [0.0,0.0]
            #            pulse_exp.config['pedestal_time'] = [0.0,0.0]
            #            pulse_exp.config['tan_pulse_up_duration'] = [[10.0,5.0],[10.0,5.0]]
            #            pulse_exp.config['tan_pulse_slopes'] = [[0.01,0.01],[0.01,0.01]]
            #            pulse_exp.config['tan_pulse_midpt'] = [-0.6,-0.0]
            #            pulse_exp.config['flux_pulse_type'] = [3,3]
            #            pulse_exp.config['trap_time'] = [110.0,80.0]
            #            pulse_exp.config['comp_pulse_height'] = [[-0.016,-0.02,-0.03],[-0.003,-0.008,-0.01]]
            #            pulse_exp.config['comp_pulse_length'] = [[200.0,1500.0],[300.0,1500.0]]
            #            pulse_exp.config['comp_pulse_height'][0][0] = 2.5*pulse_exp.config['comp_pulse_height'][0][0]
            #            pulse_exp.config['comp_pulse_height'][0][1:3] = 2.5*array(pulse_exp.config['comp_pulse_height'][0][1:3])

            #pulse_exp.config['comp_pulse_height'] = 0*array(pulse_exp.config['comp_pulse_height'])

            if after_pulse:

                num_freq_steps = 24 + 1
                pulse_exp.config['second_pulse_height'] = [0.3, 0.3]
                pulse_exp.config['second_pulse_length'] = [10.0, 10.0]
                pulse_exp.xdata = linspace(-10, 500, 25)
                for j in range(2):
                    if num_freq_steps == 1:
                        #Use this to tweak the compensation
                        pulse_exp.config['start_agilent'] = True
                        pulse_exp.config['num_avgs'] = 10000
                        pulse_exp.xdata = linspace(-10, 100, 50)
                        pulse_exp.config['drive_frequencies'][j] = [pulse_exp.config['drive_freq'][j] + 0 * 500e6]
                    else:
                        pulse_exp.config['drive_frequencies'][j] = pulse_exp.config['drive_freq'][j] + array(
                            linspace(-400e6, 150e6, num_freq_steps))

            elif before_pulse:

                num_freq_steps = 50
                pulse_exp.config['second_pulse_height'] = [0.2, 0.2]
                pulse_exp.xdata = linspace(-100, 10, 50)
                for j in range(2):
                    pulse_exp.config['drive_frequencies'][j] = pulse_exp.config['drive_freq'][j] + array(
                        linspace(-50e6, 1000e6, num_freq_steps))


            else:

                #do spectroscopy *during* the flux pulse            
                #set the frequencies for the spectroscopy

                num_freq_steps = 100
                pulse_exp.config['second_pulse_height'] = [0.05, 0.0]

                pulse_exp.config['flux_pulse_type'][0] = 2
                pulse_exp.config['trap_slope_up'] = [-0.99 / 1.0, -0.99 / 1.0]
                pulse_exp.config['trap_slope_down'] = pulse_exp.config['trap_slope_up']
                pulse_exp.config['probe_pulse_time'] = -85.0

                pulse_exp.config['flux_height'][1] = 0.0
                pulse_exp.config['flux_pulse_type'][1] = 2

                pulse_exp.config['flux_vary'] = 'time'

                if pulse_exp.config['flux_vary'] == 'flux1':
                    pulse_exp.config['pedestal_height'][0] = 0.0

                if pulse_exp.config['flux_vary'] == 'time':
                    pulse_exp.xdata = linspace(-150, 10, 50)
                    #pulse_exp.xdata = linspace(-70,-50,10) 
                else:

                    pulse_exp.xdata = linspace(-0.99, -0.02, 50)

                for j in range(2):
                    pulse_exp.config['drive_frequencies'][j] = pulse_exp.config['drive_freq'][j] + array(
                        linspace(-50e6, 50e6, num_freq_steps))


                    #pulse_exp.config['drive_frequencies'][1] = pulse_exp.config['drive_freq'][1] + array(linspace(0e6,0e6,num_freq_steps))





    elif pulse_type == 7:

        seq_file = seq_file + "pulsetest.awg"

        #pulse_exp.config['pulsed_drive'] = False

        #pulse_exp.config['data_window'] = [[200,850],[200,900]]

        #optimized!
        pulse_exp.config['pulse_length'] = 10.0
        pulse_exp.config['pulse_phase'] = pi / 2
        pulse_exp.config['pulse_gap'] = 10.0
        pulse_exp.config['num_pulses'] = 20
        pulse_exp.config['drag_pulses'] = [True, True]
        #pulse_exp.config['pulse_heights']= linspace(0.4,0.5,10)

        pulse_exp.config['pulse_heights'] = array([-0.008, -0.004, 0, 0.004, 0.008]) - 0.51
        #pulse_exp.config['pulse_heights']= array([-0.006,-0.004,-0.002,0.0,0.002,0.004,0.006])+0.156
        pulse_exp.config['pulse_factor'] = 1
        #pulse_exp.config['pulse_heights']= array([0.447])
        pulse_exp.config['pulse_heights'] = list(pulse_exp.config['pulse_heights'])
        #pulse_exp.config['pulse_heights'].insert(0,0.312)
        #pulse_exp.config['pulse_heights'].insert(0,0.0)

        pulse_exp.config['drag_prefactor'] = [-2.0 / 4.0 / (2 * pi * 0.2), -2.0 / 4.0 / (2 * pi * 0.2)]

        #pulse_exp.config['pulse_heights'] = array(pulse_exp.config['pulse_heights'])/2.0


        pulse_exp.config['monitor_pulses'] = False
        pulse_exp.config['total_length'] = 8192
        #pulse_exp.config['awg_amp'][0] = [0.8, 0.9, 0.8, 0.9]


        pulse_exp.config['num_avgs'] = 10000
        #pulse_exp.config['num_avgs'] = 100000
        #pulse_exp.config['acq_length'] = 2048

        pulse_exp.config['start_agilent'] = False

    elif pulse_type == 8:

        seq_file = seq_file + "rabi2.awg"

        pulse_exp.config['rabi_vary'] = "width"

        #pulse_exp.config['drive_sideband'] = [True,True]
        #pulse_exp.config['drive_sideband_freq'] = [50.0e6,-50.0e6]

        pulse_exp.config['PHASE'] = ['LBPHASE1', '']
        pulse_exp.config['read_phase'] = [60.0, 0.0]

        pulse_exp.config['rabi_height1'] = [1.0, 1.0]
        pulse_exp.config['rabi_height2'] = [1.0, 1.0]
        pulse_exp.config['rabi_length1'] = [25.0, 25.0]
        pulse_exp.config['rabi_length2'] = [50.0, 50.0]
        pulse_exp.config['rabi_sideband_freq1'] = [-233e6, 0.0e6]
        pulse_exp.config['rabi_sideband_freq2'] = [0.0e6, 230.0e6]
        pulse_exp.config['pulse1_2_gap'] = [0.0, 0.0]
        pulse_exp.config['read_sweep'] = False
        pulse_exp.config['sideband_sweep'] = False
        pulse_exp.config['read_freq'] = [read_freq[0], read_freq[1]]
        pulse_exp.config['awg_delay'] = [0, -10.0 / 2, -8.0 / 2, -6.0 / 2, 320, 340]

        if pulse_exp.config['read_sweep']:
            pulse_exp.config['read_sweep_end'] = [pulse_exp.config['read_freq'][0] + 15e6,
                                                  pulse_exp.config['read_freq'][1]]
            pulse_exp.config['read_sweep_points'] = [26, 26]
            pulse_exp.config['SB_read_points'] = [
                linspace(pulse_exp.config['read_freq'][0], pulse_exp.config['read_sweep_end'][0],
                         pulse_exp.config['read_sweep_points'][0]),
                linspace(pulse_exp.config['read_freq'][1], pulse_exp.config['read_sweep_end'][1],
                         pulse_exp.config['read_sweep_points'][1])]

        if pulse_exp.config['sideband_sweep']:
            pulse_exp.config['sideband_sweep_end'] = [pulse_exp.config['rabi_sideband_freq2'][0] + 15e6,
                                                      pulse_exp.config['rabi_sideband_freq2'][1]]
            pulse_exp.config['sideband_sweep_points'] = [16, 16]
            pulse_exp.config['SB_points'] = [
                linspace(pulse_exp.config['rabi_sideband_freq2'][0], pulse_exp.config['sideband_sweep_end'][0],
                         pulse_exp.config['sideband_sweep_points'][0]),
                linspace(pulse_exp.config['rabi_sideband_freq2'][1], pulse_exp.config['sideband_sweep_end'][1],
                         pulse_exp.config['sideband_sweep_points'][1])]

        #different data window if trying to see |1>-->|2> oscillations
        pulse_exp.config['data_window'] = [[200, 850], [200, 950]]

        for j in range(2):
            pulse_exp.config['drive_freq'][j] -= pulse_exp.config['rabi_sideband_freq1'][j]

            #pulse_exp.config['rabi_vary'] = 'height' #h

        if pulse_exp.config['rabi_vary'] == "height":
            pulse_exp.xdata = linspace(0.0, 1.0, 101)
        elif pulse_exp.config['rabi_vary'] == "width":
            pulse_exp.xdata = linspace(5, 300, 296)
        else:
            raise NameError("Invalid rabi vary")

        pulse_exp.config['monitor_pulses'] = False

        pulse_exp.config['num_avgs'] = 100000

    elif pulse_type == 11:

        seq_file = seq_file + "filtertest.awg"

        exp_path = exp_path + "filter_adiabaticity\\"

        #specify the type of filter test
        pulse_exp.config['filter_vary'] = 'adiabatic_calib2'

        pulse_exp.config['first_pulse_length'] = [15.0, 15.0]
        pulse_exp.config['second_pulse_length'] = [15.0, 15.0]

        #this does a pi pulse on qubit 1
        pulse_exp.config['first_pulse_height'][0] = 0.245
        pulse_exp.config['second_pulse_height'][0] = 0.245

        #this does the probe pulse
        pulse_exp.config['second_pulse_height'][1] = 0.0

        #this is the qubit 2 probe frequencies
        #pulse_exp.config['drive_freq'][1] = 6.683e9

        #pulse_exp.config['load_agilent'] = False
        pulse_exp.config['awg_amp'][1] = [2.0, 2.0]

        pulse_exp.config['flux_height'] = [-0.99, -0.0]
        pulse_exp.config['pedestal_height'] = [-0.087 * 0, 0.0]
        pulse_exp.config['pedestal_time'] = [5.0 * 0, 0.0]

        pulse_exp.config['tan_pulse_up_duration'] = [[22.0, 5.0], [22.0, 5.0]]
        pulse_exp.config['tan_pulse_slopes'] = [[0.2, 0.2], [0.2, 0.2]]
        pulse_exp.config['tan_pulse_midpt'] = [-0.6, -0.6]

        #only trapezoid pulses are allowed
        pulse_exp.config['flux_pulse_type'] = [3, 3]

        pulse_exp.config['load_agilent'] = True

        #time for the flux pulse
        #note: when testing adiabaticity this has to be at least twice the slowest ramp time
        pulse_exp.config['trap_time'] = [110.0, 80.0]
        #pulse_exp.config['trap_time'] = [72.5,65.8]

        pulse_exp.config['pulse2_delay'] = (pulse_exp.config['trap_time'][0] - pulse_exp.config['trap_time'][1]) / 2.0

        pulse_exp.config['pulse_wait'] = pulse_exp.config['trap_time'][0] + 2 * pulse_exp.config['pedestal_time'][
            0] + 100.0


        #pulse_exp.config['comp_pulse_height'] = [[-0.005,-0.037,-0.042],[-0.003,-0.008,-0.01]]
        pulse_exp.config['comp_pulse_height'] = [[-0.016, -0.02, -0.03], [-0.003, -0.008, -0.01]]
        pulse_exp.config['comp_pulse_length'] = [[200.0, 1500.0], [300.0, 1500.0]]

        #pulse_exp.config['comp_pulse_height'] = [[-0.005,-0.037,-0.042],[-0.003,-0.017,-0.02]]

        num_pts = 150

        #trapezoid flux options
        pulse_exp.config['trap_slope_down'] = (array(pulse_exp.config['flux_height']) - array(
            pulse_exp.config['pedestal_height'])) / 40
        pulse_exp.config['trap_slope_up'][1] = (
                                               pulse_exp.config['flux_height'][1] - pulse_exp.config['pedestal_height'][
                                                   1]) / 5
        pulse_exp.config['trap_slope_up'][0] = (
                                               pulse_exp.config['flux_height'][0] - pulse_exp.config['pedestal_height'][
                                                   0]) / 40

        if pulse_exp.config['filter_vary'] == 'adiabatic2' or pulse_exp.config['filter_vary'] == 'adiabatic1':
            #pulse_exp.xdata = linspace(1,15,50)*pulse_exp.config['trap_slope_up'][0]
            pulse_exp.xdata = pulse_exp.config['flux_height'][0] / array(linspace(0.5, 54, num_pts))
            pulse_exp.xdata = array(linspace(0.5, 49.0, num_pts))

            pulse_exp.config['comp_scaling'] = (1.0 + (49.0 - pulse_exp.xdata) * 1.5 / 39.0) * 1.00
            #pulse_exp.config['comp_scaling'] = 1.0+(25.0-pulse_exp.xdata)*1.5/50.0

        elif pulse_exp.config['filter_vary'] == 'adiabatic_calib' or pulse_exp.config[
            'filter_vary'] == 'adiabatic_calib2':
            pulse_exp.xdata = linspace(0, 0.4, 30)
            pulse_exp.config['comp_pulse_height'][0] = array(pulse_exp.config['comp_pulse_height'][0]) * (
            1.0 + (25.0 - pulse_exp.config['tan_pulse_up_duration'][0][0]) * 1.5 / 50.0)

        elif pulse_exp.config['filter_vary'] == 'lifetime':

            pulse_exp.config['first_pulse_height'][0] = 0.25
            pulse_exp.config['second_pulse_height'][1] = 0.25
            pulse_exp.config['pedestal_height'] = [0.0, 0.0]
            pulse_exp.config['pedestal_time'] = [0.0, 0.0]
            pulse_exp.config['trap_time'] = [50.0, 50.0]
            pulse_exp.config['flux_pulse_type'] = [2, 2]
            pulse_exp.config['total_length'] = 8192 * 3
            pulse_exp.config['comp_pulse_height'] = [[-0.000, -0.00, -0.00], [-0.003, -0.017, -0.032]]
            pulse_exp.config['comp_pulse_height'] = [[-0.000, -0.00, -0.00], [-0.005, -0.037, -0.042]]

            pulse_exp.xdata = linspace(100, 7000, 10)

    elif pulse_type == 12:

        seq_file = seq_file + "mixer.awg"

        pulse_exp.config['pulse1_height'] = [0.115, 0.139]
        pulse_exp.config['pulse2_height'] = pulse_exp.config['pulse1_height']

        pulse_exp.config['mixer_vary'] = 'qpulse'

        if pulse_exp.config['mixer_vary'] == 'pulse1':
            pulse_exp.xdata = linspace(-0.4, 0.4, 50)
        elif pulse_exp.config['mixer_vary'] == 'qpulse':
            pulse_exp.xdata = linspace(-pi / 2, pi / 2, 50)

        pulse_exp.config['start_agilent'] = False

    elif pulse_type == 14:

        seq_file = seq_file + "sideband_rabi.awg"

        #pulse_exp.config['drive_pwr'][1] = 7.0



        pulse_exp.config['sideband_vary'] = "flux_width"
        #'pi_width','pi_height','delay',flux_height','flux_width'

        pulse_exp.config['pi_pulse_height'] = [0.45 * 0, 0.45 * 0]
        pulse_exp.config['pi_pulse_width'] = [10.0, 10.0]
        pulse_exp.config['pi_to_flux_delay'] = [10.0, 10.0]
        pulse_exp.config['flux_pulse_height'] = [0.1 * 0, 0.2]
        pulse_exp.config['flux_pulse_width'] = [50.0, 50.0]
        pulse_exp.config['sideband_freq'] = [0.0, 0.0]  #0.615,0.788,0.9438 #GHz
        pulse_exp.config['drag_pulse'] = [False, False]
        pulse_exp.config['freq_sweep'] = False
        pulse_exp.config['height_sweep'] = False
        #parameter options: 'pi_width','pi_height','delay',flux_height','flux_width'

        if pulse_exp.config['freq_sweep']:
            pulse_exp.config['sideband_sweep_end'] = [1.7 - 1.213, 0.0]
            pulse_exp.config['sideband_sweep_points'] = [11, 11]
            pulse_exp.config['SB_freq_points'] = [
                linspace(pulse_exp.config['sideband_freq'][0], pulse_exp.config['sideband_sweep_end'][0],
                         pulse_exp.config['sideband_sweep_points'][0]),
                linspace(pulse_exp.config['sideband_freq'][1], pulse_exp.config['sideband_sweep_end'][1],
                         pulse_exp.config['sideband_sweep_points'][1])]
        if pulse_exp.config['height_sweep']:
            pulse_exp.config['height_sweep_end'] = [0.0, 0.5]
            pulse_exp.config['height_sweep_points'] = [11, 11]
            pulse_exp.config['SB_height_points'] = [
                linspace(pulse_exp.config['flux_pulse_height'][0], pulse_exp.config['height_sweep_end'][0],
                         pulse_exp.config['height_sweep_points'][0]),
                linspace(pulse_exp.config['flux_pulse_height'][1], pulse_exp.config['height_sweep_end'][1],
                         pulse_exp.config['height_sweep_points'][1])]

        if pulse_exp.config['sideband_vary'] == "pi_height":
            pulse_exp.xdata = linspace(-0.7, 0.7, 100)
        elif pulse_exp.config['sideband_vary'] == "pi_width":
            pulse_exp.xdata = linspace(5, 101, 100)
        elif pulse_exp.config['sideband_vary'] == "flux_height":
            pulse_exp.xdata = linspace(-0.99, 0.99, 50)
        elif pulse_exp.config['sideband_vary'] == "flux_width":
            pulse_exp.xdata = linspace(5, 200, 79)
        elif pulse_exp.config['sideband_vary'] == "delay":
            pulse_exp.xdata = linspace(0, 100, 100)
        else:
            raise NameError("Invalid sideband vary")

        #uncomment to test generator null
        #pulse_exp.config['generator_enable_delay'] = 3000

        pulse_exp.config['monitor_pulses'] = False

        pulse_exp.config['drive_angle'] = pi / 2 * 0

        pulse_exp.config['start_agilent'] = True

        #enable this to run for a fixed number of averages        
        pulse_exp.config['num_avgs'] = 20000

    elif pulse_type == 15:

        seq_file = seq_file + "blue_sideband.awg"

        #pulse_exp.config['drive_pwr'][1] = 7.0

        pulse_exp.config['read_pwr'][1] = -30.0

        pulse_exp.config['sideband_vary'] = "pi_width"
        #'pi_width','pi_height','delay',flux_height','flux_width'
        pulse_exp.config['total_length'] = 8192
        pulse_exp.config['pi_pulse_height'] = [0.1, 0.45 * 0]
        pulse_exp.config['pi_pulse_width'] = [25.0, 25.0]
        pulse_exp.config['pi_to_flux_delay'] = [0.0, 0.0]
        pulse_exp.config['flux_pulse_height'] = [0.1 * 0, 0.25 * 0]
        pulse_exp.config['flux_pulse_width'] = [1.0, 1.0]
        pulse_exp.config['sideband_freq'] = [0.0, 0.0]  #0.615,0.788,0.9438 #GHz
        pulse_exp.config['drag_pulse'] = [False, False]
        pulse_exp.config['freq_sweep'] = True
        pulse_exp.config['read_sweep'] = False
        pulse_exp.config['height_sweep'] = False
        pulse_exp.config['measure_e'] = [False, False]
        pulse_exp.config['drive_freq'] = [qubit_freq[0] - 10e6, qubit_freq[1]]
        pulse_exp.config['read_freq'] = [read_freq[0], read_freq[1]]
        pulse_exp.config['PHASE'] = ['LBPHASE1', '']
        pulse_exp.config['read_phase'] = [40.0, 0.0]
        #parameter options: 'pi_width','pi_height','delay',flux_height','flux_width'

        if pulse_exp.config['freq_sweep']:
            pulse_exp.config['freq_sweep_end'] = [pulse_exp.config['drive_freq'][0] + 20e6,
                                                  pulse_exp.config['drive_freq'][1] + 10e6]
            pulse_exp.config['freq_sweep_points'] = [11, 11]
            pulse_exp.config['SB_freq_points'] = [
                linspace(pulse_exp.config['drive_freq'][0], pulse_exp.config['freq_sweep_end'][0],
                         pulse_exp.config['freq_sweep_points'][0]),
                linspace(pulse_exp.config['drive_freq'][1], pulse_exp.config['freq_sweep_end'][1],
                         pulse_exp.config['freq_sweep_points'][1])]
        if pulse_exp.config['read_sweep']:
            pulse_exp.config['read_sweep_end'] = [pulse_exp.config['read_freq'][0] + 2e6,
                                                  pulse_exp.config['read_freq'][1]]
            pulse_exp.config['read_sweep_points'] = [11, 11]
            pulse_exp.config['SB_read_points'] = [
                linspace(pulse_exp.config['read_freq'][0], pulse_exp.config['read_sweep_end'][0],
                         pulse_exp.config['read_sweep_points'][0]),
                linspace(pulse_exp.config['read_freq'][1], pulse_exp.config['read_sweep_end'][1],
                         pulse_exp.config['read_sweep_points'][1])]
        if pulse_exp.config['height_sweep']:
            pulse_exp.config['height_sweep_end'] = [0.0, 0.5]
            pulse_exp.config['height_sweep_points'] = [25, 25]
            pulse_exp.config['SB_height_points'] = [
                linspace(pulse_exp.config['flux_pulse_height'][0], pulse_exp.config['height_sweep_end'][0],
                         pulse_exp.config['height_sweep_points'][0]),
                linspace(pulse_exp.config['flux_pulse_height'][1], pulse_exp.config['height_sweep_end'][1],
                         pulse_exp.config['height_sweep_points'][1])]

        if pulse_exp.config['sideband_vary'] == "pi_height":
            pulse_exp.xdata = linspace(-0.7, 0.7, 100)
        elif pulse_exp.config['sideband_vary'] == "pi_width":
            pulse_exp.xdata = linspace(5, 300, 60)
        elif pulse_exp.config['sideband_vary'] == "flux_height":
            pulse_exp.xdata = linspace(0.01, 1.0, 100)
        elif pulse_exp.config['sideband_vary'] == "flux_width":
            pulse_exp.xdata = linspace(5, 100, 96)
        elif pulse_exp.config['sideband_vary'] == "delay":
            pulse_exp.xdata = linspace(0, 1000, 101)
        else:
            raise NameError("Invalid sideband vary")

        #uncomment to test generator null
        #pulse_exp.config['generator_enable_delay'] = 3000

        pulse_exp.config['monitor_pulses'] = False

        pulse_exp.config['drive_angle'] = pi / 2 * 0

        pulse_exp.config['start_agilent'] = False
        pulse_exp.config['load_agilent'] = False

        #enable this to run for a fixed number of averages        
        pulse_exp.config['num_avgs'] = 40000
        pulse_exp.config['num_calib_seqs'] = 0

    elif pulse_type == 16:

        seq_file = seq_file + "awgtest.awg"

        pulse_exp.config['AWG'][1] = 'TEK2'

        pulse_exp.config['front_buffer'] = 7000

        pulse_exp.config['awg_delay'][4:5] = [0., 0.]

        pulse_exp.config['pulse_vary'] = "width"

        #pulse_exp.config['drive_sideband'] = [True,True]
        #pulse_exp.config['drive_sideband_freq'] = [50.0e6,-50.0e6]
        pulse_exp.config['pi_pulse_height'] = [0.15, 1.0]
        pulse_exp.config['pi_pulse_width'] = [75.0, 75.0]
        pulse_exp.config['pi_pulse_freq'] = 5.528e9 + 4.1918e9
        pulse_exp.config['pulse_height'] = [0.15, 0.99]
        pulse_exp.config['pulse_length'] = [10.0, 10.0]

        pulse_exp.config['height_scaling'] = [1.0, 1.0]

        pulse_exp.config['pulse_freq'] = [5.528e9 - 4.1918e9 + 4e6]

        pulse_exp.config['total_length'] = 16384

        #pulse_exp.config['pulse_freq'] = [pulse_exp.config['drive_freq'][0]+20e6,pulse_exp.config['drive_freq'][0]-20e6]

        #pulse_exp.config['pulse_freq'] = array([pulse_exp.config['drive_freq'][0],pulse_exp.config['drive_freq'][0]-218e6,pulse_exp.config['drive_freq'][0]-425e6])

        pulse_exp.config['pulse_freq_offset'] = 0e6

        if pulse_exp.config['pulse_vary'] == "height":
            pulse_exp.xdata = linspace(0, 1.0, 51)
        elif pulse_exp.config['pulse_vary'] == "width":
            pulse_exp.xdata = linspace(5, 1000, 200)
        elif pulse_exp.config['pulse_vary'] == "delay":
            pulse_exp.xdata = linspace(50, 200, 100)
        else:
            raise NameError("Invalid rabi vary")

        pulse_exp.config['start_agilent'] = True

        pulse_exp.config['num_avgs'] = 100000

    elif pulse_type == 17:

        seq_file = seq_file + "rabi_ef.awg"

        pulse_exp.config['rabi_vary'] = "height"

        #pulse_exp.config['drive_sideband'] = [True,True]
        #pulse_exp.config['drive_sideband_freq'] = [50.0e6,-50.0e6]

        pulse_exp.config['pi_pulse_height'] = [0.15, 1.0]
        pulse_exp.config['pi_pulse_width'] = [75.0, 75.0]
        pulse_exp.config['pi_pulse_freq'] = 5.528e9 + 4.1918e9
        pulse_exp.config['rabi_height'] = [0.4, 0.4]
        pulse_exp.config['rabi_length'] = [10.0, 10.0]
        pulse_exp.config['drag_pulse'] = [False, False]
        pulse_exp.config['height_scaling'] = [1.0, 1.0]

        if pulse_exp.config['rabi_vary'] == "height":
            pulse_exp.xdata = linspace(0, 1.0, 101)
        elif pulse_exp.config['rabi_vary'] == "width":
            pulse_exp.xdata = linspace(5, 200, 79)
        else:
            raise NameError("Invalid rabi vary")

        #uncomment to test generator null
        #pulse_exp.config['generator_enable_delay'] = 3000

        pulse_exp.config['monitor_pulses'] = False
        pulse_exp.config['shaped_pulse'] = True

        pulse_exp.config['drive_angle'] = pi / 2 * 0

        pulse_exp.config['start_agilent'] = False

        pulse_exp.config['num_avgs'] = 20000

    elif pulse_type == 18:

        seq_file = seq_file + "rabi_enhanced_readout.awg"

        pulse_exp.config['PHASE'] = ['LBPHASE1', '']
        pulse_exp.config['read_phase'] = [40.0, 0.0]
        pulse_exp.config['read_freq'] = [read_freq[0], read_freq[1]]
        pulse_exp.config['rabi_vary'] = "height"

        pulse_exp.config['total_length'] = 8192
        #pulse_exp.config['drive_sideband'] = [True,True]
        #pulse_exp.config['drive_sideband_freq'] = [50.0e6,-50.0e6]

        pulse_exp.config['rabi_height'] = [0.1, 0.5]
        pulse_exp.config['rabi_length'] = [15.0, 15.0]
        pulse_exp.config['drag_pulse'] = [False, False]
        pulse_exp.config['height_scaling'] = [1.0, 1.0]
        pulse_exp.config['sideband_freq'] = [1000e6, 0.0]
        pulse_exp.config['sideband_height'] = [0.5, 0.0]

        if pulse_exp.config['rabi_vary'] == "height":
            pulse_exp.xdata = linspace(0, 0.6, 61)
        elif pulse_exp.config['rabi_vary'] == "width":
            pulse_exp.xdata = linspace(5, 300, 60)
        else:
            raise NameError("Invalid rabi vary")

        #uncomment to test generator null
        #pulse_exp.config['generator_enable_delay'] = 3000

        pulse_exp.config['monitor_pulses'] = False
        pulse_exp.config['shaped_pulse'] = True

        pulse_exp.config['drive_angle'] = pi / 2 * 0

        pulse_exp.config['load_agilent'] = True
        pulse_exp.config['start_agilent'] = True

        pulse_exp.config['num_avgs'] = 100000

    #-----------------------------

    pulse_exp.config['seq_file'] = seq_file

    if not pulse_type == 13:

        #load into awg
        if not do_sequential:
            s = raw_input("Load Pulses into AWG (Y/N)?")

        if not do_sequential and s.upper() == "Y":
            print "Loading Data"
            pulse_exp.load_exp()

        if not do_sequential:
            if pulse_type == 9:
                pulse_exp.run_qubit_exp(ScriptPlotter())
            else:
                pulse_exp.run_qubit_exp(ScriptPlotter())
        else:
            pulse_exp.run_qubit_exp_sequential(1024 * 10, ScriptPlotter())

        s = raw_input("Save Data (Y/N)?")

        if s.upper() != "N":
            print "Saving Data"
            pulse_exp.save_data(exp_path, pulse_exp.config['exp_id'], overwrite=False)

    else:

        #pulse type 13 (process tomography does autosave)
        pulse_exp.run_qubit_exp(ScriptPlotter())


#this is a specialized vacuum rabi for optimizing the read frequency/phase/power
def optimize_read_script(exp_path, prefix="vac_rabi_pwr7"):
    plotter = ScriptPlotter()

    seq_file = exp_path + "sequence_files\\vac_rabi.awg"

    a = vacuum_rabi()
    a.config['qubit_enable'] = [True, True]

    #flux is fixed!!
    a.config['flux_start'] = [0.062, 0.549]
    a.config['flux_end'] = a.config['flux_start']
    a.config['flux_pts'] = 1

    #set this based on the current attenuation!
    a.config['read_pwr2'] = [-27.0, -30.0]

    #best initial guess
    a.config['read_freq'] = [4.19338e9, 4.65316e9]

    #qubit frequency
    a.config['drive_freq'] = [6.170e9 + 28.0e6, 6.413e9 - 11.5e6 + 1.0e6]

    #set frequency sweep
    a.config['start_freq'] = array(a.config['read_freq']) - 3.0e6
    a.config['end_freq'] = array(a.config['read_freq']) + 3.0e6
    a.config['num_fpts'] = 100

    a.config['read_delay'] = -100 / 1.2 * 0  #100
    a.config['card_delay'] = 300 / 1.2 + a.config['read_delay']  #300+
    a.config['read_length'] = 2048
    a.config['acq_length'] = 1024

    a.config['awg_amp'][0] = [0.8, 0.8, 0.8, 0.8]
    a.config['awg_offset'][0][0] = 0.008  #-0.081 #0.008 #6GHz
    a.config['awg_offset'][0][1] = 0.016 - 0.001  #0.098 #0.016 #6GHz
    a.config['awg_offset'][0][2] = -0.001 + 0.000  #-0.002 #6GHz
    a.config['awg_offset'][0][3] = 0.021  #0.021 #6GHz

    #set pulse properties
    a.config['pulse_height'] = [0.24, 0.29]

    phase_array = [0, 22.5, 45.0, 67.5, 90.0, 111.5]
    phase_array = linspace(0, 130.0, 20)

    a.config['num_avgs'] = 1024 * 10
    a.config['seq_file'] = seq_file

    take_window_data = True
    if take_window_data:
        a.config['save_fulldata'] = True
        a.config['start_freq'] = array(a.config['read_freq'])
        a.config['end_freq'] = a.config['start_freq']
        a.config['num_fpts'] = 1
        phase_array = [70]
        prefix = "vac_rabi_fulldata"

    for i in range(2 * len(phase_array)):

        if mod(i, 2):
            a.config['do_pulse'] = [True, True]
        else:
            a.config['do_pulse'] = [False, False]

        a.config['read_phase'] = [phase_array[int(floor(i / 2))], phase_array[int(floor(i / 2))]]

        a.take_read_data(plotter)

        a.save_data(exp_path + "optimize_read\\", prefix, False)


def take_vacuum_rabi_data_script(exp_path, prefix="vac_rabi"):
    plotter = ScriptPlotter()

    seq_file = exp_path + "sequence_files\\vac_rabi.awg"

    a = vacuum_rabi()
    a.config['qubit_enable'] = [True, False]
    '''Flux set by voltage in Volts!'''
    a.config['flux_start'] = [0.07, 0.355]
    a.config['flux_end'] = [0.07, 0.355]
    a.config['flux_pts'] = 10

    #best initial guess
    a.config['read_freq'] = [4.193e9, 4.65316e9]
    a.config['PHASE'] = ['LBPHASE1', '']
    a.config['read_phase'] = [0.0, 0.0]

    #set frequency sweep
    a.config['start_freq'] = array(a.config['read_freq']) - 3e6
    a.config['end_freq'] = array(a.config['read_freq']) + 3e6
    a.config['num_fpts'] = 61

    a.config['data_window'] = [[200, 800], [0, 600]]
    a.config['read_delay'] = -100 / 1.2 * 0  #100
    a.config['card_delay'] = 300 / 1.2 + a.config['read_delay']  #300+
    a.config['read_length'] = 2048
    a.config['acq_length'] = 1024
    a.config['meas_range'] = [0.1, 0.1]

    a.config['awg_marker_amps'][0][4] = 2.0
    a.config['awg_marker_amps'][0][2] = 2.0

    a.config['awg_amp'][0] = [0.7, 0.7, 0.8, 0.8]
    #settings for blue sideband, LB3(qubit 1) at 5.0GHz, LB1(sideband) at 9.19025GHz
    a.config['awg_offset'][0][0] = -0.003
    a.config['awg_offset'][0][1] = -0.029
    a.config['awg_offset'][0][2] = 0.015
    a.config['awg_offset'][0][3] = 0.012

    a.config['do_pulse'] = [True, False]
    a.config['drive'][0] = 'LB1'
    a.config['drive_freq'][0] = 6.02e9
    a.config['drive_pwr'] = [10.0, 10.0]
    a.config['pulse_length'] = [10.0, 10.0]
    a.config['pulse_height'] = [0.4, 0.5]
    a.config['num_avgs'] = 1024 * 20
    a.config['seq_file'] = seq_file

    a.config['cw_meas'] = False

    a.take_read_data(plotter)

    s = raw_input("Save Data (Y/N)?")

    if s.upper() != "N":
        print "Saving Data"
        a.save_data(exp_path, prefix, False)


def take_pulse_probe_data_script(exp_path, prefix="pulse_probe"):
    plotter = ScriptPlotter()

    seq_file = exp_path + "sequence_files\\pulse_probe.awg"

    single_slice = False

    do_full_acq_trace = False  #single_slice also needs to be set to true

    b = pulse_probe_spectroscopy()

    b.config['seq_file'] = seq_file
    b.config['drive_freq'][0] = 5.4e9
    b.config['drive'][0] = 'LB1'
    #flux_ratio = 4.4

    b.config['card_delay'] = -1000 * 0 + 265
    b.config['read_delay'] = b.config['card_delay'] - 200
    b.config['read_length'] = 2000
    b.config['acq_length'] = 1024 * 1
    b.config['PHASE'] = ['LBPHASE1', '']
    b.config['read_phase'] = [110.0, 0.0]

    b.config['qubit_enable'] = [True, False]

    b.config['awg_offset'][0][0] = 0.008 + 0.000  #-0.081 #0.008 #6GHz
    b.config['awg_offset'][0][1] = 0.015 - 0.000  #0.098 #0.016 #6GHz
    b.config['awg_offset'][0][2] = -0.001 + 0.000  #-0.002 #6GHz
    b.config['awg_offset'][0][3] = 0.020  #0.021 #6GHz

    #settings for blue sideband, LB3(qubit 1) at 5.0GHz, LB1(sideband) at 9.19025GHz
    b.config['awg_offset'][0][0] = -0.003
    b.config['awg_offset'][0][1] = -0.029
    b.config['awg_offset'][0][2] = 0.015
    b.config['awg_offset'][0][3] = 0.011
    b.config['awg_amp'][0] = [0.7, 0.7, 0.8, 0.8]
    #b.config['awg_amp'][1.0] = [2,2.0]

    #for now with Lab brick, need to set marker channel 3, marker 1 to 2.0V
    b.config['awg_marker_amps'][0][4] = 2.0
    b.config['awg_marker_amps'][0][2] = 2.0

    b.config['read_pulse_wait'] = 300

    b.config['do_fast_flux'] = False
    b.config['cw_meas'] = False

    #need to add a lot of attenuation to the flux line for this to work!
    b.config['add_flux_noise'] = False
    b.config['awg_amp'][1] = [0.08, 0.08]

    #b.config['drive'][1] = 'RF2'

    if single_slice:

        '''Flux set by current in Amps!'''
        b.config['flux_start'] = [0.1, 0.6865]
        b.config['flux_end'] = b.config['flux_start']
        b.config['flux_pts'] = 1
        b.config['num_fpts'] = 401
        b.config['freq_start'] = [5.57e9, 6.24e9]
        b.config['freq_end'] = [5.61e9, 6.3e9]
        b.config['drive_pulse_height'] = [0.5, 0.5]
        b.config['drive_pulse_length'] = [3000, 3000]
        b.config['generator_enable_delay'] = 250
        b.config['drive_pwr'] = [-30.0, -30.0]
        b.config['total_length'] = 1024 * 10

        if b.config['do_fast_flux']:
            b.config['dc_flux'] = [0.062, 0.549]
            b.config['flux_start'] = [-0.5, -0.5]
            b.config['flux_pulse_length'] = [100, 100]
            b.config['flux_pulse_delay'] = [-30, -30]
            b.config['drive_pulse_length'] = [10, 10]
            b.config['card_delay'] = 200
            b.config['read_delay'] = -100

        if do_full_acq_trace:
            b.config['front_buffer'] = 1000
            b.config['acq_length'] = 1024 * 10
            b.config['total_length'] = b.config['acq_length'] + 2 * 1024
            b.config['card_delay'] = -(b.config['drive_pulse_length'] + b.config['front_buffer'])
            b.config['read_delay'] = -(b.config['drive_pulse_length'] + b.config['front_buffer'] / 2)
            b.config['read_length'] = b.config['acq_length']


    else:

        #        b.config['flux_start'] = [0.035, 0.5036]
        #        b.config['flux_end'] = [0.025, 0.506]
        #        b.config['flux_end'] = b.config['flux_start']

        '''Flux set by current in Amps!'''
        b.config['flux_start'] = [0.07, 0.355]
        b.config['flux_end'] = [0.07, 0.355]
        b.config['flux_pts'] = 1

        #-0.0436683417085 0.52248040201
        #b.config['flux_start'] = [0.033216080402, 0.504028140704]
        #b.config['flux_end'] = [0.0241708542714, 0.506198994975]
        #b.config['flux_start'] = [0.0286934673367, 0.505113567839]
        #b.config['flux_end'] = b.config['flux_start']
        b.config['num_fpts'] = 101
        b.config['freq_start'] = [6.0e9, 6.24e9]
        b.config['freq_end'] = [6.05e9, 6.3e9]
        b.config['drive_pulse_height'] = [0.5, 0.5]
        b.config['drive_pulse_length'] = [3000, 3000]
        #b.config['drive_pulse_length'] = [100,100]
        #b.config['num_drive_pulses'] = 15
        #b.config['wait_time'] = 200.0
        b.config['generator_enable_delay'] = 250
        b.config['drive_pwr'] = [-20.0, -20.0]
        b.config['drive_end_pwr'] = [-20.0, -20.0]
        b.config['total_length'] = 1024 * 10
        #b.config['read_pwr'][0] = 12.0

    #b.config['num_avgs'] = 1024*80
    b.config['num_avgs'] = 1024 * 20

    #b.config['read_pwr'][0] = 16.0

    #load vacuum rabi data to get read peak values
    if (0):
        a = vacuum_rabi()
        a.load_data(exp_path, "vac_rabi_data_qubit1.h5", False)

        flux1_pts = linspace(b.config['flux_start'][0], b.config['flux_end'][0], b.config['flux_pts'])
        print flux1_pts
        flux2_pts = linspace(b.config['flux_start'][1], b.config['flux_end'][1], b.config['flux_pts'])
        print flux2_pts
        read_freqs = a.cavity_peak(flux1_pts) + 0.2e6  #flux=0 1.0001
    else:
        read_freqs = zeros((2, b.config['flux_pts']))
        for i in range(b.config['flux_pts']):
            #read_freqs[0][i] = 4.194e9-0.15e6+.4e6/(b.config['flux_pts']-1)*i
            #read_freqs[0][i] = 4.194e9+0.7e6+.45e6/(b.config['flux_pts']-1)*i
            read_freqs[0][i] = 4.1933e9  #+ 0.5e6*i/(b.config['flux_pts']-1)#4.19192e9
            read_freqs[1][i] = 4.653e9  #+ 2.0e6*i/(b.config['flux_pts']-1)

    do_windowed_scan = False

    if do_windowed_scan:

        b.config['freq_range'] = 500e6

        #Use these qubit parameters
        bb = cQED()
        maxfreq = 9.5e9
        bb.config['EC'] = 160e6
        bb.config['EJ'] = maxfreq ** 2 / 8 / bb.config['EC']
        bb.config['g'] = 70e6
        bb.config['f0'] = 4.1955e9
        bb.config['flux_to_phi'] = pi / (0.2 + 0.49)
        bb.config['flux_volt_offset'] = -(0.49 - 0.2) / 2

        qfreq = around(bb.get_qubit_freq(b.config['flux_pts'][1]), -6)
        print qfreq

        drive_pts = qfreq[0]

    else:

        drive_pts = [None, None]

    im = InstrumentManager()
    for j in range(2):
        im[b.config['flux_instr'][j][0]].set_mode('voltage')
        if b.config['do_fast_flux']:
            im[b.config['flux_instr'][j][0]].set_volt(b.config['dc_flux'][j], b.config['flux_instr'][j][1])
        else:
            im[b.config['flux_instr'][j][0]].set_volt(b.config['flux_start'][j], b.config['flux_instr'][j][1])

    #time.sleep(5) #flux ramp up time
    b.take_spec_data(plotter, read_freqs, drive_pts)

    #im[flux1_instr].set_volt(0) #reset flux
    #im[flux2_instr].set_volt(0,1)
    s = raw_input("Save Data (Y/N)?")

    if s.upper() != "N":
        print "Saving Data"
        b.save_data(exp_path, prefix, False)


def lab_brick_test():
    a = random.uniform(low=5.0e9, high=10.0e9)

    im = InstrumentManager()

    lab_brick = im['LB1']

    lab_brick.set_frequency(a)
    print a, lab_brick.get_frequency()


def compute_qubit_params_script(exp_path):
    a = pulse_probe_spectroscopy()
    a.load_data(exp_path, "pulse_probe.h5", False)
    b = a
    a = vacuum_rabi()
    a.load_data(exp_path, "vac_rabi_data.h5", False)
    a.fit_data(None, True)
    c = a
    a = cQED()
    a.load_qubit_char_data(c, b, 9.85e9)
    a.compute_qubit_parameters(ScriptPlotter(), False)
    print a.config


def test_agilent_awg(blank=False):
    load_full_seq = True

    num_seqs = 20

    awg_seq = pulse_sequence(num_seqs, 8192, 8192 * 3)

    #create a pulse sequence and load into 

    if load_full_seq:

        for i in range(num_seqs):

            for j in range(2):

                if not blank:
                    ch_ind = awg_seq.channel_index['flux'][j]
                    ind = awg_seq.add_analog_pulse(ch_ind)
                    awg_seq.analog_seq_table[ch_ind][i] = ind
                    awg_seq.analogwf[ch_ind][ind].trapezoid(300 + i * 0 * 100, 1 / 100., 1 / 100., 500.0, 1.0)
                    #awg_seq.analogwf[ch_ind][ind].offset(0.001)
                    #awg_seq.analogwf[ch_ind][ind].pulse[0] = 0.0
                    #awg_seq.analogwf[ch_ind][ind].pulse[-1] = 0.0

        awg_seq.load_full_into_agilent('AGILENT811')

    im = InstrumentManager()
    awg = im['AGILENT811']

    awg.prep_experiment()
    awg.run()


def test_tek_awg(exp_path):
    awg_seq = pulse_sequence(10, 8192, 8192 * 3)

    #create a pulse sequence and load into 

    for i in range(10):

        for j in range(2):
            ch_ind = awg_seq.channel_index['drive_I'][j]
            ind = awg_seq.add_analog_pulse(ch_ind)
            awg_seq.analog_seq_table[ch_ind][i] = ind
            awg_seq.analogwf[ch_ind][ind].trapezoid(300 + i * 0 * 100, 1 / 100., 500.0, 1.0)

    awg_seq.load_into_awg(exp_path + "sequence_files\\test.awg", "TEK")

    im = InstrumentManager()
    awg = im['TEK']

    awg.prep_experiment()
    awg.run()


def tek_blank(exp_path):
    awg_seq = pulse_sequence(1, 8192, 8192 * 3)

    #create a pulse sequence and load into 

    awg_seq.analogwf[awg_seq.channel_index['drive_I'][0]][0].offset(0 * 0.001)
    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][0]][0].offset(0 * -0.001)
    awg_seq.analogwf[awg_seq.channel_index['drive_I'][1]][0].offset(0 * -0.001)
    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][1]][0].offset(0 * -0.001)

    awg_seq.load_into_awg(exp_path + "sequence_files\\test.awg", "TEK")

    im = InstrumentManager()

    drive1 = im['LB2']
    drive2 = im['LB1']

    drive1.set_pulse_ext(True)
    drive2.set_pulse_ext(True)

    awg = im['TEK']

    awg.prep_experiment()
    #awg.set_amps_offsets([0.7, 0.7, 0.8, 0.8],[-0.014,-0.046,-0.002,0.020],[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    awg.set_amps_offsets([0.7, 0.7, 0.8, 0.8], [-0.003, -0.029, 0.015, 0.011], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    awg.run()


#Test the different AWG's using a Ramsey z-gate
def generator_test(exp_path):
    #Do a Ramsey on qubit 1 with a flux pulse to test the various AWG's
    qubit_enable = [False, True]
    flux = [0.062, 0.549]  #with bias tee #0.598

    seq_file = exp_path + "sequence_files\\"

    #if this is -1 then use the predicted qubit freq
    qubit_freq = [6.170e9 + 31.0e6, 6.413e9 - 8.5e6]
    read_freq = [4.193e9 + 0.0e6, 4.653e9 + 0.25e6]

    do_homodyne = True
    do_sequential = False
    pulse_exp = ramsey()

    pulse_exp.config['read_freq'] = read_freq
    pulse_exp.config['drive_freq'] = qubit_freq

    pulse_exp.config['flux'] = flux
    pulse_exp.config['qubit_enable'] = qubit_enable
    pulse_exp.config['meas_range'] = [1.0, 1.0]
    pulse_exp.config['drive_pwr'] = [13, 13]

    if do_homodyne:
        pulse_exp.config['LO'] = pulse_exp.config['RF']
        pulse_exp.config['read_pwr'] = [16.0, 16.0]

    pulse_exp.config['save_each_run'] = False
    pulse_exp.config['save_each_run_fulldata'] = False
    pulse_exp.config['num_avgs'] = 10000

    if pulse_exp.config['save_each_run'] or pulse_exp.config['save_each_run_fulldata']:
        pulse_exp.config['num_avgs'] = 5000

    #for now with Lab brick, need to set marker channel 3, marker 1 to 2.0V
    pulse_exp.config['awg_marker_amps'][0][4] = 2.0
    pulse_exp.config['drive'][1] = 'LB1'
    pulse_exp.config['data_window'] = [[300, 900], [200, 900]]

    pulse_exp.config['read_delay'] = -100 / 1.2  #100
    pulse_exp.config['card_delay'] = 300 / 1.2 + pulse_exp.config['read_delay']  #300+
    pulse_exp.config['read_length'] = 2048
    pulse_exp.config['acq_length'] = 1024

    pulse_exp.config['load_agilent'] = True

    #--------------------


    seq_file = seq_file + "ramsey.awg"
    pulse_exp.config['seq_file'] = seq_file

    #optimized!
    pulse_exp.config['ramsey_pulse_height'] = [0.6, 0.6]
    pulse_exp.config['ramsey_pulse_length'] = [10.0, 10.0]
    pulse_exp.config['phase_advance_rate'] = [0.00, 0.00]
    pulse_exp.config['ramsey_time'] = [200.0, 200.0]

    pulse_exp.config['ramsey_vary'] = 'flux'
    pulse_exp.xdata = linspace(-0.9, -0.99, 100)

    pulse_exp.config['monitor_pulses'] = False
    pulse_exp.config['total_length'] = 8192
    pulse_exp.config['awg_amp'][0] = [0.8, 0.9, 0.8, 0.9]
    pulse_exp.config['awg_amp'][0][1] = 2.0
    pulse_exp.config['awg_offset'][0][0] = 0.004  #6GHz
    pulse_exp.config['awg_offset'][0][1] = 0.015 * 0  #6GHz
    pulse_exp.config['awg_offset'][0][2] = 0.005  #6GHz
    pulse_exp.config['awg_offset'][0][3] = 0.014  #6GHz
    #ram_obj.config['awg_offset'][0] = 0.011 #8GHz
    #ram_obj.config['awg_offset'][1] = 0.003 #8GHz
    #ram_obj.config['awg_offset'][2] = 0.11

    pulse_exp.config['do_fast_flux'] = True
    pulse_exp.config['ramsey_fast_flux'] = [pulse_exp.xdata[0], pulse_exp.xdata[0]]
    pulse_exp.config['flux_wait_time_percent'] = [0.5,
                                                  0.5]  #fraction of the wait time that that the flux pulse is applied
    pulse_exp.config['flux_up_time'] = [25,
                                        25]  #time to the fast flux pulse height (for smooth square...ignored if doing trap pulses)
    pulse_exp.config['trap_flux_pulse'] = [True, True]  #Do trapezoidal pulse
    pulse_exp.config['trap_slope'] = [pulse_exp.xdata[0] / 10, pulse_exp.xdata[0] / 10]  #Slope of the trap pulse
    pulse_exp.config['awg_amp'][1] = [0.9, 0.9]

    s = raw_input("Load Sequence (Y/N)?")

    if not do_sequential and s.upper() == "Y":
        print "Loading Data"
        pulse_exp.load_exp()

    pulse_exp.run_qubit_exp(ScriptPlotter())

    s = raw_input("Save Data (Y/N)?")

    if s.upper() != "N":
        print "Saving Data"
        pulse_exp.save_data(exp_path, pulse_exp.config['exp_id'], overwrite=False)


#Test the adiabaticity of crossing the filter
#Use ramsey sequence, but more or less this is a rabi with a flux pulse
def filter_adiabaticity_test(exp_path):
    #Do a Ramsey on qubit 1 with a flux pulse to test the various AWG's
    qubit_enable = [True, False]
    flux = [0.062, 0.549]  #with bias tee #0.598

    seq_file = exp_path + "sequence_files\\"

    #if this is -1 then use the predicted qubit freq
    qubit_freq = [6.170e9 + 31.0e6, 6.413e9 - 8.5e6]
    read_freq = [4.193e9 + 0.0e6, 4.653e9 + 0.25e6]

    do_homodyne = True
    do_sequential = False
    pulse_exp = ramsey()

    pulse_exp.config['read_freq'] = read_freq
    pulse_exp.config['drive_freq'] = qubit_freq

    pulse_exp.config['flux'] = flux
    pulse_exp.config['qubit_enable'] = qubit_enable
    pulse_exp.config['meas_range'] = [1.0, 1.0]
    pulse_exp.config['drive_pwr'] = [13, 13]

    if do_homodyne:
        pulse_exp.config['LO'] = pulse_exp.config['RF']
        pulse_exp.config['read_pwr'] = [16.0, 16.0]

    pulse_exp.config['save_each_run'] = False
    pulse_exp.config['save_each_run_fulldata'] = False
    pulse_exp.config['num_avgs'] = 10000

    if pulse_exp.config['save_each_run'] or pulse_exp.config['save_each_run_fulldata']:
        pulse_exp.config['num_avgs'] = 5000

    #for now with Lab brick, need to set marker channel 3, marker 1 to 2.0V
    pulse_exp.config['awg_marker_amps'][0][4] = 2.0
    pulse_exp.config['drive'][1] = 'LB1'
    pulse_exp.config['data_window'] = [[300, 900], [200, 900]]

    pulse_exp.config['read_delay'] = -100 / 1.2  #100
    pulse_exp.config['card_delay'] = 300 / 1.2 + pulse_exp.config['read_delay']  #300+
    pulse_exp.config['read_length'] = 2048
    pulse_exp.config['acq_length'] = 1024

    pulse_exp.config['load_agilent'] = True

    #--------------------


    seq_file = seq_file + "ramsey.awg"
    pulse_exp.config['seq_file'] = seq_file

    #optimized!
    pulse_exp.config['ramsey_pulse_height'] = [0.6, 0.6]
    pulse_exp.config['ramsey_pulse_height2'] = [0.0, 0.0]
    pulse_exp.config['ramsey_pulse_length'] = [20.0, 20.0]
    pulse_exp.config['phase_advance_rate'] = [0.00, 0.00]
    pulse_exp.config['ramsey_time'] = [400.0, 400.0]

    pulse_exp.config['ramsey_vary'] = 'flux2'
    pulse_exp.xdata = -0.9 / linspace(50, 2, 100)

    pulse_exp.config['monitor_pulses'] = False
    pulse_exp.config['total_length'] = 8192
    pulse_exp.config['awg_amp'][0] = [0.8, 0.9, 0.8, 0.9]
    pulse_exp.config['awg_amp'][0][1] = 2.0
    pulse_exp.config['awg_offset'][0][0] = 0.004  #6GHz
    pulse_exp.config['awg_offset'][0][1] = 0.015 * 0  #6GHz
    pulse_exp.config['awg_offset'][0][2] = 0.005  #6GHz
    pulse_exp.config['awg_offset'][0][3] = 0.014  #6GHz
    #ram_obj.config['awg_offset'][0] = 0.011 #8GHz
    #ram_obj.config['awg_offset'][1] = 0.003 #8GHz
    #ram_obj.config['awg_offset'][2] = 0.11

    pulse_exp.config['do_fast_flux'] = True
    pulse_exp.config['ramsey_fast_flux'] = [-0.9, -0.9]
    pulse_exp.config['flux_wait_time_percent'] = [0.5,
                                                  0.5]  #fraction of the wait time that that the flux pulse is applied
    pulse_exp.config['flux_up_time'] = [25,
                                        25]  #time to the fast flux pulse height (for smooth square...ignored if doing trap pulses)
    pulse_exp.config['trap_flux_pulse'] = [True, True]  #Do trapezoidal pulse
    pulse_exp.config['trap_slope'] = array(pulse_exp.config['ramsey_fast_flux']) / 100  #Slope of the trap pulse
    pulse_exp.config['awg_amp'][1] = [2.0, 2.0]

    s = raw_input("Load Sequence (Y/N)?")

    if not do_sequential and s.upper() == "Y":
        print "Loading Data"
        pulse_exp.load_exp()

    pulse_exp.run_qubit_exp(ScriptPlotter())

    s = raw_input("Save Data (Y/N)?")

    if s.upper() != "N":
        print "Saving Data"
        pulse_exp.save_data(exp_path, pulse_exp.config['exp_id'], overwrite=False)


def take_esr_data(exp_path=''):
    #Script to run Anthony's experiment

    esr_exp1 = esr_exp()

    seq_file = exp_path + "sequence_files\\"
    seq_file = seq_file + "esr.awg"
    esr_exp1.config['seq_file'] = seq_file

    esr_exp1.config['num_avgs'] = 1000
    esr_exp1.config['tau_start'] = 10
    esr_exp1.config['tau_end'] = 10
    esr_exp1.config['tau_pts'] = 1

    esr_exp1.config['drive'][0] = 'RF_1'
    esr_exp1.config['AWG'][1] = 'AWG'

    esr_exp1.take_esr_data(ScriptPlotter())

    s = raw_input("Save Data (Y/N)?")

    if s.upper() != "N":
        print "Saving Data"
        esr_exp1.save_data(exp_path, esr_exp1.config['exp_id'], overwrite=False)


if __name__ == "__main__":

    #qubit script to run
    qubit_script = 1

    #path to save data
    #expt_path="S:\\_Data\\130405 - Updated Qubit Characterization\\fluxsideband\\"
    #expt_path = "S:\\_Data\\131024 - chip T4 - gate data\\"
    expt_path = "S:\\_Data\\140707 - Anthony ESR\\"

    if qubit_script == 1:
        take_esr_data(expt_path)

    #1: Vacuum rabi
    #2: Pulse probe
    #3: Rabi/Ramsey/T1 (manual)
    #4: lab brick test
    #5: Compute the qubit parameters
    #6: Test the agilent AWG
    #7: Test the TEK
    #8: Clear the TEK
    #9: Ramsey to test the AWG
    #10: Filter Adiabaticity test
    #11: Read frequency optimization

    # if qubit_script == 1:
    #     take_vacuum_rabi_data_script(expt_path)
    #
    # elif qubit_script == 2:
    #     take_pulse_probe_data_script(expt_path)
    #
    # elif qubit_script == 3:
    #     qubit_pulses_man(expt_path)
    #
    # elif qubit_script == 4:
    #     lab_brick_test()
    #
    # elif qubit_script == 5:
    #     compute_qubit_params_script(expt_path)
    #
    # elif qubit_script == 6:
    #     test_agilent_awg(True)
    #
    # elif qubit_script == 7:
    #     test_tek_awg(expt_path)
    #
    # elif qubit_script == 8:
    #     tek_blank(expt_path)
    #
    # elif qubit_script == 9:
    #     generator_test(expt_path)
    #
    # elif qubit_script == 10:
    #     filter_adiabaticity_test(expt_path)
    #
    # elif qubit_script == 11:
    #     optimize_read_script(expt_path)


    else:
        print "Not a valid selection!"