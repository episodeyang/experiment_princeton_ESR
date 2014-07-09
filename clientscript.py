# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 23:29:53 2013

@author: Ge
"""
from time import sleep
from slab.instruments import InstrumentManager
from numpy import *
from Experiment import *
from winsound import Beep

def rinse_n_fire(exp0):
    exp0.note("unbias the trap for a second")
    exp0.res.set_volt(-3)
    exp0.trap.set_volt(-3)
    exp0.srs.set_output(1, True)
    exp0.srs.set_output(2, True)
    exp0.note("make sure the probe is off before the baseline")
    exp0.lb.set_output(False)
    time.sleep(1)
    
    exp0.note('firing the filament')
    exp0.res.set_volt(1.5)    
    exp0.trap.set_volt(1.5)
    exp0.fil.fire_filament(400,0.01)
    
    exp0.note("Now wait for cooldown while taking traces")
    while exp0.fridge.get_mc_temperature() > 60e-3 or (time.time() - exp0.t0) < 360:
        exp0.run(None, None)
    exp0.note("fridge's cold, start sweeping...")
    exp0.note("sweep probe frequency and trap electrode")    
    
im=InstrumentManager()
heman=im['heman']
srs=im['SRS']
#rf1 = im['RFSwitch1']
#awg= im['AWG']
na = im['NWA']
#res = im['res']
#trap = im['trap']
lb = im['labbrick']
fil = im['fil']
rfsrc = im['RFSRC1']
#awg.setup_driver(10e-3,0.0,1.25,'on')
#awg.set_autorange('off')
#awg.pulse_voltage(5, 1.0)

#for i in arange(-4,4,0.1):
#    awg.set_volt(i)
#awg.set_volt(0)
##for i in range(1):
##    heman.set_gas(True,True)
##    sleep(1)
##    heman.set_gas(False,True)
##    heman.set_pump(True,True)
##    sleep(1)
##    heman.set_pump(False,True)
#
#heman.set_gas(False,True)
#heman.set_cryostat(True,True)
#heman.set_pump(True,True)
#
#while True:
#    heman.set_pump(False,True)
#    heman.set_gas(True,True)
#    heman.set_gas(False,True)
#    heman.set_pump(True,True)
#    heman.get_manifold_status()

#heman.set_cryostat(False,True)

def initLB(power=1, freq=7e9):
    lb.set_pulse_ext(False)
    lb.set_mod(False)
    lb.set_power(power)#10dBm max
    lb.set_frequency(freq)
    lb.set_output(True)

def FireFilament(delay=10):
    while True:
        exp0.fil.fire_filament(2000, 0.01)
        sleep(delay)
        Beep(2000,200)

na.params = {
    'power': -20.,
    'center': 6.8734133e9,
    'span' :20e6,
    'sweep_pts' :200,
    'ifbw' :10e3,
    'delay' :0.0,
    'avg' :1}
    
def updateNWA():
    #### NWA parameters
    #na.set_default_state()
    print na.params
    na.set_averages(na.params['avg'])
    na.set_power(na.params['power'])
    na.set_center_frequency(na.params['center'])
    na.set_span(na.params['span'])
    na.set_ifbw(na.params['ifbw'])
    na.set_sweep_points(na.params['sweep_pts'])
    na.set_trigger_source('INTERNAL')
    na.set_trigger_average_mode(True)
    na.set_timeout(10000)  
    na.set_format('lmag')#('slog')    

if __name__=="__main__":
     #### Fridge params
    fridgeParams = {
        'wait_for_temp': 0.040,
        'min_temp_wait_time' : 60 #11 minutes
        }    
    filamentParams = {
        "fil_amp":4.2,
        "fil_off":-0.5,
        "fil_freq":113e3,
        "fil_duration":40e-3,  
        "fil_delay":.01,
        "pulses":150}
    naParams = {
        'BigPeak': {
            'power': -25.,
            'center': 6.88025e9,
            'span' :50e6,
            'sweep_pts' :1601,
            'ifbw' :150e3,
            'avg' :1},
        'SmallPeak': {
            'power':-15.,
            'avg': 1,
            'center': 8008626000.0,
            'ifbw': 50000.0,
            'power': -15.0,
            'span': 2500000.0,
            'sweep_pts': 1601},
        "PowerSettings": {
            "-50": {
                'power': -50.,
                'center': 8.0375e9-0.015e9,#8.008126e9,
                'span' :40e6,
                'sweep_pts' :1601,
                'ifbw' :10,#50e3,
                'avg' :1},
            "-15": {
                'power':-15.,
                'avg': 1,
                'center': 8.0375e9-0.015e9,
                'ifbw': 50000.0,
                'span': 40e6,
                'sweep_pts': 1601},
            "-5": {
                'power': -5.,
                'avg': 1,
                'center': 8.0375e9-0.015e9,
                'ifbw': 10000.0,
                'span': 40e6,
                'sweep_pts': 1601},
            "5": {
                'power':5.,
                'avg': 1,
                'center': 8.0375e9-0.015e9,
                'ifbw': 10000.0,
                'span': 40e6,
                'sweep_pts': 1601},
            }
        }
    labbrickParams={}
    
    exp0=experiment('000_puffs', 'M007v5Trident', fridgeParams,filamentParams,naParams,labbrickParams)
    exp0.res = lambda:0
    exp0.trap = lambda:0
    
    def set_volt_res(volt):
        exp0.srs.set_volt(volt, channel=1)
    def set_volt_trap(volt):
        exp0.srs.set_volt(volt, channel=2)
    def get_volt_res():
        exp0.srs.get_volt(channel=1)
    def get_volt_trap():
        exp0.srs.get_volt(channel=2)
    
    exp0.res.set_volt = set_volt_res
    exp0.trap.set_volt = set_volt_trap
    exp0.res.get_volt = get_volt_res
    exp0.trap.get_volt = get_volt_trap

    exp0.na.update(exp0.na.params['PowerSettings']['-5'])
    exp0.fil.update(exp0.fil.params)
    
    #exp0.fil.fire_filament(filamentParams['pulses']*10,filamentParams['fil_delay'])
    #exp0.fil.fire_filament(2000,0.01)
