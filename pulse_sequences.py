# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 22:13:59 2012

@author: Dr Dave
"""

# This script does all the pulse sequences (T1, Rabi, Ramsey) in a single AWG sequence

from __future__ import division
from slab import *
from slab.instruments import InstrumentManager
from slab.instruments import Alazar, AlazarConfig
from slab.awgpulses import *
from scipy.signal import decimate
from numpy import *
import cProfile
import pstats
from slab.awgpulses import *
import ctypes
from slab.datamanagement import load_array
from slab.instruments.awg.TekPattern import *
from slab.instruments.awg.TekPattern2 import *
from collections import namedtuple
import sys


def expfunc2(p, x, x0):
    """p[0]+p[1]*exp(-(x-p[2])/p[3])"""
    return p[0] + p[1] * math.e ** (-(x - x0) / p[2])


#class that makes pulses for the AWG 
class pulse():
    def __init__(self, pulse_length=8192, delay=0, clk_factor=1.0):

        self.pulse_length = pulse_length
        self.pulse_range = range(pulse_length)
        self.pulse = zeros(pulse_length)
        self.delay = delay  #The delay of the AWG *in units of the AWG clock*
        self.clk_factor = clk_factor  #The ratio between the specified time units and the AWG clock

    def clear(self):

        self.pulse[0:-1] = 0.0

    def square_pulse(self, start_pulse, pulse_width, pulse_height=1, scale_times=True):

        if scale_times:
            start_pulse = start_pulse * self.clk_factor - self.delay
            pulse_width *= self.clk_factor

        if start_pulse < 0 or (start_pulse + pulse_width) > self.pulse_length:
            raise NameError(
                "Invalid pulse parameters." + str(start_pulse) + " " + str(pulse_width) + " " + str(self.pulse_length))

        self.pulse[floor(start_pulse):ceil(start_pulse + pulse_width)] = pulse_height
        return ceil(start_pulse + pulse_width)

    def square_pulse2(self, pulse_center, pulse_width, pulse_height=1, scale_times=True):

        return self.square_pulse(pulse_center - pulse_width / 2, pulse_width, pulse_height, scale_times)

    def gauss_pulse(self, gauss_center, gauss_sigma, gauss_height, start_of_pulse=False):

        gauss_center = gauss_center * self.clk_factor - self.delay
        gauss_sigma *= self.clk_factor

        #print gauss_center

        #pulse is truncated to +- 2sigma of center     
        if start_of_pulse:
            gauss_center2 = gauss_center + 2 * gauss_sigma
        else:
            gauss_center2 = gauss_center

        self.pulse += gauss_height * (
        exp(-1.0 * (array(self.pulse_range) - gauss_center2) ** 2.0 / (2.0 * gauss_sigma ** 2.0)) - exp(-2.0)) / (
                      1 - exp(-2.0)) * (array(self.pulse_range) >= (gauss_center2 - 2 * gauss_sigma)) * (
                      array(self.pulse_range) <= (gauss_center2 + 2 * gauss_sigma))

        return 2 * gauss_sigma + gauss_center

    def smooth_square(self, pulse_center, smooth_time, flat_time, pulse_height):

        pulse_center = pulse_center * self.clk_factor - self.delay
        smooth_time *= self.clk_factor
        flat_time *= self.clk_factor

        self.pulse += pulse_height * (exp(-1.0 * (array(self.pulse_range) - (pulse_center - flat_time / 2)) ** 2.0 / (
        2.0 * smooth_time ** 2.0)) - exp(-2.0)) / (1.0 - exp(-2.0)) * (
                      array(self.pulse_range) <= (pulse_center - flat_time / 2)) * (
                      array(self.pulse_range) >= (pulse_center - flat_time / 2 - 2 * smooth_time))
        self.pulse += pulse_height * (exp(-1.0 * (array(self.pulse_range) - (pulse_center + flat_time / 2)) ** 2.0 / (
        2.0 * smooth_time ** 2.0)) - exp(-2.0)) / (1.0 - exp(-2.0)) * (
                      array(self.pulse_range) >= (pulse_center + flat_time / 2)) * (
                      array(self.pulse_range) <= (pulse_center + flat_time / 2 + 2 * smooth_time))
        self.square_pulse2(pulse_center, flat_time, pulse_height, False)
        return flat_time / 2 + pulse_center + smooth_time

    def smooth_square_w_chirp(self, pulse_center, smooth_time, flat_time, pulse_height, f1, f2):

        ret_val = self.smooth_square(pulse_center, smooth_time, flat_time, pulse_height)

        f1 = f1 / self.clk_factor
        f2 = f2 / self.clk_factor
        pulse_center = pulse_center * self.clk_factor - self.delay
        flat_time *= self.clk_factor
        smooth_time *= self.clk_factor

        t = array(self.pulse_range) - (pulse_center - flat_time / 2.0 - 2 * smooth_time)
        farray = f1 + (f2 - f1) / (flat_time + 4 * smooth_time) * t
        #        print farray[floor(pulse_center)]
        #        print farray[floor(pulse_center-flat_time/2.0)]
        #        print farray[floor(pulse_center+flat_time/2.0)]
        self.pulse *= cos(2 * pi * farray * t)
        return ret_val

    def trapezoid(self, pulse_center, slope1, slope2, total_time, pulse_height, pulse_height2=None):

        #slope1 is the slope up and slope2 is the slope down        

        #a trapezoid up to pulse_height2 and a square pulse up to pulse_height       
        if pulse_height2 is None:
            pulse_height2 = pulse_height

        pulse_center = pulse_center * self.clk_factor - self.delay
        slope1 /= self.clk_factor
        slope2 /= self.clk_factor
        total_time *= self.clk_factor

        #print pulse_center, slope, total_time, pulse_height

        if sign(pulse_height) != sign(slope1) or sign(pulse_height) != sign(slope2):
            raise NameError("Slope and Pulse Height Need the Same Sign for Trapezoid")

        flat_time = total_time - abs(pulse_height2 / slope1) - abs(pulse_height2 / slope2)

        flat_center = pulse_center + (abs(pulse_height2 / slope1) - abs(pulse_height2 / slope2)) / 2

        if flat_time < 0:
            raise NameError("Slope/Total_time Error in Trapesoid Pulse")

        self.pulse += slope1 * (array(self.pulse_range) - (pulse_center - total_time / 2.0)) * (
        array(self.pulse_range) <= (flat_center - flat_time / 2)) * (
                      array(self.pulse_range) >= (pulse_center - total_time / 2))
        self.pulse += -1 * slope2 * (array(self.pulse_range) - (pulse_center + total_time / 2.0)) * (
        array(self.pulse_range) >= (flat_center + flat_time / 2)) * (
                      array(self.pulse_range) <= (pulse_center + total_time / 2))
        self.square_pulse2(flat_center, flat_time, pulse_height, False)
        return total_time / 2 + pulse_center

    #this is a trapezoid that has a pedestal at the bottom
    def trapezoid_w_pedestal(self, pulse_center, slope1, slope2, total_time, total_trap_time, pedestal_height,
                             pulse_height, pulse_height2=None):

        #slope1 is the slope up and slope2 is the slope down        

        #a trapezoid up to pulse_height2 and a square pulse up to pulse_height       
        if pulse_height2 is None:
            pulse_height2 = pulse_height

        pulse_center = pulse_center * self.clk_factor - self.delay
        #print pulse_center        
        slope1 /= self.clk_factor
        slope2 /= self.clk_factor
        total_time *= self.clk_factor
        total_trap_time *= self.clk_factor

        delta_height = pulse_height2 - pedestal_height

        #print pulse_center, slope, total_time, pulse_height

        if sign(delta_height) != sign(slope1) or sign(pulse_height) != sign(slope2):
            raise NameError("Slope and Pulse Height Need the Same Sign for Trapezoid")

        if total_trap_time > total_time:
            raise NameError("Total time needs to be larger than the trapezoid time")

        flat_time = total_trap_time - abs(delta_height / slope1) - abs(delta_height / slope2)

        flat_center = pulse_center + (abs(delta_height / slope1) - abs(delta_height / slope2)) / 2

        if flat_time < 0:
            raise NameError(
                "Slope/Total_time Error in Trapesoid Pulse." + str(total_trap_time) + " " + str(delta_height))

        #add pedestal
        self.square_pulse2(pulse_center, total_time, pedestal_height, False)

        self.pulse += slope1 * (array(self.pulse_range) - (pulse_center - total_trap_time / 2.0)) * (
        array(self.pulse_range) <= (flat_center - flat_time / 2)) * (
                      array(self.pulse_range) >= (pulse_center - total_trap_time / 2))
        self.pulse += -1 * slope2 * (array(self.pulse_range) - (pulse_center + total_trap_time / 2.0)) * (
        array(self.pulse_range) >= (flat_center + flat_time / 2)) * (
                      array(self.pulse_range) <= (pulse_center + total_trap_time / 2))
        self.square_pulse2(flat_center, flat_time, pulse_height, False)
        return total_time / 2 + pulse_center

    def gauss_pulse_with_freq(self, gauss_center, gauss_sigma, gauss_height, pulse_freq, phase=0.0,
                              start_of_pulse=False):

        gauss_center = gauss_center * self.clk_factor - self.delay
        gauss_sigma *= self.clk_factor
        pulse_freq /= self.clk_factor

        #pulse is truncated to +- 2sigma of center     
        if start_of_pulse:
            gauss_center2 = gauss_center + 2 * gauss_sigma
        else:
            gauss_center2 = gauss_center

        self.pulse += cos(2.0 * pi * pulse_freq * (array(self.pulse_range) - gauss_center2) + phase) * gauss_height * (
        exp(-1.0 * (array(self.pulse_range) - gauss_center2) ** 2.0 / (2.0 * gauss_sigma ** 2.0)) - exp(-2.0)) / (
                      1 - exp(-2.0)) * (array(self.pulse_range) >= (gauss_center2 - 2 * gauss_sigma)) * (
                      array(self.pulse_range) <= (gauss_center2 + 2 * gauss_sigma))

        return 2 * gauss_sigma + gauss_center

    def gauss_pulse_with_chirp(self, gauss_center, gauss_sigma, gauss_height, pulse_freq_start, pulse_freq_end,
                               phase=0.0, start_of_pulse=False):

        gauss_center = gauss_center * self.clk_factor - self.delay
        gauss_sigma *= self.clk_factor
        pulse_freq_start /= self.clk_factor
        pulse_freq_end /= self.clk_factor

        #pulse is truncated to +- 2sigma of center     
        if start_of_pulse:
            gauss_center2 = gauss_center + 2 * gauss_sigma
        else:
            gauss_center2 = gauss_center

        xdata = array(self.pulse_range) - gauss_center2

        self.pulse += cos(2.0 * pi * (pulse_freq_start + (xdata - gauss_sigma * 2) / (4 * gauss_sigma) * (
        pulse_freq_end - pulse_freq_start)) * xdata + phase) * gauss_height * (
                      exp(-1.0 * (array(self.pulse_range) - gauss_center2) ** 2.0 / (2.0 * gauss_sigma ** 2.0)) - exp(
                          -2.0)) / (1 - exp(-2.0)) * (array(self.pulse_range) >= (gauss_center2 - 2 * gauss_sigma)) * (
                      array(self.pulse_range) <= (gauss_center2 + 2 * gauss_sigma))

        return 2 * gauss_sigma + gauss_center

    def gauss_pulse_with_multiple_freq(self, gauss_center, gauss_sigma, gauss_height, pulse_freq, phase=0.0,
                                       start_of_pulse=False):

        gauss_center = gauss_center * self.clk_factor - self.delay
        gauss_sigma *= self.clk_factor
        pulse_freq = divide(pulse_freq, self.clk_factor)

        #pulse is truncated to +- 2sigma of center     
        if start_of_pulse:
            gauss_center2 = gauss_center + 2 * gauss_sigma
        else:
            gauss_center2 = gauss_center

        for f in pulse_freq:
            self.pulse += cos(2.0 * pi * f * (array(self.pulse_range) - gauss_center2) + phase)

        self.pulse *= gauss_height * (
        exp(-1.0 * (array(self.pulse_range) - gauss_center2) ** 2.0 / (2.0 * gauss_sigma ** 2.0)) - exp(-2.0)) / (
                      1 - exp(-2.0)) * (array(self.pulse_range) >= (gauss_center2 - 2 * gauss_sigma)) * (
                      array(self.pulse_range) <= (gauss_center2 + 2 * gauss_sigma))

        return 2 * gauss_sigma + gauss_center

    def deriv_gauss_pulse(self, gauss_center, gauss_sigma, gauss_height, start_of_pulse=False):

        gauss_center = gauss_center * self.clk_factor - self.delay
        gauss_sigma *= self.clk_factor

        #print gauss_center

        #pulse is truncated to +- 2sigma of center     
        if start_of_pulse:
            gauss_center2 = gauss_center + 2 * gauss_sigma
        else:
            gauss_center2 = gauss_center

        self.pulse += -1.0 * (array(self.pulse_range) - gauss_center2) / (gauss_sigma ** 2.0) * gauss_height * (
        exp(-1.0 * (array(self.pulse_range) - gauss_center2) ** 2.0 / (2.0 * gauss_sigma ** 2.0)) - exp(-2.0)) / (
                      1 - exp(-2.0)) * (array(self.pulse_range) >= (gauss_center2 - 2 * gauss_sigma)) * (
                      array(self.pulse_range) <= (gauss_center2 + 2 * gauss_sigma))

        return 2 * gauss_sigma + gauss_center

    def offset(self, offset_val, start_time=-1, end_time=-1):

        if start_time == -1:
            self.pulse += offset_val
        else:

            start_time = start_time * self.clk_factor - self.delay
            end_time = end_time * self.clk_factor - self.delay

            self.pulse[floor(start_time):floor(end_time)] = self.pulse[floor(start_time):floor(end_time)] + offset_val

    #sets first and last point to offset_val    
    def offset2(self, offset_val):

        self.pulse[0] = offset_val
        self.pulse[-1] = offset_val

    def line(self, start_time, line_duration, start_val, stop_val):

        start_time = int(floor(start_time * self.clk_factor - self.delay))
        line_duration *= self.clk_factor
        stop_time = int(floor(start_time + line_duration))

        self.pulse[start_time:stop_time] = (stop_val - start_val) / (stop_time - start_time) * (
        array(range(start_time, stop_time)) - start_time) + start_val


    def tan_pulse_w_pedestal(self, pulse_center, duration1, duration2, total_time, mid_val, end_val, pedestal_height,
                             pedestal_time, slope1_factor=0.9, slope2_factor=0.9):
        #slope_factor is a number between 0->1 which sets how steep
        #slope is through the middle region

        wait_time = total_time - 2 * (duration1 + duration2)

        if wait_time < 0:
            raise NameError("Total pulse time must be longer")

        pulse_center = int(floor(pulse_center * self.clk_factor - self.delay))
        duration1 = int(floor(duration1 * self.clk_factor))
        duration2 = int(floor(duration2 * self.clk_factor))
        wait_time = int(floor(wait_time * self.clk_factor))
        pedestal_time = int(floor(pedestal_time * self.clk_factor))

        #up slope
        times = zeros(8)

        times[0] = pulse_center - (wait_time / 2 + duration1 + duration2 + pedestal_time)
        times[1] = times[0] + pedestal_time
        times[2] = times[1] + duration1
        times[3] = times[2] + duration2
        times[4] = times[3] + wait_time
        times[5] = times[4] + duration2
        times[6] = times[5] + duration1
        times[7] = times[6] + pedestal_time

        #print mid_val, duration1, duration2, total_time

        if slope1_factor < 0 or slope1_factor >= 1 or slope2_factor < 0 or slope2_factor >= 1:
            raise NameError("Slope factors must be greater than 0 and less than 1")

        #add pedestals
        if pedestal_time > 0:
            self.square_pulse(times[0], pedestal_time, pedestal_height, False)
            self.square_pulse(times[6], pedestal_time, pedestal_height, False)

        start_val = pedestal_height

        #main tan pulse
        self.pulse[times[1]:times[2]] = mid_val + (start_val - mid_val) / tan(pi / 2 * slope1_factor) * tan(
            pi / 2 * slope1_factor * (duration1 * 1.0 - array(range(duration1))) / duration1)
        self.pulse[times[2]:times[3]] = mid_val + (end_val - mid_val) / tan(pi / 2 * slope2_factor) * tan(
            pi / 2 * slope2_factor * array(range(duration2)) / duration2)
        self.pulse[times[3]:times[4]] = end_val
        self.pulse[times[4]:times[5]] = mid_val + (end_val - mid_val) / tan(pi / 2 * slope2_factor) * tan(
            pi / 2 * slope2_factor * (duration2 * 1.0 - array(range(duration2))) / duration2)
        self.pulse[times[5]:times[6]] = mid_val + (start_val - mid_val) / tan(pi / 2 * slope1_factor) * tan(
            pi / 2 * slope1_factor * array(range(duration1)) / duration1)


    #class for loading pulse sequences (eg. T1, Rabi, Ramsey)


#DM: Aug 31, rewrite for TEK
#For now all TEK pulses have to be the same length
class pulse_sequence():
    def __init__(self, seq_length=1, pulse_length=8192, agilent_length=8192, TEK_clk=1.2, Agilent_clk=4.2,
                 delays=[0, 0, 0, 0, 0, 0]):

        #agilent can have a different pulse length than the TEK because the clock is different
        self.pulse_length = int(pulse_length)

        self.TEK_clk = TEK_clk
        self.Agilent_clk = Agilent_clk
        self.delays = delays

        #note that agilent length needs to be divisible by 32
        agilent_length2 = int(floor(agilent_length / 32) * 32)

        if agilent_length2 != agilent_length:
            print "Original agilent sequence length of " + str(agilent_length) + " changed to " + str(
                agilent_length2) + " to be divisible by 32."
            agilent_length = agilent_length2

        self.agilent_length = agilent_length

        #add a least one pulse to each channel

        #there are 4 TEK and 2 Agilent analog channels
        self.analogwf = list()
        for i in range(6):
            self.analogwf.append(list())
            self.add_analog_pulse(i)


        #there are 8 marker channels
        self.marker = list()
        for i in range(12):
            self.marker.append(list())
            self.add_marker_pulse(i)


        #this is the sequencing length information
        #by default the sequence will point to first pulse 
        self.num_seqs = seq_length
        self.marker_seq_table = zeros((12, seq_length))
        self.analog_seq_table = zeros((6, seq_length))

        #this index can be used by other classes to indicate which channels
        #are connected to which devices
        self.channel_index = dict()
        self.channel_index['card_trig'] = 0
        self.channel_index['read_trig'] = [1, 3]
        self.channel_index['drive_trig'] = [2, 4]
        self.channel_index['flux_trig'] = 5
        self.channel_index['drive_I'] = [0, 2]
        self.channel_index['drive_Q'] = [1, 3]
        self.channel_index['flux'] = [4, 5]

        #Uncomment to use the TEK for flux pulses
        #self.channel_index['drive_Q'] = [4,5]
        #self.channel_index['flux'] = [1,3]

    #ALWAYS add pulses through this method to ensure proper timing   
    def add_marker_pulse(self, ch_num):

        if ch_num < 8:
            self.marker[ch_num].append(pulse(self.pulse_length, self.delays[int(floor(ch_num / 2))], self.TEK_clk))
        else:
            self.marker[ch_num].append(
                pulse(self.agilent_length, self.delays[int(floor(ch_num / 2))], self.Agilent_clk))
        return len(self.marker[ch_num]) - 1

    #ALWAYS add pulses through this method to ensure proper timing    
    def add_analog_pulse(self, ch_num):

        if ch_num <= 3:
            self.analogwf[ch_num].append(pulse(self.pulse_length, self.delays[ch_num], self.TEK_clk))
        else:
            self.analogwf[ch_num].append(pulse(self.agilent_length, self.delays[ch_num], self.Agilent_clk))
        return len(self.analogwf[ch_num]) - 1

    #this adds the I and Q part of the gauss pulse
    def gauss_pulse(self, gauss_center, gauss_sigma, gauss_height, gauss_phase, qubit_num, indI, indQ,
                    start_of_pulse=False):

        self.analogwf[self.channel_index['drive_I'][qubit_num]][indI].gauss_pulse(gauss_center, gauss_sigma,
                                                                                  cos(gauss_phase) * gauss_height,
                                                                                  start_of_pulse)
        self.analogwf[self.channel_index['drive_Q'][qubit_num]][indQ].gauss_pulse(gauss_center, gauss_sigma,
                                                                                  sin(gauss_phase) * gauss_height,
                                                                                  start_of_pulse)


    #this adds the I and Q part of the DRAG pulse
    def DRAG_pulse(self, gauss_center, gauss_sigma, gauss_height, gauss_phase, qubit_num, indI, indQ, drag_prefactor,
                   start_of_pulse=False):

        #main Gaussian
        self.analogwf[self.channel_index['drive_I'][qubit_num]][indI].gauss_pulse(gauss_center, gauss_sigma,
                                                                                  cos(gauss_phase) * gauss_height,
                                                                                  start_of_pulse)
        self.analogwf[self.channel_index['drive_Q'][qubit_num]][indQ].gauss_pulse(gauss_center, gauss_sigma,
                                                                                  sin(gauss_phase) * gauss_height,
                                                                                  start_of_pulse)

        #derivative in the other axis
        self.analogwf[self.channel_index['drive_I'][qubit_num]][indI].deriv_gauss_pulse(gauss_center, gauss_sigma, -sin(
            gauss_phase) * gauss_height * drag_prefactor, start_of_pulse)
        self.analogwf[self.channel_index['drive_Q'][qubit_num]][indQ].deriv_gauss_pulse(gauss_center, gauss_sigma, cos(
            gauss_phase) * gauss_height * drag_prefactor, start_of_pulse)


    #this loads into the TEK
    def load_into_awg(self, filename, awg_name):

        #filename is where to build the sequence file

        #create AWG data file
        awgdata = dict()
        key_names = ['ch12', 'ch34']
        for i in range(1, 5):
            for j in range(1, 3):
                key_names.append('ch{0}m{1}'.format(i, j))

        for i in range(len(key_names)):
            awgdata[key_names[i]] = dict()
            awgdata[key_names[i]]['wfLib'] = dict()
            awgdata[key_names[i]]['linkList'] = list()



        #add wf's to awgdata

        #analog
        #combine channels 1 and 2, 3 and 4 together
        for i in range(2):
            data_key = 'ch{0}{1}'.format(2 * i + 1, 2 * i + 2)
            for j in range(max((len(self.analogwf[2 * i]), len(self.analogwf[2 * i + 1])))):
                awgdata[data_key]['wfLib'][str(j)] = self.analogwf[2 * i][j % len(self.analogwf[2 * i])].pulse + (
                                                                                                                 self.analogwf[
                                                                                                                     2 * i + 1][
                                                                                                                     j % len(
                                                                                                                         self.analogwf[
                                                                                                                             2 * i + 1])].pulse) * 1j

            #create sequence
            for j in range(self.num_seqs):
                awgdata[data_key]['linkList'].append(list())
                awgdata[data_key]['linkList'][j].append(namedtuple('a', 'key isTimeAmp length repeat'))
                awgdata[data_key]['linkList'][j][0].key = "{0:g}".format(self.analog_seq_table[2 * i][j])
                awgdata[data_key]['linkList'][j][0].isTimeAmp = False

        for i in range(4):
            for j in range(2):
                marker_key = 'ch{0}m{1}'.format(i + 1, j + 1)
                for k in range(len(self.marker[2 * i + j])):
                    awgdata[marker_key]['wfLib'][str(k)] = self.marker[2 * i + j][k].pulse

                    #create sequence
                for k in range(self.num_seqs):
                    awgdata[marker_key]['linkList'].append(list())
                    awgdata[marker_key]['linkList'][k].append(namedtuple('a', 'key isTimeAmp length repeat'))
                    awgdata[marker_key]['linkList'][k][0].key = "{0:g}".format(self.marker_seq_table[2 * i + j][k])
                    awgdata[marker_key]['linkList'][k][0].isTimeAmp = False

        write_Tek_file(awgdata, filename, 'seq1', None, False)

        #load the file into the TEK
        im = InstrumentManager()
        awg = im[awg_name]
        awg.pre_load()
        awg.load_sequence_file(filename)


    #load the full waveform sequence into the Agilent 
    def load_full_into_agilent(self, awg_str, do_seq=True):


        im = InstrumentManager()
        awg = im[awg_str]

        awg.reset()
        awg.presetup_for_sequences()

        for i in range(2):

            awg.select_channel(i + 1)
            #ch = self.channel_index['flux'][i]
            ch = i + 4

            for j in range(len(self.analogwf[ch])):

                if len(self.marker[i * 2 + 1 + 7]) <= j:
                    marker1data = zeros(len(self.analogwf[ch][j].pulse) / 4)
                else:
                    marker1data = self.marker[i * 2 + 1 + 7][j].pulse[0:len(self.marker[i * 2 + 1 + 7][j].pulse):4]

                #print  len(self.marker[i*2+2+7]),len(self.marker[i*2+2+7][j].pulse)                

                if len(self.marker[i * 2 + 2 + 7]) <= j:
                    marker2data = zeros(len(self.analogwf[ch][j].pulse) / 4)
                else:
                    marker2data = self.marker[i * 2 + 2 + 7][j].pulse[0:len(self.marker[i * 2 + 2 + 7][j].pulse):4]

                #print len(self.analogwf[ch][j].pulse), len(marker1data), len(marker2data)

                awg.add_floatsegment(self.analogwf[ch][j].pulse, segnum=j + 1, marker1=marker1data, marker2=marker2data,
                                     float_range=(-1., 1.))
                print "Loading Waveform" + str(j) + " for Channel " + str(i)


        #define the sequence
        if do_seq:
            self.load_sequence_into_agilent(awg_str)


    #just rewrite the sequence information for the Agilent
    def load_sequence_into_agilent(self, awg_str, preset=False):

        im = InstrumentManager()
        awg = im[awg_str]

        if preset:
            awg.reset(True)  #deletes the sequences only!
            awg.preset_for_sequences()

        for i in range(2):

            awg.select_channel(i + 1)
            #ch = self.channel_index['flux'][i]
            ch = i + 4

            for j in range(len(self.analog_seq_table[ch])):
                awg.define_sequence_step(step=j + 1, seg_num=self.analog_seq_table[ch][j] + 1, loops=1,
                                         jump_flag=False)  #will advance to the next segment on trigger
                print "Defining Sequence Step " + str(j) + " for Channel " + str(i)

    #load into TEK700001
    #Will *only* load analog channel #4 (first agilent channel)
    def load_into_tek2(self, folder, awg_str):

        im = InstrumentManager()
        awg = im[awg_str]

        awg.pre_load()

        #make waveform files
        for j in range(len(self.analogwf[4])):
            filename = os.path.join(folder, 'A' + str(j) + '.wfmx')

            #code in 'TekPattern2'
            create_waveform_file(filename, self.analogwf[4][j].pulse)

            print "Loading Waveform File" + filename + " into TEK70001"

            #add file
            awg.load_waveform_file(filename)
            awg.operation_complete()

        #add new sequence
        awg.new_sequence(num_steps=len(self.analog_seq_table[4]))

        for j in range(len(self.analog_seq_table[4])):
            #assign waveform
            print "A{:g}".format(self.analog_seq_table[4][j])
            awg.assign_seq_waveform(step=j + 1, waveform="A{:g}".format(self.analog_seq_table[4][j]),
                                    last_step=((j + 1) == len(self.analog_seq_table[4])))


#general class for loading, running and saving a qubit experiment
class qubit_exp():
    def __init__(self):

        #xdata is the array which is being varied in the experiment
        self.xdata = []

        #result of the experiment (2 x len(xdata) array)
        self.ydata_avg = []

        #full averaged data (2 x len(xdata) x len(acq_length) array)
        self.full_trace_data_avg = []

        #These next two are only optionally saved because they can be quite large        
        #result of each run (2 x len(xdata) x len(num_avgs) array)
        self.ydata = []

        #full data from each run (2 x len(xdata) x len(acq_length) x len(num_avgs) array)
        self.full_trace_data = []

        #define defaults for the config settings
        #for settings that depend on the qubit, these will be a list of two options
        self.config = dict(
            monitor_pulses=False,
            read_delay=50,
            card_delay=256,  # 256*2+100
            read_length=2048,
            acq_length=1024,
            meas_range=[0.1, 0.1],
            total_length=8192,
            front_buffer=500,
            generator_enable_delay=200,  #150

            #Qubit drive properties
            drive_freq=[7.1405e9, 7.1405e9],
            drive_sideband=[False, False],
            drive_sideband_freq=[50.0e6, 50.0e6],
            drive_pwr=[13.0, 13.0],  #note that this doesn't do anything when driving with the BNC
            drive=['RF3', 'RF4'],
            shaped_pulse=True,
            pulsed_drive=[True, True],
            drag_pulse=[False, False],
            drag_prefactor=[-2.0 / 4.0 / (2 * pi * 0.2), -2.0 / 4.0 / (2 * pi * 0.2)],

            #Qubit read properties 
            read_freq=[5.6e9, 5e9],
            lo_power=[16.0, 16.0],
            read_pwr=[-5.0, -5.0],
            RF=['RF2', 'RF1'],
            LO=['RF2', 'RF1'],
            PHASE=['LBPHASE1', 'LBPHASE2'],
            read_phase=[70.0, 70.0],
            IFreq=[10.0e6, 10.0e6],  #for heterodyne

            #flux properties
            flux_instr=[['YOKO2', 1], ['YOKO1', 3]],
            flux=[0.0, 0.0],
            flux_type='voltage',

            qubit_enable=[True, True],  #1 or 2 qubit experiment

            #awg properties
            #first awg is the main "trigger" AWG
            seq_file="S:\\_Data\\sequence_files\temp.awg",
            AWG=['TEK', 'AGILENT811'],
            awg_amp=[[2.0, 2.0, 2.0, 2.0], [1.0, 1.0]],
            awg_offset=[[0.0, 0.0, 0.0, 0.0], [0.0, 0.0]],
            awg_marker_amps=[[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0]],
            awg_clocks=[1.2, 4.2],  #GSa/s (Note: these need to be set manually)
            awg_delay=[0, -6.0 * 0, -0.0, -1.0, 320 + 10 + 10, 345 - 5],
            #These are the delays between the various channels (note: there is an arbitrary overall delay and the delays are defined w.r.t. the LOCAL clock)
            load_agilent=True,  #if this is false then the agilent will not load under any condition!
            start_agilent=True,  #if this is false then the agilent outputs will be set to off
            awg_pulse_offset=[0.0, 0.0, 0.0, 0.0],  #this is the offset when the drive pulse is on

            exp_id='',  #filled out by the specified experiment (rabi, ramsey, etc.)

            #acquisition properties
            num_avgs=-1,  #-1 is run continuously,
            data_window=[[0, 1024], [0, 1024]],  #part of the acqusition to use for determining the qubit state,
            save_each_run=False,  #this can only be selected if num_avgs is a finite number
            save_each_run_fulldata=False,
            #this can only be selected if num_avgs is a finite number (and should be small!)

            num_calib_seqs=0
        )

    def save_data(self, exp_path, file_prefix, overwrite=False):

        #print expt_path+fname
        try:
            if overwrite:
                f = SlabFile(exp_path + file_prefix + '.h5', 'w')
            else:
                fname = get_next_filename(exp_path, prefix=file_prefix, suffix='.h5')
                print fname
                f = SlabFile(exp_path + fname, 'w')
        except NameError:
            print "Error opening file for saving rabi data"
            return

        #save
        f.create_dataset('xdata', data=self.xdata)
        f.create_dataset('ydata', data=self.ydata)
        f.create_dataset('ydata_avg', data=self.ydata_avg)
        f.create_dataset('full_trace_data_avg', data=self.full_trace_data_avg)
        f.create_dataset('full_trace_data', data=self.full_trace_data)

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

    def load_data(self, exp_path, filename, settings_only=True, exp_id=''):

        try:
            f = SlabFile(exp_path + filename, 'r')
        except:
            print "Error opening file, will use default settings"
            return

        #check that this is a vacuum rabi file
        if f.attrs['exp_id'] != 'exp_id':
            print "Wrong data file! Not " + exp_id
            f.close()
            return

        #load
        if not settings_only:
            self.xdata = load_array(f, 'xdata')
            self.ydata = load_array(f, 'ydata')
            self.ydata_avg = load_array(f, 'ydata_avg')
            self.full_trace_data_avg = load_array(f, 'full_trace_data_avg')
            self.full_trace_data = load_array(f, 'full_trace_data')


        #load keys (allow for new keys to be added later)
        config2 = f.load_settings()

        for i, key in enumerate(config2.keys()):
            self.config[key] = config2[key]

        f.close()

    #run the experiment
    def run_qubit_exp(self, plotter=None):

        #define instrument manager    

        im = InstrumentManager()



        #stop the awg
        for i in range(len(self.config['AWG'])):
            awg = im[self.config['AWG'][i]]
            awg.stop()

        if self.config['LO'][0] == self.config['RF'][0]:
            homodyne_meas = True
            print "Homodyne measurement"
        else:
            homodyne_meas = False

        for i in range(2):


            #setup drive
            drive = im[self.config['drive'][i]]

            drive_freq = self.config['drive_freq'][i] - self.config['drive_sideband_freq'][i] * \
                                                        self.config['drive_sideband'][i]

            if not self.config['qubit_enable'][i]:
                drive.set_output(False)

            else:

                #setup drive
                if self.config['drive'][i][0:2] == "LB":
                    drive.set_output(True)
                    drive.set_pulse_ext(mod=self.config['pulsed_drive'][i])

                    #Internal pulse must be set to false!
                    drive.set_mod(False)
                    drive.set_power(self.config['drive_pwr'][i])  #-28 when direct driving
                    print drive_freq
                    drive.set_frequency(drive_freq)
                    time.sleep(0.1)

                    for j in range(5):
                        if abs(drive.get_frequency() - drive_freq) > 100:
                            drive.set_frequency(drive_freq)
                        else:
                            break
                    if abs(drive.get_frequency() - drive_freq) > 100:
                        print "Lab brick frequency error!"
                    time.sleep(.2)

                else:
                    drive.set_power(self.config['drive_pwr'][i])
                    drive.set_frequency(drive_freq)
                    drive.set_ext_pulse(self.config['pulsed_drive'][i])
                    drive.set_output(True)

                print self.config['drive_freq'][i], drive.get_frequency()



            #setup LO and RF


            RF = im[self.config['RF'][i]]
            if not self.config['qubit_enable'][i]:
                #RF.set_output(False)
                pass
            else:
                RF.set_output(True)
                RF.set_mod(True)
                RF.set_power(self.config['read_pwr'][i])
                RF.set_frequency(self.config['read_freq'][i])
                RF.set_ext_pulse()

            if not homodyne_meas:

                LO = im[self.config['LO'][i]]
                if not self.config['qubit_enable'][i]:
                    LO.set_output(False)
                else:
                    LO.set_power(self.config['lo_power'][i])
                    LO.set_mod(True)
                    LO.set_frequency(self.config['read_freq'][i] + self.config['IFreq'][i])
                    LO.set_output(True)

            #set phase
            if self.config['PHASE'][i] != '':
                phase_instr = im[self.config['PHASE'][i]]
                phase_instr.set_working_frequency(self.config['read_freq'][i])
                phase_instr.set_phase(self.config['read_phase'][i])

                #setup the flux
            #note: even if the qubit is not participating in the experiment, setup the flux!
            flux = im[self.config['flux_instr'][i][0]]

            if self.config['flux_instr'][i][0] == 'flux2' or self.config['flux_instr'][i][0] == 'flux1':
                flux.set_function("DC")
                flux.set_offset(self.config['flux'][i])
            else:
                if self.config['flux_type'] == 'current':
                    flux.set_mode('current')
                    flux.set_current(self.config['flux'][i])
                    print 'Flux %s: %s A' % (i, flux.get_current())
                else:
                    flux.set_mode('voltage')
                    flux.set_volt(self.config['flux'][i])

        #setup the Alazar card
        card = self.__setup_Alazar_card(-1)

        #set to sequence mode and set and the channel amplitudes and offsets
        for i in range(len(self.config['AWG'])):
            awg = im[self.config['AWG'][i]]
            awg.prep_experiment()
            awg.set_amps_offsets(self.config['awg_amp'][i], self.config['awg_offset'][i],
                                 self.config['awg_marker_amps'][i])



        #do the acquision
        self.__card_get_data(card, plotter, homodyne_meas)

    def run_qubit_exp_sequentially():
        #To do
        pass

    def load_exp():
        #This is specific to the child class
        pass

    #wrapper class for setting up the card
    def __setup_Alazar_card(self, records=-1):


        #records = 0x7FFFFFFF is for continuous streaming  
        if records == -1:
            records = 0x7FFFFFFF

        if self.config['monitor_pulses']:
            records = total_seqs
            ch1_range = 4
            ch2_range = 4
            ch1_coupling = 'DC'
            ch2_coupling = 'DC'
            ch1_enable = True
            ch2_enable = True
        else:
            ch1_range = self.config['meas_range'][0]  #0.1
            ch2_range = self.config['meas_range'][1]
            ch1_coupling = 'DC'
            ch2_coupling = 'DC'
            ch1_enable = self.config['qubit_enable'][0]
            ch2_enable = self.config['qubit_enable'][1]

        print "Configuring card"
        config = {'clock_edge': 'rising', 'clock_source': 'reference',
                  'trigger_coupling': 'DC', 'trigger_operation': 'or',
                  'trigger_source2': 'disabled', 'trigger_level2': 1.0, 'trigger_edge2': 'rising',
                  'trigger_level1': 0.6, 'trigger_edge1': 'rising', 'trigger_source1': 'external',
                  'ch1_enabled': ch1_enable, 'ch1_coupling': ch1_coupling, 'ch1_range': ch1_range, 'ch1_filter': False,
                  'ch2_enabled': ch2_enable, 'ch2_coupling': ch2_coupling, 'ch2_range': ch2_range, 'ch2_filter': False,
                  'bufferCount': 20, 'recordsPerBuffer': len(self.xdata) + self.config['num_calib_seqs'],
                  'trigger_delay': 0, 'timeout': 1000,
                  'samplesPerRecord': self.config['acq_length'], 'recordsPerAcquisition': records,
                  'sample_rate': 1000000
        }


        # 0x7FFFFFFF
        scope_settings = AlazarConfig(config)

        card = Alazar(scope_settings)
        card.configure(scope_settings)

        if records == 0x7FFFFFFF or self.config['monitor_pulses']:
            card.Az.AlazarAbortAsyncRead(card.get_handle())
            card.post_buffers()

        return card


    def __card_get_data(self, card, plotter, homodyne_meas=False):

        buffersCompleted = 0
        buffersPerAcquisition = card.config.recordsPerAcquisition / card.config.recordsPerBuffer
        recordsPerBuffer = card.config.recordsPerBuffer
        samplesPerRecord = card.config.samplesPerRecord
        #ch1_sum = np.zeros((total_seqs, samplesPerRecord),dtype=float)
        ch_sum = [zeros(recordsPerBuffer * samplesPerRecord, dtype=float),
                  zeros(recordsPerBuffer * samplesPerRecord, dtype=float)]
        ch_inst = zeros((2, recordsPerBuffer * samplesPerRecord), dtype=float)
        timepts = linspace(0, (samplesPerRecord * 1.0) / (card.config.sample_rate * 1.0e3), samplesPerRecord)

        num_pts = recordsPerBuffer
        num_data_pts = len(self.xdata)
        num_calib_pts = self.config['num_calib_seqs']
        total_length = samplesPerRecord * recordsPerBuffer

        #print num_pts,num_data_pts,num_calib_pts

        #map config settings to local variables
        num_avgs = self.config['num_avgs']
        channels = self.config['qubit_enable']
        data_window = self.config['data_window']
        get_fulldata = self.config['save_each_run']
        get_fulldata2d = self.config['save_each_run_fulldata']

        print data_window

        if plotter is not None:
            for i in range(2):
                if channels[i]:
                    plotter.init_plot("FullData" + str(i), rank=2, accum=False)
                    plotter.init_plot("Scope" + str(i), accum=False)
                    plotter.init_plot("Data" + str(i), accum=False)
                    if not homodyne_meas:
                        plotter.init_plot("Data2" + str(i), accum=False)
                    if num_calib_pts > 0:
                        plotter.init_plot("Calib" + str(i), accum=False)

        self.full_trace_data_avg = zeros((2, num_pts, samplesPerRecord))
        self.ydata_avg = zeros((2, num_pts))

        if (get_fulldata or get_fulldata2d) and num_avgs == -1:
            raise NameError("Can only collect full data if the number of averages is finite")

        if get_fulldata:
            self.ydata = zeros((2, num_avgs, num_pts))
        else:
            self.ydata = []

        if get_fulldata2d:
            self.full_trace_data = zeros((2, num_avgs, num_pts, samplesPerRecord))
        else:
            self.full_trace_data = []

        ret = card.Az.AlazarStartCapture(card.handle)
        print "Start Acquisition", ret

        #turn on the awg's
        #NOTE: TURN ON THE FIRST AWG LAST!!!
        im = InstrumentManager()
        for i in range(len(self.config['AWG']) - 1, -1, -1):
            awg = im[self.config['AWG'][i]]
            if not (i == 1 and (not self.config['start_agilent'])):
                awg.run()
                while not awg.query("*OPC?"):
                    pass

        print num_avgs

        if self.config['monitor_pulses']:
            card.capture_buffer_async()
            plotter.init_plot("Channel 1", rank=2)
            plotter.init_plot("Channel 2", rank=2)
            for ch1_rec, ch2_rec in card.get_records():
                plotter.plot(ch1_rec, 'Channel 1')
                plotter.plot(ch2_rec, 'Channel 2')
            #stop acquisition
            card.Az.AlazarAbortAsyncRead(card.handle)
        else:

            if channels[0]:
                start_index = [0, total_length]
            else:
                start_index = [0, 0]
            end_index = array(start_index) + total_length

            ch_range = array([card.config.ch1_range, card.config.ch2_range])

            while 1:
                try:
                    #Current buffer
                    buf_idx = buffersCompleted % card.config.bufferCount

                    #Get Data from current buffer
                    ret = card.capture_buffer_async(buf_idx)

                    if ret != 512:
                        #print Alazar.dsalazar.ret_to_str(ret, card.Az)
                        break

                    #ret = card.repost_buffer(buf_idx)
                    buffersCompleted += 1

                    prefact = (buffersCompleted - 1.0) / buffersCompleted * 1.0

                    #read all the data out *ONCE* into a preallocated array
                    #Otherwise the processing is too slow and the buffers overflow

                    if not (get_fulldata or get_fulldata2d):
                        for j in range(2):
                            if channels[j]:
                                #for some reason the commented out line is slow
                                #ch_sum[j] = ch_sum[j]*prefact +(card.arrs[buf_idx][start_index[j]:end_index[j]] - 128.0) * ch_range[j] / 128.0/buffersCompleted
                                ch_sum[j] = ch_sum[j] * prefact + card.arrs[buf_idx][start_index[j]:end_index[j]] * (
                                ch_range[j] / 128.0 / buffersCompleted) - (ch_range[j] / buffersCompleted)

                    else:

                        #get the data from this buffer!
                        for j in range(2):
                            if channels[j]:
                                ch_inst[j] = card.arrs[buf_idx][start_index[j]:end_index[j]] * (ch_range[j] / 128.0) - \
                                             ch_range[j]
                                ch_sum[j] = ch_sum[j] * prefact + ch_inst[j] / buffersCompleted

                                for i in range(recordsPerBuffer):
                                    if get_fulldata:
                                        if not homodyne_meas:
                                            amp, phase, _, _ = heterodyne(timepts, ch_pts=ch_inst[j][(
                                            i * samplesPerRecord + data_window[j][0]):(
                                            i * samplesPerRecord + data_window[j][1])], IFfreq=self.config['IFreq'][j],
                                                                          anti_alias=False)
                                            self.ydata[j, buffersCompleted - 1, i] = amp + phase * 1j
                                        else:
                                            self.ydata[j, buffersCompleted - 1, i] = mean(ch_inst[j][(
                                            i * samplesPerRecord + data_window[j][0]):(
                                            i * samplesPerRecord + data_window[j][1])])

                                    if get_fulldata2d:
                                        self.full_trace_data[j, buffersCompleted - 1, i] = ch_inst[j][
                                                                                           (i * samplesPerRecord):(
                                                                                           (i + 1) * samplesPerRecord)]

                    #Plot every 100 acquisitions
                    if buffersCompleted % 100 == 0:
                        print buffersCompleted

                        for i in range(recordsPerBuffer):

                            for j in range(2):
                                if channels[j]:
                                    if not homodyne_meas:
                                        amp, phase, _, _ = heterodyne(timepts, ch_pts=ch_sum[j, (
                                        i * samplesPerRecord + data_window[j][0]):(
                                        i * samplesPerRecord + data_window[j][1])], IFfreq=self.config['IFreq'][j],
                                                                      anti_alias=False)
                                        self.ydata_avg[j, i] = amp + phase * 1j
                                    else:
                                        self.ydata_avg[j, i] = mean(ch_sum[j][
                                                                    (i * samplesPerRecord + data_window[j][0]):(
                                                                    i * samplesPerRecord + data_window[j][1])])

                        for j in range(2):

                            if channels[j] and plotter is not None:
                                if not homodyne_meas:
                                    plotter.plot((self.xdata, re(self.ydata_avg[j, 0:num_data_pts])), "Data" + str(j))
                                    plotter.plot((self.xdata, im(self.ydata_avg[j, 0:num_data_pts])), "Data2" + str(j))
                                else:
                                    plotter.plot((self.xdata, self.ydata_avg[j, 0:num_data_pts]), "Data" + str(j))

                                if num_calib_pts > 0:
                                    plotter.plot((range(num_calib_pts), self.ydata_avg[j, num_data_pts:num_pts]),
                                                 "Calib" + str(j))

                                run_plot = 10
                                plotter.plot((timepts, ch_sum[j][(run_plot * samplesPerRecord):(
                                (run_plot + 1) * samplesPerRecord)]), "Scope" + str(j))

                                for jj in range(num_pts):
                                    self.full_trace_data_avg[j, jj] = ch_sum[j][(jj * samplesPerRecord):(
                                    (jj + 1) * samplesPerRecord)]
                                plotter.plot(self.full_trace_data_avg[j], "FullData" + str(j))

                        if num_avgs != -1 and buffersCompleted >= num_avgs:
                            print "Finished Averaging"
                            card.Az.AlazarAbortAsyncRead(card.get_handle())
                            break


                except (KeyboardInterrupt, SystemExit):

                    print "Exiting"
                    card.Az.AlazarAbortAsyncRead(card.get_handle())
                    self.config['num_avgs'] = buffersCompleted

                except:

                    print "Error", sys.exc_info()
                    card.Az.AlazarAbortAsyncRead(card.get_handle())
                    print str(buffersCompleted)


class rabi(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'RABIOSC'

        #custom Rabi experiment settings
        self.config['rabi_height'] = [0.3, 0.3]
        self.config['rabi_length'] = [10, 10]
        self.config['rabi_vary'] = 'height'  #height or width ..then the xdata is the vary parameter
        self.config['drive_angle'] = 0.0
        self.config['height_scaling'] = [1.0, 1.0]


    #load the rabi oscillation into the AWG
    def load_exp(self, plotter=None, rabi_index=-1):

        if len(self.xdata) == 0:
            raise NameError("Rabi sequence not initialized")

        if self.config['rabi_vary'] == 'height':
            self._rabi_vs_width = False
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (
            max(self.config['rabi_length']))
        elif self.config['rabi_vary'] == 'width':
            self._rabi_vs_width = True
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (self.xdata[-1])

        if rabi_index == -1:
            rabi_vary = self.xdata
        else:
            rabi_vary = zeros(1)
            rabi_vary[0] = self.xdata[rabi_index]


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #allocate pulse array
        #self.rabi_fullpulsesequence = zeros((max(len(self.rabi_height),len(self.rabi_length)),6,self.config['total_length']*awg_modifier))

        awg_seq = pulse_sequence(len(rabi_vary), self.config['total_length'] * self.config['awg_clocks'][0],
                                 TEK_clk=self.config['awg_clocks'][0], Agilent_clk=self.config['awg_clocks'][1],
                                 delays=self.config['awg_delay'])


        #define card trigger and read trigger once (these are the same waveforms each time)
        for j in range(2):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (max_rabi + self.config['read_delay']), self.config['read_length'])

        if trigger_card_early or self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((max_rabi + self.config['card_delay']),
                                                                               self.config['read_length'])

        for i in range(len(rabi_vary)):


            print rabi_vary[i]

            for j in range(2):

                if not self.config['qubit_enable'][j]:
                    continue

                if self._rabi_vs_width:
                    rabi_pulse_i = floor(rabi_vary[i])
                    pulse_height_i = self.config['rabi_height'][j]
                else:
                    rabi_pulse_i = floor(self.config['rabi_length'][j])
                    pulse_height_i = rabi_vary[i]
                pulse_height_i *= self.config['height_scaling'][j]
                if self.config['shaped_pulse']:

                    if len(rabi_vary) > 1:
                        ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind1
                        ind2a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind2a
                        ind2b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind2b
                    else:
                        ind1 = 0
                        ind2a = 0
                        ind2b = 0

                    if self.config['pulsed_drive'][j]:
                        awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(
                            (max_rabi - rabi_pulse_i), (2 * rabi_pulse_i + 2 * self.config['generator_enable_delay']))

                    if self.config['drive_sideband'][j]:
                        awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse_with_freq(
                            (max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                            self.config['drive_sideband_freq'][j] / 1.0e9, 0.0, False)
                        awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse_with_freq(
                            (max_rabi - rabi_pulse_i), rabi_pulse_i / 2,
                            -1 * sign(self.config['drive_sideband_freq'][j]) * pulse_height_i,
                            self.config['drive_sideband_freq'][j] / 1.0e9, pi / 2.0, False)
                    else:
                        if self.config['drag_pulse'][j]:
                            awg_seq.DRAG_pulse((max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                                               self.config['drive_angle'], j, ind2a, ind2b,
                                               self.config['drag_prefactor'][j], False)
                        else:
                            awg_seq.gauss_pulse((max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                                                self.config['drive_angle'], j, ind2a, ind2b, False)

                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
                    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])


                else:

                    if len(rabi_vary) > 1:
                        ind = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind
                    else:
                        ind = 0

                    awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind].square_pulse2((max_rabi - rabi_pulse_i),
                                                                                              rabi_pulse_i)

            if plotter is not None:
                cur_pulse_seq.plot_pulses(plotter)

        #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #self.rabi_fullpulsesequence[i] = cur_pulse_seq.output_pulses()


    #run the rabi oscillation    
    def run_rabi_sequential(self, num_avgs=1024, plotter=None):

        #define instrument manager    
        im = InstrumentManager()

        #setup awg
        awg = im[self.config['AWG']]

        if self.config['LO'] == self.config['RF']:
            homodyne_meas = True
            print "Homodyne measurement"
        else:
            homodyne_meas = False

        #setup drive

        drive = im[self.config['drive']]

        #setup drive
        if self.config['drive'][0:2] == "LB":
            drive.set_output(True)
            drive.set_pulse_ext(mod=True)
            #Internal pulse must be set to false!
            drive.set_mod(False)
            drive.set_power(self.config['rabi_drive_pwr'])  #-28 when direct driving
            drive.set_frequency(self.config['rabi_drive_freq'])
        else:
            drive.set_power(self.config['rabi_drive_pwr'])
            drive.set_mod(True)
            drive.set_frequency(self.config['rabi_drive_freq'])
            drive.set_ext_pulse()
            drive.set_output(True)

        for i in range(5):
            if abs(drive.get_frequency() - self.config['rabi_drive_freq']) > 100:
                drive.set_frequency(self.config['rabi_drive_freq'])
            else:
                break
        if abs(drive.get_frequency() - self.config['rabi_drive_freq']) > 100:
            print "Lab brick frequency error!"
        time.sleep(.2)

        print self.config['rabi_drive_freq'], drive.get_frequency()

        #setup initial drive (if applicable)
        if self.config['rabi_b_drive_pulse_length'] > 0:
            drive2 = im[self.config['drive_b']]
            #setup drive2 (assuming it's a lab brick)
            drive2.set_output(True)
            drive2.set_pulse_ext(mod=True)
            #Internal pulse must be set to false!
            drive2.set_mod(False)
            drive2.set_power(self.config['rabi_b_drive_pwr'])  #-28 when direct driving
            drive2.set_frequency(self.config['rabi_b_drive_freq'])
            print drive2.get_frequency()

        #setup LO and RF

        RF = im[self.config['RF']]
        RF.set_output(True)
        RF.set_mod(True)
        RF.set_power(self.config['rabi_read_pwr'])
        RF.set_frequency(self.config['rabi_read_freq'])
        RF.set_ext_pulse()

        if not homodyne_meas:
            LO = im[self.config['LO']]
            LO.set_power(self.config['rabi_lo_power'])
            LO.set_mod(True)
            LO.set_frequency(self.config['rabi_read_freq'] + self.config['IFreq'])
            LO.set_output(True)



        #setup the flux
        flux = im[self.config['flux']]

        if self.config['flux'] == 'flux2' or self.config['flux'] == 'flux1':
            flux.set_function("DC")
            flux.set_offset(self.config['rabi_flux'])
        else:
            #flux.set_mode('current')
            flux.set_volt(self.config['rabi_flux'])

        #setup the Alazar card
        card = setup_Alazar_card(512, self.config['acq_length'], self.config['meas_range'],
                                 self.config['monitor_pulses'], True, num_avgs)

        if self._rabi_vs_width:
            xdata = self.rabi_length
        else:
            xdata = self.rabi_height

            #setup plots
        if plotter is not None:
            plotter.init_plot("FullData", rank=2, accum=False)
            plotter.init_plot("Scope", accum=False)
            plotter.init_plot("Amp", accum=False)
            plotter.init_plot("Phase", accum=False)

        self.rabi_amp = zeros(len(xdata))
        self.rabi_phase = zeros(len(xdata))
        self.rabi_fulldata = zeros((len(xdata), self.config['acq_length']))

        for i in range(len(xdata)):

            #program sequence
            awg.pre_experiment()
            self.load_rabi(False, None, i)
            awg.run_experiment(self.config['seq_file'])

            tpts, ch1_pts, ch2_pts = card.acquire_avg_data()
            if not homodyne_meas:
                amp, phase, _, _ = heterodyne(tpts, ch1_pts, ch2_pts, self.config['IFreq'])
            else:
                amp = mean(ch1_pts)
                phase = amp

            self.rabi_phase[i] = phase
            self.rabi_amp[i] = amp
            self.rabi_fulldata[i] = ch1_pts

            if plotter is not None:
                plotter.plot((xdata[0:(i + 1)], self.rabi_phase[0:(i + 1)]), "Phase")
                plotter.plot((xdata[0:(i + 1)], self.rabi_amp[0:(i + 1)]), "Amp")
                plotter.plot((range(self.config['acq_length']), ch1_pts), "Scope")
                plotter.plot(self.rabi_fulldata, "FullData")


    def fit_rabi(self):

        #start the fit after fit_start
        if len(self.rabi_height) > 1:
            xdata = self.rabi_height
        else:
            xdata = self.rabi_length

        start_index = -1
        for i in range(len(xdata)):
            if xdata[i] >= self.config['fit_start']:
                start_index = i
                break

        if start_index == -1:
            start_index = 0

        #fit the t1 data to an decaying sinusoid 
        if self.config['fit_phase']:
            ydata = self.rabi_phase[start_index:-1]
            xdata = xdata[start_index:-1]

        else:
            ydata = self.rabi_amp[start_index:-1]
            xdata = xdata[start_index:-1]

        self.config['fit_start'] = xdata[0]

        fitvals = fitdecaysin(xdata, ydata)

        self.config['tau'] = fitvals[3]
        self.config['A'] = fitvals[0]
        self.config['y0'] = fitvals[4]
        self.config['omega'] = fitvals[1]
        self.config['phase'] = fitvals[2]

        print fitvals

    def display_rabi(self, plotter, pulse_sequence=False, display_fit=True):


        plotter.init_plot("FullData", rank=2, accum=False)
        plotter.init_plot("Amp", rank=1, accum=False)
        plotter.init_plot("Phase", rank=1, accum=False)

        plotter.plot((self.rabi_fulldata), "FullData")

        p = [self.config['A'], self.config['omega'], self.config['phase'], self.config['tau'], self.config['y0'],
             self.config['fit_start']]

        if len(self.rabi_height) > 1:
            xdata = self.rabi_height
        else:
            xdata = self.rabi_length

        if self.config['fit_phase'] and display_fit:
            plotter.plot((concatenate((xdata, xdata)), concatenate((self.rabi_phase, decaysin(p, xdata)))), "Phase")
        else:
            plotter.plot((xdata, self.rabi_phase), "Phase")

        time.sleep(.2)

        if (not self.config['fit_phase']) and display_fit:
            plotter.plot((concatenate((xdata, xdata)), concatenate((self.rabi_amp, decaysin(p, xdata)))), "Amp")
        else:
            plotter.plot((xdata, self.rabi_amp * 1.0), "Amp")


#drive Rabi |1>-->|2> using sideband 
class rabi2(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'RABI2'

        #custom Rabi experiment settings
        self.config['rabi_height1'] = [0.3, 0.3]
        self.config['rabi_height2'] = [0.3, 0.3]
        self.config['rabi_length1'] = [10, 10]
        self.config['rabi_length2'] = [10, 10]
        self.config['rabi_sideband_freq1'] = [120e6, 120e6]
        self.config['rabi_sideband_freq2'] = [-120e6, -120e6]
        self.config['pulse1_2_gap'] = [15, 15]
        self.config['rabi_vary'] = 'height'  #height or width ..then the xdata is the vary parameter
        self.config['read_sweep'] = False
        self.config['sideband_sweep'] = False
        self.config['read_freq'] = [5.5e9, 5.5e9]
        #self.config['drive_angle'] = 0.0
        if self.config['read_sweep']:
            self.config['read_sweep_end'] = [0.7, 0.7]
            self.config['read_sweep_points'] = [50, 50]
            self.config['SB_read_points'] = [linspace(self.config['read_freq'][0], self.config['read_sweep_end'][0],
                                                      self.config['read_sweep_points'][0]),
                                             linspace(self.config['read_freq'][1], self.config['read_sweep_end'][1],
                                                      self.config['read_sweep_points'][1])]
        if self.config['sideband_sweep']:
            self.config['sideband_sweep_end'] = [self.config['rabi_sideband_freq2'][0] + 25e6,
                                                 self.config['rabi_sideband_freq2'][1]]
            self.config['sideband_sweep_points'] = [26, 26]
            self.config['SB_points'] = [
                linspace(self.config['rabi_sideband_freq2'][0], self.config['sideband_sweep_end'][0],
                         self.config['sideband_sweep_points'][0]),
                linspace(self.config['rabi_sideband_freq2'][1], self.config['sideband_sweep_end'][1],
                         self.config['sideband_sweep_points'][1])]


        #don't use the agilent for this!
        self.config['start_agilent'] = False

    def run_qubit_exp(self, plotter=None):

        if self.config['read_sweep']:

            if self.config['num_avgs'] == -1:
                raise NameError('Must run a finite number of averages')

            #XXXX
            temp_ydata = zeros((2, len(self.config['SB_read_points'][0]), len(self.xdata)))

            num_avgs = self.config['num_avgs']

            if plotter is not None:
                for i in range(2):
                    plotter.init_plot("FullSpectData" + str(i), rank=2, accum=False)

            #XXX               
            for i in range(len(self.config['SB_read_points'][0])):


                #set the drive frequencies  
                #XXX

                self.config['read_freq'] = [self.config['SB_read_points'][0][i], self.config['SB_read_points'][1][i]]

                #XXX
                print self.config['read_freq']

                rabi2.load_exp(self)

                #run!
                qubit_exp.run_qubit_exp(self, plotter)

                if self.config['num_avgs'] != num_avgs:
                    #was interrupted prematurely
                    print("Sideband sweep spectroscopy interrupted")
                    break

                for j in range(2):
                    temp_ydata[j, i] = self.ydata_avg[j]

                if plotter is not None:
                    for j in range(2):
                        plotter.plot(temp_ydata[j], "FullSpectData" + str(j))

            self.ydata = temp_ydata

        elif self.config['sideband_sweep']:

            if self.config['num_avgs'] == -1:
                raise NameError('Must run a finite number of averages')

            #XXXX
            temp_ydata = zeros((2, len(self.config['SB_points'][0]), len(self.xdata)))

            num_avgs = self.config['num_avgs']

            if plotter is not None:
                for i in range(2):
                    plotter.init_plot("FullSpectData" + str(i), rank=2, accum=False)

            #XXX               
            for i in range(len(self.config['SB_points'][0])):


                #set the drive frequencies  
                #XXX

                self.config['rabi_sideband_freq2'] = [self.config['SB_points'][0][i], self.config['SB_points'][1][i]]

                #XXX
                print self.config['rabi_sideband_freq2']

                rabi2.load_exp(self)

                #run!
                qubit_exp.run_qubit_exp(self, plotter)

                if self.config['num_avgs'] != num_avgs:
                    #was interrupted prematurely
                    print("Sideband sweep spectroscopy interrupted")
                    break

                for j in range(2):
                    temp_ydata[j, i] = self.ydata_avg[j]

                if plotter is not None:
                    for j in range(2):
                        plotter.plot(temp_ydata[j], "FullSpectData" + str(j))

            self.ydata = temp_ydata

        else:
            #run normally!
            qubit_exp.run_qubit_exp(self, plotter)

    #load the rabi oscillation into the AWG
    def load_exp(self, plotter=None, rabi_index=-1):

        if len(self.xdata) == 0:
            raise NameError("Rabi sequence not initialized")

        #calculate the maximum length of the rabi
        max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (
        max(self.config['rabi_length1'])) + max(self.config['pulse1_2_gap'])
        if self.config['rabi_vary'] == 'height':
            self._rabi_vs_width = False
            max_rabi += 2 * (max(self.config['rabi_length2']))
        elif self.config['rabi_vary'] == 'width':
            self._rabi_vs_width = True
            max_rabi += 2 * (self.xdata[-1])

        if rabi_index == -1:
            rabi_vary = self.xdata
        else:
            rabi_vary = zeros(1)
            rabi_vary[0] = self.xdata[rabi_index]


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #allocate pulse array
        #self.rabi_fullpulsesequence = zeros((max(len(self.rabi_height),len(self.rabi_length)),6,self.config['total_length']*awg_modifier))

        awg_seq = pulse_sequence(len(rabi_vary), self.config['total_length'] * self.config['awg_clocks'][0],
                                 TEK_clk=self.config['awg_clocks'][0], Agilent_clk=self.config['awg_clocks'][1],
                                 delays=self.config['awg_delay'])


        #define card trigger and read trigger once (these are the same waveforms each time)
        for j in range(2):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (max_rabi + self.config['read_delay']), self.config['read_length'])

        if trigger_card_early or self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((max_rabi + self.config['card_delay']),
                                                                               self.config['read_length'])

        for i in range(len(rabi_vary)):


            print rabi_vary[i]

            for j in range(2):

                if not self.config['qubit_enable'][j]:
                    continue

                if self._rabi_vs_width:
                    rabi_pulse_i = floor(rabi_vary[i])
                    pulse_height_i = self.config['rabi_height2'][j]
                else:
                    rabi_pulse_i = floor(self.config['rabi_length2'][j])
                    pulse_height_i = rabi_vary[i]

                if len(rabi_vary) > 1:
                    ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                    awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind1
                    ind2a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind2a
                    ind2b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind2b
                else:
                    ind1 = 0
                    ind2a = 0
                    ind2b = 0

                #always drive sideband
                drive_freq1 = self.config['rabi_sideband_freq1'][j]
                drive_freq2 = self.config['rabi_sideband_freq2'][j]

                #first pulse
                pulse1_len = self.config['rabi_length1'][j]
                pulse1_height = self.config['rabi_height1'][j]

                #second pulse
                pulse2_len = self.config['rabi_length2'][j]
                pulse2_height = self.config['rabi_height2'][j]

                pulse1_center = max_rabi - 2 * rabi_pulse_i - self.config['pulse1_2_gap'][j] - 3 * pulse1_len
                pulse2_center = max_rabi - rabi_pulse_i - 2 * pulse1_len

                #drive trigger
                awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(
                    (pulse1_center + max_rabi - pulse1_len) / 2,
                    max_rabi - pulse1_len - pulse1_center + 2 * self.config['generator_enable_delay'])

                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse_with_freq(pulse1_center,
                                                                                                   pulse1_len / 2,
                                                                                                   pulse1_height, abs(
                        drive_freq1) / 1.0e9, 0.0, False)
                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse_with_freq(pulse1_center,
                                                                                                   pulse1_len / 2,
                                                                                                   1 * sign(
                                                                                                       drive_freq1) * pulse1_height,
                                                                                                   abs(
                                                                                                       drive_freq1) / 1.0e9,
                                                                                                   pi / 2.0, False)

                #second pulse
                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse_with_freq(pulse2_center,
                                                                                                   rabi_pulse_i / 2,
                                                                                                   pulse_height_i, abs(
                        drive_freq2) / 1.0e9, 0.0, False)
                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse_with_freq(pulse2_center,
                                                                                                   rabi_pulse_i / 2,
                                                                                                   1 * sign(
                                                                                                       drive_freq2) * pulse_height_i,
                                                                                                   abs(
                                                                                                       drive_freq2) / 1.0e9,
                                                                                                   pi / 2.0, False)

                #pulse for swapping |e> and |f> population
                #                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse_with_freq(max_rabi-2*pulse1_len-pulse2_len, pulse2_len/2, pulse2_height, abs(drive_freq2)/1.0e9, 0.0, False)
                #                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse_with_freq(max_rabi-2*pulse1_len-pulse2_len, pulse2_len/2, 1*sign(drive_freq2)*pulse2_height, abs(drive_freq2)/1.0e9, pi/2.0, False)

                #pulse for swapping |g> and |e> population
                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse_with_freq(
                    max_rabi - pulse1_len, pulse1_len / 2, pulse1_height, abs(drive_freq1) / 1.0e9, 0.0, False)
                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse_with_freq(
                    max_rabi - pulse1_len, pulse1_len / 2, 1 * sign(drive_freq1) * pulse1_height,
                    abs(drive_freq1) / 1.0e9, pi / 2.0, False)

            if plotter is not None:
                cur_pulse_seq.plot_pulses(plotter)

        #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #self.rabi_fullpulsesequence[i] = cur_pulse_seq.output_pulses()


class t1(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'T1'

        self.config['t1_pulse_height'] = [0.3, 0.3]
        self.config['t1_pulse_length'] = [15, 15]

        self.config['pulse_type'] = [1, 1]  #1: gaussian, 2: smooth square (useful if T1 is short)


    #load the rabi oscillation into the AWG
    def load_exp(self, plotter=None, wait_index=-1):

        if len(self.xdata) == 0:
            raise NameError("T1 sequence not initialied")

        if wait_index == -1:
            wait_times = self.xdata
        else:
            wait_times = zeros(1)
            wait_times[0] = self.xdata[wait_index]
            print wait_times

        max_t1_wait = 2 * ceil(self.xdata[-1])


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False


        #where to start the card and read pulses
        max_t1_wait = self.config['front_buffer'] + self.config['generator_enable_delay'] + self.xdata[-1] + 1 * max(
            self.config['t1_pulse_length'])

        awg_seq = pulse_sequence(len(wait_times), self.config['total_length'] * self.config['awg_clocks'][0],
                                 TEK_clk=self.config['awg_clocks'][0], Agilent_clk=self.config['awg_clocks'][1],
                                 delays=self.config['awg_delay'])

        #read_triggers
        for j in range(2):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (max_t1_wait + self.config['read_delay']), self.config['read_length'])

        #card trigger
        if self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            #Read output triggers the card now (not the marker bit)
            if trigger_card_early:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
            else:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (max_t1_wait + self.config['card_delay']), self.config['read_length'])

        for i in range(len(wait_times)):

            print i

            for j in range(2):

                if self.config['qubit_enable'][j]:

                    t1_pulse_center = max_t1_wait - wait_times[i] - self.config['t1_pulse_length'][j] / 2

                    if self.config['shaped_pulse']:

                        if len(wait_times) > 1:
                            ind1 = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                            awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind1
                            ind2 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                            awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind2
                        else:
                            ind1 = 0
                            ind2 = 0

                        awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse2(t1_pulse_center, (
                        1 * self.config['t1_pulse_length'][j] + 2 * self.config['generator_enable_delay']))

                        if self.config['pulse_type'][j] == 1:
                            awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1].gauss_pulse(t1_pulse_center,
                                                                                                    self.config[
                                                                                                        't1_pulse_length'][
                                                                                                        j] / 2,
                                                                                                    self.config[
                                                                                                        't1_pulse_height'][
                                                                                                        j], False)
                        elif self.config['pulse_type'][j] == 2:
                            awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1].smooth_square(t1_pulse_center,
                                                                                                      25, self.config[
                                    't1_pulse_length'][j], self.config['t1_pulse_height'][j])

                    else:

                        if len(wait_times) > 1:
                            ind = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                            awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind
                        else:
                            ind = 0

                        awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind].square_pulse(
                            (1 + max_t1_wait - wait_times[i] - self.config['t1_pulse_length']) * awg_modifier,
                            self.config['t1_pulse_length'] * awg_modifier)


        #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

    #run the t1 scan, but one point at a time
    def run_t1_sequential(self, num_avgs=1024, plotter=None):

        #define instrument manager    
        im = InstrumentManager()

        #setup awg
        awg = im[self.config['AWG']]

        if self.config['LO'] == self.config['RF']:
            homodyne_meas = True
            print "Homodyne measurement"
        else:
            homodyne_meas = False

            #setup drive,RF,LO
        drive = im[self.config['drive']]

        #setup drive
        if self.config['drive'][0:2] == "LB":
            drive.set_output(True)
            drive.set_pulse_ext(mod=True)
            #Internal pulse must be set to false!
            drive.set_mod(False)
            drive.set_power(self.config['t1_drive_pwr'])  #-28 when direct driving
            drive.set_frequency(self.config['t1_drive_freq'])
        else:
            drive.set_output(True)
            drive.set_ext_pulse()
            drive.set_mod(True)
            drive.set_power(self.config['t1_drive_pwr'])  #-28 when direct driving
            drive.set_frequency(self.config['t1_drive_freq'])

        print self.config['t1_drive_freq']

        #setup LO and RF 

        RF = im[self.config['RF']]
        RF.set_output(True)
        RF.set_ext_pulse()
        RF.set_mod(True)
        RF.set_power(self.config['t1_read_pwr'])
        RF.set_frequency(self.config['t1_read_freq'])

        if not homodyne_meas:
            LO = im[self.config['LO']]
            LO.set_output(True)
            LO.set_internal_pulse(40e-6)
            LO.set_mod(True)
            LO.set_power(self.config['t1_lo_power'])
            LO.set_frequency(self.config['t1_read_freq'] + self.config['IFreq'])




        #setup the flux
        flux = im[self.config['flux']]
        flux.set_mode('voltage')
        if self.config['flux'] == 'flux2' or self.config['flux'] == 'flux1':
            flux.set_function("DC")
            flux.set_offset(self.config['t1_flux'])
        else:
            flux.set_volt(self.config['t1_flux'])

        #setup the Alazar card
        card = setup_Alazar_card(512, self.config['acq_length'], self.config['meas_range'],
                                 self.config['monitor_pulses'], True, num_avgs)

        #setup plots
        if plotter is not None:
            plotter.init_plot("FullData", rank=2, accum=False)
            plotter.init_plot("Scope", accum=False)
            plotter.init_plot("Amp", accum=False)
            plotter.init_plot("Phase", accum=False)

        self.t1_amp = zeros(len(self.t1_wait))
        self.t1_phase = zeros(len(self.t1_wait))
        self.t1_fulldata = zeros((len(self.t1_wait), self.config['acq_length']))

        #program sequence
        awg.pre_experiment()
        for i in range(len(self.t1_wait)):

            awg.stop()
            self.load_t1(False, None, i)
            awg.prep_experiment(self.config['seq_file'])
            awg.set_amps_offsets(self.config['awg_amp'], self.config['awg_offset'])
            awg.run()

            tpts, ch1_pts, ch2_pts = card.acquire_avg_data()
            if not homodyne_meas:
                amp, phase, _, _ = heterodyne(tpts, ch1_pts, ch2_pts, self.config['IFreq'])
            else:
                amp = mean(ch1_pts)
                phase = amp

            self.t1_phase[i] = phase
            self.t1_amp[i] = amp
            self.t1_fulldata[i] = ch1_pts

            if plotter is not None:
                plotter.plot((self.t1_wait[0:(i + 1)], self.t1_phase[0:(i + 1)]), "Phase")
                plotter.plot((self.t1_wait[0:(i + 1)], self.t1_amp[0:(i + 1)]), "Amp")
                plotter.plot((range(self.config['acq_length']), ch1_pts), "Scope")
                plotter.plot(self.t1_fulldata, "FullData")


    def fit_t1(self):

        #start the fit after fit_start
        start_index = -1
        for i in range(len(self.t1_wait)):
            if self.t1_wait[i] > self.config['fit_start']:
                start_index = i
                break

        if start_index == -1:
            start_index = 0

        #fit the t1 data to an exponential
        if self.config['fit_phase']:
            ydata = self.t1_phase[start_index:-1]
            xdata = self.t1_wait[start_index:-1]

        else:
            ydata = self.t1_amp[start_index:-1]
            xdata = self.t1_wait[start_index:-1]

        fitvals = fitgeneral(xdata, ydata, lambda p, x: expfunc2(p, x, self.config['fit_start']),
                             [ydata[-1], (ydata[0] - ydata[-1]), xdata[-1] / 2])

        self.config['tau'] = fitvals[2]
        self.config['A'] = fitvals[1]
        self.config['y0'] = fitvals[0]

        print fitvals

    def display_t1(self, plotter, pulse_sequence=False, display_fit=True):


        plotter.init_plot("FullData", rank=2, accum=False)
        plotter.init_plot("Amp", accum=False)
        plotter.init_plot("Phase", accum=False)

        plotter.plot((self.t1_fulldata), "FullData")

        if self.config['fit_phase'] and display_fit:
            plotter.plot((concatenate((self.t1_wait, self.t1_wait)), concatenate((self.t1_phase, expfunc2(
                [self.config['y0'], self.config['A'], self.config['tau']], self.t1_wait, self.config['fit_start'])))),
                         "Phase")
        else:
            plotter.plot((self.t1_wait, self.t1_phase), "Phase")

        if not self.config['fit_phase'] and display_fit:
            plotter.plot((concatenate((self.t1_wait, self.t1_wait)), concatenate((self.t1_amp, expfunc2(
                [self.config['y0'], self.config['A'], self.config['tau']], self.t1_wait, self.config['fit_start'])))),
                         "Amp")
        else:
            plotter.plot((self.t1_wait, self.t1_amp), "Amp")


class ramsey(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'RAMSEY'

        #custom Ramsey experiment settings
        self.config['ramsey_pulse_height'] = [0.3, 0.3]  #first ramsey pulse height
        self.config['ramsey_pulse_height2'] = [0.3, 0.3]  #second ramsey pulse height
        self.config['ramsey_pulse_length'] = [10, 10]
        self.config['ramsey_time'] = [100, 100]
        self.config['ramsey_vary'] = 'time'  #time, time2 or flux or flux2
        self.config['phase_advance_rate'] = [0.0, 0.0]

        #spin echo
        self.config['spin_echo'] = False
        self.config['echo_pulse_height'] = [0.3, 0.3]
        self.config['echo_pulse_length'] = [20, 20]

        #Properties if we are doing a flux pulse
        self.config['do_fast_flux'] = False  #Note if ramsey_vary is flux then this will also do the fast flux
        self.config['ramsey_fast_flux'] = [0.0, 0.0]
        self.config['flux_wait_time_percent'] = [0.75,
                                                 0.75]  #fraction of the wait time that that the flux pulse is applied
        self.config['flux_up_time'] = [25, 25]  #time to the fast flux pulse height
        self.config['trap_flux_pulse'] = [False, False]
        self.config['trap_slope'] = [0.1, 0.1]


    #load the ramsey oscillation into the AWG
    def load_exp(self, plotter=None, ramsey_index=-1):


        if len(self.xdata) == 0:
            raise NameError("Ramsey sequence not initialied - Ramsey Variation Array is Null")



        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False


        #where to start the card and read pulses
        max_ramsey_wait = self.config['front_buffer'] + self.config['generator_enable_delay'] + 5 * max(
            self.config['ramsey_pulse_length'])


        #determine if this is a Ramsey where we vary the wait time OR a z Ramsey where we vary the flux pulse height
        if self.config['ramsey_vary'] == 'flux' or self.config['ramsey_vary'] == 'flux2':
            vary_array = self.xdata
            vary_flux = True
            max_ramsey_wait += max(self.config['ramsey_time'])
        else:
            vary_array = self.xdata
            vary_flux = False
            max_ramsey_wait += self.xdata[-1]

        if ramsey_index == -1:
            #wait_times = self.ramsey_time
            pass
        else:
            vary_array = self.xdata[ramsey_index]
            print vary_array

        if self.config['spin_echo']:
            max_ramsey_wait += 2 * max(self.config['echo_pulse_length'])

        awg_seq = pulse_sequence(len(vary_array), int(ceil(self.config['total_length'] * self.config['awg_clocks'][0])),
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][1])),
                                 self.config['awg_clocks'][0], self.config['awg_clocks'][1], self.config['awg_delay'])

        for i in range(2):
            if self.config['qubit_enable'][i]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][i]][0].square_pulse(
                    (max_ramsey_wait + self.config['read_delay']), self.config['read_length'])


        #set up the card trigger
        if self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            if trigger_card_early:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (1 + self.config['generator_enable_delay']), self.config['read_length'])
            else:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (max_ramsey_wait + self.config['card_delay']), self.config['read_length'])

        flux_pulsing = False

        for i in range(len(vary_array)):

            for j in range(2):

                if self.config['qubit_enable'][j]:

                    if vary_flux:
                        wait_time = self.config['ramsey_time'][j]
                    else:
                        wait_time = vary_array[i]

                    ramsey_pulse_center1 = max_ramsey_wait - wait_time
                    ramsey_pulse_center2 = max_ramsey_wait
                    if self.config['ramsey_vary'] == 'time2':
                        if j == 0:
                            ramsey_pulse_center1 = max_ramsey_wait - self.config['ramsey_time'][1] / 2. - wait_time / 2.
                            ramsey_pulse_center2 = max_ramsey_wait - self.config['ramsey_time'][1] / 2. + wait_time / 2.
                        else:
                            ramsey_pulse_center1 = max_ramsey_wait - self.config['ramsey_time'][1]

                    spin_echo_center = (ramsey_pulse_center1 + ramsey_pulse_center2) / 2.0 - \
                                       self.config['echo_pulse_length'][j]

                    if self.config['spin_echo']:
                        ramsey_pulse_center1 -= 2 * self.config['echo_pulse_length'][j]

                    #set the fast flux pluse
                    if vary_flux or self.config['ramsey_fast_flux'][j] != 0 or self.config['do_fast_flux']:
                        flux_pulsing = True
                        if vary_flux and self.config['ramsey_vary'] == 'flux':
                            flux_val = vary_array[i]
                        else:
                            flux_val = self.config['ramsey_fast_flux'][j]

                        if self.config['ramsey_vary'] == 'flux2':
                            slope_up = vary_array[i]
                        else:
                            slope_up = self.config['trap_slope'][j]

                        indflux = awg_seq.add_analog_pulse(awg_seq.channel_index['flux'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['flux'][j]][i] = indflux

                        if flux_val != 0:

                            if self.config['trap_flux_pulse'][j]:
                                awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].trapezoid(
                                    (ramsey_pulse_center1 + ramsey_pulse_center2) / 2.0, slope_up,
                                    self.config['trap_slope'][j], wait_time * self.config['flux_wait_time_percent'][j],
                                    flux_val, self.config['ramsey_fast_flux'][j])

                                #do negative pulse after experiment is done
                                awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].trapezoid((self.config[
                                                                                                           'total_length'] - 2 * (
                                                                                                       wait_time *
                                                                                                       self.config[
                                                                                                           'flux_wait_time_percent'][
                                                                                                           j])),
                                                                                                      -slope_up, -
                                    self.config['trap_slope'][j], wait_time * self.config['flux_wait_time_percent'][j],
                                                                                                      -flux_val, -
                                    self.config['ramsey_fast_flux'][j])
                            else:
                                awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].smooth_square(
                                    (ramsey_pulse_center1 + ramsey_pulse_center2) / 2.0 * awg_clk_ratio,
                                    self.config['flux_up_time'][j],
                                    wait_time * self.config['flux_wait_time_percent'][j], flux_val)

                                #do negative pulse after experiment is done
                                awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].smooth_square((self.config[
                                                                                                               'total_length'] - 2 * (
                                                                                                           3 *
                                                                                                           self.config[
                                                                                                               'flux_up_time'][
                                                                                                               j] + wait_time *
                                                                                                           self.config[
                                                                                                               'flux_wait_time_percent'][
                                                                                                               j])),
                                                                                                          self.config[
                                                                                                              'flux_up_time'][
                                                                                                              j],
                                                                                                          wait_time *
                                                                                                          self.config[
                                                                                                              'flux_wait_time_percent'][
                                                                                                              j],
                                                                                                          -flux_val)


                                #setup the Ramsey pulses
                    if (not vary_flux and len(vary_array) > 1) or awg_seq.channel_index['flux'][j] != (j + 4):
                        if self.config['shaped_pulse']:
                            ind1a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                            awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind1a

                            ind1b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                            awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind1b

                        ind2 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind2
                    else:
                        ind1a = 0
                        ind1b = 0
                        ind2 = 0

                    if ind1a > 0 or (ind1a == 0 and i == 0):

                        if self.config['shaped_pulse']:


                            if self.config['drive_sideband'][j]:
                                #pulse 1
                                #cur_pulse_seq.drive_trigger.square_pulse2(ramsey_pulse_center1*awg_modifier, (3*self.config['ramsey_pulse_length'] + 2*self.config['generator_enable_delay'])*awg_modifier) 
                                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse_with_freq(
                                    ramsey_pulse_center1, self.config['ramsey_pulse_length'][j] / 2,
                                    self.config['ramsey_pulse_height'][j], self.config['drive_sideband_freq'][j] / 1e9,
                                    False)

                                #pulse 2
                                #cur_pulse_seq.drive_trigger.square_pulse2(ramsey_pulse_center2*awg_modifier, (3*self.config['ramsey_pulse_length'] + 2*self.config['generator_enable_delay'])*awg_modifier) 
                                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse_with_freq(
                                    ramsey_pulse_center2, self.config['ramsey_pulse_length'][j] / 2,
                                    self.config['ramsey_pulse_height2'][j], self.config['drive_sideband_freq'][j] / 1e9,
                                    False)
                                #awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].gauss_pulse(ramsey_pulse_center2, self.config['ramsey_pulse_length'][j]/2, sin(self.config['phase_advance_rate'][j]*wait_time)*self.config['ramsey_pulse_height2'][j], False)

                            else:
                                #pulse 1
                                #cur_pulse_seq.drive_trigger.square_pulse2(ramsey_pulse_center1*awg_modifier, (3*self.config['ramsey_pulse_length'] + 2*self.config['generator_enable_delay'])*awg_modifier) 
                                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse(
                                    ramsey_pulse_center1, self.config['ramsey_pulse_length'][j] / 2,
                                    self.config['ramsey_pulse_height'][j], False)

                                #pulse 2
                                #cur_pulse_seq.drive_trigger.square_pulse2(ramsey_pulse_center2*awg_modifier, (3*self.config['ramsey_pulse_length'] + 2*self.config['generator_enable_delay'])*awg_modifier) 
                                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse(
                                    ramsey_pulse_center2, self.config['ramsey_pulse_length'][j] / 2,
                                    cos(self.config['phase_advance_rate'][j] * wait_time) *
                                    self.config['ramsey_pulse_height2'][j], False)
                                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].gauss_pulse(
                                    ramsey_pulse_center2, self.config['ramsey_pulse_length'][j] / 2,
                                    sin(self.config['phase_advance_rate'][j] * wait_time) *
                                    self.config['ramsey_pulse_height2'][j], False)

                            #keep the drive switch on for the whole time 
                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse2(
                                (ramsey_pulse_center2 + ramsey_pulse_center1) / 2, (
                                ramsey_pulse_center2 - ramsey_pulse_center1 + 3 * self.config['ramsey_pulse_length'][
                                    j] + 2 * self.config['generator_enable_delay']))

                            #spin echo
                            if self.config['spin_echo']:
                                awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse2(
                                    spin_echo_center, (3 * self.config['echo_pulse_length'][j] + 2 * self.config[
                                        'generator_enable_delay']))
                                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse(
                                    spin_echo_center, self.config['echo_pulse_length'][j] / 2,
                                    self.config['echo_pulse_height'][j], False)


                        else:

                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse2(
                                ramsey_pulse_center1, self.config['ramsey_pulse_length'][j])
                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse2(
                                ramsey_pulse_center2, self.config['ramsey_pulse_length'][j])
                            if self.config['spin_echo']:
                                awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse2(
                                    spin_echo_center, self.config['echo_pulse_length'][j])

                    if plotter is not None:
                        cur_pulse_seq.plot_pulses(plotter)


        #load pulses into the AGILENT AWG
        if flux_pulsing and self.config['load_agilent']:
            #trigger the flux AWG
            awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)

            awg_seq.load_full_into_agilent(self.config['AWG'][1])

        #load pulses into the awg (TEK)
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])


    #run the ramsey pulse sequence    
    def run_ramsey_sequential(self, num_avgs=1024, plotter=None):

        #define instrument manager    
        im = InstrumentManager()

        if self.config['LO'] == self.config['RF']:
            homodyne_meas = True
            print "Homodyne measurement"
        else:
            homodyne_meas = False

            #setup awg
        awg = im[self.config['AWG']]


        #setup drive,RF,LO

        drive = im[self.config['drive']]

        #setup drive
        if self.config['drive'][0:2] == "LB":
            drive.set_output(True)
            drive.set_pulse_ext(mod=True)
            #Internal pulse must be set to false!
            drive.set_mod(False)
            drive.set_power(self.config['ramsey_drive_pwr'])  #-28 when direct driving
            drive.set_frequency(self.config['ramsey_drive_freq'])
        else:
            drive.set_output(True)
            drive.set_ext_pulse()
            #Internal pulse must be set to false!
            drive.set_mod(True)
            drive.set_power(self.config['ramsey_drive_pwr'])  #-28 when direct driving
            drive.set_frequency(self.config['ramsey_drive_freq'])

        #        for i in range(5):
        #            if abs(drive.get_frequency()-self.config['ramsey_drive_freq'])>100:
        #                drive.set_frequency(self.config['ramsey_drive_freq'])
        #            else:
        #                pass
        #                #break
        #        if abs(drive.get_frequency()-self.config['ramsey_drive_freq'])>100:
        #            print "Lab brick frequency error!"
        #        time.sleep(.2)

        print self.config['ramsey_drive_freq'], drive.get_frequency()

        #setup LO and RF

        RF = im[self.config['RF']]
        RF.set_output(True)
        RF.set_ext_pulse()
        RF.set_mod(True)
        RF.set_power(self.config['ramsey_read_pwr'])
        RF.set_frequency(self.config['ramsey_read_freq'])

        if not homodyne_meas:
            LO = im[self.config['LO']]
            LO.set_output(True)
            LO.set_internal_pulse(40e-6)
            LO.set_mod(True)
            LO.set_power(self.config['ramsey_lo_power'])
            LO.set_frequency(self.config['ramsey_read_freq'] + self.config['IFreq'])

        #setup the flux
        flux = im[self.config['flux']]
        flux.set_mode('voltage')
        if self.config['flux'] == 'flux2' or self.config['flux'] == 'flux1':
            flux.set_function("DC")
            flux.set_offset(self.config['ramsey_flux'])
        else:
            flux.set_volt(self.config['ramsey_flux'])

        #setup the Alazar card
        card = setup_Alazar_card(512, self.config['acq_length'], self.config['meas_range'],
                                 self.config['monitor_pulses'], True, num_avgs)

        #setup plots
        if plotter is not None:
            plotter.init_plot("FullData", rank=2, accum=False)
            plotter.init_plot("Scope", accum=False)
            plotter.init_plot("Amp", accum=False)
            plotter.init_plot("Phase", accum=False)

        self.ramsey_amp = zeros(len(self.ramsey_time))
        self.ramsey_phase = zeros(len(self.ramsey_time))
        self.ramsey_fulldata = zeros((len(self.ramsey_time), self.config['acq_length']))

        for i in range(len(self.ramsey_time)):

            #program sequence
            awg.pre_experiment()
            self.load_ramsey(False, None, i)
            awg.prep_experiment(self.config['seq_file'])
            awg.run()

            tpts, ch1_pts, ch2_pts = card.acquire_avg_data()
            if not homodyne_meas:
                amp, phase, _, _ = heterodyne(tpts, ch1_pts, ch2_pts, self.config['IFreq'])
            else:
                amp = mean(ch1_pts)
                phase = amp

            self.ramsey_phase[i] = phase
            self.ramsey_amp[i] = amp
            self.ramsey_fulldata[i] = ch1_pts

            if plotter is not None:
                plotter.plot((self.ramsey_time[0:(i + 1)], self.ramsey_phase[0:(i + 1)]), "Phase")
                plotter.plot((self.ramsey_time[0:(i + 1)], self.ramsey_amp[0:(i + 1)]), "Amp")
                plotter.plot((range(self.config['acq_length']), ch1_pts), "Scope")
                plotter.plot(self.ramsey_fulldata, "FullData")


    def fit_ramsey(self, fit_exp=False):

        #start the fit after fit_start
        start_index = -1
        for i in range(len(self.ramsey_time)):
            if self.ramsey_time[i] >= self.config['fit_start']:
                start_index = i
                break

        if start_index == -1:
            start_index = 0

        #fit the t1 data to an decaying sinusoid 
        if self.config['fit_phase']:
            ydata = self.ramsey_phase[start_index:-1]
            xdata = self.ramsey_time[start_index:-1]

        else:
            ydata = self.ramsey_amp[start_index:-1]
            xdata = self.ramsey_time[start_index:-1]

        self.config['fit_start'] = xdata[0]

        if not fit_exp:
            fitvals = fitdecaysin(xdata, ydata)
            self.config['tau'] = fitvals[3]
            self.config['A'] = fitvals[0]
            self.config['y0'] = fitvals[4]
            self.config['omega'] = fitvals[1]
            self.config['phase'] = fitvals[2]
        else:
            fitvals = fitgeneral(xdata, ydata, lambda p, x: expfunc2(p, x, self.config['fit_start']),
                                 [ydata[-1], (ydata[0] - ydata[-1]), xdata[-1] / 2])
            self.config['tau'] = fitvals[2]
            self.config['A'] = fitvals[1]
            self.config['y0'] = fitvals[0]
            self.config['omega'] = 0.
            self.config['phase'] = 90.

        print fitvals

    def display_ramsey(self, plotter, pulse_sequence=False, display_fit=True):


        plotter.init_plot("FullData", rank=2, accum=False)
        plotter.init_plot("Amp", rank=1, accum=False)
        plotter.init_plot("Phase", rank=1, accum=False)

        plotter.plot((self.ramsey_fulldata), "FullData")

        p = [self.config['A'], self.config['omega'], self.config['phase'], self.config['tau'], self.config['y0'],
             self.config['fit_start']]

        if self.config['fit_phase'] and display_fit:
            plotter.plot((concatenate((self.ramsey_time, self.ramsey_time)),
                          concatenate((self.ramsey_phase, decaysin(p, self.ramsey_time)))), "Phase")
        else:
            plotter.plot((self.ramsey_time, self.ramsey_phase), "Phase")

        time.sleep(.2)

        if (not self.config['fit_phase']) and display_fit:
            plotter.plot((concatenate((self.ramsey_time, self.ramsey_time)),
                          concatenate((self.ramsey_amp, decaysin(p, self.ramsey_time)))), "Amp")
        else:
            plotter.plot((self.ramsey_time, self.ramsey_amp * 1.0), "Amp")


#Do a ramsey on 1 qubit and a rabi on the other qubit with optional flux pulses
class ramsey1rabi2(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        #vary parameter is the height of the rabi pulse!  
        #NOTE: the qubit_enable parameter DOES NOT APPLY!!! Always a two qubit experiment!

        self.config['exp_id'] = 'RAMSEY1RABI2'

        #Ramsey parameters for Qubit 1
        self.config['ramsey_height'] = 0.3
        self.config['ramsey_pulse_length'] = 10  #note: the rabi pulse will be the same length!
        self.config['ramsey_time'] = 100
        self.config['phase_advance_rate'] = 0.0


        #Properties if we are doing a flux pulse
        self.config['do_fast_flux'] = False  #Note if ramsey_vary is flux then this will also do the fast flux
        self.config['flux_height'] = [0.0, 0.0]

        self.config['flux_time'] = [20, 30]

        #shape of the flux pulse
        self.config['flux_up_time'] = [25, 25]  #time to the fast flux pulse height
        self.config['trap_flux_pulse'] = [False, False]
        self.config['trap_slope'] = [0.1, 0.1]


    #load the ramsey oscillation into the AWG
    def load_exp(self, plotter=None, ramsey_index=-1):


        if len(self.xdata) == 0:
            raise NameError("Sequence not initialized")

        #ALWAYS TWO QUBIT
        self.config['qubit_enable'] = [True, True]

        #ALWAYS Shaped pulse
        self.config['shaped_pulse'] = True


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #where to start the card and read pulses
        max_ramsey_wait = self.config['front_buffer'] + self.config['generator_enable_delay'] + self.config[
            'ramsey_time'] + 5 * self.config['ramsey_pulse_length']

        if ramsey_index == -1:
            vary_array = self.xdata
            pass
        else:
            vary_array = self.xdata[ramsey_index]
            print vary_array

        awg_clk_ratio = self.config['awg_clocks'][1] / self.config['awg_clocks'][0]

        awg_seq = pulse_sequence(len(vary_array), self.config['total_length'],
                                 int(ceil(self.config['total_length'] * awg_clk_ratio)))

        for i in range(2):
            if self.config['qubit_enable'][i]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][i]][0].square_pulse(
                    (max_ramsey_wait + self.config['read_delay']), self.config['read_length'])

        #trigger the flux AWG
        awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)

        #set up the card trigger
        if self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            if trigger_card_early:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (1 + self.config['generator_enable_delay']), self.config['read_length'])
            else:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (max_ramsey_wait + self.config['card_delay']), self.config['read_length'])

        flux_pulsing = False

        #setup qubit 1 (always the same!)
        #---------------------------------
        wait_time = self.config['ramsey_time']
        ramsey_pulse_center1 = max_ramsey_wait - wait_time
        ramsey_pulse_center2 = max_ramsey_wait

        #pulse 1
        awg_seq.analogwf[awg_seq.channel_index['drive_I'][0]][0].gauss_pulse(ramsey_pulse_center1,
                                                                             self.config['ramsey_pulse_length'] / 2,
                                                                             self.config['ramsey_pulse_height'], False)

        #pulse 2
        awg_seq.analogwf[awg_seq.channel_index['drive_I'][0]][0].gauss_pulse(ramsey_pulse_center2,
                                                                             self.config['ramsey_pulse_length'] / 2,
                                                                             cos(self.config[
                                                                                     'phase_advance_rate'] * wait_time) *
                                                                             self.config['ramsey_pulse_height'], False)
        awg_seq.analogwf[awg_seq.channel_index['drive_Q'][0]][0].gauss_pulse(ramsey_pulse_center2,
                                                                             self.config['ramsey_pulse_length'] / 2,
                                                                             sin(self.config[
                                                                                     'phase_advance_rate'] * wait_time) *
                                                                             self.config['ramsey_pulse_height'], False)

        #keep the drive switch on for the whole time 
        awg_seq.marker[awg_seq.channel_index['drive_trig'][0]][0].square_pulse2(
            (ramsey_pulse_center2 + ramsey_pulse_center1) / 2, (
            ramsey_pulse_center2 - ramsey_pulse_center1 + 3 * self.config['ramsey_pulse_length'] + 2 * self.config[
                'generator_enable_delay']))

        #set the fast flux pluse
        for j in range(2):

            if self.config['do_fast_flux']:

                flux_pulsing = True
                flux_val = self.config['flux_height'][j]

                indflux = 0

                if flux_val != 0:

                    flux_time = self.config['flux_time'][j] * awg_clk_ratio
                    flux_cntr = (ramsey_pulse_center1 + ramsey_pulse_center2) / 2.0 * awg_clk_ratio - \
                                self.config['awg_delay'][j]

                    if self.config['trap_flux_pulse'][j]:
                        awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].trapezoid(flux_cntr,
                                                                                              self.config['trap_slope'][
                                                                                                  j] / awg_clk_ratio,
                                                                                              flux_time, flux_val)

                        #do negative pulse after experiment is done
                        awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].trapezoid((self.config[
                                                                                                   'total_length'] - 2 * (
                                                                                               3 * self.config[
                                                                                                   'flux_up_time'][
                                                                                                   j] + flux_time / awg_clk_ratio)) * awg_clk_ratio,
                                                                                              -
                                                                                              self.config['trap_slope'][
                                                                                                  j] / awg_clk_ratio,
                                                                                              flux_time, -flux_val)
                    else:
                        awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].smooth_square(flux_cntr,
                                                                                                  self.config[
                                                                                                      'flux_up_time'][
                                                                                                      j] * awg_clk_ratio,
                                                                                                  flux_time, flux_val)

                        #do negative pulse after experiment is done
                        awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].smooth_square((self.config[
                                                                                                       'total_length'] - 2 * (
                                                                                                   3 * self.config[
                                                                                                       'flux_up_time'][
                                                                                                       j] + flux_time / awg_clk_ratio)) * awg_clk_ratio,
                                                                                                  self.config[
                                                                                                      'flux_up_time'][
                                                                                                      j] * awg_clk_ratio,
                                                                                                  flux_time, -flux_val)


                        #setup qubit 2 Rabi oscillations
        for i in range(len(vary_array)):


            #setup the Rabi pulse for qubit 2
            ind1a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][1])
            awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][1]][i] = ind1a

            ind2 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][1])
            awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][1]][i] = ind2

            #Rabi Pulse
            awg_seq.analogwf[awg_seq.channel_index['drive_I'][1]][ind1a].gauss_pulse(ramsey_pulse_center1, self.config[
                'ramsey_pulse_length'] / 2, vary_array[i], False)

            #Trigger
            awg_seq.marker[awg_seq.channel_index['drive_trig'][1]][ind2].square_pulse2(ramsey_pulse_center1,
                                                                                       3 * self.config[
                                                                                           'ramsey_pulse_length'] + 2 *
                                                                                       self.config[
                                                                                           'generator_enable_delay'])

            if plotter is not None:
                cur_pulse_seq.plot_pulses(plotter)

        #load pulses into the awg (TEK)
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #load pulses into the AGILENT AWG
        if flux_pulsing:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])


#Do a pulsed flux experiment with 2 qubits
#This class lets one vary a number of the gate parameters
class qubit_fastflux(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        #vary parameter is either the pedestal height of qubit 1, the length of 
        #the flux pulse for qubit 2 or the pedestal height for qubit 2
        self.config['fastflux_vary'] = "pedestal1"  #pedestal1, length2, pedestal2
        #NOTE: the qubit_enable parameter DOES NOT APPLY!!! 
        #Always a two qubit experiment!

        #do a Ramsey on qubit if 'fastflux_vary' is the second qubit length
        self.config['qubit2_ramsey'] = False

        self.config['exp_id'] = 'QUBITFASTFLUX'

        #total time between pulses
        self.config['pulse_wait'] = 150

        #pulse 1 parameters
        self.config['first_pulse_height'] = [0.3, 0.3]
        self.config['first_pulse_length'] = [10, 10]
        self.config['first_pulse_phase'] = [0.0, 0.0]

        #pulse 2 parameters
        self.config['second_pulse_height'] = [0.3, 0.3]
        self.config['second_pulse_length'] = [10, 10]
        self.config['second_pulse_phase'] = [0.0, 0.0]


        #Properties if we are doing a flux pulse
        self.config['do_fast_flux'] = False
        self.config['flux_height'] = [0.0, 0.0]
        #flux pulse duration (total up plus wait time)
        self.config['trap_time'] = [20, 20]



        #shape of the flux pulse
        self.config['flux_pulse_type'] = [2, 2]  #1: smooth square, 2: trapezoid, 3: tan
        self.config['flux_up_time'] = [25, 25]  #time to the fast flux pulse height

        #parameters for a trapezoidal flux pulse
        self.config['trap_slope_up'] = [0.1, 0.1]
        self.config['trap_slope_down'] = [0.1, 0.1]

        #parameters for the tan pulse
        self.config['tan_pulse_up_duration'] = [[5, 5], [5, 5]]
        self.config['tan_pulse_slopes'] = [[0.9, 0.9], [0.9, 0.9]]
        self.config['tan_pulse_midpt'] = [-0.4, -0.4]

        #pedestal before the flux pulse
        self.config['pedestal_height'] = [0, 0]
        self.config['pedestal_time'] = [20, 30]

        self.config['comp_pulse_height'] = [[0, 0, 0], [0, 0, 0]]
        self.config['comp_pulse_length'] = [[100, 1000], [100, 1000]]

        self.config['num_calib_seqs'] = 20

        #tomography options
        #note that both qubit 1 and 2 have the same length pulse
        self.config['tomo_pulse_length'] = 15

        #first list is the X(I),Y(Q) pulse for qubit 1...second list for qubit 2
        self.config['tomo_pulse_height'] = [[0.4, 0.4], [0.4, 0.4]]
        self.config['tomo_pulse_height_neg'] = [[0.4, 0.4], [0.4,
                                                             0.4]]  #pulse heights for negative rotations (NOTE: these are the absolute values!!!)     
        self.config['tomo_pi_pulse'] = [0.2, 0.2]  #these are pi pulses on *I* for qubit 1 and qubit 2
        self.config['tomo_wait'] = 20

        #mixer calibration
        #self.config['q_phase'] = [0.172,0]
        self.config['q_phase'] = [-0.172 * 0, 0]


    #load the ramsey oscillation into the AWG
    def load_exp(self, plotter=None, ramsey_index=-1):


        if len(self.xdata) == 0:
            raise NameError("Sequence not initialized")

        #ALWAYS TWO QUBIT
        self.config['qubit_enable'] = [True, True]

        #ALWAYS Shaped pulse
        self.config['shaped_pulse'] = True


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #where to start the card and read pulses
        max_ramsey_wait = self.config['front_buffer'] + self.config['generator_enable_delay'] + self.config[
            'pulse_wait'] + 5 * max(self.config['first_pulse_length']) + 2 * self.config['tomo_pulse_length']

        if ramsey_index == -1:
            vary_array = self.xdata
        else:
            vary_array = self.xdata[ramsey_index]
            print vary_array

        awg_seq = pulse_sequence(len(vary_array) + self.config['num_calib_seqs'],
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][0])),
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][1])),
                                 self.config['awg_clocks'][0], self.config['awg_clocks'][1], self.config['awg_delay'])

        for i in range(2):
            if self.config['qubit_enable'][i]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][i]][0].square_pulse(
                    (max_ramsey_wait + self.config['read_delay']), self.config['read_length'])

        #trigger the flux AWG
        awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)

        #set up the card trigger
        if self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            if trigger_card_early:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (1 + self.config['generator_enable_delay']), self.config['read_length'])
            else:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (max_ramsey_wait + self.config['card_delay']), self.config['read_length'])

        flux_pulsing = False

        #setup qubit 1 (always the same!)
        #---------------------------------
        wait_time = self.config['pulse_wait']
        tomo_wait = self.config['tomo_wait']

        pulse_center3 = max_ramsey_wait
        pulse_center2 = max_ramsey_wait - (
        max(self.config['second_pulse_length']) + self.config['tomo_pulse_length'] + tomo_wait)  #last gate pulse
        pulse_center1 = pulse_center2 - wait_time



        #set the fast flux pluse
        for i in range(len(vary_array) + self.config['num_calib_seqs']):

            for j in range(2):

                if self.config['do_fast_flux']:

                    indflux = awg_seq.add_analog_pulse(awg_seq.channel_index['flux'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['flux'][j]][i] = indflux

                    flux_pulsing = True
                    flux_val = self.config['flux_height'][j]

                    #indflux = 0

                    if flux_val != 0 and (i < len(vary_array) or i > (len(vary_array) + 6)):

                        pedestal_time = self.config['pedestal_time'][j]
                        trap_time = self.config['trap_time'][j]
                        pedestal_height = self.config['pedestal_height'][j]

                        if i < len(vary_array):

                            if self.config['fastflux_vary'] == "pedestal1" and j == 0:
                                pedestal_height = vary_array[i]
                            elif self.config['fastflux_vary'] == "pedestal2" and j == 1:
                                pedestal_time = vary_array[i]
                            elif self.config['fastflux_vary'] == "length2" and j == 1:
                                trap_time = vary_array[i]

                        flux_cntr = (pulse_center1 + pulse_center2) / 2.0

                        self.__flux_pulses__(awg_seq, indflux, flux_cntr, flux_val, trap_time, pedestal_time,
                                             pedestal_height, j, max_ramsey_wait)



                        #setup pulses
        for i in range(len(vary_array) + self.config['num_calib_seqs']):

            for j in range(2):


                pulse1_height = self.config['first_pulse_height'][j]
                pulse2_height = self.config['second_pulse_height'][j]

                pulse1_len = self.config['first_pulse_length'][j]
                pulse2_len = self.config['second_pulse_length'][j]

                if self.config['fastflux_vary'] == "pulse_height2" and j == 1 and i < len(vary_array):
                    pulse1_height = vary_array[i]
                elif self.config['fastflux_vary'] == "pulse_height1" and j == 0 and i < len(vary_array):
                    pulse1_height = vary_array[i]
                    pulse2_height = vary_array[i]

                if (self.config['fastflux_vary'] == "length2" and j == 1 and self.config['qubit2_ramsey'] and i < len(
                        vary_array)):
                    pulse_center1a = flux_cntr - vary_array[i] / 2 - pulse1_len - 10.0
                    pulse_center2a = flux_cntr + vary_array[i] / 2 + pulse2_len + 10.0
                else:
                    pulse_center1a = pulse_center1
                    pulse_center2a = pulse_center2

                num_test_pulses1 = 6
                num_test_pulses2 = 13

                if i > len(vary_array):
                    pulse1_height = 0

                if i > (len(vary_array)) and i < (len(vary_array) + num_test_pulses1 + 1):
                    pulse2_height *= (i - len(vary_array) - 1) / (num_test_pulses1 - 1) * 2.5
                elif i > (len(vary_array)):
                    pulse2_height *= (i - len(vary_array) - num_test_pulses1 - 1) / (num_test_pulses2 - 1) * 3.0

                self.__add_pulses__(awg_seq, pulse_center1a, pulse_center2a, pulse1_height, pulse2_height, pulse1_len,
                                    pulse2_len, i, j)

                if plotter is not None:
                    cur_pulse_seq.plot_pulses(plotter)


        #load pulses into the awg (TEK)
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])


        #load pulses into the AGILENT AWG
        if flux_pulsing and self.config['load_agilent']:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])

    def __add_pulses__(self, awg_seq, pulse_center1, pulse_center2, pulse1_height, pulse2_height, pulse1_len,
                       pulse2_len, i, j):

        #Add the qubit pulses to the awg
        ind1a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
        awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind1a

        ind1b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
        awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind1b

        ind2 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind2

        phase1 = self.config['first_pulse_phase'][j]
        phase2 = self.config['second_pulse_phase'][j]


        #ind1a = 0
        #ind1b = 0
        #ind2 = 0     

        if self.config['drag_pulse'][j]:

            awg_seq.DRAG_pulse(pulse_center1, pulse1_len / 2.0, pulse1_height, phase1, j, ind1a, ind1b,
                               self.config['drag_prefactor'][j], False)

            awg_seq.DRAG_pulse(pulse_center2, pulse2_len / 2.0, pulse2_height, phase2, j, ind1a, ind1b,
                               self.config['drag_prefactor'][j], False)

        else:

            #pulse 1
            awg_seq.gauss_pulse(pulse_center1, pulse1_len / 2.0, pulse1_height, phase1, j, ind1a, ind1b, False)
            #awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse(pulse_center1, pulse1_len/2, (cos(phase1)+sin(phase1)*tan(self.config['q_phase'][j]))*pulse1_height, False)
            #awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].gauss_pulse(pulse_center1, pulse1_len/2, sin(phase1)*pulse1_height, False)

            #pulse 2
            awg_seq.gauss_pulse(pulse_center2, pulse2_len / 2.0, pulse2_height, phase2, j, ind1a, ind1b, False)
            #awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse(pulse_center2, pulse2_len/2, (cos(phase2)+sin(phase2)*tan(self.config['q_phase'][j]))*pulse2_height, False)
            #awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].gauss_pulse(pulse_center2, pulse2_len/2, sin(phase2)*pulse2_height, False)

        #add offsets
        awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].offset(
            self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
        awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].offset(
            self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])


        #keep the drive switch on for the whole time 
        #awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse2((pulse_center2+pulse_center1)/2, (pulse_center2-pulse_center1+3*self.config['first_pulse_length'][j] + 2*self.config['generator_enable_delay'])) 
        if self.config['pulsed_drive'][j]:
            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse(
                pulse_center1 - self.config['generator_enable_delay'] - self.config['first_pulse_length'][j],
                pulse_center2 - pulse_center1 + self.config['first_pulse_length'][j] +
                self.config['second_pulse_length'][j] + 1.2 * self.config['generator_enable_delay'])

        return ind1a, ind1b, ind2

    def __flux_pulses__(self, awg_seq, indflux, flux_cntr, flux_val, trap_time, pedestal_time, pedestal_height, j,
                        max_ramsey_wait, comp_pulse_scaling=1.0):

        flux_cntr2 = max_ramsey_wait + 3000

        if self.config['flux_pulse_type'][j] == 2:  #trapezoidal pulse

            awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].trapezoid_w_pedestal(flux_cntr, self.config[
                'trap_slope_up'][j], self.config['trap_slope_down'][j], trap_time + 2 * pedestal_time, trap_time,
                                                                                             pedestal_height, flux_val,
                                                                                             flux_val)

            #do negative pulse after experiment is done

            awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].trapezoid_w_pedestal(flux_cntr2, -
            self.config['trap_slope_up'][j], -self.config['trap_slope_down'][j], trap_time + 2 * pedestal_time,
                                                                                             trap_time,
                                                                                             -pedestal_height,
                                                                                             -flux_val, -flux_val)

        elif self.config['flux_pulse_type'][j] == 3:  #tan pulse

            print self.config['tan_pulse_up_duration'][j][0]

            awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].tan_pulse_w_pedestal(flux_cntr, self.config[
                'tan_pulse_up_duration'][j][0], self.config['tan_pulse_up_duration'][j][1], trap_time, self.config[
                                                                                                 'tan_pulse_midpt'][j],
                                                                                             flux_val, pedestal_height,
                                                                                             pedestal_time, self.config[
                    'tan_pulse_slopes'][j][0], self.config['tan_pulse_slopes'][j][1])

            awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].tan_pulse_w_pedestal(flux_cntr2, self.config[
                'tan_pulse_up_duration'][j][0], self.config['tan_pulse_up_duration'][j][1], trap_time, -self.config[
                'tan_pulse_midpt'][j], -flux_val, -pedestal_height, pedestal_time, self.config['tan_pulse_slopes'][j][
                                                                                                 0], self.config[
                                                                                                 'tan_pulse_slopes'][j][
                                                                                                 1])


        elif self.config['flux_pulse_type'][j] == 1:  #smooth square

            awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].smooth_square(flux_cntr,
                                                                                      self.config['flux_up_time'][j],
                                                                                      trap_time, flux_val)

            #do negative pulse after experiment is done
            awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].smooth_square(
                (self.config['total_length'] - 2 * (3 * self.config['flux_up_time'][j] + flux_time)),
                self.config['flux_up_time'][j], flux_time, -flux_val)

        else:

            raise NameError('Invalid flux pulse type')

        #positive and negative compensation pulses

        start_pt = 0
        for i in range(len(self.config['comp_pulse_length'][j])):
            awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].line(
                flux_cntr + (trap_time + 2 * pedestal_time) / 2 + start_pt, self.config['comp_pulse_length'][j][i],
                self.config['comp_pulse_height'][j][i] * comp_pulse_scaling,
                self.config['comp_pulse_height'][j][i + 1] * comp_pulse_scaling)

            awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].line(
                flux_cntr2 + (trap_time + 2 * pedestal_time) / 2 + start_pt, self.config['comp_pulse_length'][j][i],
                -self.config['comp_pulse_height'][j][i] * comp_pulse_scaling,
                -self.config['comp_pulse_height'][j][i + 1] * comp_pulse_scaling)

            start_pt += self.config['comp_pulse_length'][j][i]


#Do a pulsed flux experiment and then do qubit spectroscopy
class qubit_fastflux_spectroscopy(qubit_fastflux):
    def __init__(self):

        qubit_fastflux.__init__(self)

        self.config['num_calib_seqs'] = 0

        #set the drive frequencies to run over
        self.config['drive_frequencies'] = [[0, 0], [0, 0]]

        self.config['probe_pulse_time'] = 0  #if varying the flux, this is when to probe

        self.config['flux_vary'] = 'time'  #time or flux1, flux2


    #override the qubit_exp
    def run_qubit_exp(self, plotter=None):

        if self.config['num_avgs'] == -1:
            raise NameError('Must run a finite number of averages')

        temp_ydata = zeros((2, len(self.config['drive_frequencies'][0]), len(self.xdata)))

        num_avgs = self.config['num_avgs']

        if plotter is not None:
            for i in range(2):
                plotter.init_plot("FullSpectData" + str(i), rank=2, accum=False)

        for i in range(len(self.config['drive_frequencies'][0])):


            #set the drive frequencies            
            self.config['drive_freq'] = [self.config['drive_frequencies'][0][i], self.config['drive_frequencies'][1][i]]

            print self.config['drive_freq']

            #run!
            qubit_fastflux.run_qubit_exp(self, plotter)

            if self.config['num_avgs'] != num_avgs:
                #was interrupted prematurely
                print("Fast flux spectroscopy interrupted")
                break

            for j in range(2):
                temp_ydata[j, i] = self.ydata_avg[j]

            if plotter is not None:
                for j in range(2):
                    plotter.plot(temp_ydata[j], "FullSpectData" + str(j))




        #save the spectroscopy data in the ydata variable
        self.ydata = temp_ydata

    #load the ramsey oscillation into the AWG
    def load_exp(self, plotter=None, ramsey_index=-1):


        if len(self.xdata) == 0:
            raise NameError("Sequence not initialized")

        #ALWAYS TWO QUBIT
        self.config['qubit_enable'] = [True, True]

        #ALWAYS Shaped pulse
        self.config['shaped_pulse'] = True


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #where to start the card and read pulses
        max_ramsey_wait = self.config['front_buffer'] + self.config['generator_enable_delay'] + self.config[
            'pulse_wait'] + 5 * max(self.config['first_pulse_length']) + self.xdata[-1]


        #xdata is the spacing from the flux pulse to the 
        #spectroscopy pulse

        if ramsey_index == -1:
            vary_array = self.xdata
        else:
            vary_array = self.xdata[ramsey_index]
            print vary_array

        awg_seq = pulse_sequence(len(vary_array) + self.config['num_calib_seqs'],
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][0])),
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][1])),
                                 self.config['awg_clocks'][0], self.config['awg_clocks'][1], self.config['awg_delay'])

        for i in range(2):
            if self.config['qubit_enable'][i]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][i]][0].square_pulse(
                    (max_ramsey_wait + self.config['read_delay']), self.config['read_length'])

        #trigger the flux AWG
        awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)

        #set up the card trigger
        if self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            if trigger_card_early:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (1 + self.config['generator_enable_delay']), self.config['read_length'])
            else:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (max_ramsey_wait + self.config['card_delay']), self.config['read_length'])

        flux_pulsing = False

        #setup qubit 1 (always the same!)
        #---------------------------------
        wait_time = self.config['pulse_wait']
        pulse_center1 = max_ramsey_wait - wait_time - self.xdata[-1]
        pulse_center2 = max_ramsey_wait - self.xdata[-1]


        #set the fast flux pulse

        for i in range(len(vary_array)):

            for j in range(2):

                if self.config['do_fast_flux']:

                    if self.config['flux_vary'] != 'time':
                        indflux = awg_seq.add_analog_pulse(awg_seq.channel_index['flux'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['flux'][j]][i] = indflux
                    else:
                        indflux = 0
                        if i > 0:
                            break

                    flux_pulsing = True
                    flux_val = self.config['flux_height'][j]

                    #indflux = 0

                    if flux_val != 0:

                        if self.config['flux_vary'] == 'flux1' and j == 0:
                            flux_val = self.xdata[i]
                        elif self.config['flux_vary'] == 'flux2' and j == 1:
                            flux_val = self.xdata[i]

                        print flux_val

                        pedestal_time = self.config['pedestal_time'][j]
                        trap_time = self.config['trap_time'][j]
                        pedestal_height = self.config['pedestal_height'][j]

                        flux_cntr = (pulse_center1 + pulse_center2) / 2.0

                        self.__flux_pulses__(awg_seq, indflux, flux_cntr, flux_val, trap_time, pedestal_time,
                                             pedestal_height, j, max_ramsey_wait)



                        #setup pulses
        for i in range(len(vary_array)):

            for j in range(2):

                #note: there is no first pulse!
                pulse1_height = self.config['first_pulse_height'][j] * 0
                pulse2_height = self.config['second_pulse_height'][j]

                pulse1_len = self.config['first_pulse_length'][j]
                pulse2_len = self.config['second_pulse_length'][j]

                if self.config['flux_vary'] == 'time':
                    pulse_center2b = pulse_center2 + self.xdata[i]
                else:
                    pulse_center2b = pulse_center2 + self.config['probe_pulse_time']

                self.__add_pulses__(awg_seq, pulse_center1, pulse_center2b, pulse1_height, pulse2_height, pulse1_len,
                                    pulse2_len, i, j)

                if plotter is not None:
                    cur_pulse_seq.plot_pulses(plotter)


        #load pulses into the awg (TEK)
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #load pulses into the AGILENT AWG
        if flux_pulsing and self.config['load_agilent']:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])


#Do state tomography on the gate
#Inherit the qubit_fastflux class to get the gate options
#note that there is some code duplication...eventually I would like to grab
class gate_w_statetomography(qubit_fastflux):
    def __init__(self):

        #note that some of the qubit_fastflux options will not apply
        qubit_fastflux.__init__(self)

        self.config['tomo_vary'] = "calib1"  #calib1,calib2,tomo,basis
        #calib1 varies the X (I) tomography pulses
        #calib2 varies the Y (Q) tomography pulses
        #tomo does the tomography
        #basis takes just the basis states |gg>,|eg>,|ge>,|ee>

        #NOTE: the qubit_enable parameter DOES NOT APPLY!!! 
        #Always a two qubit experiment!

        self.config['exp_id'] = 'GATE_W_STATETOMO'

        self.config['num_calib_seqs'] = 0

        #do an expanded tomography set (i.e. including -pi/2 rotation)
        self.config['expanded_tomo'] = False

    #load the ramsey oscillation into the AWG
    def load_exp(self, plotter=None, ramsey_index=-1):


        if len(self.xdata) == 0 and self.config['tomo_vary'][0:5] == "calib":
            raise NameError("Calibration sequence not initialized")

        tomo_index = [[[0, 0], [0, 0]],  #Z Z
                      [[1, 0], [0, 0]],  #X Z
                      [[0, 1], [0, 0]],  #Y Z
                      [[0, 0], [1, 0]],  #Z X
                      [[0, 0], [0, 1]],  #Z Y
                      [[1, 0], [1, 0]],  #X X
                      [[1, 0], [0, 1]],  #X Y
                      [[0, 1], [1, 0]],  #Y X
                      [[0, 1], [0, 1]]]  #Y Y

        if self.config['expanded_tomo']:
            tomo_index.append([[-1, 0], [0, 0]])  #-X Z
            tomo_index.append([[0, -1], [0, 0]])  #-Y Z
            tomo_index.append([[0, 0], [-1, 0]])  #Z -X
            tomo_index.append([[0, 0], [0, -1]])  #Z -Y
            tomo_index.append([[-1, 0], [-1, 0]])  #-X -X
            tomo_index.append([[-1, 0], [0, -1]])  #-X -Y
            tomo_index.append([[0, -1], [-1, 0]])  #-Y -X
            tomo_index.append([[0, -1], [0, -1]])  #-Y -Y

        #only do the basis set pulses
        if self.config['tomo_vary'] == 'basis':
            tomo_index = []

        tomo_index.append([[0, 0], [0, 0]])  #|gg>
        tomo_index.append(
            [[self.config['tomo_pi_pulse'][0] / self.config['tomo_pulse_height'][0][0], 0], [0, 0]])  #|eg>
        tomo_index.append(
            [[0, 0], [self.config['tomo_pi_pulse'][1] / self.config['tomo_pulse_height'][1][0], 0]])  #|ge>
        tomo_index.append([[self.config['tomo_pi_pulse'][0] / self.config['tomo_pulse_height'][0][0], 0],
                           [self.config['tomo_pi_pulse'][1] / self.config['tomo_pulse_height'][1][0], 0]])  #|ee>  

        num_tomo_pulses = len(tomo_index) - 4

        #For the actual tomography there are 4 measurement calibrations plus
        #9 tomography pulses
        if self.config['tomo_vary'] == 'tomo' or self.config['tomo_vary'] == 'basis':
            self.xdata = range(len(tomo_index))

            #ALWAYS TWO QUBIT
        self.config['qubit_enable'] = [True, True]

        #ALWAYS Shaped pulse
        self.config['shaped_pulse'] = True


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #where to start the card and read pulses
        max_ramsey_wait = self.config['front_buffer'] + self.config['generator_enable_delay'] + self.config[
            'pulse_wait'] + 5 * max(self.config['first_pulse_length']) + 2 * self.config['tomo_pulse_length']

        if ramsey_index == -1:
            vary_array = self.xdata
            if self.config['tomo_vary'] == 'calib3c':
                vary_array = zeros(3 * len(self.xdata))
                for i in range(len(self.xdata)):
                    vary_array[(3 * i):(3 * i + 3)] = self.xdata[i]
                self.xdata = vary_array.copy()
                for i in range(int(len(self.xdata) / 3)):
                    self.xdata[3 * i + 1] += (vary_array[3] - vary_array[0]) / 3.0
                    self.xdata[3 * i + 2] += 2 * (vary_array[3] - vary_array[0]) / 3.0
        else:
            vary_array = self.xdata[ramsey_index]
            print vary_array

        awg_seq = pulse_sequence(len(vary_array), int(ceil(self.config['total_length'] * self.config['awg_clocks'][0])),
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][1])),
                                 self.config['awg_clocks'][0], self.config['awg_clocks'][1], self.config['awg_delay'])

        for i in range(2):
            if self.config['qubit_enable'][i]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][i]][0].square_pulse(
                    (max_ramsey_wait + self.config['read_delay']), self.config['read_length'])

        #trigger the flux AWG
        awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)

        #set up the card trigger
        if self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            if trigger_card_early:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (1 + self.config['generator_enable_delay']), self.config['read_length'])
            else:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (max_ramsey_wait + self.config['card_delay']), self.config['read_length'])

        flux_pulsing = False

        #setup qubit 1 (always the same!)
        #---------------------------------
        tomo_wait = self.config['tomo_wait']
        wait_time = self.config['pulse_wait']
        pulse_center3 = max_ramsey_wait  #tomography pulse
        pulse_center2 = max_ramsey_wait - (
        max(self.config['second_pulse_length']) + self.config['tomo_pulse_length'] + tomo_wait)  #last gate pulse
        pulse_center1 = pulse_center2 - wait_time

        for i in range(len(vary_array)):

            for j in range(2):

                #set the fast flux pluse
                #------------
                if self.config['do_fast_flux']:

                    if i == 0:

                        #indflux = awg_seq.add_analog_pulse(awg_seq.channel_index['flux'][j])
                        #awg_seq.analog_seq_table[awg_seq.channel_index['flux'][j]][i] = indflux
                        indflux = 0

                        flux_pulsing = True
                        flux_val = self.config['flux_height'][j]

                        #indflux = 0

                        if flux_val != 0:


                            pedestal_time = self.config['pedestal_time'][j]
                            trap_time = self.config['trap_time'][j]
                            pedestal_height = self.config['pedestal_height'][j]

                            if self.config['tomo_vary'] == 'calib4b':
                                flux_val = 0
                                pedestal_height = 0

                            print j, flux_val, pedestal_height

                            flux_cntr = (pulse_center1 + pulse_center2) / 2.0

                            self.__flux_pulses__(awg_seq, indflux, flux_cntr, flux_val, trap_time, pedestal_time,
                                                 pedestal_height, j, max_ramsey_wait)

                #------------

                #Do the gate pulses!!
                #-------------
                pulse1_height = self.config['first_pulse_height'][j]
                pulse2_height = self.config['second_pulse_height'][j]

                pulse1_len = self.config['first_pulse_length'][j]
                pulse2_len = self.config['second_pulse_length'][j]

                pulse_center2b = pulse_center2

                #calibration pulses of the 4 measurement states
                if (self.config['tomo_vary'] == 'tomo' or self.config['tomo_vary'] == 'basis') and i >= (
                num_tomo_pulses):
                    pulse1_height = 0
                    pulse2_height = 0

                #if calibrating the tomography pulses then don't do these pulses
                if self.config['tomo_vary'][0:5] == 'calib':
                    pulse1_height = 0
                    pulse2_height = 0

                if self.config['tomo_vary'][0:6] == 'calib3':
                    pulse2_height = vary_array[i]
                    if self.config['tomo_vary'] == 'calib3b':
                        self.config['second_pulse_phase'][j] = pi / 2
                    elif self.config['tomo_vary'] == 'calib3c':
                        pulse_center2b = pulse_center2 + vary_array[i]
                        if mod(i, 3) == 0:
                            pulse2_height = 0
                        elif mod(i, 3) == 1:
                            pulse2_height = self.config['tomo_pulse_height'][j][0]
                            pulse2_height = self.config['tomo_pi_pulse'][j]
                        else:
                            pulse2_height = -self.config['tomo_pulse_height_neg'][j][0]
                            pulse2_height = -self.config['tomo_pi_pulse'][j]

                if self.config['tomo_vary'][0:6] == 'calib4':
                    pulse1_height = vary_array[i]
                    if self.config['tomo_vary'] == 'calib4c':
                        self.config['first_pulse_phase'][j] = pi / 2

                ind1a, ind1b, ind2 = self.__add_pulses__(awg_seq, pulse_center1, pulse_center2b, pulse1_height,
                                                         pulse2_height, pulse1_len, pulse2_len, i, j)

                #-------------  

                #Do the tomography pulses!!
                #-------------


                tomo_pulse_len = self.config['tomo_pulse_length']

                tomo_pulseI_height = 0
                tomo_pulseQ_height = 0

                if self.config['tomo_vary'] == 'calib1':
                    tomo_pulseI_height = vary_array[i]


                elif self.config['tomo_vary'] == 'calib2':
                    tomo_pulseQ_height = vary_array[i]

                elif self.config['tomo_vary'] == 'tomo' or self.config['tomo_vary'] == 'basis':

                    if tomo_index[i][j][0] < 0:
                        tomo_pulseI_height = self.config['tomo_pulse_height_neg'][j][0]
                    else:
                        tomo_pulseI_height = self.config['tomo_pulse_height'][j][0]

                    if tomo_index[i][j][1] < 0:
                        tomo_pulseQ_height = self.config['tomo_pulse_height_neg'][j][1]
                    else:
                        tomo_pulseQ_height = self.config['tomo_pulse_height'][j][1]

                    tomo_pulseI_height *= tomo_index[i][j][0]
                    tomo_pulseQ_height *= tomo_index[i][j][1]


                #tomography pulse
                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse(pulse_center3,
                                                                                         tomo_pulse_len / 2,
                                                                                         tomo_pulseI_height + tan(
                                                                                             self.config['q_phase'][
                                                                                                 j]) * tomo_pulseQ_height,
                                                                                         False)
                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].gauss_pulse(pulse_center3,
                                                                                         tomo_pulse_len / 2,
                                                                                         tomo_pulseQ_height, False)

                #keep the drive switch on for the whole time 
                drive_trig_start = pulse_center1 - self.config['generator_enable_delay'] - \
                                   self.config['first_pulse_length'][j]
                drive_trig_length = 1 * (pulse_center3 - pulse_center1) + tomo_pulse_len + \
                                    self.config['second_pulse_length'][j] + 1.2 * self.config[
                    'generator_enable_delay'] + max(abs(array(self.config['awg_delay'][0:4])))
                #drive_trig_start = 0
                #drive_trig_length = 8e3        
                awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].clear()
                if self.config['pulsed_drive'][j]:
                    awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse(drive_trig_start,
                                                                                              drive_trig_length)



                    #-------------  

                if plotter is not None:
                    cur_pulse_seq.plot_pulses(plotter)


        #load pulses into the awg (TEK)
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #load pulses into the AGILENT AWG
        if flux_pulsing and self.config['load_agilent']:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])


#Filter tests
#similar to the gate, but some differences
#no tomography
class filter_tests(qubit_fastflux):
    def __init__(self):

        #note that some of the qubit_fastflux options will not apply
        qubit_fastflux.__init__(self)

        self.config['filter_vary'] = "adiabatic"  #adiabatic1, adiabatic2, adiabatic_calib, adiabatic_calib2, lifetime
        #adiabatic1: varies the flux pulse slope and does spectroscopy on qubit 2
        #adiabatic2: varies the flux pulse slop without raising qubit 2
        #adiabatic_calib: does a rabi oscillation after the flux pulse from adiabatic2 for calibrating |g> and |e>
        #adiabatic_calib2: does a rabi oscillation before the flux pulse to calibrate pi pulse
        #lifetime varies the time in the filter

        #NOTE: the qubit_enable parameter DOES NOT APPLY!!! 
        #Always a two qubit experiment!

        self.config['comp_scaling'] = []

        self.config['exp_id'] = 'FILTER_TEST'

        self.config['num_calib_seqs'] = 2

        #this is the delay between the qubit 1 flux ramp and the qubit 2 spectroscopy pulse
        self.config['pulse2_delay'] = 20

        #this is for the fitler adiabaticity/lifetime tests
        #if true then include calibrations of the population
        self.config['include_calibs'] = True


    #load the ramsey oscillation into the AWG
    def load_exp(self, plotter=None, ramsey_index=-1):


        if len(self.xdata) == 0:
            raise NameError("Calibration sequence not initialized")


        #ALWAYS TWO QUBIT
        self.config['qubit_enable'] = [True, True]

        #ALWAYS FAST FLUX
        self.config['do_fast_flux'] = [True, True]

        #ALWAYS Shaped pulse
        self.config['shaped_pulse'] = True

        if self.config['include_calibs'] and (
                self.config['filter_vary'] == 'adiabatic2' or self.config['filter_vary'] == 'lifetime'):
            self.config['num_calib_seqs'] = 2 * len(self.xdata)


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #where to start the card and read pulses
        max_ramsey_wait = self.config['front_buffer'] + self.config['generator_enable_delay'] + self.config[
            'pulse_wait'] + 5 * max(self.config['first_pulse_length'])

        if self.config['filter_vary'] == 'lifetime':
            max_ramsey_wait += self.xdata[-1]
            #self.config['num_calib_seqs'] = 0

        if ramsey_index == -1:
            vary_array = self.xdata
        else:
            vary_array = self.xdata[ramsey_index]
            print vary_array

        awg_seq = pulse_sequence(len(vary_array) + self.config['num_calib_seqs'],
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][0])),
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][1])),
                                 self.config['awg_clocks'][0], self.config['awg_clocks'][1], self.config['awg_delay'])

        for i in range(2):
            if self.config['qubit_enable'][i]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][i]][0].square_pulse(
                    (max_ramsey_wait + self.config['read_delay']), self.config['read_length'])

        #trigger the flux AWG
        awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)

        #set up the card trigger
        if self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            if trigger_card_early:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (1 + self.config['generator_enable_delay']), self.config['read_length'])
            else:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (max_ramsey_wait + self.config['card_delay']), self.config['read_length'])

        flux_pulsing = False

        #setup qubit 1 (always the same!)
        #---------------------------------

        wait_time = self.config['pulse_wait']


        #the 2nd pulse and flux are centered and at a fixed time

        if self.config['filter_vary'][0:9] == 'adiabatic':
            if self.config['filter_vary'] == 'adiabatic1':
                pulse_center2 = max_ramsey_wait - wait_time / 2.0 + self.config['pulse2_delay']
            else:
                pulse_center2 = max_ramsey_wait

            pulse_center1 = max_ramsey_wait - wait_time
            if self.config['trap_time'][1] > self.config['trap_time'][0]:
                raise NameError(
                    "For the adiabatic filter test the qubit 2 flux time must be less than the qubit 1 flux time")

        elif self.config['filter_vary'] == 'lifetime':
            pulse_center2 = max_ramsey_wait
            #center of the second flux pulse
            pulse_center2b = max_ramsey_wait - wait_time / 2.0

        for i in range(len(vary_array) + self.config['num_calib_seqs']):

            for j in range(2):

                #set the fast flux pluse
                #------------
                if self.config['do_fast_flux']:


                    indflux = awg_seq.add_analog_pulse(awg_seq.channel_index['flux'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['flux'][j]][i] = indflux

                    flux_pulsing = True
                    flux_val = self.config['flux_height'][j]

                    #indflux = 0

                    if flux_val != 0:

                        #no pedestals for the filter tests
                        pedestal_time = self.config['pedestal_time'][j]
                        pedestal_height = self.config['pedestal_height'][j]

                        trap_time = self.config['trap_time'][j]

                        if self.config['filter_vary'] == 'adiabatic2' or self.config['filter_vary'] == 'lifetime':
                            vary_i = vary_array[mod(i, len(vary_array))]
                        else:


                            if i < len(vary_array):
                                vary_i = vary_array[i]
                            else:
                                vary_i = vary_array[0]

                        #print j,flux_val,pedestal_height
                        if j == 1:
                            flux_cntr = max_ramsey_wait - wait_time / 2.0

                            if self.config['filter_vary'] == 'adiabatic1':
                                flux_cntr = pulse_center2
                            elif self.config['filter_vary'] == 'adiabatic2':
                                self.config['tan_pulse_up_duration'][1][0] = vary_i
                            elif self.config['filter_vary'] == 'lifetime':
                                flux_cntr = pulse_center2b
                                self.config['trap_slope_up'][1] = flux_val / 1.0

                                #flux_val = 0
                        else:
                            if self.config['filter_vary'][0:9] == 'adiabatic':
                                #for the adaibatic test both qubit flux pulses are centered on each other
                                flux_cntr = max_ramsey_wait - wait_time / 2.0

                                if self.config['flux_pulse_type'][j] == 2:

                                    #if using adiabatic_calib then just use the specificed values
                                    #for the slope                                    

                                    if self.config['filter_vary'] == 'adiabatic1':
                                        self.config['trap_slope_up'][0] = vary_i
                                        self.config['trap_slope_down'][0] = flux_val / 2.0
                                    elif self.config['filter_vary'] == 'adiabatic2':
                                        self.config['trap_slope_up'][0] = vary_i
                                        self.config['trap_slope_down'][0] = vary_i


                                elif self.config['flux_pulse_type'][j] == 3:

                                    if self.config['filter_vary'] == 'adiabatic1':
                                        self.config['tan_pulse_up_duration'][j][0] = vary_i

                                    elif self.config['filter_vary'] == 'adiabatic2':
                                        self.config['tan_pulse_up_duration'][j][0] = vary_i



                                else:
                                    raise NameError("Adiabaticity tests can only be run with the trapezoid pulses")

                            elif self.config['filter_vary'] == 'lifetime':

                                if self.config['flux_pulse_type'][j] != 2:
                                    raise NameError("Lifetime Can only be run with trap pulses")

                                flux_cntr = pulse_center2b - vary_i
                                pulse_center1 = flux_cntr - wait_time / 2.0

                                #down slope is very large (so that the pulse is diabatic)
                                self.config['trap_slope_down'][0] = flux_val / 1.0

                        if len(self.config['comp_scaling']) == 0:
                            comp_scale = 1.0
                        else:
                            comp_scale = self.config['comp_scaling'][mod(i, len(vary_array))]

                        #print comp_scale

                        self.__flux_pulses__(awg_seq, indflux, flux_cntr, flux_val, trap_time, pedestal_time,
                                             pedestal_height, j, max_ramsey_wait, comp_scale)

                #------------

                #Do the gate pulses!!
                #-------------
                pulse1_height = 0.0
                pulse2_height = 0.0

                if i >= len(vary_array):


                    pulse1_height = 0.0
                    pulse2_height = 0.0

                    if self.config['filter_vary'] == 'adiabatic2':
                        if i >= (2 * len(vary_array)) and j == 0:
                            pulse2_height = self.config['second_pulse_height'][j]
                    elif self.config['filter_vary'] == 'lifetime':
                        if i >= (2 * len(vary_array)) and j == 1:
                            pulse2_height = self.config['second_pulse_height'][j]
                    else:
                        if i == (len(vary_array) + 1) and j == 0:
                            pulse2_height = self.config['first_pulse_height'][j]


                else:

                    if self.config['filter_vary'] == 'lifetime':

                        if j == 0:
                            pulse1_height = self.config['first_pulse_height'][0]

                    else:

                        pulse1_height = self.config['first_pulse_height'][j]
                        if j == 1:
                            #there is no first pulse for qubit 2
                            pulse1_height = 0.0

                            #no qubit 2 pulsing for adiabatic2
                            if self.config['filter_vary'] == 'adiabatic2':
                                pulse2_height = 0

                        pulse2_height = self.config['second_pulse_height'][j]
                        if j == 0:
                            #there is no second pulse for qubit 1
                            pulse2_height = 0.0

                            #if calibrating the adiabatic pulses
                            if self.config['filter_vary'] == 'adiabatic_calib':
                                pulse2_height = vary_array[i]
                                pulse1_height = 0.0
                            elif self.config['filter_vary'] == 'adiabatic_calib2':
                                pulse1_height = vary_array[i]
                                pulse2_height = 0.0

                pulse1_len = self.config['first_pulse_length'][j]
                pulse2_len = self.config['second_pulse_length'][j]

                ind1a, ind1b, ind2 = self.__add_pulses__(awg_seq, pulse_center1, pulse_center2, pulse1_height,
                                                         pulse2_height, pulse1_len, pulse2_len, i, j)

                #-------------  

                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].offset(
                    self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].offset(
                    self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])

                #-------------  

                if plotter is not None:
                    cur_pulse_seq.plot_pulses(plotter)


        #load pulses into the awg (TEK)
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #load pulses into the AGILENT AWG
        if flux_pulsing and self.config['load_agilent']:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])


#Do single qubit pulse tests
class pulse_tests(qubit_exp):
    def __init__(self):

        #note that some of the qubit_fastflux options will not apply
        qubit_exp.__init__(self)

        #IGNORE XDATA, the vary parameters is the pulse height

        #max number of times to repeat the pulse
        self.config['num_pulses'] = 10

        self.config['exp_id'] = 'PULSE_TEST'

        #the length of the pulse being tested
        self.config['pulse_length'] = 15.0
        self.config['pulse_phase'] = 0.0

        #gap between pulses
        self.config['pulse_gap'] = 5.0

        #different pulse heights
        self.config['pulse_heights'] = [0.1, 0.2]

        #number of pulses to do together
        #NEEDS TO BE AN INT!
        self.config['pulse_factor'] = 1

    #load the ramsey oscillation into the AWG
    def load_exp(self, plotter=None, ramsey_index=-1):


        self.xdata = range(self.config['num_pulses'] * len(self.config['pulse_heights']))


        #ALWAYS Shaped pulse
        self.config['shaped_pulse'] = True


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        pulse_factor = self.config['pulse_factor']

        total_pulse_time = (2 * self.config['pulse_length'] + self.config['pulse_gap']) * self.config[
            'num_pulses'] * pulse_factor

        #where to start the card and read pulses
        max_ramsey_wait = self.config['front_buffer'] + self.config['generator_enable_delay'] + total_pulse_time

        if ramsey_index == -1:
            vary_array = self.xdata
        else:
            vary_array = self.xdata[ramsey_index]
            print vary_array

        awg_seq = pulse_sequence(len(vary_array), int(ceil(self.config['total_length'] * self.config['awg_clocks'][0])),
                                 int(ceil(self.config['total_length'] * self.config['awg_clocks'][1])),
                                 self.config['awg_clocks'][0], self.config['awg_clocks'][1], self.config['awg_delay'])

        for i in range(2):
            if self.config['qubit_enable'][i]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][i]][0].square_pulse(
                    (max_ramsey_wait + self.config['read_delay']), self.config['read_length'])

        #trigger the flux AWG
        #awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1,100)        

        #set up the card trigger
        if self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            if trigger_card_early:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (1 + self.config['generator_enable_delay']), self.config['read_length'])
            else:
                awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(
                    (max_ramsey_wait + self.config['card_delay']), self.config['read_length'])

        flux_pulsing = False

        for i in range(len(vary_array)):

            print i

            for j in range(2):

                if self.config['qubit_enable'][j]:

                    num_pulses = mod(i, self.config['num_pulses']) + 1 * 0
                    pulse_height = self.config['pulse_heights'][int(floor(i / self.config['num_pulses']))]
                    phase = self.config['pulse_phase']
                    pulse_len = self.config['pulse_length']


                    #Add the qubit pulses to the awg
                    ind1a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind1a

                    ind1b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind1b

                    ind2 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                    awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind2


                    #ind1a = 0
                    #ind1b = 0
                    #ind2 = 0
                    cur_loc = max_ramsey_wait

                    for k in range(num_pulses * pulse_factor):

                        #add pulse
                        if self.config['drive_sideband'][j]:

                            awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].gauss_pulse_with_freq(cur_loc,
                                                                                                               pulse_len / 2.0,
                                                                                                               pulse_height,
                                                                                                               self.config[
                                                                                                                   'drive_sideband_freq'][
                                                                                                                   j] / 1.0e9,
                                                                                                               0.0,
                                                                                                               False)
                            awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].gauss_pulse_with_freq(cur_loc,
                                                                                                               pulse_len / 2.0,
                                                                                                               -1 * sign(
                                                                                                                   self.config[
                                                                                                                       'drive_sideband_freq'][
                                                                                                                       j]) * pulse_height,
                                                                                                               self.config[
                                                                                                                   'drive_sideband_freq'][
                                                                                                                   j] / 1.0e9,
                                                                                                               pi / 2.0,
                                                                                                               False)

                        else:
                            if self.config['drag_pulse'][j]:
                                awg_seq.DRAG_pulse(cur_loc, pulse_len / 2.0, pulse_height, phase, j, ind1a, ind1b,
                                                   self.config['drag_prefactor'][j], False)
                            else:
                                awg_seq.gauss_pulse(cur_loc, pulse_len / 2.0, pulse_height, phase, j, ind1a, ind1b,
                                                    False)

                        cur_loc -= self.config['pulse_gap'] + 2 * pulse_len

                    #keep the drive switch on for the whole time 
                    if self.config['pulsed_drive'][j]:
                        awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind2].square_pulse(
                            max_ramsey_wait - total_pulse_time - self.config['generator_enable_delay'],
                            total_pulse_time + 1.2 * self.config['generator_enable_delay'])

                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind1a].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
                    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind1b].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])

                    #-------------  

                if plotter is not None:
                    cur_pulse_seq.plot_pulses(plotter)


        #load pulses into the awg (TEK)
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #load pulses into the AGILENT AWG
        if flux_pulsing and self.config['load_agilent']:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])


class mixer_test(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'MIX_TEST'

        #custom experiment settings
        #These define the "I" channel pi/2 pulses
        self.config['pulse1_height'] = [0.3, 0.3]
        self.config['pulse1_length'] = [15.0, 15.0]

        self.config['pulse2_height'] = [0.3, 0.3]
        self.config['pulse2_length'] = [15.0, 15.0]

        self.config['pulse_gap'] = 10.0

        self.config[
            'mixer_vary'] = 'pulse1'  #pulse1: vary pulse 1 height to find the pi/2 pulse, qpulse: vary the q pulse phase to determine the orthogonal axis

        self.config['drive_angle'] = [0.0, 0.0]


    #load the mixer test into the AWG
    def load_exp(self, plotter=None, rabi_index=-1):

        if len(self.xdata) == 0:
            raise NameError("Mixer test sequence not initialized")

        max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (
        max(self.config['pulse1_length'])) + 2 * (max(self.config['pulse2_length'])) + self.config['pulse_gap']

        if rabi_index == -1:
            rabi_vary = self.xdata
        else:
            rabi_vary = zeros(1)
            rabi_vary[0] = self.xdata[rabi_index]


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #allocate pulse array
        #self.rabi_fullpulsesequence = zeros((max(len(self.rabi_height),len(self.rabi_length)),6,self.config['total_length']*awg_modifier))

        if self.config['mixer_vary'] == 'pulse1':
            self.config['num_calib_seqs'] = 0
        elif self.config['mixer_vary'] == 'qpulse':
            self.config['num_calib_seqs'] = 2
        else:
            raise NameError("Invalid mixer vary")

        awg_seq = pulse_sequence(len(rabi_vary) + self.config['num_calib_seqs'],
                                 self.config['total_length'] * self.config['awg_clocks'][0],
                                 TEK_clk=self.config['awg_clocks'][0], Agilent_clk=self.config['awg_clocks'][1],
                                 delays=self.config['awg_delay'])


        #define card trigger and read trigger once (these are the same waveforms each time)
        for j in range(2):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (max_rabi + self.config['read_delay']), self.config['read_length'])

        if trigger_card_early or self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((max_rabi + self.config['card_delay']),
                                                                               self.config['read_length'])

        for i in range(len(rabi_vary) + self.config['num_calib_seqs']):


            print i

            for j in range(2):

                if not self.config['qubit_enable'][j]:
                    continue

                pulse1_height = self.config['pulse1_height'][j]
                pulse1_length = self.config['pulse1_length'][j]
                pulse2_height = self.config['pulse2_height'][j]
                pulse2_length = self.config['pulse2_length'][j]
                pulse2_phase = self.config['drive_angle'][j]

                if i >= len(rabi_vary):
                    pulse2_height = 0
                else:
                    if self.config['mixer_vary'] == 'pulse1':
                        pulse1_height = rabi_vary[i]
                        pulse2_height = 0
                    elif self.config['mixer_vary'] == 'qpulse':
                        pulse2_phase = rabi_vary[i]

                if self.config['shaped_pulse']:

                    if len(rabi_vary) > 1:
                        ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind1
                        ind2a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind2a
                        ind2b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind2b
                    else:
                        ind1 = 0
                        ind2a = 0
                        ind2b = 0

                    pulse1_center = max_rabi - 2 * pulse2_length - pulse1_length - self.config['pulse_gap']
                    pulse2_center = max_rabi - pulse2_length

                    #drive drigger
                    awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(
                        (pulse1_center + pulse2_center) / 2, (
                        pulse2_center - pulse1_center + pulse1_length + pulse2_length + 2 * self.config[
                            'generator_enable_delay']))

                    #first pulse
                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse(pulse1_center,
                                                                                             pulse1_length / 2,
                                                                                             pulse1_height, False)

                    #second pulse
                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse(pulse2_center,
                                                                                             pulse2_length / 2, sin(
                            pulse2_phase) * pulse2_height, False)
                    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse(pulse2_center,
                                                                                             pulse2_length / 2, cos(
                            pulse2_phase) * pulse2_height, False)

                    #add extra offset to try and null mixer
                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
                    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])


                else:

                    pass

            if plotter is not None:
                cur_pulse_seq.plot_pulses(plotter)

        #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #self.rabi_fullpulsesequence[i] = cur_pulse_seq.output_pulses()


class gate_w_process_tomography(gate_w_statetomography):
    def __init__(self):

        gate_w_statetomography.__init__(self)

        self.config['exp_id'] = 'PROCESSTOMO'

        self.config['first_pulse_pi'] = [0.1, 0.1]  #pi pulse for the first pulse (I channel)
        self.config['first_pulse_Q'] = [0.1, 0.1]  #pi/2 pulse for the first pulse (Q channel)
        self.config['exp_path'] = ''

    def run_qubit_exp(self, plotter=None):

        #run through the different process tomography steps
        self.config['tomo_vary'] = 'tomo'
        self.config['save_each_run'] = True

        #set the wait time between the gate pulse and the tomography pulse to negative
        #so that the gate pulse and the tomography pulse will be identical
        self.config['tomo_wait'] = -(max(self.config['second_pulse_length']) + self.config['tomo_pulse_length'])

        #No second pulse for the process tomography
        self.config['second_pulse_height'] = [0, 0]

        process_index = [[0.0, 0.0], self.config['first_pulse_height'], self.config['first_pulse_Q'],
                         self.config['first_pulse_pi']]
        process_phase = [[0.0, 0.0], [0.0, 0.0], [pi / 2, pi / 2], [0.0, 0.0]]

        #run the 16 input states through the filter!
        for i in range(4):
            for j in range(4):
                self.config['first_pulse_height'] = [process_index[i][0], process_index[j][1]]
                self.config['first_pulse_phase'] = [process_phase[i][0], process_phase[j][1]]

                #load new settings                
                self.load_exp()

                #run the experiment
                gate_w_statetomography.run_qubit_exp(self, plotter)

                #save data
                self.save_data(self.config['exp_path'], self.config['exp_id'] + "_" + str(i) + str(j), overwrite=False)


class sideband_rabi(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'SIDEBAND'

        self.config['pi_pulse_height'] = [1.0, 1.0]
        self.config['pi_pulse_width'] = [1.0, 1.0]
        self.config['drive_angle'] = 0.0
        self.config['pi_to_flux_delay'] = [1.0, 1.0]
        self.config['flux_pulse_height'] = [1.0, 1.0]
        self.config['flux_pulse_width'] = [1.0, 1.0]
        self.config['sideband_freq'] = [0.5, 0.5]
        self.config['freq_sweep'] = False
        self.config['height_sweep'] = False

        #parameter options: 'pi_width','pi_height','delay',flux_height','flux_width'
        self.config['sideband_vary'] = 'pi_width'
        if self.config['freq_sweep']:
            self.config['sideband_sweep_end'] = [0.7, 0.7]
            self.config['sideband_sweep_points'] = [50, 50]
            self.config['SB_freq_points'] = [
                linspace(self.config['sideband_freq'][0], self.config['sideband_sweep_end'][0],
                         self.config['sideband_sweep_points'][0]),
                linspace(self.config['sideband_freq'][1], self.config['sideband_sweep_end'][1],
                         self.config['sideband_sweep_points'][1])]
        if self.config['height_sweep']:
            self.config['height_sweep_end'] = [1.0, 1.0]
            self.config['height_sweep_points'] = [10, 10]
            self.config['SB_height_points'] = [
                linspace(self.config['flux_pulse_height'][0], self.config['height_sweep_end'][0],
                         self.config['height_sweep_points'][0]),
                linspace(self.config['flux_pulse_height'][1], self.config['height_sweep_end'][1],
                         self.config['height_sweep_points'][1])]

    def run_qubit_exp(self, plotter=None):

        if self.config['freq_sweep']:

            if self.config['height_sweep']:
                raise NameError('Cannot run both frequency and height sweep')

            if self.config['num_avgs'] == -1:
                raise NameError('Must run a finite number of averages')

            #XXXX
            temp_ydata = zeros((2, len(self.config['SB_freq_points'][0]), len(self.xdata)))

            num_avgs = self.config['num_avgs']

            if plotter is not None:
                for i in range(2):
                    plotter.init_plot("FullSpectData" + str(i), rank=2, accum=False)

            #XXX               
            for i in range(len(self.config['SB_freq_points'][0])):


                #set the drive frequencies  
                #XXX

                self.config['sideband_freq'] = [self.config['SB_freq_points'][0][i],
                                                self.config['SB_freq_points'][1][i]]

                #XXX
                print self.config['sideband_freq']

                sideband_rabi.load_exp(self)

                #run!
                qubit_exp.run_qubit_exp(self, plotter)

                if self.config['num_avgs'] != num_avgs:
                    #was interrupted prematurely
                    print("Sideband sweep spectroscopy interrupted")
                    break

                for j in range(2):
                    temp_ydata[j, i] = self.ydata_avg[j]

                if plotter is not None:
                    for j in range(2):
                        plotter.plot(temp_ydata[j], "FullSpectData" + str(j))

            self.ydata = temp_ydata

        elif self.config['height_sweep']:

            if self.config['num_avgs'] == -1:
                raise NameError('Must run a finite number of averages')

            #XXXX
            temp_ydata = zeros((2, len(self.config['SB_height_points'][0]), len(self.xdata)))

            num_avgs = self.config['num_avgs']

            if plotter is not None:
                for i in range(2):
                    plotter.init_plot("FullSpectData" + str(i), rank=2, accum=False)

            #XXX               
            for i in range(len(self.config['SB_height_points'][0])):


                #set the drive frequencies  
                #XXX

                self.config['flux_pulse_height'] = [self.config['SB_height_points'][0][i],
                                                    self.config['SB_height_points'][1][i]]

                #XXX
                print self.config['flux_pulse_height']

                sideband_rabi.load_exp(self)

                #run!
                qubit_exp.run_qubit_exp(self, plotter)

                if self.config['num_avgs'] != num_avgs:
                    #was interrupted prematurely
                    print("Sideband sweep spectroscopy interrupted")
                    break

                for j in range(2):
                    temp_ydata[j, i] = self.ydata_avg[j]

                if plotter is not None:
                    for j in range(2):
                        plotter.plot(temp_ydata[j], "FullSpectData" + str(j))

            self.ydata = temp_ydata

        else:
            #run normally!
            qubit_exp.run_qubit_exp(self, plotter)

    def load_exp(self):

        if len(self.xdata) == 0:
            raise NameError("Sideband Rabi sequence not initialized")

        if self.config['sideband_vary'] == 'pi_height':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + max(
                self.config['pi_to_flux_delay']) + 2 * (
            max(self.config['pi_pulse_width']) + max(self.config['flux_pulse_width']))
        elif self.config['sideband_vary'] == 'pi_width':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + max(
                self.config['pi_to_flux_delay']) + 2 * (self.xdata[-1] + max(self.config['flux_pulse_width']))
        elif self.config['sideband_vary'] == 'flux_height':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + max(
                self.config['pi_to_flux_delay']) + 2 * (
            max(self.config['pi_pulse_width']) + max(self.config['flux_pulse_width']))
        elif self.config['sideband_vary'] == 'flux_width':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + max(
                self.config['pi_to_flux_delay']) + 2 * (max(self.config['pi_pulse_width']) + self.xdata[-1])
        elif self.config['sideband_vary'] == 'delay':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + self.xdata[-1] + 2 * (
            max(self.config['pi_pulse_width']) + max(self.config['flux_pulse_width']))

        sideband_vary = self.xdata


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #allocate pulse array
        #self.rabi_fullpulsesequence = zeros((max(len(self.rabi_height),len(self.rabi_length)),6,self.config['total_length']*awg_modifier))

        awg_seq = pulse_sequence(len(sideband_vary), self.config['total_length'] * self.config['awg_clocks'][0],
                                 TEK_clk=self.config['awg_clocks'][0], Agilent_clk=self.config['awg_clocks'][1],
                                 delays=self.config['awg_delay'])


        #define card trigger and read trigger once (these are the same waveforms each time)
        for j in range(2):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (end_time + self.config['read_delay']), self.config['read_length'])

        #trigger the flux AWG
        awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)

        if trigger_card_early or self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((end_time + self.config['card_delay']),
                                                                               self.config['read_length'])

        for i in range(len(sideband_vary)):


            print sideband_vary[i]

            for j in range(2):

                if not self.config['qubit_enable'][j]:
                    continue

                vary_rabi = False
                vary_flux = False

                if self.config['sideband_vary'] == 'pi_height':
                    pi_height_i = sideband_vary[i]
                    pi_width_i = self.config['pi_pulse_width'][j]
                    flux_height_i = self.config['flux_pulse_height'][j]
                    flux_width_i = self.config['flux_pulse_width'][j]
                    delay_i = self.config['pi_to_flux_delay'][j]
                    vary_rabi = True
                elif self.config['sideband_vary'] == 'pi_width':
                    pi_height_i = self.config['pi_pulse_height'][j]
                    pi_width_i = sideband_vary[i]
                    flux_height_i = self.config['flux_pulse_height'][j]
                    flux_width_i = self.config['flux_pulse_width'][j]
                    delay_i = self.config['pi_to_flux_delay'][j]
                    vary_rabi = True
                elif self.config['sideband_vary'] == 'flux_height':
                    pi_height_i = self.config['pi_pulse_height'][j]
                    pi_width_i = self.config['pi_pulse_width'][j]
                    flux_height_i = sideband_vary[i]
                    flux_width_i = self.config['flux_pulse_width'][j]
                    delay_i = self.config['pi_to_flux_delay'][j]
                    vary_flux = True
                elif self.config['sideband_vary'] == 'flux_width':
                    pi_height_i = self.config['pi_pulse_height'][j]
                    pi_width_i = self.config['pi_pulse_width'][j]
                    flux_height_i = self.config['flux_pulse_height'][j]
                    flux_width_i = sideband_vary[i]
                    delay_i = self.config['pi_to_flux_delay'][j]
                    vary_flux = True
                    vary_rabi = True
                elif self.config['sideband_vary'] == 'delay':
                    pi_height_i = self.config['pi_pulse_height'][j]
                    pi_width_i = self.config['pi_pulse_width'][j]
                    flux_height_i = self.config['flux_pulse_height'][j]
                    flux_width_i = self.config['flux_pulse_width'][j]
                    delay_i = sideband_vary[i]
                    vary_rabi = True

                pi_center_i = end_time - 2 * flux_width_i - delay_i - pi_width_i
                flux_center_i = end_time - flux_width_i

                if vary_rabi:

                    ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                    awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind1
                    ind2a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind2a
                    ind2b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind2b

                else:
                    ind1 = 0
                    ind2a = 0
                    ind2b = 0

                if vary_flux:

                    indflux = awg_seq.add_analog_pulse(awg_seq.channel_index['flux'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['flux'][j]][i] = indflux

                else:
                    indflux = 0

                if vary_rabi or i == 0:

                    if self.config['pulsed_drive'][j]:
                        pass
                    #                        if j==0:
                    #                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(pi_center_i, (2*pi_width_i + 2*self.config['generator_enable_delay']))

                    if self.config['drag_pulse'][j]:
                        awg_seq.DRAG_pulse(pi_center_i, pi_width_i / 2, pi_height_i, self.config['drive_angle'], j,
                                           ind2a, ind2b, self.config['drag_prefactor'][j], False)
                    else:
                        awg_seq.gauss_pulse(pi_center_i, pi_width_i / 2, pi_height_i, self.config['drive_angle'], j,
                                            ind2a, ind2b, False)

                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
                    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])

                if vary_flux or i == 0:

                    if size(self.config['sideband_freq'][j]) == 1:
                        #awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].gauss_pulse_with_freq(flux_center_i, flux_width_i/2, flux_height_i, self.config['sideband_freq'][j])

                        if j == 1:
                            awg_seq.gauss_pulse(flux_center_i, flux_width_i / 2, flux_height_i, 0.0, j, ind2a, ind2b)
                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(flux_center_i,
                                                                                                       2 * (
                                                                                                       flux_width_i +
                                                                                                       self.config[
                                                                                                           'generator_enable_delay']))
                    else:
                        awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].gauss_pulse_with_multiple_freq(
                            flux_center_i, flux_width_i / 2, flux_height_i, self.config['sideband_freq'][j])


        #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #load pulses into the AGILENT AWG
        if self.config['load_agilent']:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])


class blue_sideband(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'BLUE_SIDEBAND'

        self.config['pi_pulse_height'] = [1.0, 1.0]
        self.config['pi_pulse_width'] = [1.0, 1.0]
        self.config['drive_angle'] = 0.0
        self.config['pi_to_flux_delay'] = [1.0, 1.0]
        self.config['flux_pulse_height'] = [1.0, 1.0]
        self.config['flux_pulse_width'] = [1.0, 1.0]
        self.config['sideband_freq'] = [0.5, 0.5]
        self.config['freq_sweep'] = False
        self.config['read_sweep'] = False
        self.config['height_sweep'] = False
        self.config['measure_e'] = [False, False]
        self.config['drive_freq'] = [6.0e9, 6.0e9]
        self.config['read_freq'] = [6.0e9, 6.0e9]

        #parameter options: 'pi_width','pi_height','delay',flux_height','flux_width'
        self.config['sideband_vary'] = 'pi_width'
        if self.config['freq_sweep']:
            self.config['freq_sweep_end'] = [0.7, 0.7]
            self.config['freq_sweep_points'] = [50, 50]
            self.config['SB_freq_points'] = [linspace(self.config['drive_freq'][0], self.config['freq_sweep_end'][0],
                                                      self.config['freq_sweep_points'][0]),
                                             linspace(self.config['drive_freq'][1], self.config['freq_sweep_end'][1],
                                                      self.config['freq_sweep_points'][1])]
        if self.config['read_sweep']:
            self.config['read_sweep_end'] = [0.7, 0.7]
            self.config['read_sweep_points'] = [50, 50]
            self.config['SB_read_points'] = [linspace(self.config['read_freq'][0], self.config['read_sweep_end'][0],
                                                      self.config['read_sweep_points'][0]),
                                             linspace(self.config['read_freq'][1], self.config['read_sweep_end'][1],
                                                      self.config['read_sweep_points'][1])]
        if self.config['height_sweep']:
            self.config['height_sweep_end'] = [1.0, 1.0]
            self.config['height_sweep_points'] = [10, 10]
            self.config['SB_height_points'] = [
                linspace(self.config['flux_pulse_height'][0], self.config['height_sweep_end'][0],
                         self.config['height_sweep_points'][0]),
                linspace(self.config['flux_pulse_height'][1], self.config['height_sweep_end'][1],
                         self.config['height_sweep_points'][1])]

        self.config['num_calib_seqs'] = 20
        self.config['rabi_height'] = [0.3, 0.3]
        self.config['rabi_length'] = [20, 20]
        self.config['height_scaling'] = [1.0, 1.0]

    def run_qubit_exp(self, plotter=None):

        if self.config['freq_sweep']:

            if self.config['height_sweep']:
                raise NameError('Cannot run both frequency and height sweep')

            if self.config['num_avgs'] == -1:
                raise NameError('Must run a finite number of averages')

            #XXXX
            temp_ydata = zeros((2, len(self.config['SB_freq_points'][0]), len(self.xdata)))

            num_avgs = self.config['num_avgs']

            if plotter is not None:
                for i in range(2):
                    plotter.init_plot("FullSpectData" + str(i), rank=2, accum=False)

            #XXX               
            for i in range(len(self.config['SB_freq_points'][0])):


                #set the drive frequencies  
                #XXX

                self.config['drive_freq'] = [self.config['SB_freq_points'][0][i], self.config['SB_freq_points'][1][i]]

                #XXX
                print self.config['drive_freq']

                blue_sideband.load_exp(self)

                #run!
                qubit_exp.run_qubit_exp(self, plotter)

                if self.config['num_avgs'] != num_avgs:
                    #was interrupted prematurely
                    print("Sideband sweep spectroscopy interrupted")
                    break

                for j in range(2):
                    temp_ydata[j, i] = self.ydata_avg[j]

                if plotter is not None:
                    for j in range(2):
                        plotter.plot(temp_ydata[j], "FullSpectData" + str(j))

            self.ydata = temp_ydata

        elif self.config['read_sweep']:

            if self.config['height_sweep']:
                raise NameError('Cannot run both frequency and height sweep')

            if self.config['num_avgs'] == -1:
                raise NameError('Must run a finite number of averages')

            #XXXX
            temp_ydata = zeros((2, len(self.config['SB_read_points'][0]), len(self.xdata)))

            num_avgs = self.config['num_avgs']

            if plotter is not None:
                for i in range(2):
                    plotter.init_plot("FullSpectData" + str(i), rank=2, accum=False)

            #XXX               
            for i in range(len(self.config['SB_read_points'][0])):


                #set the drive frequencies  
                #XXX

                self.config['read_freq'] = [self.config['SB_read_points'][0][i], self.config['SB_read_points'][1][i]]
                #self.config['read_phase'][0] = i*180.0/len(self.config['SB_read_points'][0])

                blue_sideband.load_exp(self)

                #run!
                qubit_exp.run_qubit_exp(self, plotter)

                if self.config['num_avgs'] != num_avgs:
                    #was interrupted prematurely
                    print("Sideband sweep spectroscopy interrupted")
                    break

                for j in range(2):
                    temp_ydata[j, i] = self.ydata_avg[j]

                if plotter is not None:
                    for j in range(2):
                        plotter.plot(temp_ydata[j], "FullSpectData" + str(j))

            self.ydata = temp_ydata

        elif self.config['height_sweep']:

            if self.config['num_avgs'] == -1:
                raise NameError('Must run a finite number of averages')

            #XXXX
            temp_ydata = zeros((2, len(self.config['SB_height_points'][0]), len(self.xdata)))

            num_avgs = self.config['num_avgs']

            if plotter is not None:
                for i in range(2):
                    plotter.init_plot("FullSpectData" + str(i), rank=2, accum=False)

            #XXX               
            for i in range(len(self.config['SB_height_points'][0])):


                #set the drive frequencies  
                #XXX

                self.config['flux_pulse_height'] = [self.config['SB_height_points'][0][i],
                                                    self.config['SB_height_points'][1][i]]

                #XXX
                print self.config['flux_pulse_height']

                blue_sideband.load_exp(self)

                #run!
                qubit_exp.run_qubit_exp(self, plotter)

                if self.config['num_avgs'] != num_avgs:
                    #was interrupted prematurely
                    print("Sideband sweep spectroscopy interrupted")
                    break

                for j in range(2):
                    temp_ydata[j, i] = self.ydata_avg[j]

                if plotter is not None:
                    for j in range(2):
                        plotter.plot(temp_ydata[j], "FullSpectData" + str(j))

            self.ydata = temp_ydata

        else:
            #run normally!
            qubit_exp.run_qubit_exp(self, plotter)

    def load_exp(self):

        if len(self.xdata) == 0:
            raise NameError("Sideband Rabi sequence not initialized")

        if self.config['sideband_vary'] == 'pi_height':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + max(
                self.config['pi_to_flux_delay']) + 2 * (
            max(self.config['pi_pulse_width']) + max(self.config['flux_pulse_width']))
        elif self.config['sideband_vary'] == 'pi_width':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + max(
                self.config['pi_to_flux_delay']) + 2 * (self.xdata[-1] + max(self.config['flux_pulse_width']))
        elif self.config['sideband_vary'] == 'flux_height':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + max(
                self.config['pi_to_flux_delay']) + 2 * (
            max(self.config['pi_pulse_width']) + max(self.config['flux_pulse_width']))
        elif self.config['sideband_vary'] == 'flux_width':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + max(
                self.config['pi_to_flux_delay']) + 2 * (max(self.config['pi_pulse_width']) + self.xdata[-1])
        elif self.config['sideband_vary'] == 'delay':
            end_time = self.config['front_buffer'] + self.config['generator_enable_delay'] + self.xdata[-1] + 2 * (
            max(self.config['pi_pulse_width']) + max(self.config['flux_pulse_width']))

        sideband_vary = self.xdata


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #allocate pulse array
        #self.rabi_fullpulsesequence = zeros((max(len(self.rabi_height),len(self.rabi_length)),6,self.config['total_length']*awg_modifier))

        awg_seq = pulse_sequence(len(sideband_vary) + self.config['num_calib_seqs'],
                                 self.config['total_length'] * self.config['awg_clocks'][0],
                                 TEK_clk=self.config['awg_clocks'][0], Agilent_clk=self.config['awg_clocks'][1],
                                 delays=self.config['awg_delay'])


        #define card trigger and read trigger once (these are the same waveforms each time)
        for j in range(1):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (end_time + self.config['read_delay']), self.config['read_length'])

        #trigger the flux AWG
        awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)

        if trigger_card_early or self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((end_time + self.config['card_delay']),
                                                                               self.config['read_length'])

        for i in range(len(sideband_vary)):


            print sideband_vary[i]

            for j in range(2):

                if not self.config['qubit_enable'][j]:
                    continue

                vary_rabi = False
                vary_flux = False

                if self.config['sideband_vary'] == 'pi_height':
                    pi_height_i = sideband_vary[i]
                    pi_width_i = self.config['pi_pulse_width'][j]
                    flux_height_i = self.config['flux_pulse_height'][j]
                    flux_width_i = self.config['flux_pulse_width'][j]
                    delay_i = self.config['pi_to_flux_delay'][j]
                    vary_rabi = True
                elif self.config['sideband_vary'] == 'pi_width':
                    pi_height_i = self.config['pi_pulse_height'][j]
                    pi_width_i = sideband_vary[i]
                    flux_height_i = self.config['flux_pulse_height'][j]
                    flux_width_i = self.config['flux_pulse_width'][j]
                    delay_i = self.config['pi_to_flux_delay'][j]
                    vary_rabi = True
                elif self.config['sideband_vary'] == 'flux_height':
                    pi_height_i = self.config['pi_pulse_height'][j]
                    pi_width_i = self.config['pi_pulse_width'][j]
                    flux_height_i = sideband_vary[i]
                    flux_width_i = self.config['flux_pulse_width'][j]
                    delay_i = self.config['pi_to_flux_delay'][j]
                    vary_flux = True
                    vary_rabi = True
                elif self.config['sideband_vary'] == 'flux_width':
                    pi_height_i = self.config['pi_pulse_height'][j]
                    pi_width_i = self.config['pi_pulse_width'][j]
                    flux_height_i = self.config['flux_pulse_height'][j]
                    flux_width_i = sideband_vary[i]
                    delay_i = self.config['pi_to_flux_delay'][j]
                    vary_flux = True
                    vary_rabi = True
                elif self.config['sideband_vary'] == 'delay':
                    pi_height_i = self.config['pi_pulse_height'][j]
                    pi_width_i = self.config['pi_pulse_width'][j]
                    flux_height_i = self.config['flux_pulse_height'][j]
                    flux_width_i = self.config['flux_pulse_width'][j]
                    delay_i = sideband_vary[i]
                    vary_rabi = True

                if self.config['measure_e'][0] or self.config['measure_e'][1]:

                    pi_center_i = end_time - 2 * flux_width_i - delay_i - pi_width_i - 2 * max(
                        self.config['pi_pulse_width'])
                    flux_center_i = end_time - flux_width_i - 2 * max(self.config['pi_pulse_width'])

                else:

                    pi_center_i = end_time - 2 * flux_width_i - delay_i - pi_width_i
                    flux_center_i = end_time - flux_width_i

                if vary_rabi:

                    ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                    awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind1
                    ind2a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind2a
                    ind2b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind2b

                else:
                    ind1 = 0
                    ind2a = 0
                    ind2b = 0

                if vary_flux:

                    indflux = awg_seq.add_analog_pulse(awg_seq.channel_index['flux'][j])
                    awg_seq.analog_seq_table[awg_seq.channel_index['flux'][j]][i] = indflux

                else:
                    indflux = 0

                if vary_rabi or i == 0:

                    if self.config['pulsed_drive'][j]:
                        if j == 0:
                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(pi_center_i, (
                            2 * pi_width_i + 2 * self.config['generator_enable_delay']))

                    if self.config['drag_pulse'][j]:
                        awg_seq.DRAG_pulse(pi_center_i, pi_width_i / 2, pi_height_i, self.config['drive_angle'], j,
                                           ind2a, ind2b, self.config['drag_prefactor'][j], False)
                    else:
                        awg_seq.gauss_pulse(pi_center_i, pi_width_i / 2, pi_height_i, self.config['drive_angle'], j,
                                            ind2a, ind2b, False)

                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
                    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])

                if vary_flux or i == 0:

                    if size(self.config['sideband_freq'][j]) == 1:
                        #awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].gauss_pulse_with_freq(flux_center_i, flux_width_i/2, flux_height_i, self.config['sideband_freq'][j])

                        if j == 1:
                            awg_seq.gauss_pulse(flux_center_i, flux_width_i / 2, flux_height_i, 0.0, j, ind2a, ind2b)
                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(flux_center_i,
                                                                                                       2 * (
                                                                                                       flux_width_i +
                                                                                                       self.config[
                                                                                                           'generator_enable_delay']))
                    else:
                        awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].gauss_pulse_with_multiple_freq(
                            flux_center_i, flux_width_i / 2, flux_height_i, self.config['sideband_freq'][j])

                if self.config['measure_e'][j]:
                    awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(end_time - pi_width_i, (
                    2 * pi_width_i + 2 * self.config['generator_enable_delay']))
                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse(end_time - pi_width_i,
                                                                                             pi_width_i / 2,
                                                                                             pi_height_i, False)

        rabi_vary = linspace(0.0, 0.5, self.config['num_calib_seqs'])
        #max_rabi = self.config['front_buffer']+self.config['generator_enable_delay']+2*(max(self.config['rabi_length']))
        #awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((max_rabi+self.config['card_delay']),self.config['read_length'])
        for i in range(self.config['num_calib_seqs']):

            index_offset = len(sideband_vary)
            #do rabi oscillation!

            j = 0

            rabi_pulse_i = floor(self.config['rabi_length'][j])
            pulse_height_i = rabi_vary[i]
            pulse_height_i *= self.config['height_scaling'][j]

            if len(rabi_vary) > 1:
                ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i + index_offset] = ind1
                ind2a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i + index_offset] = ind2a
                ind2b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i + index_offset] = ind2b
            else:
                ind1 = 0
                ind2a = 0
                ind2b = 0

            if self.config['pulsed_drive'][j]:
                awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2((end_time - rabi_pulse_i), (
                2 * rabi_pulse_i + 2 * self.config['generator_enable_delay']))

            if self.config['drive_sideband'][j]:
                awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse_with_freq(
                    (end_time - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                    self.config['drive_sideband_freq'][j] / 1.0e9, 0.0, False)
                awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse_with_freq(
                    (end_time - rabi_pulse_i), rabi_pulse_i / 2,
                    -1 * sign(self.config['drive_sideband_freq'][j]) * pulse_height_i,
                    self.config['drive_sideband_freq'][j] / 1.0e9, pi / 2.0, False)
            else:
                if self.config['drag_pulse'][j]:
                    awg_seq.DRAG_pulse((end_time - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                                       self.config['drive_angle'], j, ind2a, ind2b, self.config['drag_prefactor'][j],
                                       False)
                else:
                    awg_seq.gauss_pulse((end_time - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                                        self.config['drive_angle'], j, ind2a, ind2b, False)

            awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].offset(
                self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
            awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].offset(
                self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])


            #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #load pulses into the AGILENT AWG
        if self.config['load_agilent']:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])


#test tek7k by direct driving the qubit using the awg
class awg_test(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'AWGTEST'

        #custom Rabi experiment settings
        self.config['pi_pulse_height'] = [1.0, 1.0]
        self.config['pi_pulse_width'] = [1.0, 1.0]
        self.config['pi_pulse_freq'] = 5e9
        self.config['pulse_height'] = [0.3, 0.3]
        self.config['pulse_length'] = [10, 10]
        self.config['pulse_vary'] = 'height'  #height,width or delay
        self.config['pulse_freq'] = 5e9
        self.config['drive_angle'] = 0.0
        self.config['height_scaling'] = [1.0, 1.0]
        self.config['pulse_freq_offset'] = 30e6


    #load the rabi oscillation into the AWG
    def load_exp(self, plotter=None, rabi_index=-1):

        if len(self.xdata) == 0:
            raise NameError("Rabi sequence not initialized")

        if self.config['pulse_vary'] == 'height':
            self._rabi_vs_width = 0
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (
            max(self.config['pulse_length']))
        elif self.config['pulse_vary'] == 'width':
            self._rabi_vs_width = 1
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (
            self.xdata[-1] + self.config['pi_pulse_width'][0])
        elif self.config['pulse_vary'] == 'delay':
            self._rabi_vs_width = 2
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (self.xdata[-1])

        if rabi_index == -1:
            rabi_vary = self.xdata
        else:
            rabi_vary = zeros(1)
            rabi_vary[0] = self.xdata[rabi_index]


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #allocate pulse array
        #self.rabi_fullpulsesequence = zeros((max(len(self.rabi_height),len(self.rabi_length)),6,self.config['total_length']*awg_modifier))

        awg_seq = pulse_sequence(len(rabi_vary), self.config['total_length'] * self.config['awg_clocks'][0],
                                 agilent_length=262144, TEK_clk=self.config['awg_clocks'][0], Agilent_clk=50.0,
                                 delays=self.config['awg_delay'])


        #define card trigger and read trigger once (these are the same waveforms each time)
        for j in range(2):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (max_rabi + self.config['read_delay']), self.config['read_length'])

        if trigger_card_early or self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((max_rabi + self.config['card_delay']),
                                                                               self.config['read_length'])

        for i in range(len(rabi_vary)):


            print rabi_vary[i]

            for j in range(1):

                if not self.config['qubit_enable'][j]:
                    continue

                if self._rabi_vs_width == 1:
                    rabi_pulse_i = floor(rabi_vary[i])
                    pulse_height_i = self.config['pulse_height'][j]
                elif self._rabi_vs_width == 0:
                    rabi_pulse_i = floor(self.config['pulse_length'][j])
                    pulse_height_i = rabi_vary[i]

                else:
                    rabi_pulse_i = floor(self.config['pulse_length'][j])
                    pulse_height_i = self.config['pulse_height'][j]

                pulse_height_i *= self.config['height_scaling'][j]

                if len(rabi_vary) > 1:
                    ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['flux_trig'])
                    awg_seq.marker_seq_table[awg_seq.channel_index['flux_trig']][i] = ind1

                    ind2a = awg_seq.add_analog_pulse(4)
                    awg_seq.analog_seq_table[4][i] = ind2a

                else:
                    ind1 = 0
                    ind2a = 0
                    ind2b = 0

                pulse_center = max_rabi - rabi_pulse_i - self.config['pi_pulse_width'][j] - 1500

                #add trigger for tek
                awg_seq.marker[awg_seq.channel_index['flux_trig']][ind1].square_pulse(pulse_center - 6 * rabi_pulse_i,
                                                                                      10)


                #pulse on tek
                if len(self.config['pulse_freq']) == 2:
                    awg_seq.analogwf[4][ind2a].gauss_pulse_with_chirp(3 * rabi_pulse_i, rabi_pulse_i / 2,
                                                                      pulse_height_i,
                                                                      self.config['pulse_freq'][0] / 1.0e9,
                                                                      self.config['pulse_freq'][1] / 1.0e9, 0.0, False)
                elif len(self.config['pulse_freq']) > 2:
                    #awg_seq.analogwf[4][ind2a].gauss_pulse_with_multiple_freq(3*rabi_pulse_i, rabi_pulse_i/2, pulse_height_i/len(self.config['pulse_freq']), self.config['pulse_freq']/1.0e9, 0.0, False)
                    awg_seq.analogwf[4][ind2a].gauss_pulse_with_freq(3 * rabi_pulse_i, rabi_pulse_i / 2,
                                                                     self.config['pulse_height'][j],
                                                                     self.config['pulse_freq'][0] / 1.0e9, 0.0, False)
                    #print self.config['pulse_freq'][1]/1.0e9                    
                    awg_seq.analogwf[4][ind2a].gauss_pulse_with_freq(6 * rabi_pulse_i, rabi_pulse_i / 2, pulse_height_i,
                                                                     self.config['pulse_freq'][1] / 1.0e9, 0.0, False)
                    #awg_seq.analogwf[4][ind2a].gauss_pulse_with_freq(9*rabi_pulse_i, rabi_pulse_i/2, pulse_height_i, self.config['pulse_freq'][2]/1.0e9, 0.0, False)

                else:
                    awg_seq.analogwf[4][ind2a].gauss_pulse_with_freq(3 * self.config['pi_pulse_width'][j],
                                                                     self.config['pi_pulse_width'][j] / 2,
                                                                     self.config['pi_pulse_height'][j],
                                                                     self.config['pi_pulse_freq'] / 1.0e9, 0.0, False)
                    awg_seq.analogwf[4][ind2a].gauss_pulse_with_freq(
                        4 * self.config['pi_pulse_width'][j] + rabi_pulse_i, rabi_pulse_i / 2, pulse_height_i,
                        self.config['pulse_freq'][0] / 1.0e9, 0.0, False)
                if self._rabi_vs_width == 2:
                    #add second pulse
                    awg_seq.analogwf[4][ind2a].gauss_pulse_with_freq(3 * rabi_pulse_i + rabi_vary[i], rabi_pulse_i / 2,
                                                                     pulse_height_i, self.config['pulse_freq'] / 1.0e9,
                                                                     rabi_vary[i] * (
                                                                     self.config['pulse_freq'] + self.config[
                                                                         'pulse_freq_offset']) / 1.0e9 * 2 * pi, False)

            if plotter is not None:
                cur_pulse_seq.plot_pulses(plotter)

        #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        awg_seq.load_into_tek2('S:\\David M\\python scripts\\testing new tek', 'TEK2')


class rabi_ef(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'RABI_EF'

        #custom Rabi experiment settings
        self.config['pi_pulse_height'] = [1.0, 1.0]
        self.config['pi_pulse_width'] = [1.0, 1.0]
        self.config['pi_pulse_freq'] = 5e9
        self.config['rabi_height'] = [0.3, 0.3]
        self.config['rabi_length'] = [10, 10]
        self.config['rabi_vary'] = 'height'  #height or width ..then the xdata is the vary parameter
        self.config['drive_angle'] = 0.0
        self.config['height_scaling'] = [1.0, 1.0]


    #load the rabi oscillation into the AWG
    def load_exp(self, plotter=None, rabi_index=-1):

        if len(self.xdata) == 0:
            raise NameError("Rabi sequence not initialized")

        if self.config['rabi_vary'] == 'height':
            self._rabi_vs_width = False
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (
            max(self.config['rabi_length']) + max(self.config['pi_pulse_width']))
        elif self.config['rabi_vary'] == 'width':
            self._rabi_vs_width = True
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (
            self.xdata[-1] + max(self.config['pi_pulse_width']))

        if rabi_index == -1:
            rabi_vary = self.xdata
        else:
            rabi_vary = zeros(1)
            rabi_vary[0] = self.xdata[rabi_index]


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #allocate pulse array
        #self.rabi_fullpulsesequence = zeros((max(len(self.rabi_height),len(self.rabi_length)),6,self.config['total_length']*awg_modifier))

        awg_seq = pulse_sequence(len(rabi_vary), self.config['total_length'] * self.config['awg_clocks'][0],
                                 TEK_clk=self.config['awg_clocks'][0], Agilent_clk=self.config['awg_clocks'][1],
                                 delays=self.config['awg_delay'])


        #define card trigger and read trigger once (these are the same waveforms each time)
        for j in range(2):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (max_rabi + self.config['read_delay']), self.config['read_length'])

        if trigger_card_early or self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((max_rabi + self.config['card_delay']),
                                                                               self.config['read_length'])

        for i in range(len(rabi_vary)):


            print rabi_vary[i]

            for j in range(2):

                if not self.config['qubit_enable'][j]:
                    continue

                if self._rabi_vs_width:
                    rabi_pulse_i = floor(rabi_vary[i])
                    pulse_height_i = self.config['rabi_height'][j]
                else:
                    rabi_pulse_i = floor(self.config['rabi_length'][j])
                    pulse_height_i = rabi_vary[i]
                pulse_height_i *= self.config['height_scaling'][j]
                if self.config['shaped_pulse']:

                    if len(rabi_vary) > 1:
                        ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind1
                        ind2a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind2a
                        ind2b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind2b
                    else:
                        ind1 = 0
                        ind2a = 0
                        ind2b = 0

                    if self.config['pulsed_drive'][j]:
                        if j == 1:  #pi
                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(
                                (max_rabi - 2 * rabi_pulse_i - self.config['pi_pulse_width']),
                                2 * self.config['pi_pulse_width'] + 2 * self.config['generator_enable_delay'])
                        else:  #rabi
                            awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(
                                (max_rabi - rabi_pulse_i),
                                (2 * rabi_pulse_i + 2 * self.config['generator_enable_delay']))

                    if self.config['drive_sideband'][j]:
                        awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse_with_freq(
                            (max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                            self.config['drive_sideband_freq'][j] / 1.0e9, 0.0, False)
                        awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse_with_freq(
                            (max_rabi - rabi_pulse_i), rabi_pulse_i / 2,
                            -1 * sign(self.config['drive_sideband_freq'][j]) * pulse_height_i,
                            self.config['drive_sideband_freq'][j] / 1.0e9, pi / 2.0, False)
                    else:
                        if self.config['drag_pulse'][j]:
                            awg_seq.DRAG_pulse((max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                                               self.config['drive_angle'], j, ind2a, ind2b,
                                               self.config['drag_prefactor'][j], False)
                        else:
                            if j == 1:  #pi
                                awg_seq.gauss_pulse((max_rabi - 2 * rabi_pulse_i - self.config['pi_pulse_width']),
                                                    self.config['pi_pulse_width'] / 2, self.config['pi_pulse_height'],
                                                    self.config['drive_angle'], j, ind2a, ind2b, False)
                            else:  #rabi
                                awg_seq.gauss_pulse((max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                                                    self.config['drive_angle'], j, ind2a, ind2b, False)

                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
                    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])


                else:

                    if len(rabi_vary) > 1:
                        ind = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind
                    else:
                        ind = 0

                    awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind].square_pulse2((max_rabi - rabi_pulse_i),
                                                                                              rabi_pulse_i)

            if plotter is not None:
                cur_pulse_seq.plot_pulses(plotter)

        #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #self.rabi_fullpulsesequence[i] = cur_pulse_seq.output_pulses()


    #run the rabi oscillation    
    def run_rabi_sequential(self, num_avgs=1024, plotter=None):

        #define instrument manager    
        im = InstrumentManager()

        #setup awg
        awg = im[self.config['AWG']]

        if self.config['LO'] == self.config['RF']:
            homodyne_meas = True
            print "Homodyne measurement"
        else:
            homodyne_meas = False

        #setup drive

        drive = im[self.config['drive']]

        #setup drive
        if self.config['drive'][0:2] == "LB":
            drive.set_output(True)
            drive.set_pulse_ext(mod=True)
            #Internal pulse must be set to false!
            drive.set_mod(False)
            drive.set_power(self.config['rabi_drive_pwr'])  #-28 when direct driving
            drive.set_frequency(self.config['rabi_drive_freq'])
        else:
            drive.set_power(self.config['rabi_drive_pwr'])
            drive.set_mod(True)
            drive.set_frequency(self.config['rabi_drive_freq'])
            drive.set_ext_pulse()
            drive.set_output(True)

        for i in range(5):
            if abs(drive.get_frequency() - self.config['rabi_drive_freq']) > 100:
                drive.set_frequency(self.config['rabi_drive_freq'])
            else:
                break
        if abs(drive.get_frequency() - self.config['rabi_drive_freq']) > 100:
            print "Lab brick frequency error!"
        time.sleep(.2)

        print self.config['rabi_drive_freq'], drive.get_frequency()

        #setup initial drive (if applicable)
        if self.config['rabi_b_drive_pulse_length'] > 0:
            drive2 = im[self.config['drive_b']]
            #setup drive2 (assuming it's a lab brick)
            drive2.set_output(True)
            drive2.set_pulse_ext(mod=True)
            #Internal pulse must be set to false!
            drive2.set_mod(False)
            drive2.set_power(self.config['rabi_b_drive_pwr'])  #-28 when direct driving
            drive2.set_frequency(self.config['rabi_b_drive_freq'])
            print drive2.get_frequency()

        #setup LO and RF

        RF = im[self.config['RF']]
        RF.set_output(True)
        RF.set_mod(True)
        RF.set_power(self.config['rabi_read_pwr'])
        RF.set_frequency(self.config['rabi_read_freq'])
        RF.set_ext_pulse()

        if not homodyne_meas:
            LO = im[self.config['LO']]
            LO.set_power(self.config['rabi_lo_power'])
            LO.set_mod(True)
            LO.set_frequency(self.config['rabi_read_freq'] + self.config['IFreq'])
            LO.set_output(True)



        #setup the flux
        flux = im[self.config['flux']]

        if self.config['flux'] == 'flux2' or self.config['flux'] == 'flux1':
            flux.set_function("DC")
            flux.set_offset(self.config['rabi_flux'])
        else:
            #flux.set_mode('current')
            flux.set_volt(self.config['rabi_flux'])

        #setup the Alazar card
        card = setup_Alazar_card(512, self.config['acq_length'], self.config['meas_range'],
                                 self.config['monitor_pulses'], True, num_avgs)

        if self._rabi_vs_width:
            xdata = self.rabi_length
        else:
            xdata = self.rabi_height

            #setup plots
        if plotter is not None:
            plotter.init_plot("FullData", rank=2, accum=False)
            plotter.init_plot("Scope", accum=False)
            plotter.init_plot("Amp", accum=False)
            plotter.init_plot("Phase", accum=False)

        self.rabi_amp = zeros(len(xdata))
        self.rabi_phase = zeros(len(xdata))
        self.rabi_fulldata = zeros((len(xdata), self.config['acq_length']))

        for i in range(len(xdata)):

            #program sequence
            awg.pre_experiment()
            self.load_rabi(False, None, i)
            awg.run_experiment(self.config['seq_file'])

            tpts, ch1_pts, ch2_pts = card.acquire_avg_data()
            if not homodyne_meas:
                amp, phase, _, _ = heterodyne(tpts, ch1_pts, ch2_pts, self.config['IFreq'])
            else:
                amp = mean(ch1_pts)
                phase = amp

            self.rabi_phase[i] = phase
            self.rabi_amp[i] = amp
            self.rabi_fulldata[i] = ch1_pts

            if plotter is not None:
                plotter.plot((xdata[0:(i + 1)], self.rabi_phase[0:(i + 1)]), "Phase")
                plotter.plot((xdata[0:(i + 1)], self.rabi_amp[0:(i + 1)]), "Amp")
                plotter.plot((range(self.config['acq_length']), ch1_pts), "Scope")
                plotter.plot(self.rabi_fulldata, "FullData")


    def fit_rabi(self):

        #start the fit after fit_start
        if len(self.rabi_height) > 1:
            xdata = self.rabi_height
        else:
            xdata = self.rabi_length

        start_index = -1
        for i in range(len(xdata)):
            if xdata[i] >= self.config['fit_start']:
                start_index = i
                break

        if start_index == -1:
            start_index = 0

        #fit the t1 data to an decaying sinusoid 
        if self.config['fit_phase']:
            ydata = self.rabi_phase[start_index:-1]
            xdata = xdata[start_index:-1]

        else:
            ydata = self.rabi_amp[start_index:-1]
            xdata = xdata[start_index:-1]

        self.config['fit_start'] = xdata[0]

        fitvals = fitdecaysin(xdata, ydata)

        self.config['tau'] = fitvals[3]
        self.config['A'] = fitvals[0]
        self.config['y0'] = fitvals[4]
        self.config['omega'] = fitvals[1]
        self.config['phase'] = fitvals[2]

        print fitvals

    def display_rabi(self, plotter, pulse_sequence=False, display_fit=True):


        plotter.init_plot("FullData", rank=2, accum=False)
        plotter.init_plot("Amp", rank=1, accum=False)
        plotter.init_plot("Phase", rank=1, accum=False)

        plotter.plot((self.rabi_fulldata), "FullData")

        p = [self.config['A'], self.config['omega'], self.config['phase'], self.config['tau'], self.config['y0'],
             self.config['fit_start']]

        if len(self.rabi_height) > 1:
            xdata = self.rabi_height
        else:
            xdata = self.rabi_length

        if self.config['fit_phase'] and display_fit:
            plotter.plot((concatenate((xdata, xdata)), concatenate((self.rabi_phase, decaysin(p, xdata)))), "Phase")
        else:
            plotter.plot((xdata, self.rabi_phase), "Phase")

        time.sleep(.2)

        if (not self.config['fit_phase']) and display_fit:
            plotter.plot((concatenate((xdata, xdata)), concatenate((self.rabi_amp, decaysin(p, xdata)))), "Amp")
        else:
            plotter.plot((xdata, self.rabi_amp * 1.0), "Amp")


class rabi_enhanced_readout(qubit_exp):
    def __init__(self):

        qubit_exp.__init__(self)

        self.config['exp_id'] = 'RABI_ENHANCED'

        #custom Rabi experiment settings
        self.config['rabi_height'] = [0.3, 0.3]
        self.config['rabi_length'] = [10, 10]
        self.config['rabi_vary'] = 'height'  #height or width ..then the xdata is the vary parameter
        self.config['drive_angle'] = 0.0
        self.config['height_scaling'] = [1.0, 1.0]
        self.config['sideband_freq'] = [0.0, 0.0]
        self.config['sideband_height'] = [0.0, 0.0]


    #load the rabi oscillation into the AWG
    def load_exp(self, plotter=None, rabi_index=-1):

        if len(self.xdata) == 0:
            raise NameError("Rabi sequence not initialized")

        if self.config['rabi_vary'] == 'height':
            self._rabi_vs_width = False
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (
            max(self.config['rabi_length']))
        elif self.config['rabi_vary'] == 'width':
            self._rabi_vs_width = True
            max_rabi = self.config['front_buffer'] + self.config['generator_enable_delay'] + 2 * (self.xdata[-1])

        if rabi_index == -1:
            rabi_vary = self.xdata
        else:
            rabi_vary = zeros(1)
            rabi_vary[0] = self.xdata[rabi_index]


        #Trigger the card early if you want to measure the drive pulses directly with the Alazar    
        trigger_card_early = False

        #allocate pulse array
        #self.rabi_fullpulsesequence = zeros((max(len(self.rabi_height),len(self.rabi_length)),6,self.config['total_length']*awg_modifier))

        awg_seq = pulse_sequence(len(rabi_vary), self.config['total_length'] * self.config['awg_clocks'][0],
                                 self.config['total_length'] * self.config['awg_clocks'][1],
                                 TEK_clk=self.config['awg_clocks'][0], Agilent_clk=self.config['awg_clocks'][1],
                                 delays=self.config['awg_delay'])

        awg_seq.marker[awg_seq.channel_index['flux_trig']][0].square_pulse(1, 100)
        #define card trigger and read trigger once (these are the same waveforms each time)
        for j in range(2):
            if self.config['qubit_enable'][j]:
                awg_seq.marker[awg_seq.channel_index['read_trig'][j]][0].square_pulse(
                    (max_rabi + self.config['read_delay']), self.config['read_length'])

        if trigger_card_early or self.config['monitor_pulses']:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse(1, self.config['read_length'])
        else:
            awg_seq.marker[awg_seq.channel_index['card_trig']][0].square_pulse((max_rabi + self.config['card_delay']),
                                                                               self.config['read_length'])

        for i in range(len(rabi_vary)):


            print rabi_vary[i]

            for j in range(2):

                if not self.config['qubit_enable'][j]:
                    continue

                if self._rabi_vs_width:
                    rabi_pulse_i = floor(rabi_vary[i])
                    pulse_height_i = self.config['rabi_height'][j]
                else:
                    rabi_pulse_i = floor(self.config['rabi_length'][j])
                    pulse_height_i = rabi_vary[i]
                pulse_height_i *= self.config['height_scaling'][j]
                if self.config['shaped_pulse']:

                    if len(rabi_vary) > 1:
                        ind1 = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind1
                        ind2a = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_I'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['drive_I'][j]][i] = ind2a
                        ind2b = awg_seq.add_analog_pulse(awg_seq.channel_index['drive_Q'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['drive_Q'][j]][i] = ind2b

                        indflux = awg_seq.add_analog_pulse(awg_seq.channel_index['flux'][j])
                        awg_seq.analog_seq_table[awg_seq.channel_index['flux'][j]][i] = indflux

                    else:
                        ind1 = 0
                        ind2a = 0
                        ind2b = 0

                        indflux = 0

                    if self.config['pulsed_drive'][j]:
                        awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind1].square_pulse2(
                            (max_rabi - rabi_pulse_i), (2 * rabi_pulse_i + 2 * self.config['generator_enable_delay']))

                    if self.config['drive_sideband'][j]:
                        awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].gauss_pulse_with_freq(
                            (max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                            self.config['drive_sideband_freq'][j] / 1.0e9, 0.0, False)
                        awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].gauss_pulse_with_freq(
                            (max_rabi - rabi_pulse_i), rabi_pulse_i / 2,
                            -1 * sign(self.config['drive_sideband_freq'][j]) * pulse_height_i,
                            self.config['drive_sideband_freq'][j] / 1.0e9, pi / 2.0, False)
                    else:
                        if self.config['drag_pulse'][j]:
                            awg_seq.DRAG_pulse((max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                                               self.config['drive_angle'], j, ind2a, ind2b,
                                               self.config['drag_prefactor'][j], False)
                        else:
                            awg_seq.gauss_pulse((max_rabi - rabi_pulse_i), rabi_pulse_i / 2, pulse_height_i,
                                                self.config['drive_angle'], j, ind2a, ind2b, False)

                    awg_seq.analogwf[awg_seq.channel_index['drive_I'][j]][ind2a].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_I'][j]])
                    awg_seq.analogwf[awg_seq.channel_index['drive_Q'][j]][ind2b].offset(
                        self.config['awg_pulse_offset'][awg_seq.channel_index['drive_Q'][j]])

                    awg_seq.analogwf[awg_seq.channel_index['flux'][j]][indflux].gauss_pulse_with_freq(
                        (max_rabi + self.config['read_delay'] + self.config['read_length'] / 2),
                        self.config['read_length'] / 4, self.config['sideband_height'][j],
                        self.config['sideband_freq'][j])



                else:

                    if len(rabi_vary) > 1:
                        ind = awg_seq.add_marker_pulse(awg_seq.channel_index['drive_trig'][j])
                        awg_seq.marker_seq_table[awg_seq.channel_index['drive_trig'][j]][i] = ind
                    else:
                        ind = 0

                    awg_seq.marker[awg_seq.channel_index['drive_trig'][j]][ind].square_pulse2((max_rabi - rabi_pulse_i),
                                                                                              rabi_pulse_i)

            if plotter is not None:
                cur_pulse_seq.plot_pulses(plotter)

        #load pulses into the awg
        awg_seq.load_into_awg(self.config['seq_file'], self.config['AWG'][0])

        #load pulses into the AGILENT AWG
        if self.config['load_agilent']:
            awg_seq.load_full_into_agilent(self.config['AWG'][1])

            #self.rabi_fullpulsesequence[i] = cur_pulse_seq.output_pulses()
