# -*- coding: utf-8 -*-
"""
Created on Mon July 08 23:34:06 2013
@author: Dave McKay, Ge Yang
"""
from qubit import *
from pulse_sequences import *

if __name__ == "__main__":
    expt_path = r'S:\_Data\140707 - Anthony ESR\data'
    prefix = 'esr'

    alazarConfig = {'clock_edge': 'rising', 'trigger_delay': 0,
                    'ch1_filter': False, 'ch1_enabled': True,
                    'samplesPerRecord': 49024,
                    'recordsPerBuffer': 200,
                    'recordsPerAcquisition': 20000,
                    # 'samples_per_buffer': 2000064,
                    # 'samples_per_record': 2000064,
                    # 'seconds_per_buffer': 1.0e-6,
                    # 'seconds_per_record': 1.0e-6,
                    # 'records_per_buffer': 1,
                    # 'records_per_acquisition': 1,
                    # 'bytes_per_buffer': 4000128,
                    'bufferCount': 1, 'trigger_edge1': 'rising', 'trigger_edge2': 'rising',
                    'ch2_range': 1.0, 'clock_source': 'reference', 'trigger_level2': .5, 'trigger_level1': 0.5,
                    'ch2_coupling': 'DC', 'trigger_coupling': 'DC', 'ch2_filter': False, 'trigger_operation': 'or',
                    'ch1_coupling': 'DC', 'trigger_source2': 'disabled', 'trigger_source1': 'external',
                    'sample_rate': 50000, 'timeout': 30000, 'ch1_range': 1.0,
                    'ch2_enabled': True}

    esr = esrExperiment(expt_path, prefix, alazarConfig)

    # seq_file = exp_path + "sequence_files\\"
    # seq_file = seq_file + "esr.awg"
    # esr._config['seq_file'] = seq_file

    esr._config['num_avgs'] = 1000
    esr._config['tau_start'] = 10
    esr._config['tau_end'] = 10
    esr._config['tau_pts'] = 1

    esr._config['drive'][0] = 'RF_1'
    esr._config['AWG'][1] = 'AWG'

    esr._config['awg_offset'][1][0] = 0.02
    esr._config['awg_offset'][1][1] = -0.02

    esr.take_esr_data(ScriptPlotter())
    #
    # s = raw_input("Save Data (Y/N)?")
    #
    # # esr_exp.rf.set_out
    #
    # if s.upper() != "N":
    #     print "Saving Data"
    #     esr.save_data(exp_path, esr._config['exp_id'], overwrite=False)
