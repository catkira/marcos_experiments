#!/usr/bin/env python3
#
# loopback test using ocra-pulseq
#

import numpy as np
import matplotlib.pyplot as plt
import pdb
import time

import external
import experiment as ex
import os
import scipy.io as sio
from shutil import copyfile
from flocra_pulseq_interpreter import PSInterpreter
st = pdb.set_trace

if __name__ == "__main__":
    lo_freq = 17.283 # MHz
    tx_t = 1 # us
    grad_t = 10 #

    gamma = 42570000 # Hz/T

    # value for tabletopMRI  gradient coil
    grad_B_per_m_per_current = 0.02 # [T/m/A], approximate value for tabletop gradient coil
    R_coil = 2

    # value for tabletopMRI hf coil
    hf_B_per_m_current = 2.483E-4 # [T/A] theoretical value

    # values for gpa fhdo
    gpa_current_per_volt = 2.5 # gpa fhdo 6A configuration
    max_dac_voltage = 2.5

    # values for red pitaya 
    hf_max_dac_voltage = 1 # +-

    # HF-PA
    hf_PA_gain = 20 # dB

    #grad_max_Hz_per_m = max_dac_voltage * gpa_current_per_volt * grad_B_per_m_per_current * gamma	
    grad_max_Hz_per_m = 15E6 # experimental value
    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))

    #hf_max_Hz_per_m = np.sqrt(1/50 * 10**(hf_PA_gain/10) / R_coil) * hf_B_per_m_current * gamma
    hf_max_Hz_per_m = 3400 # experimental value, needs to take into account coil tuning! (do flip-angle calibration)
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    #tx_warmup = 0 # already handled by delay in RF block
    seq_file = "tabletop_radial_v2_2d_pulseq.seq"

    psi = PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=rf_amp_max,
                        tx_t=tx_t,
                        grad_t=grad_t,
                        grad_max=grad_max)
    od, pd = psi.interpret(seq_file)   

    TR = pd['TR']
    Nx = pd['Nx']
    Nspokes = int(pd['Nspokes'])
    sliceThickness = int(pd['SliceThickness'])
    angle = int(pd['angle'])

    from datetime import datetime
    now = datetime.now()
    current_time = now.strftime("%y-%m-%d %H_%M_%S")
    filename = f"data2d radial v2 {current_time} Nspokes {Nspokes} Nx {Nx} TR {TR} angle {angle} SliceThickness {sliceThickness}"
    copyfile(seq_file,filename+".seq")        

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_t/3) 
    expt.add_flodict(od)

    expt.gradb.calibrate(channels=[0,1,2], max_current=6, num_calibration_points=30, averages=5, poly_degree=5)

    rxd, msgs = expt.run()
    nSamples = pd['readout_number']
    
    np.save(filename + ".npy",rxd['rx0'])
    sio.savemat(filename + ".mat",rxd)
    plt.close()
    plt.ioff()


  

 
