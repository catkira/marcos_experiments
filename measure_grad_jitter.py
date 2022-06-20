#!/usr/bin/env python3
#
# loopback test using ocra-pulseq
#

import numpy as np
import matplotlib.pyplot as plt
import pdb

import external
import experiment as ex
import os
import time
import flocra_pulseq.interpreter
st = pdb.set_trace
from mri_config import lo_freq, grad_max_Hz_per_m, grad_max_x_Hz_per_m, grad_max_y_Hz_per_m, grad_max_z_Hz_per_m, hf_max_Hz_per_m, gamma, shim, max_grad_current 

if __name__ == "__main__":
    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))    
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    tx_warmup = 200 # already handled by delay in RF block
    psi = flocra_pulseq.interpreter.PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=rf_amp_max,
                        gx_max=grad_max_x_Hz_per_m,
                        gy_max=grad_max_y_Hz_per_m,
                        gz_max=grad_max_z_Hz_per_m,
                        tx_warmup=tx_warmup)
    od, pd = psi.interpret("measure_grad_jitter.seq")
    grad_interval = pd['grad_t']

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_interval/3,
                         flush_old_rx=True,
                         halt_and_reset=True) 
    for i in range(100):
        expt.add_flodict(od)
        #expt.gradb.calibrate(channels=[0,1,2,3], max_current=max_grad_current, num_calibration_points=30, averages=5, poly_degree=5)  
        rxd, msgs = expt.run()
        print("run " + str(i))

    expt.close_server(True) 
    # st()    

  

 
