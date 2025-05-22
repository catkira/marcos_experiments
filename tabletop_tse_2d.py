#!/usr/bin/env python3
#
# loopback test using ocra-pulseq
#

import numpy as np
import matplotlib.pyplot as plt
import pdb
import time
import scipy.signal as sig

import external
import experiment as ex
import os
import flocra_pulseq.interpreter
from mri_config import lo_freq, grad_max_Hz_per_m, hf_max_Hz_per_m, gamma, shim, max_grad_current
st = pdb.set_trace

if __name__ == "__main__":
    tx_t = 1 # us
    tx_warmup = 20 # already handled by delay in RF block
    adc_pad = 85 # padding to prevent junk in rx buffer        
    num_grad_channels = 3

    gamma = 42570000 # Hz/T

    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    #tx_warmup = 0 # already handled by delay in RF block
    psi = flocra_pulseq.interpreter.PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=hf_max_Hz_per_m,
                        tx_t=tx_t,
                        tx_warmup = tx_warmup,
                        grad_max=grad_max_Hz_per_m)
    od, pd = psi.interpret("tse_2d.seq")   

    # grad_interval = pd['grad_t']

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                        #  gpa_fhdo_offset_time=grad_interval/3,
                         flush_old_rx=True,                         
                         halt_and_reset=True) 
    expt.add_flodict(od)

    # expt.gradb.calibrate(channels=[0, 1, 2], max_current=max_grad_current, num_calibration_points=30, averages=5, poly_degree=5)

    expt.plot_sequence()
    plt.show()

    rxd, msgs = expt.run()
    data = rxd['rx0']

    # st()    

  

 
