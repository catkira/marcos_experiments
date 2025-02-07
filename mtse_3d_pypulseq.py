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
    adc_pad = 85 # padding to prevent junk in rx buffer        
    num_grad_channels = 3

    gamma = 42570000 # Hz/T

    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    tx_warmup = 0 # already handled by delay in RF block
    psi = flocra_pulseq.interpreter.PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=hf_max_Hz_per_m,
                        tx_t=tx_t,
                        tx_warmup = tx_warmup,
                        grad_max=grad_max_Hz_per_m)
    # od, pd = psi.interpret("tabletop_se_2d_pulseq.seq")
    od, pd = psi.interpret("2d_tse_ptb.seq")
    # od, pd = psi.interpret("mtse_3d_pypulseq.seq")

    grad_interval = pd['grad_t']
    # TR = pd['TR']
    # Nx = int(pd['Nx'])
    # Ny = int(pd['Ny'])
    # delayTR = pd['delayTR']
    # oversampling_factor = pd['oversampling_factor']

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                        #  gpa_fhdo_offset_time=grad_interval/3,
                         flush_old_rx=True,                         
                         halt_and_reset=True) 
    expt.add_flodict(od)

    # expt.plot_sequence()
    # plt.show()

    rxd, msgs = expt.run()
    # data = rxd['rx0'][adc_pad:]


    # reconstruction
    if False:
        oversampling_factor = int(np.round(data2d.shape[1]/Nx))
        data2d = sig.decimate(data2d, oversampling_factor, ftype='iir', axis=1)    
        plt.figure(1)
        plt.subplot(1, 3, 1)
        plt.imshow(10*np.log(np.abs(data2d)),aspect='auto',interpolation='none', origin='lower')
        plt.subplot(1, 3, 2)
        plt.imshow(np.angle(data2d),aspect='auto',interpolation='none')
        plt.subplot(1, 3, 3)
        img = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(data2d))))
        plt.imshow(img, aspect='auto',cmap='gray',interpolation='none')
        plt.show()
    
    # st()    

  

 
