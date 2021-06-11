#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pdb

import external
import experiment as ex
import os
import time
#from pulseq_assembler import PSAssembler
import flocra_pulseq.interpreter
st = pdb.set_trace
from mri_config import lo_freq, grad_max_Hz_per_m, hf_max_Hz_per_m, gamma, shim

if __name__ == "__main__":  
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    psi = flocra_pulseq.interpreter.PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=rf_amp_max,
                        grad_max=grad_max) # very large, just for testing
    od, pd = psi.interpret("tabletop_noise_monitor.seq")         
    grad_interval = pd['grad_t']


    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_interval/3,
                         flush_old_rx=True,
                         halt_and_reset=True) 
    expt.add_flodict(od)
    expt.gradb.calibrate(channels=[0,1,2], max_current=6, num_calibration_points=30, averages=5, poly_degree=5)

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    while True:
        rxd, msgs = expt.run()
        data = rxd['rx0']
        data = data[6:]
        nSamples_orig = pd['readout_number']   
        nSamples = len(data)
        dt = pd['rx_t']
        t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
        noise = np.abs(np.std(np.real(np.fft.fft(data))))
        fig.suptitle('Noise Monitor [n={:d}, lo_freq={:f} Mhz]\nNoise={:f}'.format(nSamples_orig,lo_freq,noise))
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax1.plot(t_axis, np.abs(data)*15)
        ax1.set_ylabel('abs [mV]')   # 500mVpp is max RP122 can receive
        ax2.set_xlabel('time [us]')
        ax2.plot(t_axis, data.real*15)
        ax2.set_ylabel('real [mV]')
        f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))
        ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))/np.sqrt(nSamples)))
        plt.ion()
        plt.show()
        plt.pause(0.001)        
        fig.tight_layout()
        time.sleep(1)
    expt.gradb.init_hw()  # set gradient currents back to zero
    expt.close_server(True) 
    # st()    

  

 
