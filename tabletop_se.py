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
import flocra_pulseq.interpreter
st = pdb.set_trace
from mri_config import lo_freq, grad_max_Hz_per_m, hf_max_Hz_per_m, gamma, shim

if __name__ == "__main__":  
    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    #tx_warmup = 0 # already handled by delay in RF block
    psi = flocra_pulseq.interpreter.PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=rf_amp_max,
                        grad_max=grad_max,
                        tx_warmup=200)
    od, pd = psi.interpret("tabletop_se_pulseq.seq")         
    grad_interval = pd['grad_t']

    if True:
        # Shim
        grads = ['grad_vx', 'grad_vy', 'grad_vz', 'grad_vz2']
        for ch in range(4):
            od[grads[ch]] = (od[grads[ch]][0], od[grads[ch]][1] + shim[ch])

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_interval/3,
                         flush_old_rx=True,
                         halt_and_reset=True) 
    expt.add_flodict(od)
    expt.gradb.calibrate(channels=[0,1,2], max_current=6, num_calibration_points=30, averages=5, poly_degree=5)

    rxd, msgs = expt.run()
    expt.gradb.init_hw()  # set gradient currents back to zero


    data = rxd['rx0']
    data = data[6:]
    nSamples_orig = pd['readout_number']   
    nSamples = len(data)
    Noise = np.abs(np.std(np.real(np.fft.fft(data))[int(data.size/2)-3:int(data.size/2)+3]))
    SNR=np.max(np.abs(np.fft.fft(data)))/Noise
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Spin Echo [n={:d}, lo_freq={:f} Mhz]\nSNR={:f}'.format(nSamples_orig,lo_freq,SNR))
    dt = pd['rx_t']
    t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
    ax1.plot(t_axis, np.abs(data)*3.3)
    ax1.set_ylabel('abs [mV]')
    ax2.set_xlabel('time [us]')
    ax2.plot(t_axis, data.real*3.3)
    ax2.set_ylabel('real [mV]')
    #f_axis = np.linspace(-1/dt*nSamples,1/dt*nSamples,nSamples)
    #nFFT_window = 127
    #f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))[int(nSamples/2)-nFFT_window:int(nSamples/2)+nFFT_window]
    #ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))[int(nSamples/2)-nFFT_window:int(nSamples/2)+nFFT_window]/np.sqrt(nSamples)))
    f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))
    ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))/np.sqrt(nSamples)))
    ax3.set_ylabel('Spectrum')
    ax3.set_xlabel('Hz')
    plt.show()
    fig.tight_layout()
    expt.close_server(True) 
    # st()    

  

 
