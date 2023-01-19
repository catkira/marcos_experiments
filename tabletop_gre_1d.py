#!/usr/bin/env python3
#
# loopback test using ocra-pulseq
#

import numpy as np
import matplotlib.pyplot as plt
import pdb
import time
from mri_config import lo_freq, grad_max_Hz_per_m, hf_max_Hz_per_m, gamma, shim, max_grad_current


import external
import experiment as ex
import os
import flocra_pulseq.interpreter
st = pdb.set_trace
from mri_config import lo_freq, grad_max_Hz_per_m, hf_max_Hz_per_m, gamma, shim

if __name__ == "__main__":
    tx_t = 1 # us
    num_grad_channels = 3

    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    adc_pad = 85 # padding to prevent junk in rx buffer
    tx_warmup = 200 # already handled by delay in RF block
    psi = flocra_pulseq.interpreter.PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=hf_max_Hz_per_m,
                        tx_t=tx_t,
                        tx_warmup = tx_warmup,
                        grad_max=grad_max_Hz_per_m)
    od, pd = psi.interpret("tabletop_gre_1d_pulseq.seq")
    grad_interval = pd['grad_t']
     

    if True:
        # Shim
        grads = ['grad_vx', 'grad_vy', 'grad_vz']
        for ch in range(3):
            od[grads[ch]] = (od[grads[ch]][0], od[grads[ch]][1] + shim[ch])

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_interval/3,
                         flush_old_rx=True,                         
                         halt_and_reset=True) 
    expt.add_flodict(od)

    expt.gradb.calibrate(channels=[0, 1, 2], max_current=6, num_calibration_points=30, averages=5, poly_degree=5)

    rxd, msgs = expt.run()
    expt.gradb.init_hw()  # set gradient currents back to zero

    # from datetime import datetime
    # now = datetime.now()
    # current_time = now.strftime("%y-%d-%m %H_%M_%S")
    # filename = f"data1d ben Nx {nSamples} {current_time}.npz"
    # if os.path.exists(filename):
    #     os.remove(filename)
    # np.savez(filename,data=data,dt=dt,nSamples=int(nSamples),lo_freq=lo_freq,data1d=data)

    data = rxd['rx0'][6:]
    nSamples = len(data)
    dt = pd['rx_t']    
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Gradient Echo [n={:d}, lo_freq={:f} Mhz]\n'.format(nSamples,lo_freq))
    t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
    ax1.plot(t_axis, np.abs(data)*3.3)
    ax1.set_ylabel('abs [mV]')
    ax2.set_xlabel('time [us]')
    ax2.plot(t_axis, data.real*3.3)
    ax2.set_ylabel('real [mV]')
    #f_axis = np.linspace(-1/dt*nSamples,1/dt*nSamples,nSamples)
    #nFFT_window = 60
    #f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))[int(nSamples/2)-nFFT_window:int(nSamples/2)+nFFT_window]
    #ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))[int(nSamples/2)-nFFT_window:int(nSamples/2)+nFFT_window]/np.sqrt(nSamples)))
    f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))
    ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))/np.sqrt(nSamples)))
    plt.show()
    fig.tight_layout()

    # st()    

  

 
