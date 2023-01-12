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
    tx_warmup = 200 # already handled by delay in RF block
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
    od, pd = psi.interpret("tabletop_se_2d_pulseq.seq")   

    grad_interval = pd['grad_t']
    TR = pd['TR']
    Nx = int(pd['Nx'])
    Ny = int(pd['Ny'])
    delayTR = pd['delayTR']
    oversampling_factor = pd['oversampling_factor']

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_interval/3,
                         flush_old_rx=True,                         
                         halt_and_reset=True) 
    expt.add_flodict(od)

    expt.gradb.calibrate(channels=[0, 1, 2], max_current=max_grad_current, num_calibration_points=30, averages=5, poly_degree=5)

    adc_pad = 0
    nSamples = int(Nx * oversampling_factor - adc_pad)
    data2d = np.zeros((Ny, nSamples),dtype=np.complex64)
    import copy as cp
    for k, phase_factor in enumerate(np.linspace(-1, 1, Ny)):
        print('phase step {:d}\n'.format(k))
        phase_x, phase_y = od['grad_vy']
        phase_y_new = phase_y * phase_factor
        expt.add_flodict({'grad_vy':(phase_x,phase_y_new)},append=False)
        rxd, msgs = expt.run()
        data = rxd['rx0'][adc_pad:]
        if False:
            fig, (ax1, ax2, ax3) = plt.subplots(3)
            Noise = np.abs(np.std(np.real(np.fft.fft(data))[int(data.size/2)-3:int(data.size/2)+3]))
            SNR=np.max(np.abs(np.fft.fft(data)))/Noise
            dt = params['rx_t']
            fig.suptitle('Spin Echo [n={:d}, lo_freq={:f} Mhz]\nSNR={:f}'.format(nSamples,lo_freq,SNR))
            t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
            ax1.plot(t_axis, np.abs(data)*3.3)
            ax1.set_ylabel('voltage [V]')
            ax2.set_xlabel('time [us]')
            ax2.plot(t_axis, data.real*3.3)
            ax2.set_ylabel('voltage [V]')
            f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))
            ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))/np.sqrt(nSamples)))
            plt.show()
            fig.tight_layout()
        data2d[k, :] = data
        if True:
            if k == 0:
                plt.figure(1)
                plt.subplot(1, 3, 1)
                im = plt.imshow(10*np.log(np.abs(data2d)),aspect='auto',interpolation='none', origin='lower')
                plt.ion()
                plt.show()
                plt.pause(0.001)
            else:
                im.set_data(10*np.log(np.abs(data2d)))
                plt.pause(0.001)

        time.sleep(delayTR)

    from datetime import datetime
    now = datetime.now()
    current_time = now.strftime("%y-%d-%m %H_%M_%S")
    filename = f"data2d se {current_time} Nx {Nx} Ny {Ny} TR {TR}.npy"
    if os.path.exists(filename):
        os.remove(filename)
    np.save(filename,data2d)
    plt.close()
    plt.ioff()

    # reconstruction
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

  

 
