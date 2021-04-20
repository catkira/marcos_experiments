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
from flocra_pulseq_interpreter import PSInterpreter
st = pdb.set_trace

if __name__ == "__main__":
    lo_freq = 17.286# MHz
    tx_t = 1 # us
    num_grad_channels = 3
    grad_t = 10 # us between [num_grad_channels] channel updates

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
    hf_max_Hz_per_m = 4200 # experimental value
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    #tx_warmup = 0 # already handled by delay in RF block
    psi = PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=rf_amp_max,
                        tx_t=tx_t,
                        grad_t=grad_t,
                        grad_max=grad_max)
    od, pd = psi.interpret("tabletop_se_v2_2d_pulseq.seq")   

    TR = pd['TR']
    Nx = pd['Nx']
    Ny = int(pd['Ny'])

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_t/3) 
    expt.add_flodict(od)

    expt.gradb.calibrate(channels=[0,1], max_current=6, num_calibration_points=30, averages=5, poly_degree=5)

    rxd, msgs = expt.run()
    data2d = rxd['rx0'].reshape(200,Ny)
    adc_pad = 0
    #data2d=data2d[adc_pad:][:]
    nSamples = pd['readout_number'] - adc_pad
    
    from datetime import datetime
    now = datetime.now()
    current_time = now.strftime("%y-%d-%m %H_%M_%S")
    filename = f"data2d se v2 Nx {nSamples} Ny {Ny} TR {TR} {current_time}.npy"
    if os.path.exists(filename):
        os.remove(filename)
    np.save(filename,data2d)
    plt.close()
    plt.ioff()
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

  

 
