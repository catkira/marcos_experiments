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
#from pulseq_assembler import PSAssembler
from flocra_pulseq_interpreter import PSInterpreter
st = pdb.set_trace

if __name__ == "__main__":
    lo_freq = 17.301 # MHz
    #lo_freq = 1 # MHz  # only for debugging on scope
    tx_t = 1.001 # us
    rx_t = 0.497
    clk_t = 0.007
    num_grad_channels = 3
    grad_interval = 10.003 # us between [num_grad_channels] channel updates

    gamma = 42570000 # Hz/T

    # value for tabletopMRI  gradient coil
    grad_B_per_m_per_current = 0.02 # [T/m/A], approximate value for tabletop gradient coil
    R_coil = 2

    # value for tabletopMRI hf coil
    hf_B_per_m_current = 2.483E-3 # [T/A] theoretical value

    # values for gpa fhdo
    gpa_current_per_volt = 2.5 # gpa fhdo 6A configuration
    max_dac_voltage = 2.5

    # values for red pitaya 
    hf_max_dac_voltage = 1 # +-

    # HF-PA
    hf_PA_gain = 20 # dB

    grad_max_Hz_per_m = max_dac_voltage * gpa_current_per_volt * grad_B_per_m_per_current * gamma	
    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))

    #hf_max_Hz_per_m = np.sqrt(1/50 * 10**(hf_PA_gain/10) / R_coil) * hf_B_per_m_current * gamma
    hf_max_Hz_per_m = 4200 # experimental value
    print(hf_max_Hz_per_m)
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    tx_warmup = 0 # already handled by delay in RF block
    adc_pad = 85 # padding to prevent junk in rx buffer
    psi = PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=rf_amp_max,
                        rf_ends_zero=True,
                        tx_t=tx_t,
                        grad_t=grad_interval,
                        grad_max=grad_max) # very large, just for testing
    od, pd = psi.interpret("tabletop_se_pulseq.seq")         
    od['grad_vy'] = od['grad_vy'][0] + grad_interval/3, od['grad_vy'][1]
    od['grad_vz'] = od['grad_vz'][0] + 2*grad_interval/3, od['grad_vz'][1]       
    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True) 
    expt.add_flodict(od)

    rxd, msgs = expt.run()


    data = rxd['rx0']
    data = data[6:]
    #nSamples = pd['readout_number']   
    nSamples = len(data)
    Noise = np.abs(np.std(np.real(np.fft.fft(data))[int(data.size/2)-3:int(data.size/2)+3]))
    SNR=np.max(np.abs(np.fft.fft(data)))/Noise
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Spin Echo [n={:d}, lo_freq={:f} Mhz]\nSNR={:f}'.format(nSamples,lo_freq,SNR))
    dt = pd['rx_t']
    t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
    ax1.plot(t_axis, np.abs(data)*3.3)
    ax1.set_ylabel('voltage [V]')
    ax2.set_xlabel('time [us]')
    ax2.plot(t_axis, data.real*3.3)
    ax2.set_ylabel('voltage [V]')
    #f_axis = np.linspace(-1/dt*nSamples,1/dt*nSamples,nSamples)
    #nFFT_window = 127
    #f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))[int(nSamples/2)-nFFT_window:int(nSamples/2)+nFFT_window]
    #ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))[int(nSamples/2)-nFFT_window:int(nSamples/2)+nFFT_window]/np.sqrt(nSamples)))
    f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))
    ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))/np.sqrt(nSamples)))
    plt.show()
    fig.tight_layout()
    expt.close_server(True) 
    # st()    

  

 
