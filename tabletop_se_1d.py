#!/usr/bin/env python3
#
# loopback test using ocra-pulseq
#

import numpy as np
import matplotlib.pyplot as plt
import pdb

import external
import experiment as ex
from pulseq_assembler import PSAssembler
st = pdb.set_trace

if __name__ == "__main__":
    lo_freq = 17.315 # MHz
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
    hf_B_per_m_current = 2.483E-4 # [T/A] theoretical value

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
    hf_max_Hz_per_m = 4166 # experimental value
    print(hf_max_Hz_per_m)
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    grad_max = grad_max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	
    rf_amp_max = hf_max_Hz_per_m # factor used to normalize RF amplitude, should be max value of system used!
    tx_warmup = 0 # already handled by delay in RF block
    adc_pad = 0 # padding to prevent junk in rx buffer
    ps = PSAssembler(rf_center=lo_freq*1e6,
        # how many Hz the max amplitude of the RF will produce; i.e. smaller causes bigger RF V to compensate
        rf_amp_max=rf_amp_max,
        grad_max=grad_max,
        clk_t=clk_t,
        tx_t=tx_t,
        grad_t=grad_interval,
        tx_warmup=tx_warmup,
        adc_pad=adc_pad,
		rf_delay_preload=True)
    tx_arr, grad_arr, cb, params = ps.assemble('tabletop_se_1d_pulseq.seq')

    # Temporary hack, until next ocra-pulseq update
    if 'rx_t' not in params:
        params['rx_t'] = rx_t    

    exp = ex.Experiment(samples=params['readout_number'], 
        lo_freq=lo_freq,
        tx_t=tx_t,
		rx_t=params['rx_t'],
        grad_channels=num_grad_channels,
        grad_t=grad_interval/num_grad_channels,
		acq_retry_limit=500000,
        assert_errors=False)
        
    exp.define_instructions(cb)
    exp.add_tx(ps.tx_arr)
    exp.add_grad(ps.grad_arr)

    # plt.plot(ps.gr_arr[0]);plt.show()

    exp.calibrate_gpa_fhdo(max_current = 4,
        num_calibration_points=10,
        gpa_current_per_volt=gpa_current_per_volt) 

    # set all channels back to 0 A
    for ch in range(num_grad_channels):
        dac_code = exp.ampere_to_dac_code(0)
        dac_code = exp.calculate_corrected_dac_code(ch,dac_code)
        exp.write_gpa_dac(ch, dac_code)      


    data = exp.run() # Comment out this line to avoid running on the hardware
    # set all channels back to 0 A
    for ch in range(num_grad_channels):
        dac_code = exp.ampere_to_dac_code(0)
        dac_code = exp.calculate_corrected_dac_code(ch,dac_code)
        exp.write_gpa_dac(ch, dac_code)  

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Spin Echo [n={:d}, lo_freq={:f} Mhz]'.format(params['readout_number'],lo_freq))
    dt = params['rx_t']
    nSamples = params['readout_number']
    t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
    ax1.plot(t_axis, np.abs(data)*3.3)
    ax1.set_ylabel('voltage [V]')
    ax2.set_xlabel('time [us]')
    ax2.plot(t_axis, data.real*3.3)
    ax2.set_ylabel('voltage [V]')
    #f_axis = np.linspace(-1/dt*nSamples,1/dt*nSamples,nSamples)
    nFFT_window = 127
    f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))[int(nSamples/2)-nFFT_window:int(nSamples/2)+nFFT_window]
    ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))[int(nSamples/2)-nFFT_window:int(nSamples/2)+nFFT_window]/np.sqrt(nSamples)))
    plt.show()
    fig.tight_layout()

    # st()    

  

 
