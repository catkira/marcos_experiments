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
    lo_freq = 17.31 # MHz
    tx_t = 1.001 # us
    rx_t = 0.497
    clk_t = 0.007
    num_grad_channels = 3
    grad_interval = 10.003 # us between [num_grad_channels] channel updates

    gamma = 42570000 # Hz/T

    # value for tabletopMRI  gradient coil
    B_per_m_per_current = 0.02 # T/m/A, approximate value for tabletop gradient coil

    # values for gpa fhdo
    gpa_current_per_volt = 2.5 # gpa fhdo 6A configuration
    max_dac_voltage = 2.5

    max_Hz_per_m = max_dac_voltage * gpa_current_per_volt * B_per_m_per_current * gamma	


    grad_max = max_Hz_per_m # factor used to normalize gradient amplitude, should be max value of the gpa used!	

    rf_amp_max = 10 # factor used to normalize RF amplitude, should be max value of system used!
    tx_warmup = 500 # us
    adc_pad = 0 # padding to prevent junk in rx buffer
    ps = PSAssembler(rf_center=lo_freq*1e6,
        # how many Hz the max amplitude of the RF will produce; i.e. smaller causes bigger RF V to compensate
        rf_amp_max=rf_amp_max,
        grad_max=grad_max,
        clk_t=clk_t,
        tx_t=tx_t,
        grad_t=grad_interval,
        tx_warmup=tx_warmup,
        adc_pad=adc_pad)
    tx_arr, grad_arr, cb, params = ps.assemble('tabletop_se.seq')

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
    x = np.linspace(0,2*np.pi, 100)
    ramp_sine = np.sin(2*x)
    exp.add_tx(ps.tx_arr)
    exp.add_grad(ps.grad_arr)

    # plt.plot(ps.gr_arr[0]);plt.show()


    exp.calibrate_gpa_fhdo(max_current = 2,
        num_calibration_points=10,
        gpa_current_per_volt=gpa_current_per_volt) 

    # set all channels back to 0 A
    for ch in range(num_grad_channels):
        dac_code = exp.ampere_to_dac_code(0)
        dac_code = exp.calculate_corrected_dac_code(ch,dac_code)
        exp.write_gpa_dac(ch, dac_code)      


    data = exp.run() # Comment out this line to avoid running on the hardware
    plt.plot((np.sqrt(data.imag**2 + data.real**2)))
    plt.show()

    # st()    

  

 
