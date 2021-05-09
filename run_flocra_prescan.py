#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pdb
import time
import argparse
from mri_config import lo_freq, grad_max_Hz_per_m, hf_max_Hz_per_m, gamma

import external
import experiment as ex
import os
import scipy.io as sio
from shutil import copyfile
from flocra_pulseq_interpreter import PSInterpreter
import scipy
st = pdb.set_trace

if __name__ == "__main__":
    #parser = argparse.ArgumentParser(description='Run a pulseq sequence on FLOCRA')
    #parser.add_argument('seq', metavar='seq', type=str, nargs='+',
    #                help='pulseq generated .seq file')
    #args = parser.parse_args()                    
    #print(F"seq file: {args.seq}")

    tx_t = 1 # us
    grad_t = 10 #

    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    seq_file = "tabletop_prescan_se.seq"
    psi = PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=hf_max_Hz_per_m,
                        tx_t=tx_t,
                        grad_t=grad_t,
                        grad_max=grad_max_Hz_per_m)
    od, pd = psi.interpret(seq_file)   
    Nx = pd['Nx']    

    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_t/3) 
    expt.add_flodict(od)

    expt.gradb.calibrate(channels=[0,1,2], max_current=max_grad_current, num_calibration_points=30, averages=5, poly_degree=5)

    rxd, msgs = expt.run()
    #nSamples = pd['readout_number']
    data = rxd['rx0'][6:]
    nSamples = len(data)    
    dt = pd['rx_t'] 

    fig, (ax1, ax2, ax3) = plt.subplots(3)    

    fig.suptitle('Spin Echo [n={:d}, lo_freq={:f} Mhz]\n'.format(nSamples,lo_freq))
    t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
    ax1.plot(t_axis, np.abs(data)*33)
    ax1.set_ylabel('voltage [mV]')
    ax2.set_xlabel('time [us]')
    ax2.plot(t_axis, data.real*33)
    ax2.set_ylabel('voltage [mV]')

    f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))
    fft_data = np.fft.fftshift(np.fft.fft(data))
    ax3.plot(f_axis,np.abs(fft_data/np.sqrt(nSamples)))
    peaks, properties = scipy.signal.find_peaks(np.abs(fft_data),prominence=1)
    ax3.plot(peaks, x[peaks], "x")
    ax3.plt.vlines(x=peaks, ymin=x[peaks] - properties["prominences"],
           ymax = x[peaks], color = "C1")
    plt.show()
    fig.tight_layout()

    # write new lo_freq to mri_config.py
    mri_config_file = 'mri_config.py'
    with open(mri_config_file) as f:
        config=f.read()
    new_content = ""
    for line in iter(config.splitlines()):
        if line.find("lo_freq") != -1:
            new_content += "lo_freq = " + str(lo_freq) + "\n"
        else:
            new_content += line + "\n"

    with open(mri_config_file, "w") as f:
        f.write(new_content)

    plt.close()
    plt.ioff()
