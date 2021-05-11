#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pdb
import time
import argparse
import os.path
from mri_config import lo_freq, grad_max_Hz_per_m, hf_max_Hz_per_m, gamma, max_grad_current, grad_max_x_Hz_per_m, grad_max_y_Hz_per_m, grad_max_z_Hz_per_m, data_path, shim

import external
import experiment as ex
import os
import scipy.io as sio
from shutil import copyfile
from flocra_pulseq_interpreter import PSInterpreter
st = pdb.set_trace

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run a pulseq sequence on FLOCRA')
    parser.add_argument('seq', metavar='seq', type=str, nargs='+',
                    help='pulseq generated .seq file')
    args = parser.parse_args()                    
    print(F"seq file: {args.seq}")

    tx_t = 1 # us
    grad_t = 10 #

    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    seq_file = args.seq[0]
    psi = PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=hf_max_Hz_per_m,
                        tx_t=tx_t,
                        grad_t=grad_t,
                        grad_max=grad_max_Hz_per_m,
                        gx_max=grad_max_x_Hz_per_m,
                        gy_max=grad_max_y_Hz_per_m,
                        gz_max=grad_max_z_Hz_per_m)
    od, pd = psi.interpret(seq_file)   

    TR = pd['TR']
    Nx = pd['Nx']
    sliceThickness = int(pd['SliceThickness'])

    from datetime import datetime
    now = datetime.now()
    current_time = now.strftime("%y-%m-%d %H_%M_%S")
    if seq_file.find("radial") != -1:
        Nspokes = int(pd['Nspokes'])
        angle = int(pd['angle'])
        filename = f"{seq_file[:-4]} {current_time} Nspokes {Nspokes} Nx {Nx} TR {TR} angle {angle} SliceThickness {sliceThickness}"
    else:
        Ny = int(pd['Ny'])
        filename = f"{seq_file[:-4]} {current_time} Ny {Ny} Nx {Nx} TR {TR} SliceThickness {sliceThickness}"
    copyfile(seq_file,os.path.join(data_path,filename+".seq"))        
    copyfile(seq_file[:-4]+".m",os.path.join(data_path,filename+".m"))

    if True:
        # Shim
        grads = ['grad_vx', 'grad_vy', 'grad_vz']
        for ch in range(3):
            od[grads[ch]] = (np.concatenate((np.array([10.0]), od[grads[ch]][0])), np.concatenate((np.array([0]), od[grads[ch]][1])))
            od[grads[ch]] = (od[grads[ch]][0], od[grads[ch]][1] + shim[ch])


    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=grad_t/3,
                         flush_old_rx=True) 
    expt.add_flodict(od)

    expt.gradb.calibrate(channels=[0,1,2], max_current=max_grad_current, num_calibration_points=30, averages=5, poly_degree=5)

    rxd, msgs = expt.run()
    expt.gradb.init_hw()  # set gradient currents back to zero

    nSamples = pd['readout_number']
    
    np.save(os.path.join(data_path,filename) + ".npy",rxd['rx0'])
    sio.savemat(os.path.join(data_path,filename) + ".mat",rxd)
    plt.close()
    plt.ioff()
