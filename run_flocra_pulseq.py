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
import mplcursors
st = pdb.set_trace

class SnappingCursor:
    """
    A cross hair cursor that snaps to the data point of a line, which is
    closest to the *x* position of the cursor.

    For simplicity, this assumes that *x* values of the data are sorted.
    """
    def __init__(self, ax, line):
        self.ax = ax
        self.horizontal_line = ax.axhline(color='k', lw=0.8, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.8, ls='--')
        self.x, self.y = line.get_data()
        self._last_index = None
        # text location in axes coords
        self.text = ax.text(0.72, 0.9, '', transform=ax.transAxes)

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if not event.inaxes:
            self._last_index = None
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            index = min(np.searchsorted(self.x, x), len(self.x) - 1)
            if index == self._last_index:
                return  # still on the same data point. Nothing to do.
            self._last_index = index
            x = self.x[index]
            y = self.y[index]
            # update the line positions
            self.horizontal_line.set_ydata(y)
            self.vertical_line.set_xdata(x)
            self.text.set_text('x=%1.2f, y=%1.2f' % (x, y))
            self.ax.figure.canvas.draw()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run a pulseq sequence on FLOCRA')
    parser.add_argument('--seq', metavar='seq', type=str, nargs='+',
                    help='pulseq generated .seq file')
    parser.add_argument('--plot-only',nargs='?', default = 0, const = 1,
                    help='only plot the sequence')                    
    parser.add_argument('--scope-csv',nargs=1, default = 0,
                    help='only plot the sequence')                    
    args = parser.parse_args()                    
    print(F"seq file: {args.seq}")

    print('gradient max_B_per_m = {:f} mT/m'.format(grad_max_Hz_per_m/gamma*1e3))	
    print('gradient max_Hz_per_m = {:f} MHz/m'.format(grad_max_Hz_per_m/1E6))
    print('HF max_Hz_per_m = {:f} kHz'.format(hf_max_Hz_per_m/1E3))

    seq_file = args.seq[0]
    psi = PSInterpreter(rf_center=lo_freq*1e6,
                        rf_amp_max=hf_max_Hz_per_m,
                        grad_max=grad_max_Hz_per_m,
                        gx_max=grad_max_x_Hz_per_m,
                        gy_max=grad_max_y_Hz_per_m,
                        gz_max=grad_max_z_Hz_per_m,
                        tx_warmup=200)
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
        Ndummy = int(pd['Ndummy'])
        if seq_file.find("tse") != -1:
            filename = f"{seq_file[:-4]} {current_time} Ny {Ny} Nx {Nx} TR {TR} ETL {pd['ETL']} Ndummy {Ndummy} SliceThickness {sliceThickness}"
        else:
            filename = f"{seq_file[:-4]} {current_time} Ny {Ny} Nx {Nx} TR {TR} SliceThickness {sliceThickness}"
    filename = os.path.basename(filename)
    copyfile(seq_file,os.path.join(data_path,filename+".seq"))        
    copyfile(seq_file[:-4]+".m",os.path.join(data_path,filename+".m"))

    if True and args.plot_only == 0:
        # Shim
        grads = ['grad_vx', 'grad_vy', 'grad_vz']
        for ch in range(3):
            od[grads[ch]] = (np.concatenate((np.array([10.0]), od[grads[ch]][0])), np.concatenate((np.array([0]), od[grads[ch]][1])))
            od[grads[ch]] = (od[grads[ch]][0], od[grads[ch]][1] + shim[ch])


    expt = ex.Experiment(lo_freq=lo_freq,
                         rx_t=pd['rx_t'],
                         init_gpa=True,
                         gpa_fhdo_offset_time=pd['grad_t']/3,
                         flush_old_rx=True,
                         halt_and_reset=True,
                         grad_max_update_rate = 0.2) # in MSPS
    expt.add_flodict(od)

    if args.plot_only == 1:
        if args.scope_csv != 0:
            ##### settings
            grad_delay = -50*8E-3 # us
            global_delay = 0.015 #us
            channels = [[1,'scope_ix'],[0,'scope_txgate'],[2,'scope_rx'],[3,'scope_vx']]
            #####
            scope_info = np.genfromtxt(args.scope_csv[0], delimiter=':',usecols=(0,1), dtype=str)
            scope_info = dict(scope_info)
            scope_wf = np.genfromtxt(args.scope_csv[0].replace(".csv",".Wfm.csv"), delimiter=',')
            fig, axes = plt.subplots(4, 1, figsize=(12,8), sharex='col')
            (txs, grads, rxs, ios) = axes
            expt.plot_sequence(axes = axes)
            t_scope = np.linspace(float(scope_info['XStart'])*1E6,float(scope_info['XStop'])*1E6,num=scope_wf.shape[0])
            t_scope = t_scope - global_delay
            lines = np.array([])
            for channel in channels:
                if "scope_v" in channel[1] or "scope_i" in channel[1]:
                    align_index = abs(int(float(scope_info['XStart']) / float(scope_info['Resolution']))) + 1000  # assumes grad amplitude is 0 at start
                    grad = (scope_wf[:,channel[0]] - 2.5 - (scope_wf[align_index,channel[0]] - 2.5))/2.5
                    if "scope_i" in channel[1]:
                        grad *= 2
                    lines = np.append(lines, grads.step(t_scope - grad_delay, grad, where='post', label=channel[1]))
                if "scope_txgate" in channel[1]:
                    tx_gate = scope_wf[:,channel[0]]/5 
                    lines = np.append(lines, ios.step(t_scope, tx_gate, where='post', label=channel[1]))
                if "scope_rx" in channel[1]:
                    rx = scope_wf[:,channel[0]]
                    lines = np.append(lines, rxs.step(t_scope, rx, where='post', label=channel[1]))
            for axis in axes:
                axis.legend()
                axis.grid(True)
            mplcursors.cursor(lines,multiple=True) # or just mplcursors.cursor()            
        else:
            expt.plot_sequence()
        plt.show()
    else:
        expt.gradb.calibrate(channels=[0,1,2], max_current=max_grad_current, num_calibration_points=30, averages=5, poly_degree=5)

        rxd, msgs = expt.run()
        expt.gradb.init_hw()  # set gradient currents back to zero

        nSamples = pd['readout_number']
        
        np.save(os.path.join(data_path,filename) + ".npy",rxd['rx0'])
        sio.savemat(os.path.join(data_path,filename) + ".mat",rxd)
    plt.close()
    plt.ioff()
