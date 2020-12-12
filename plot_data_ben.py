import numpy as np
import matplotlib.pyplot as plt
import os

file0d = 'data ben Nx 128 20-12-12 16_07_20.npz'
file1d = 'data1d ben Nx 256 20-12-12 16_05_28.npz'
file2d = 'data2d ben Nx 166 Ny 66 TR 5 20-12-12 15_47_17.npy'

def plot0d():
   data0d_dict = np.load(file0d)
   dt = data0d_dict['dt']
   nSamples = int(data0d_dict['nSamples'])
   lo_freq = data0d_dict['lo_freq']
   data = data0d_dict['data1d']

   Noise = np.abs(np.std(np.real(np.fft.fft(data))[int(data.size/2)-3:int(data.size/2)+3]))
   SNR=np.max(np.abs(np.fft.fft(data)))/Noise
   fig, (ax1, ax2, ax3) = plt.subplots(3)
   fig.suptitle('Spin Echo [n={:d}, lo_freq={:f} Mhz]\nSNR={:f}'.format(nSamples,lo_freq,SNR))
   t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
   ax1.plot(t_axis, np.abs(data)*3.3)
   ax1.set_ylabel('voltage [V]')
   ax2.set_xlabel('time [us]')
   ax2.plot(t_axis, data.real*3.3)
   ax2.set_ylabel('voltage [V]')
   f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))
   ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data))/np.sqrt(nSamples)))
   fig.tight_layout()


def plot1d():
   data1d_dict = np.load(file1d)

   dt = data1d_dict['dt']
   nSamples = int(data1d_dict['nSamples'])
   lo_freq = data1d_dict['lo_freq']
   data1d = data1d_dict['data1d']

   fig, (ax1, ax2, ax3) = plt.subplots(3)
   fig.suptitle('Spin Echo [n={:d}, lo_freq={:f} Mhz]\n'.format(nSamples,lo_freq))
   t_axis = np.linspace(0, dt * nSamples, nSamples)  # us    
   ax1.plot(t_axis, np.abs(data1d)*3.3)
   ax1.set_ylabel('voltage [V]')
   ax2.set_xlabel('time [us]')
   ax2.plot(t_axis, data1d.real*3.3)
   ax2.set_ylabel('voltage [V]')
   f_axis = np.fft.fftshift(np.fft.fftfreq(nSamples,dt*1E-6))
   ax3.plot(f_axis,np.abs(np.fft.fftshift(np.fft.fft(data1d))/np.sqrt(nSamples)))
   fig.tight_layout()

def plot2d():
   data2d = np.load(file2d)
   plt.figure(3)
   plt.subplot(1, 3, 1)
   plt.imshow(10*np.log(np.abs(data2d)),aspect='auto',interpolation='none')
   plt.subplot(1, 3, 2)
   plt.imshow(np.angle(data2d),aspect='auto',interpolation='none')
   plt.subplot(1, 3, 3)
   img = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(data2d))))
   plt.imshow(img, aspect='auto',cmap='gray',interpolation='none')



plot0d()
plot1d()
plot2d()
plt.show()



