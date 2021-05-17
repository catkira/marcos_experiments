import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import scipy.signal as sig
from skimage.transform import iradon
from scipy.interpolate import griddata

if __name__ == "__main__":

    files = listdir("./")
    files2 = [f for f in files if f.find("data2d se Nx") != -1]

    # 2 cylinder decent
    current_file = files2[-1]
    print(current_file)
    data2d = np.load(current_file)

    plt.figure(1)
    plt.subplot(1, 3, 1)
    plt.imshow(10*np.log(np.abs(data2d)),aspect='auto',interpolation='none', origin='lower')
    plt.subplot(1, 3, 2)
    plt.imshow(np.angle(data2d),aspect='auto',interpolation='none')
    plt.subplot(1, 3, 3)
    Nspokes = data2d.shape[0]
    if False:
      alpha = 0.5
      window_x = sig.tukey(data2d.shape[0],alpha)[...,None]*np.ones([1,data2d.shape[1]])
      window_y = np.transpose(sig.tukey(data2d.shape[1],alpha)[...,None]*np.ones([1,data2d.shape[0]]))
      data2d = np.multiply(data2d,window_x)
      data2d = np.multiply(data2d,window_y)

    oversampling_factor = int(np.round(data2d.shape[1]/194))
    data2d = sig.decimate(data2d, oversampling_factor, ftype='iir', axis=1)
    Nx = data2d.shape[1]
    img = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(data2d))))

    plt.imshow(img, aspect='auto',cmap='gray',interpolation='none')
    plt.show()