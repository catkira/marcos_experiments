import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import scipy.signal as sig
from skimage.transform import iradon
from scipy.interpolate import griddata

if __name__ == "__main__":

    files = listdir("./")
    files2 = [f for f in files if f.find("data2d radial") != -1]

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

    r_len = 1
    oversampling_factor = 4
    if True:
        data2d = sig.decimate(data2d, oversampling_factor, ftype='iir', axis=1)
        Nx = data2d.shape[1]
        points2d = np.zeros( (Nspokes * Nx, 2) )
        for k, alpha in enumerate(np.linspace(0,np.pi,Nspokes,endpoint=False)):        
            rpoints = np.linspace(-r_len, r_len, Nx)
            points2d[Nx * k : Nx * (k + 1), 0] = rpoints * np.cos(alpha)
            points2d[Nx * k : Nx * (k + 1), 1] = rpoints * np.sin(alpha)

        grid_points = int(Nx)
        grid_x, grid_y = np.mgrid[-r_len:r_len:grid_points*(1j), -r_len:r_len:grid_points*(1j)]
        data2dgridded = griddata(points2d, data2d.flatten(), (grid_x, grid_y),fill_value=0,method='linear')
        data2dgridded = data2dgridded.reshape((grid_points,grid_points))
        img = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(data2dgridded))))
    else:
        data2d = sig.decimate(data2d, oversampling_factor, ftype='iir', axis=1)
        theta = np.linspace(0., 180., Nspokes, endpoint=False)        
        img = iradon(np.swapaxes(np.abs(data2d),0,1), theta=theta, filter_name='ramp')
    plt.imshow(img, aspect='auto',cmap='gray',interpolation='none')
    plt.show()