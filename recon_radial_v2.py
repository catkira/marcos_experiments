import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import scipy.signal as sig
from skimage.transform import iradon
from scipy.interpolate import griddata

if __name__ == "__main__":

    files = listdir("./")
    files2 = [f for f in files if f.find("data2d radial v2 ") != -1]

    # 2: decent star
    # 6: decent star with oversampling
    current_file = files2[-1]
    print(current_file)
    data2d = np.load(current_file)

    Nx = 200 # default value
    Nr = 80 # default value
    tokens = current_file.split(" ")
    for n in np.arange(1,len(tokens)):
        if tokens[n] == "Nx":
            Nx = int(float(tokens[n+1]))
        if tokens[n] == "Nspokes":
            Nr = int(float(tokens[n+1]))
    oversampling_factor = int(np.round(data2d.shape[0]/Nx/Nr))
    data2d = data2d.reshape(Nr,Nx*oversampling_factor)    
    adc_pad = 6
    data2d=data2d[:,adc_pad:]  

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

    if False:
        data2d = sig.decimate(data2d, oversampling_factor, ftype='iir', axis=1)
        Nx_truncated = data2d.shape[1]
        points2d = np.zeros( (Nspokes * Nx_truncated, 2) )
        r_len = 1
        for k, alpha in enumerate(np.linspace(0,np.pi,Nspokes,endpoint=False)):        
            rpoints = np.linspace(-r_len, r_len, Nx_truncated)
            points2d[Nx_truncated * k : Nx_truncated * (k + 1), 0] = rpoints * np.cos(alpha)
            points2d[Nx_truncated * k : Nx_truncated * (k + 1), 1] = rpoints * np.sin(alpha)
        grid_oversampling_factor = 1
        grid_points = int(Nx_truncated)*grid_oversampling_factor
        grid_x, grid_y = np.mgrid[-r_len:r_len:grid_points*(1j), -r_len:r_len:grid_points*(1j)]
        data2dgridded = griddata(points2d, data2d.flatten(), (grid_x, grid_y),fill_value=0,method='linear')
        data2dgridded = data2dgridded.reshape((grid_points,grid_points))
        img = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(data2dgridded))))
    else:
        data2d = sig.decimate(data2d, oversampling_factor, ftype='iir', axis=1)
        sinogram = np.abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(data2d), axis=1))) 
        theta = np.linspace(0., 180., Nspokes, endpoint=False)    
        img = iradon(np.swapaxes(sinogram,0,1), theta=theta, filter_name='hann')
    plt.imshow(img, aspect='auto',cmap='gray',interpolation='none')
    plt.show()