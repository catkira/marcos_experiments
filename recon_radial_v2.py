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
    Nspokes = 80 # default value
    tokens = current_file.split(" ")
    for n in np.arange(1,len(tokens)):
        if tokens[n] == "Nx":
            Nx = int(float(tokens[n+1]))
        if tokens[n] == "Nspokes":
            Nspokes = int(float(tokens[n+1]))
    oversampling_factor = int(np.round(data2d.shape[0]/Nx/Nspokes))
    data2d = data2d.reshape(Nspokes,Nx*oversampling_factor)    
    adc_pad = 6
    data2d=data2d[:,adc_pad:]  
    if False:
        # correct phase drift 
        center_index = int(data2d.shape[1]/2)
        phase = np.angle(data2d[0,center_index])
        for k in np.arange(1,Nspokes):
            delta = np.angle(data2d[k,center_index]) - phase
            data2d[k,:] *= np.exp(-1j*delta)
    if True:
        # plot peak distribution
        plt.figure(3)
        peaks = np.argmax(np.abs(data2d), axis=1)
        plt.plot(peaks)

    data2d = sig.decimate(data2d, oversampling_factor, ftype='fir', axis=1)   

    if False:
        # apodization, usually not needed
        alpha = 0.2
        window_x = sig.tukey(data2d.shape[0],alpha)[...,None]*np.ones([1,data2d.shape[1]])
        window_y = np.transpose(sig.tukey(data2d.shape[1],alpha)[...,None]*np.ones([1,data2d.shape[0]]))
        data2d = np.multiply(data2d,window_x)
        data2d = np.multiply(data2d,window_y)
    plt.figure(1)
    plt.subplot(1, 3, 1)
    plt.imshow(10*np.log(np.abs(data2d)),aspect='auto',interpolation='none', origin='lower')
    plt.subplot(1, 3, 2)
    plt.imshow(np.angle(data2d),aspect='auto',interpolation='none')
    if True:
        Nx_truncated = data2d.shape[1]
        #print(np.angle(data2d[:,int(Nx_truncated/2)]))        
        points2d = np.zeros((Nspokes * Nx_truncated, 2))
        r_len = 1
        rpoints = np.linspace(-r_len, r_len, Nx_truncated)
        for k, alpha in enumerate(np.linspace(0,np.pi,Nspokes,endpoint=False)):        
            points2d[Nx_truncated * k : Nx_truncated * (k + 1), 0] = rpoints * np.cos(alpha)
            points2d[Nx_truncated * k : Nx_truncated * (k + 1), 1] = rpoints * np.sin(alpha)

        grid_oversampling_factor = 2
        grid_points = int(Nx_truncated*grid_oversampling_factor)
        grid_x, grid_y = np.mgrid[-r_len:r_len:grid_points*(1j), -r_len:r_len:grid_points*(1j)]
        data2dgridded = griddata(points2d, data2d.flatten(), (grid_x, grid_y),fill_value=0,method='cubic')
        data2dgridded = data2dgridded.reshape((grid_points,grid_points))
        if True:
            # plot circular k-space
            plt.figure(2)     
            plt.imshow(10*np.log(np.abs(data2dgridded)),aspect='auto',interpolation='none', origin='lower')        
            img = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(data2dgridded))))
        if False:
            # density compensation: abs(r)*1/Nspokes (M-shape)
            offset = 0.5
            for k in np.arange(img.shape[0]):
                img[k] = img[k][:]*(np.abs((np.arange(img.shape[1])-img.shape[1]/2) + 1j*(k-img.shape[0]/2))/Nspokes + offset)
    else:
        sinogram = np.abs(np.fft.fftshift(np.fft.fft(np.fft.fftshift(data2d), axis=1))) 
        theta = np.linspace(0., 180., Nspokes, endpoint=False)    
        img = iradon(np.swapaxes(sinogram,0,1), theta=theta, filter_name='hann')
    
    plt.figure(1)
    plt.subplot(1, 3, 3)        
    plt.imshow(img, aspect='auto',cmap='gray',interpolation='none')
    plt.show()