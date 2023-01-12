import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import scipy.signal as sig
from scipy.interpolate import griddata
from mri_config import data_path
import os.path

if __name__ == "__main__":

    #data_path = 'D:\\Sciebo\\LEMB\\TabletopMRI\\tse analysis\\'
    #data_path = 'D:\\Sciebo\\LEMB\\TabletopMRI\\measurements\\'
    data_path = '/mnt/d/sciebo/LEMB/TabletopMRI/measurements/'
    files = listdir(os.path.join(data_path,""))
    #files2 = [f for f in files if f.find("se_v2_2d") != -1 and f.find(".npy") != -1]
    files2 = [f for f in files if f.find("tse_") != -1 and f.find(".npy") != -1]
    #files2 = [f for f in files if f.find("gre_v2_2d") != -1 and f.find(".npy") != -1]

    # 2 cylinder decent
    current_file = files2[-1]
    print(current_file)
    data2d = np.load(os.path.join(data_path,current_file))
    Nx = 200 
    Ny = 80 
    tokens = current_file.split(" ")
    ETL = 1
    Ndummy = 0
    for n in np.arange(1,len(tokens)):
        if tokens[n] == "Nx":
            Nx = int(float(tokens[n+1]))
        if tokens[n] == "Ny":
            Ny = int(float(tokens[n+1]))    
        if tokens[n] == "ETL":
            ETL = int(float(tokens[n+1]))    
        if tokens[n] == "Ndummy":
            Ndummy = int(float(tokens[n+1]))    
    oversampling_factor = int(np.round(data2d.shape[0]/Nx/(Ny+Ndummy*ETL)))
    data2d = data2d.reshape(Ny+Ndummy*ETL,Nx*oversampling_factor)    

    adc_pad = 6
    data2d=data2d[:,adc_pad:]   

    if Ndummy > 0:
        data2d = np.delete(data2d,np.arange(Ndummy*ETL),axis=0)

    if True and "tse" in current_file:
        nex = int(Ny / ETL)
        pe_steps=np.arange(Ny)
        pe_steps_interleaved = np.reshape(pe_steps, (ETL, nex)).astype('int'); 
        data2d_sorted = np.empty_like(data2d)
        for n in range(nex):
            for m in range(ETL):
                data2d_sorted[pe_steps_interleaved[m,n]] = data2d[n*ETL + m,:]
        data2d = data2d_sorted

    plt.figure(1)
    plt.subplot(1, 3, 1)
    plt.imshow(10*np.log(np.abs(data2d)),aspect='auto',interpolation='none', origin='lower')
    plt.subplot(1, 3, 2)
    plt.imshow(np.angle(data2d),aspect='auto',interpolation='none',origin='lower')
    plt.subplot(1, 3, 3)
    Nspokes = data2d.shape[0]
    if False:
      alpha = 0.5
      window_x = sig.tukey(data2d.shape[0],alpha)[...,None]*np.ones([1,data2d.shape[1]])
      window_y = np.transpose(sig.tukey(data2d.shape[1],alpha)[...,None]*np.ones([1,data2d.shape[0]]))
      data2d = np.multiply(data2d,window_x)
      data2d = np.multiply(data2d,window_y)

    data2d = sig.decimate(data2d, oversampling_factor, ftype='iir', axis=1)
    Nx = data2d.shape[1]
    img = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(data2d))))

    plt.imshow(img, aspect='auto',cmap='gray',interpolation='none')
    plt.show()