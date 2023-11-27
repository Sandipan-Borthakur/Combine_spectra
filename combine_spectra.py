import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import os

def raw_spec_combine(stardirpath,unnorm_combine_specpath,median_combine=False,wave_in_angstrom=True,plot_comb_spec=False,savefile=True,**kwargs):
    '''

    :param stardirpath: Directory with all the spectrum
    :param unnorm_combine_specpath: Filename with complete path to save the combine spectrum
    :param median_combine: If True, combine spectra by taking median. By default combine spectrum by taking mean.
    :param wave_in_angstrom: If True, convert any wavelength to Angstrom. By default it is true.
    :param plot_comb_spec: Plot the combine spectrum. By default, False.
    :param savefile: Save the combined spectrum. By default, True.
    :param kwargs: skiprows = n (int) to skip the first "n" rows of the spectra.
    :return: No return
    '''
    # List the all spectrum in stardirpath
    specfilelist = os.listdir(stardirpath)
    masterflux = []
    fluxerr_final = []
    wave_precisionlist = []
    flux_precisionlist = []
    # reading each file from the list. Also considers the flux error if it exists in the file.
    for ispecfile,specfile in enumerate(specfilelist):
        fluxerr_exists = False
        specfilepath = os.path.join(stardirpath, specfile)

        # for fits format spectrum
        if specfile.endswith(".fits"):
            data = fits.getdata(specfilepath)

            # checking if flux error exists
            if not np.isnan(data[0][2][0]):
                fluxerr = data[0][2]
                fluxerr_exists = True
            wave, flux = data[0][0], data[0][1]

        # for text format spectrum
        else:
            if "skiprows" in kwargs.keys():
                skiprows = kwargs["skiprows"]
            else:
                skiprows = 0
            data = np.loadtxt(specfilepath,skiprows=skiprows)

            # checking if flux error exists
            if data.shape[1]==3:
                fluxerr = data[:,2]
                fluxerr_exists = True
            else:
                print(specfile)
            wave,flux = data[:,0],data[:,1]

        if ispecfile==0:
            wave_final = wave
            masterflux.append(flux)
        else:
            # interpolating other spectrum to the wavelengths of the first spectrum
            flux_temp = np.interp(wave_final,wave,flux)
            masterflux.append(flux_temp)

        if fluxerr_exists:
            fluxerr_temp = np.interp(wave_final,wave,fluxerr)
            fluxerr_final.append(fluxerr_temp**2)

        wave_precision = np.median(np.array([len(str(i).split(".")[-1]) for i in wave],dtype='i'))
        flux_precision = np.median(np.array([len(str(i).split(".")[-1]) for i in flux],dtype='i'))
        wave_precisionlist.append(wave_precision)
        flux_precisionlist.append(flux_precision)

    wave_precision_final = np.mean(np.array(wave_precisionlist,dtype="i"),dtype="i")
    flux_precision_final = np.mean(np.array(flux_precisionlist,dtype="i"),dtype="i")
    print(wave_precision_final,flux_precision_final,type(wave_precision_final),type(flux_precision_final))
    masterflux = np.array(masterflux)
    if median_combine:
        flux_final = np.median(masterflux,axis=0)
    else:
        flux_final = np.mean(masterflux,axis=0)

    if wave_in_angstrom:
        numsize = len(str(wave_final[0]).split(".")[0])
        wavemultfactor = 4 - numsize
        wave_final = wave_final * 10 ** wavemultfactor
        wave_precision_final = wave_precision_final-1

    # plot the combined unnormalised spectrum, if plot_comb_spec is True
    if plot_comb_spec:
        plt.plot(wave_final,flux_final)
        plt.show()

    # Calculating the flux error for the combined spectrum, if it exists
    if fluxerr_exists and len(fluxerr_final)==len(masterflux):
        fluxerr_final = np.array(fluxerr_final)

        # Combined flux error is the root-mean-square of the individual errors
        fluxerr_final = np.sqrt(np.sum(fluxerr_final,axis=0))/len(fluxerr_final)

        #error propagation if the spectrum is median combined
        if median_combine:
            fluxerr_final = fluxerr_final*np.sqrt(np.pi*(1/2+1/(4*len(fluxerr_final))))
        savelist = wave_final,flux_final,fluxerr_final
        fmt = ["%.{}f".format(wave_precision_final),"%.{}f".format(flux_precision_final),"%.{}f".format(flux_precision_final)]
    else:
        savelist = wave_final,flux_final
        fmt = ["%.{}f".format(wave_precision_final),"%.{}f".format(flux_precision_final)]

    # save the combined unnormalised spectrum, if savefile is True
    if savefile:
        np.savetxt(unnorm_combine_specpath,np.c_[savelist],delimiter=" ",fmt=fmt)

stardirpath = "/home/sandipan/Documents/E/PhD_Tartu1/Ariel/Heleri_paper/Data/KELT-17"
unnorm_combine_specpath = "/home/sandipan/Documents/E/PhD_Tartu1/Ariel/Heleri_paper/Data/KELT-17/trial.txt"
# raw_spec_combine(stardirpath,unnorm_combine_specpath,plot_comb_spec=True,savefile=True)
data = np.loadtxt(unnorm_combine_specpath)
plt.errorbar(data[:,0],data[:,1],yerr=data[:,2])
plt.show()