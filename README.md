# Combine_spectra
This code combines multiple spectra together and also considers the error propagation in case there is a column on flux errors. The code reads either the fits files or text files. Currently, the code has been tested on ESO and polarbase spectrum. The structure of the code requires the first column to be wavelength, the second to be flux and the third to be flux errors. The parameters of the function are -

**stardirpath**: Directory with all the spectrum

**unnorm_combine_specpath**: Filename with complete path to save the combine spectrum

**median_combine**: If True, combine spectra by taking median. By default combine spectrum by taking mean.

**wave_in_angstrom**: If True, convert any wavelength to Angstrom. By default it is true.

**plot_comb_spec**: Plot the combine spectrum. By default, False.

**savefile**: Save the combined spectrum. By default, True.

**kwargs**: skiprows = n (int) to skip the first "n" rows of the spectra.

**return**: No return
