# Spectrum

This folder contains the programs to calculate the power spectrum plots in Chapter 5.

The program in this folder, 'spectrum.f90', calculates the first order correlation of the system and outputs the data to './data_files/spectrum.txt'. The plotting is made by running 'plot_spectrum.py', where the FFT is also calculated.

The contour plots of Chapter 5 are made in from two different programs, both in './Contour_Plots'. One program is for the spectrum as a function of increasing drive strength, \Omega, and for drive detuning, \delta. The respective data for each program is written to the 'Contour_Plots/data_files/map_(delta/omega)' folder.
