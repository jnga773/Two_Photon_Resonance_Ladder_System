# Two Photon Resonance Fluorescence In A Three-Level Ladder System

This Git Repository contains the code for the main results of my MSc Physics thesis.

There are two main folders for the code
  - Spectrum
  - Second Order Correlations

In the Resonance Fluorescence folder, the program is a Lindblad Master Equation solver using Quantum Regression Forumalae to solve for the first order correlation. The power spectrum is then calculated using a Fast Fourier Transform in the plotting python script.

The Second Order Correlations folder there are two main programs. The first calculates the second order correlation of the unfiltered fluorescence using a similiar method used for the power spectrum calculations. The second program introduces the Lorentzian filters to filter the fluorescence and calculate the second order coherence. This can be done in two ways: using a Lindblad Master Equations solve, and using Quantum Trajectory Theory.

The main programs are written in Fortran and were compiled using the Intel iFort compiler with

```shell
$ ifort -O3 -o NAME ./PROGRAM.f90
$ ./NAME
```

The programs only need to be compiled once as all the parameters can be changed by editing the local 'ParamsList.nml' file, which the program will pull all the relevant parameters from. All of the programs write the data into a './data_files' folder for neatness.

### Abstract for Thesis

In this thesis we consider a three-level ladder-type atom driven by a coherent laser. When driven on two-photon resonance, the atom is excited into its highest state by absorbing two photons simultaneously, followed by a cascaded decay.

Employing techniques derived from the Lindblad master equation, we solve for the first-order correlation function, from which we can calculate the fluorescence spectrum. We find that under a strong drive field, the fluorescence spectrum contains up to seven different peaks. We aim to explain the emergence of these peaks by looking at transitions amongst the atom's dressed states.

We then aim to characterise the nature of the emitted light by investigating the second-order correlation function. In order to obtain a more precise picture of the photon correlations, we introduce a frequency filtering technique that allows us to isolate individual transitions. Measuring the photon correlations of these transitions provides a more complete picture of the role of specific dressed states of the system. We provide mostly numerical surveys of the fluorescence spectrum and correlations, with some analytic expressions to verify results.
