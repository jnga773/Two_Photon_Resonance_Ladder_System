# Second Order Correlations

This folder contains the programs to calculate the second order correlation plots in Chapter 6.

There are two main programs here:

  - Filtered
  - Unfiltered

The Unfiltered program does not use the `ParamList.nml` convention that the other programs use but, like the spectrum programs, solve the Lindblad Master Equation and uses the Quantum Regression Formulae to solve for the second order correlation.

In the `Filtered` folder, there is a Master Equation solver and a Quantum Trajectory solver. Both program give the same results but the Quantum Trajectory solve must be averaged many times to give the exact same as the Master Equation solver. Both of these programs use the `ParamList.nml` convention and so only need to be compiled once. 
