# -*- coding: utf-8 -*-
"""
Created on Wed May  2 16:08:51 2018

@author: Jacob
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.close('all')

###############################################################################
#                                 Parameters                                  #
###############################################################################
# filename
direc = "./data_files/"
name = "corr_full"
#name = "omega40_xi05_delshift"
ext = ".txt"
filename = direc + name + ext

# Reading the values of each parameter from the txt file
omega = float(open(filename).readlines()[0])
#omega = round(omega, 2)
delta = float(open(filename).readlines()[1])
xi = float(open(filename).readlines()[2])
alpha = float(open(filename).readlines()[3])

# List of times from the txt file
tau = np.loadtxt(fname=filename, dtype='float', usecols=(0,), skiprows=5)
# List of correlations from txt file.
corr_full = np.loadtxt(fname=filename, dtype='float', usecols=(1,), skiprows=5)

###############################################################################
#                                 Plot                                        #
###############################################################################
## Plotting the correlation function against tau
#plt.figure(figsize=[6, 2.5])
#plt.plot(tau, corr_full, color='k')
#
#plt.xlim(-0.5, 20.5)
#
#plt.xlabel(r'$ \gamma \tau $', fontsize=11)
#plt.ylabel(r'$ g^{(2)}(\tau) $', fontsize=11)
#
##plt.title(r'$ \Omega / \gamma =%s, \delta / \gamma = %s$'%(omega, delta) + \
##          r'$ \alpha / \gamma = %s, \xi = %s$'%(alpha, xi))
#
##plt.grid()
#plt.tight_layout()
#plt.show()
#
#plt.savefig("./Images/" + name + ".pdf")

###############################################################################
#                          Plot for thesis                                    #
###############################################################################
## filename
#direc = "./data_files/"
#name = "2level_omega0001"
#ext = ".txt"
#filename = direc + name + ext
#
## List of times from the txt file
#tau1 = np.loadtxt(fname=filename, dtype='float', usecols=(0,), skiprows=5)
## List of correlations from txt file.
#corr1 = np.loadtxt(fname=filename, dtype='float', usecols=(1,), skiprows=5)
#
#name = "2level_omega275"
#filename = direc + name + ext
## List of times from the txt file
#tau2 = np.loadtxt(fname=filename, dtype='float', usecols=(0,), skiprows=5)
## List of correlations from txt file.
#corr2 = np.loadtxt(fname=filename, dtype='float', usecols=(1,), skiprows=5)
#
#fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True,
#                       figsize=[6, 3])
#
## omega = 5
#ax[0].plot(tau1, corr1, color='black')
#ax[0].axhline(y=1.0, color='k', alpha=0.25, ls='dashed')
#
## omega = 40
#ax[1].plot(tau2, corr2, color='black')
#ax[1].axhline(y=1.0, color='k', alpha=0.25, ls='dashed')
#
#ax[0].set_xlim(min(tau1), max(tau1))
#ax[1].set_xlim(min(tau2), max(tau2))
#
#ax[0].set_ylabel(r'$ g^{(2)}(\tau) $')
#ax[0].set_xlabel(r'$ \gamma \tau $')
#ax[1].set_xlabel(r'$ \gamma \tau $')
#
#
#fig.tight_layout()
#fig.savefig("./Images/2level_antibunching.pdf")

###############################################################################
#                          Multiple Omega                                     #
###############################################################################
#omega_list = ["0.1", "1.0", "2.5", "5.0"]
omega_list = ["0.1", "0.25", "0.5", "0.75"]
direc = "./data_files/low_drive/"
ext = ".txt"

tau = np.loadtxt(fname=direc + "0.1" + ext, dtype='float',
                 usecols=(0,), skiprows=5)

plt.figure(figsize=[6, 3.5])
for i in omega_list:
    filename = direc + i + ext
    corr = np.loadtxt(fname=filename, dtype='float', usecols=(1,), skiprows=5)
    plt.plot(tau, corr, label=r'$\Omega/\gamma = %s$'%(str(i)))

plt.xlabel(r'$ \tau / \gamma $', fontsize=11)
plt.ylabel(r'$ g^{(2)}_{\mathrm{Unfiltered}}(\tau) $', fontsize=11)

plt.xlim(0, 10)

#plt.grid()
plt.legend(fontsize=11)
plt.tight_layout()
plt.show()

plt.savefig("./Images/omegalow.pdf")
