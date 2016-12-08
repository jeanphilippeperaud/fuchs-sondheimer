# Copyright (c) Jean-Philippe Peraud 2016
# All rights reserved

from scipy import special
import numpy as np

# load the phonon data for silicon
# column 0 = frequencies
# column 1 = density of states
# column 2 = velocity
# column 3 = Delta frequencies
# column 4 = relaxation times
# column 5 = polarization
phonondat=np.loadtxt('dataSi.txt')
N = len(phonondat)   # number of frequency cells
F   = phonondat[:,0]
DOS = phonondat[:,1]
V   = phonondat[:,2]
DF  = phonondat[:,3]
tau = phonondat[:,4]
pol = phonondat[:,5]

# list of film thicknesses (in meters)
LList = [10e-9, 20e-9, 40e-9, 60e-9, 80e-9, 1.2e-7, 1.6e-7, 2e-7, 2.5e-7, 3e-7, 3.5e-7, 4e-7, 4.5e-7, 5e-7]


# temperature (K)
T = 300
# reduced Planck constant
hbar = 1.054517e-34
#Boltzmann constant
boltz = 1.38065e-23

# construct the temperature derivative of (energy-based) Bose Einstein distribution and store it
# construct the product DOS*V^2*tau*d(BE)/dT*Delta_freq and store it
dBE = np.zeros(N)
dC  = np.zeros(N)
for i in range(N):
    dBE[i] = (hbar*F[i]/T)**2/boltz*np.exp(hbar*F[i]/boltz/T)/(np.exp(hbar*F[i]/boltz/T)-1)**2
    dC[i] = DOS[i]*dBE[i]*V[i]**2*tau[i]*DF[i]

# calculate the bulk thermal conductivity
kth = np.sum(dC)*1/3

kfilm = []
for L in LList:
    # calculate the thin film correction
    corr = 0
    for i in range(N):
        a= L/V[i]/tau[i]
        tmp = special.expn(3,a)-special.expn(5,a)
        corr = corr + V[i]**3*tau[i]**2*dBE[i]*DOS[i]/2*DF[i]*tmp/L
        corr = corr - V[i]**3*tau[i]**2*dBE[i]*DOS[i]/8*DF[i]/L
    kfilm.append(kth+corr)

print "\n The bulk thermal conductivity is: " + str(kth) + " W/m/K \n"
print "\n The thin film effective thermal conductivities for the thicknesses"
print LList
print "are (in W/m/K):"
print kfilm
data = np.concatenate((np.array(LList,ndmin=2),np.array(kfilm,ndmin=2)),0)
np.savetxt("kfilm.txt", np.transpose(data))
