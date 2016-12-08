# Fuchs-Sondheimer

This is a simple python script that calculates the (theoretical) effective thermal conductivity of several thin films, based on
the Fuchs-Sondheimer theory.
Film thicknesses are entered within the script, in the list named "LList".
The materials properties must be entered as a txt file with 6 columns corresponding to, respectively, the frequencies,
the density of states, the velocities, the size of the frequency cells, the relaxation times and the polarization (the latter is not really used). See dataSi.txt for an example. Note that in dataSi.txt, the density of states for the transverse acoustic phonons is doubled in order to account for the two TA branches.
The thermal conductivities are output both in the terminal and in a .txt file.
