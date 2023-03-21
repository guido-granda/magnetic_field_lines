# Magnetic field lines
The goal of this code is to obtain magnetic field lines from three- and two-dimensional grid data coming from simulations. This code main fetaure is that allows to compute physical quantities along the magnetic field lines, which can be usefull to understand physical processes related to magnetohydrodynamics.

Currently, the physical quantitties that can be interpolated along the magnetic field lines are the following:
- Density,
- Magnetic field magnitude
- Alfénic match number
- Match number
- Velocity


Requiretments: 
It requires a uniform grid which can be obtained with codes like yt-project.

Description: 
It store the outputs in hdf5 format and these files contain information of the position of the magnetic field lines and the derived quantitites along then. 
