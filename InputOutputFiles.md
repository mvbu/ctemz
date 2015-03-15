# Input file `temzinp.txt` #

This file is read in right at the start of execution and contains the runtime values for various simulation parameters. Here is the file as currently checked in (Feb 2014), with some comments added:

```
# Input Parameters no. cells radius (<= 5), zred, dgpc, alpha, p, bavg, -PSD slope, E(e)/E(B), rsize, gmaxmn, gmrat, gmin, betau
p(laminar), betaup(turb), zeta (shockangle), thlos, opang, dust temp., obs. dust luminosity in units of 1e45 erg/s, dust torus d
istance, torus x-section rad., zdist0, vdm (vol. of Mach disk/cell vol.)
5  # number of cells in radial direction (defines density of cells viewed axially)
0.069 # zred: redshift of the source/blazar
0.295 # dgpc: distance to the source in Gpc
0.55 # alpha = (s-1)/2 is the underlying spectral index, where s = slope of electron E dist
0.80 
0.06
1.7 # negative of the slope of the power spectrum distribution
0.25e0
0.004
3000.0 # these next three parameters define the range of electron energies (gamma)
10.0
300.0
0.99d0
0.577d0
10.0d0
7.7d0
1.9d0
1200.0 # (blackbody) temperature of the dusty torus, used to calculate radiation received from dusty torus
0.1 
1.5
0.3
1.0
0.4e0 # defines volume of Mach disk relative to cell volume, which is calculated from the above params
```

(TODO: need to describe **all** the parameters above)

# Output Files #

The output files are as follows:
  * `ctemzlc.txt` Light curve at each time step
  * `ctemzspec.txt` Spectrum at each time step, for 68 frequencies
  * `ctemzpol.txt` Polarization values for each time step
  * `maps/ctemzmap[xxxx].txt` "Map" files used to create diagrams (TODO: not sure exactly what they do. Prof. Marscher would know)