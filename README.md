# MafSSO
Preliminary development of SSO capabilities for sims_maf.

Requirements:
-sims_maf (preferably v2.0+)
-oorb

Install sims_maf however you would prefer (I'll assume you are familiar with the LSST software stack).

Install oorb. Two choices: git clone the original/master repo (https://github.com/oorb/oorb) and follow all instructions there.
Or (preferred) git clone the eups-ified version created by Mario (https://github.com/EUPSForge/oorb) and follow the original
installation instructions up to the point of setting environment variables. At that point, stop and use eups (this fork has a ups directory).

Note when installing oorb:
* You need a fortran compiler. gfortran will be fine. gfortran 4.9 works, 5.0 did not last I checked (maybe it does now?).
* from the main directory: ./configure gfortran -opt ;  cd main ;  make
* You need to get the JPL binaries: cd data/JPL_ephemeris; make ; cd .. ; ./updateOBSCODE
* You need to build the python interfaces:  cd python ; make pyoorb
Test the installation by doing 'python test.py'.
Don't forget to eups declare / setup oorb.

Git clone this repository.

* setup sims_maf
* setup oorb

You're ready to go - try the ipython notebooks ExampleMoObs.ipynb and ExampleMoMetrics.ipynb.
