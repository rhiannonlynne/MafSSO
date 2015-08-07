# MafSSO
Preliminary development of SSO capabilities for sims_maf.

Welcome beta testers! As a head's up: this software is currently beta-version.

Very shortly (2-3 weeks) it will be merged from here into the sims_maf package and some of the resulting awkwardness of this
package will go away. Other improvements anticipated at that time include improving the documentation, improving how the metrics run
(you shouldn't have to calculate the apparent magnitude in your metric, we'll do that just before), and there may be some name changes for variables.
In general, things should look pretty much the same, but hopefully just easier to use.

Requirements:
-sims_maf (preferably v2.0+)
-oorb

Install sims_maf however you would prefer (the easiest thing to do is to use anaconda + conda to install the binaries - should be about a 30 minute job).
Here are more [instructions](https://github.com/LSST-nonproject/sims_maf_contrib/blob/master/tutorials/Index.ipynb) on installing sims_maf.

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
