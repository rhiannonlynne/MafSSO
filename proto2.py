import os
import numpy as np
from numpy.lib.recfunctions import merge_arrays
import matplotlib.pyplot as plt
from itertools import repeat
import scipy
from scipy import interpolate
import pandas
import pyoorb as oo

from lsst.sims.maf.db import OpsimDatabase
import lsst.sims.photUtils.Bandpass as Bandpass
import lsst.sims.photUtils.Sed as Sed

def calcColors(sedname='C.dat'):
    # Calculate SSO colors.
    filterdir = os.getenv('LSST_THROUGHPUTS_BASELINE')
    filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
    lsst ={}
    for f in filterlist:
        lsst[f] = Bandpass()
        lsst[f].readThroughput(os.path.join(filterdir, 'total_'+f+'.dat'))
    vband = Bandpass()
    vband.readThroughput('harris_V.dat')
    csed = Sed()
    csed.readSED_flambda(sedname)
    vmag = csed.calcMag(vband)
    dmags = {}
    for f in filterlist:
        dmags[f] = csed.calcMag(lsst[f]) - vmag
    return dmags

def joinObs(ephs, simdata, sedname='C.dat', tol=1e-8):
    dmagDict = calcColors(sedname)
    idxs = np.zeros(len(ephs), int)
    dmagColor = np.zeros(len(ephs), float)
    for i, obs in enumerate(ephs):
        match = np.where(np.abs(simdata['expMJD'] - obs['time']) < tol) [0]
        if len(match) != 1:
            print 'Found %d simdata matches - expected only one' %(len(match))
        idxs[i] = match
        # Find adjustment due to color of sso.
        dmagColor[i] = dmagDict[simdata[match]['filter'][0]]
    # Find losses due to trailing, using adjustment to exp time.
    vel = np.sqrt(ephs['dradt']**2 + ephs['ddecdt']**2)  #deg/day
    vel = vel / 24.0  # "/s
    # See https://listserv.lsstcorp.org/mailman/private/lsst-solarsystem/2006-December/24.html
    # should grab simObs['expTime'] .. 
    t = 30.0 # seconds
    teff = t/(1+1.6*vel*t/simdata[idxs]['finSeeing'])
    dmagTrailing = 1.25*np.log10(t/teff)
    # Convert to recarray.
    dmags = np.rec.fromarrays([dmagColor, dmagTrailing], names=['dmagColor', 'dmagTrailing'])
    # Join arrays.
    ssoObs = merge_arrays([ephs, simdata[idxs], dmags], flatten=True, asrecarray=True)
    return ssoObs


def nObsMetric(ssoObs, magH=20., Hrange=np.arange(12, 27, 0.25), snrLimit=5):
    # Given the observations for a particular object and the opsim metadata (join using joinObs)
    # Return the number of observations above SNR (in any band) as a function of H
    countObs = np.zeros(len(Hrange), int)
    # Calculate the magnitude of this object in this filter, each H magnitude. 
    for i, H in enumerate(Hrange):
        magObs = H - magH + ssoObs['magV'] + ssoObs['dmagColor']
        magLimitWithTrailing = ssoObs['fiveSigmaDepth'] - ssoObs['dmagTrailing']
        snr = 5.0 * np.power(10., 0.4*(magLimitWithTrailing - magObs))
        countObs[i] = np.where(snr >= snrLimit)[0].size
    return countObs


def nObsMetric(ssoObs, magH=20., Hrange=np.arange(12, 27, 0.25), snrLimit=5):
    # Given the observations for a particular object and the opsim metadata (join using joinObs)
    # Return the number of observations above SNR (in any band) as a function of H
    countObs = np.zeros(len(Hrange), int)
    # Calculate the magnitude of this object in this filter, each H magnitude. 
    for i, H in enumerate(Hrange):
        magObs = H - magH + ssoObs['magV'] + ssoObs['dmagColor']
        magLimitWithTrailing = ssoObs['fiveSigmaDepth'] - ssoObs['dmagTrailing']
        snr = 5.0 * np.power(10., 0.4*(magLimitWithTrailing - magObs))
        countObs[i] = np.where(snr >= snrLimit)[0].size
    return countObs



####

dbAddress = 'sqlite:///enigma_1189_sqlite.db'
ops = OpsimDatabase(dbAddress)

dbcols = ['expMJD', 'night', 'fieldRA', 'fieldDec', 'rotSkyPos', 'filter', 'finSeeing', 'fiveSigmaDepth']
simdata = ops.fetchMetricData(dbcols, sqlconstraint='')

orbits = pandas.read_table('pha20141031.des', sep=' ')
orbits = orbits.to_records()

rephs = pandas.read_table('phas_obs.txt', sep=' ')
ssoObs = joinObs(rephs, simdata)

ssos = np.unique(ssoObs['!!ObjID'])
nobsSsos = np.zeros([len(ssos), len(Hrange)], int)
outfile = open('nObs.txt', 'w')
for i, sso in enumerate(ssos):
    obs = ssoObs[np.where(ssoObs['!!ObjID'] == sso)]
    nobsSsos[i] = nObsMetric(obs)
    writestring = '%d' %(sso)
    for nobs in nobsSsos:
        writestring += ' %d' %(nobs)
    print >>outfile, writestring
outfile.close()


