import numpy as np
from itertools import repeat
import scipy
from scipy import interpolate
import pandas
import pyoorb as oo

import lsst.sims.maf.db as db
from lsst.sims.utils import haversine
import lsst.sims.maf.utils as utils


def packOorbElem(sso):
    oorbelems = [sso['!!ObjID'], sso['q'], sso['e'], np.radians(sso['i']), np.radians(sso['Omega/node']), 
             np.radians(sso['omega/argperi']), sso['t_p'],  2, sso['t_0'], 3, sso['magHv'], 0.15]
    oorbelems = np.column_stack(oorbelems)
    return oorbelems

def unpackEphs(oorbephem):
    # oorbephem = single ephem - i.e. for single object. 
    # then swap time/ephemeris element, so that we can get easy arrays of things like 'ra', 'dec'.
    oorbephem = np.swapaxes(oorbephem, 0, 1)
    dist = oorbephem[0]
    ra = oorbephem[1]
    dec = oorbephem[2]
    magV = oorbephem[3]
    time = oorbephem[4]
    dradt = oorbephem[6]
    ddecdt = oorbephem[7]
    phaseangle = oorbephem[8]
    solarelon = oorbephem[9]
    ephs = np.rec.fromarrays([time, ra, dec, dradt, ddecdt, dist, magV, phaseangle, solarelon], 
                             names=['time', 'ra', 'dec', 'dradt', 'ddecdt', 'dist', 'magV', 'phaseangle', 'solarelon'])
    return ephs

def interpolateEphs(ephs):
    interpfuncs = {}
    for n in ephs.dtype.names:
        if n == 'time':
            continue
        interpfuncs[n] = interpolate.interp1d(ephs['time'], ephs[n], kind='linear', assume_sorted=True, copy=False)
    return interpfuncs

def ssoInFov(interpfuncs, simdata, rFov=np.radians(1.75), raCol='fieldRA', decCol='fieldDec'):
    raSso = np.radians(interpfuncs['ra'](simdata['expMJD']))
    decSso = np.radians(interpfuncs['dec'](simdata['expMJD']))
    sep = haversine(raSso, decSso, simdata[raCol], simdata[decCol])
    tObs = simdata['expMJD'][np.where(sep<rFov)]
    return tObs


#####

#dbAddress = 'sqlite:///enigma_1189_sqlite.db'
#ops = db.OpsimDatabase(dbAddress)
dbAddress = 'sqlite:///opsim3_61_sqlite.db'
ops = db.OpsimDatabase(dbAddress, defaultdbTables=None, dbTables={'Summary':['Summary','obsHistID']})

dbcols = ['expMJD', 'night', 'fieldRA', 'fieldDec', 'rotSkyPos', 'filter', 'finSeeing', 'fiveSigmaDepth']
simdata = ops.fetchMetricData(dbcols, sqlconstraint='')

orbits = pandas.read_table('pha20141031.des', sep=' ')
orbits = orbits.to_records()

timestep = 2.0 / 24.0  # in days
timestart = simdata['expMJD'][0]
testlong = 10.0 # years
timeend = timestart + 365 * testlong + 1.0
times = np.arange(timestart, timeend + timestep/2.0, timestep)
# For pyoorb, we need to tag times with timescales;
# 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
ephTimes = np.array(zip(times, repeat(4, len(times))), dtype='double', order='F')

oo.pyoorb.oorb_init(ephemeris_fname="")

outfile = open('phas_obs_all.txt', 'w')
# Add Header.
writestring = '!!ObjID'
names=['time', 'ra', 'dec', 'dradt', 'ddecdt', 'dist', 'magV', 'phaseangle', 'solarelon']
for n in names:
    writestring += ' %s' %(n)
print >> outfile, writestring
for sso in orbits:
    oorbelems = packOorbElem(sso)
    oorbephems, err = oo.pyoorb.oorb_ephemeris(in_orbits = oorbelems, in_obscode=807, in_date_ephems=ephTimes)
    ephs = unpackEphs(oorbephems[0])
    interpfuncs = interpolateEphs(ephs)
    tvis = ssoInFov(interpfuncs, simdata)
    obs = np.recarray([len(tvis)], dtype=ephs.dtype)
    for n in interpfuncs:
        obs[n] = interpfuncs[n](tvis)
    obs['time'] = tvis
    for i in range(len(obs['time'])):
        writestring = '%s' %(sso['!!ObjID'])
        for n in names:
            writestring += ' %f' %(obs[n][i])
        print >> outfile, writestring
outfile.close()
