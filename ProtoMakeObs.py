import os, sys
import numpy as np
import matplotlib.pyplot as plt

from itertools import repeat
from scipy import interpolate

import pandas as pd
import pyoorb as oo

import lsst.sims.photUtils.Bandpass as Bandpass
import lsst.sims.photUtils.Sed as Sed

from lsst.sims.maf.db import OpsimDatabase
from lsst.sims.utils import haversine

# For the footprint generation and conversion between galactic/equatorial coordinates.
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.coordUtils import findChipName, observedFromICRS
from lsst.sims.catalogs.generation.db.ObservationMetaData import ObservationMetaData

mapper = LsstSimMapper()
camera = mapper.camera
epoch = 2000.0


def _updateColMap(colMap, outCol, alternativeNames, ssoCols):
    for a in alternativeNames:
        if a in ssoCols:
            colMap[outCol] = a
    return colMap

def readOrbits(orbitfile='pha20141031.des'):
    """
    Read the orbits from file 'orbitfile', returning a numpy structured array with the columns:
    'objID q e inc node argPeri tPeri epoch H g a meanAnom'
    """
    orbits = pd.read_table(orbitfile, sep='\s*', engine='python')
    # Normalize the column names, as different inputs tend to have some commonly-different names.
    ssoCols = orbits.columns.values.tolist()
    nSso = len(orbits)
    outCols = ['objId', 'q', 'e', 'inc', 'node', 'argPeri', 'tPeri', 'epoch', 'H', 'g', 'a', 'meanAnom']
    # Create mapping between column names read from disk and desired column names.
    colMap = {}
    for o in outCols:
        if o in ssoCols:
            colMap[o] = o
        else:
            # Try to find corresponding value
            if o == 'objId':
                alternatives = ['!!ObjID','objid', '!!OID']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'inc':
                alternatives = ['i']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'node':
                alternatives = ['BigOmega', 'Omega/node', 'Omega']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'argPeri':
                alternatives = ['argperi', 'omega/argperi']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'tPeri':
                alternatives = ['t_p', 'timeperi']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'epoch':
                alternatives = ['t_0']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'H':
                alternatives = ['magH', 'magHv', 'Hv']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'g':
                alternatives = ['phaseV', 'phase', 'gV']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'a':
                alternatives = ['semimajor']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)
            elif o == 'meanAnom':
                alternative = ['meanAnomaly', 'meananomaly', 'M']
                colMap = _updateColMap(colMap, o, alternatives, ssoCols)

    # Add the columns we can generate.
    if 'objId' not in colMap:
        orbids = np.arange(0, nSso, 1)
    else:
        orbids = orbits[colMap['objId']]
    if 'H' not in colMap:
        Hval = np.zeros(nSso) + 20.0
    else:
        Hval = orbits[colMap['H']]
    if 'g' not in colMap:
        gval = np.zeros(nSso) + 0.15
    else:
        gval = orbits[colMap['g']]

    # And some columns that can be converted from different orbit types.
    if 'a' not in colMap:
        aval = orbits[colMap['q']] / (1 + orbits[colMap['e']])
    else:
        aval = orbits[colMap['a']]
    period = np.sqrt(aval**3)
    if 'meanAnom' not in colMap:
        meanAnomval = 360.0*(orbits[colMap['epoch']] - orbits[colMap['tPeri']]) / (period*365.25)
    else:
        meanAnomval = orbits[colMap['meanAnom']]
    if 'q' not in colMap:
        qval = orbits[colMap['a']] * (1 - orbits[colMap['e']])
    else:
        qval = orbits[colMap['q']]
    if 'tPeri' not in colMap:
        tPerival = orbits[colMap['epoch']] - (orbits[colMap['meanAnom']]/360.0) * (period*365.25)
    else:
        tPerival = orbits[colMap['tPeri']]

    # Put it all together and turn it into a numpy array.
    orbits = np.rec.fromarrays([orbids, qval, orbits[colMap['e']], orbits[colMap['inc']],
                                 orbits[colMap['node']], orbits[colMap['argPeri']], tPerival,
                                 orbits[colMap['epoch']], Hval, gval, aval, meanAnomval],
                                names = outCols)
    return orbits



def setTimes(timestep=1., nyears=1., timestart=49353.):
    # Extend times beyond first/last observation, so that interpolation doesn't fail
    timestart = timestart - timestep
    timeend = timestart + 365 * nyears + timestep 
    times = np.arange(timestart, timeend + timestep/2.0, timestep)
    # For pyoorb, we need to tag times with timescales;
    # 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
    ephTimes = np.array(zip(times, repeat(4, len(times))), dtype='double', order='F')
    return ephTimes


def packOorbElem(sso):
    """
    Convert numpy structured array of orbital elements into the array OpenOrb needs as input.
    'sso' can be the orbital elements of a single object or of multiple objects.
    To normalize the column names to those expected here, read in data using 'readOrbits'.
    """
    if len(sso.shape) == 0:
        # Passed a single SSO
        nSso = 0
        sso0 = sso
        elem_type = 2
        epoch_type = 3
    else:
        nSso = len(sso)
        sso0 = sso[0]
        elem_type = np.zeros(nSso) + 2
        epoch_type = np.zeros(nSso) + 3
    # Check on orbit id type - pyoorb needs numbers not strings.
    if (isinstance(sso0['objId'], float) or isinstance(sso0['objId'], int)):
        orbids = sso['objId']
    else:
        if nSso == 0:
            orbids = 0
        else:
            orbids = np.arange(0, nSso, 1)
    # Convert to format for pyoorb, INCLUDING converting inclination, node, argperi to RADIANS
    oorbArray = np.column_stack((orbids, sso['q'], sso['e'], np.radians(sso['inc']),
                                 np.radians(sso['node']), np.radians(sso['argPeri']),
                                 sso['tPeri'], elem_type, sso['epoch'], epoch_type,
                                 sso['H'], sso['g']))
    return oorbArray


def unpackEphs(oorbephems, byObject=True):
    """
    Given oorb ephemeris array (shape = object / times / eph@time),
    Return an array aranged with
     columns = ['delta', 'ra', 'dec', 'mag', 'time', 'timescale', 'dradt', 'ddecdt', 'phase', 'solarelon']
     as the second
    grouped either by object (i.e. length of ra array == length of times) (default)
    or grouped by time (i.e. length of ra array == number of objects) (if byObject not true).
    """
    ephs = np.swapaxes(oorbephems, 2, 0)
    # oorbcols=['delta', 'ra', 'dec', 'magV', 'time', 'timescale', 'dradt', 'ddecdt', 'phase', 'solarelon']
    velocity = np.sqrt(ephs[6]**2 + ephs[7]**2)
    if byObject:
        ephs = np.swapaxes(ephs, 2, 1)
        velocity = np.swapaxes(velocity, 1, 0)
    # Create numpy structured array
    ephs = np.rec.fromarrays([ephs[0], ephs[1], ephs[2], ephs[3], ephs[4],
                              ephs[6], ephs[7], ephs[8], ephs[9], velocity],
                             names=['delta', 'ra', 'dec', 'magV', 'time', 'dradt',
                                    'ddecdt', 'phase', 'solarelon','velocity'])
    return ephs


# Linear interpolation
def interpolateEphs(ephs, i=0):
    interpfuncs = {}
    for n in ephs.dtype.names:
        if n == 'time':
            continue
        interpfuncs[n] = interpolate.interp1d(ephs['time'][i], ephs[n][i], kind='linear', 
                                              assume_sorted=True, copy=False)
    return interpfuncs

def ssoInFov(interpfuncs, simdata, rFov=np.radians(1.75), raCol='fieldRA', decCol='fieldDec'):
    """
    Return the indexes of the observations where the object could be seen.
    """
    raSso = np.radians(interpfuncs['ra'](simdata['expMJD']))
    decSso = np.radians(interpfuncs['dec'](simdata['expMJD']))
    sep = haversine(raSso, decSso, simdata[raCol], simdata[decCol])
    idxObs = np.where(sep<rFov)[0]
    return idxObs

def ssoInFovChip(interpfuncs, simdata, rFov=np.radians(1.75), raCol='fieldRA', decCol='fieldDec'):
    """
    Return the indexes of the observations where the object could be seen.
    """
    raSso = np.radians(interpfuncs['ra'](simdata['expMJD']))
    decSso = np.radians(interpfuncs['dec'](simdata['expMJD']))
    sep = haversine(raSso, decSso, simdata[raCol], simdata[decCol])
    idxObsRough = np.where(sep<rFov)[0]
    idxObs = []
    for idx in idxObsRough:
        mjd = simdata[idx]['expMJD']
        obs_metadata = ObservationMetaData(unrefractedRA=np.degrees(simdata[idx][raCol]),
                                           unrefractedDec=np.degrees(simdata[idx][decCol]),
                                           rotSkyPos=np.degrees(simdata[idx]['rotSkyPos']),
                                           mjd=simdata[idx]['expMJD'])
        raObj = np.radians(np.array([interpfuncs['ra'](simdata[idx]['expMJD'])]))
        decObj = np.radians(np.array([interpfuncs['dec'](simdata[idx]['expMJD'])]))
        raObj, decObj = observedFromICRS(raObj, decObj, obs_metadata=obs_metadata, epoch=epoch)
        chipNames = findChipName(ra=raObj,dec=decObj, epoch=epoch, camera=camera, obs_metadata=obs_metadata)
        if chipNames != [None]:
            idxObs.append(idx)
    idxObs = np.array(idxObs)
    return idxObs


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

def joinObs(objId, ephs, simdata, idxObs, sedname='C.dat', tol=1e-8, outfileName='out.txt',
            seeingcol='finSeeing', expTime='visitExpTime'):
    """
    Call for each object; write out the observations of each object.
    """
    # calculate extra values (color, trailing loss, and detection loss terms)
    # Calculate color term.
    dmagColor = np.zeros(len(idxObs), float)
    dmagDict = calcColors(sedname)
    filterlist = np.unique(simdata[idxObs]['filter'])
    for f in filterlist:
        if f not in dmagDict:
            raise UserWarning('Could not find filter %s in calculated colors!' %(f))
        match = np.where(simdata[idxObs]['filter'] == f)[0]
        dmagColor[match] = dmagDict[f]
    # Calculate trailing and detection loses.
    pm_LHS = np.array([  4.04959668e-04,   1.54946862e-03,  -7.32311940e-02,
                         4.13138884e-01,  -2.13021674e-01,   1.03493047e+00], float)
    pt_LHS = np.array([ -4.65513418e-04,   8.86629066e-03,  -6.20199786e-02,
                        1.81920172e-01, 1.86156205e-02,   9.97466212e-01], float)
    x = ephs['velocity'] * simdata['visitExpTime'][idxObs] / simdata['finSeeing'][idxObs] / 24.0
    trail = 2.5*np.log10(np.polyval(pt_LHS, x))
    detect = 2.5*np.log10(np.polyval(pm_LHS, x))
    dmagTrailing = np.where(trail<0, 0, trail)
    dmagDetect = np.where(detect<0, 0, detect)

    # Turn into a recarray so it's easier below.
    dmags = np.rec.fromarrays([dmagColor, dmagTrailing, dmagDetect],
                              names=['dmagColor', 'dmagTrailing', 'dmagDetect'])

    outCols = ['objId',] + list(ephs.dtype.names) + list(simdata.dtype.names) + list(dmags.dtype.names)
    if not os.path.isfile(outfileName):
        outfile = open(outfileName, 'w')
        writestring = '#'
        for col in outCols:
            writestring += '%s ' %(col)
        print >> outfile, writestring
    else:
        outfile = open(outfileName, 'a')

    # Write results.
    for eph, simdat, dm in zip(ephs, simdata[idxObs], dmags):
        writestring = '%s ' %(objId)
        for col in ephs.dtype.names:
            writestring += '%s ' %(eph[col])
        for col in simdat.dtype.names:
            writestring += '%s ' %(simdat[col])
        for col in dm.dtype.names:
            writestring += '%s ' %(dm[col])
        print >> outfile, writestring
    outfile.close()



if __name__ == '__main__':

    orbitfile = sys.argv[1]
    print 'working on %s' %orbitfile
    outfileName = orbitfile.replace('.des', '') + '_allObs.txt'


    orbits = readOrbits(orbitfile)

    nyears = 10
    ephTimes = setTimes(timestep=2./24., nyears=nyears)
    oo.pyoorb.oorb_init(ephemeris_fname="")

    dbAddress = 'enigma_1189_sqlite.db'
    ops = OpsimDatabase(dbAddress)

    dbcols = ['expMJD', 'night', 'fieldRA', 'fieldDec', 'rotSkyPos', 'filter',
              'finSeeing', 'fiveSigmaDepth', 'visitExpTime', 'solarElong']
    simdata = ops.fetchMetricData(dbcols, sqlconstraint='')

    for sso in orbits:
        oorbelems = packOorbElem(sso)
        oorbephems, err = oo.pyoorb.oorb_ephemeris(in_orbits = oorbelems, in_obscode=807, in_date_ephems=ephTimes)
        ephs = unpackEphs(oorbephems)
        interpfuncs = interpolateEphs(ephs)
        #idxObs = ssoInFov(interpfuncs, simdata)
        idxObs = ssoInFovChip(interpfuncs, simdata)
        if len(idxObs) > 0:
            tvis = simdata['expMJD'][idxObs]
            obs = np.recarray([len(tvis)], dtype=ephs.dtype)
            for n in interpfuncs:
                obs[n] = interpfuncs[n](tvis)
            obs['time'] = tvis
            joinObs(sso['objId'], obs, simdata, idxObs, outfileName=outfileName)
