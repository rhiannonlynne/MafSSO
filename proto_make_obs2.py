import sys
import numpy as np
from numpy.lib.recfunctions import merge_arrays
import matplotlib.pyplot as plt
from itertools import repeat
import scipy
from scipy import interpolate
import pandas
import pyoorb as oo


from lsst.sims.maf.db import OpsimDatabase
from lsst.sims.utils import haversine


# For the footprint generation and conversion between galactic/equatorial coordinates.
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.coordUtils import CameraCoords, AstrometryBase
from lsst.sims.catalogs.generation.db.ObservationMetaData import ObservationMetaData

mapper = LsstSimMapper()
camera = mapper.camera
myCamCoords = CameraCoords()
epoch = 2000.0
astrometryObject = AstrometryBase()

orbitfile = sys.argv[1]
print 'working on %s' %orbitfile


def packOorbElem(sso):
    ssoID = 0
    oorbelems = [ssoID, sso['q'], sso['e'], np.radians(sso['i']), np.radians(sso['Omega']),
             np.radians(sso['argperi']), sso['t_p'],  2, sso['t_0'], 3, 20., 0.15]
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

def ssoInFovChip(interpfuncs, simdata, rFov=np.radians(1.75), raCol='fieldRA', decCol='fieldDec'):
    raSso = np.radians(interpfuncs['ra'](simdata['expMJD']))
    decSso = np.radians(interpfuncs['dec'](simdata['expMJD']))
    sep = haversine(raSso, decSso, simdata[raCol], simdata[decCol])
    idxObs = np.where(sep<rFov)[0]
    tObs = []
    for idx in idxObs:
        mjd = simdata[idx]['expMJD']
        obs_metadata = ObservationMetaData(unrefractedRA=np.degrees(simdata[idx][raCol]),
                                           unrefractedDec=np.degrees(simdata[idx][decCol]),
                                           rotSkyPos=np.degrees(simdata[idx]['rotSkyPos']),
                                           mjd=simdata[idx]['expMJD'])
        raObj = np.radians(np.array([interpfuncs['ra'](simdata[idx]['expMJD'])]))
        decObj = np.radians(np.array([interpfuncs['dec'](simdata[idx]['expMJD'])]))
        raObj, decObj = astrometryObject.correctCoordinates(raObj, decObj, obs_metadata=obs_metadata, epoch=epoch)
        chipNames = myCamCoords.findChipName(ra=raObj,
                                             dec=decObj,
                                             epoch=epoch, camera=camera, obs_metadata=obs_metadata)
        if chipNames != [None]:
            tObs.append(simdata[idx]['expMJD'])
    tObs = np.array(tObs)
    return tObs


def calcColors(sedname='C.dat'):
    # Calculate SSO colors.
    import os
    import lsst.sims.photUtils.Bandpass as Bandpass
    import lsst.sims.photUtils.Sed as Sed
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

def joinObs(ephs, simdata, sedname='C.dat', tol=1e-8, outfile='out.txt'):
    outfile = open(jointObsfile, 'w')
    outnames = ['!!ObjID', 'time', 'ra', 'dec', 'dradt', 'ddecdt', 'dist', 'magV', 'phaseangle', 'solarelon',
                'expMJD', 'night', 'fieldRA', 'fieldDec', 'rotSkyPos', 'filter', 'finSeeing', 'fiveSigmaDepth',
                'dmagColor', 'dmagTrailing']
    writestring = ''
    for n in outnames:
        writestring += '%s ' %(n)
    print >> outfile, writestring

    dmagDict = calcColors(sedname)
    idxs = np.zeros(len(ephs), int)
    dmagColor = np.zeros(len(ephs), float)
    ### FIX THIS - search sorted? (?)
    for i, obs in enumerate(ephs):
        match = np.where(np.abs(simdata['expMJD'] - obs['time']) < tol) 
        if len(match) == 0 :
            continue
        match = match[0]
        if len(match) != 1:
            print 'Found %d simdata matches - expected only one' %(len(match))
            match = match[0:]
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

    # Write these out to disk.
    for eph, simdat, dmag in zip(ephs, simdata[idxs], dmags):
        writestring = ''
        for n in outnames:
            if n in ephs.dtype.names:
                writestring += '%s ' %(eph[n])
            elif n in simdat.dtype.names:
                writestring += '%s ' %(simdat[n])
            elif n in dmag.dtype.names:
                writestring += '%s ' %(dmag[n])
            else:
                print 'Could not find ', n
        print >> outfile, writestring
    outfile.close()

######

dbAddress = 'sqlite:///enigma_1189_sqlite.db'
ops = OpsimDatabase(dbAddress)
dbcols = ['expMJD', 'night', 'fieldRA', 'fieldDec', 'rotSkyPos', 'filter', 'finSeeing', 'fiveSigmaDepth']
simdata = ops.fetchMetricData(dbcols, sqlconstraint='')

ssoObsfile = orbitfile.replace('.des', '') + '_obs.txt'
jointObsfile = ssoObsfile.replace('_obs.txt', '_allObs.txt')

doVis = True
if doVis:

    orbits = pandas.read_table(orbitfile, sep='\s*')
    orbits = orbits.to_records()

    timestep =  3.0/24.0  # in days
    timestart = simdata['expMJD'].min()
    nyears = 10.0 # years
    timeend = timestart + 365 * nyears + 1.0
    timeend = simdata['expMJD'].max() + 1.0
    times = np.arange(timestart, timeend + timestep/2.0, timestep)
    print 'timestart', timestart
    print 'timeend', timeend
    # For pyoorb, we need to tag times with timescales;
    # 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
    ephTimes = np.array(zip(times, repeat(4, len(times))), dtype='double', order='F')

    oo.pyoorb.oorb_init(ephemeris_fname="")

    outfile = open(ssoObsfile, 'w')
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
        tvis = ssoInFovChip(interpfuncs, simdata)
        obs = np.recarray([len(tvis)], dtype=ephs.dtype)
        for n in interpfuncs:
            obs[n] = interpfuncs[n](tvis)
        obs['time'] = tvis
        for i in range(len(obs['time'])):
            writestring = '%s' %(sso['!!OID'])
            for n in names:
                writestring += ' %f' %(obs[n][i])
            print >> outfile, writestring
    outfile.close()

doJoin = True
if doJoin:

    rephs = pandas.read_table(ssoObsfile, sep='\s*')
    rephs = rephs.to_records()

    joinObs(rephs, simdata, outfile=jointObsfile)
