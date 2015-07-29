import os
import numpy as np
import pandas as pd

from itertools import repeat
import pyoorb as oo
from scipy import interpolate

import lsst.sims.photUtils.Bandpass as Bandpass
import lsst.sims.photUtils.Sed as Sed

from lsst.sims.utils import haversine, ObservationMetaData
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.coordUtils import findChipName, observedFromICRS

__all__ = ['MoOrbits', 'MoObs', 'runMoObs']

class MoOrbits(object):
    """
    This class stores the orbits for a set of moving objects (and reads them from disk).
    """
    def __init__(self):
        self.orbits = None
        self.ssoIds = None
        self.nSso = 0


    def _updateColMap(self, colMap, outCol, alternativeNames, ssoCols):
        """
        Private method used by readOrbits to update the map between data file and final orbit columns.

        This map is needed because data files can use different column headers for the same quantity.
        """
        for a in alternativeNames:
            if a in ssoCols:
                colMap[outCol] = a
        return colMap

    def readOrbits(self, orbitfile):
        """
        Read the orbits from file 'orbitfile', generating a numpy structured array with the columns:
        'objID q e inc node argPeri tPeri epoch H g a meanAnom'
        """
        self.orbitfile = orbitfile
        orbits = pd.read_table(orbitfile, delim_whitespace=True)
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
                    alternatives = ['!!ObjID', 'objid', '!!OID']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'inc':
                    alternatives = ['i']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'node':
                    alternatives = ['BigOmega', 'Omega/node', 'Omega']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'argPeri':
                    alternatives = ['argperi', 'omega/argperi']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'tPeri':
                    alternatives = ['t_p', 'timeperi']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'epoch':
                    alternatives = ['t_0']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'H':
                    alternatives = ['magH', 'magHv', 'Hv']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'g':
                    alternatives = ['phaseV', 'phase', 'gV']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'a':
                    alternatives = ['semimajor']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)
                elif o == 'meanAnom':
                    alternatives = ['M', 'meanAnomaly']
                    colMap = self._updateColMap(colMap, o, alternatives, ssoCols)

        # Add the columns we can generate with some guesses.
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

        # And some columns that can be generated from the input data we do have.
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

        # Put it all together into a dataframe.
        self.orbits = pd.DataFrame({'objId':orbids,
                                    'q':qval,
                                    'e':orbits[colMap['e']],
                                    'inc':orbits[colMap['inc']],
                                    'node':orbits[colMap['node']],
                                    'argPeri':orbits[colMap['argPeri']],
                                    'tPeri':tPerival,
                                    'epoch':orbits[colMap['epoch']],
                                    'H':Hval,
                                    'g':gval,
                                    'a':aval,
                                    'M':meanAnomval})
        self.ssoIds = np.unique(self.orbits['objId'])
        self.nSso = len(self.ssoIds)


class MoObs(MoOrbits):
    """
    Class to generate observations of a set of moving objects.
    Inherits from moOrbits (in order to read and save orbits).
    """

    def setTimes(self, timestep=1., ndays=1., timestart=49353.):
        """
        Generate and return an array for oorb of the ephemeris times desired.
        @ timestep : timestep for ephemeris generation (days)
        @ ndays    : number of days to generate ephemerides for (days)
        @ timestart : starting time of ephemerides (MJD)
        """
        # Extend times beyond first/last observation, so that interpolation doesn't fail
        timestart = timestart - timestep
        timeend = timestart + ndays + timestep
        times = np.arange(timestart, timeend + timestep/2.0, timestep)
        # For pyoorb, we need to tag times with timescales;
        # 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
        self.ephTimes = np.array(zip(times, repeat(4, len(times))), dtype='double', order='F')

    def _packOorbElem(self, sso=None):
        """
        Convert row from pandas dataframe (or numpy recarray) of orbital elements into the array OpenOrb needs as input.
        'sso' can be the orbital elements of a single object or of multiple objects.
        To normalize the column names to those expected here, read in data using 'readOrbits'.
        """
        # Python oorb element format:
        #
        # 0: orbitId  (cannot be a string)
        # 1 - 6: orbital elements, using radians for angles
        # 7: element type code
        #       2 = cometary - means timescale is TT, too
        #       3 = keplerians - timescale?
        # 8: epoch
        # 9: timescale for the epoch; 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
        # 10: magHv
        # 11: G
        #
        #  so we have to do a little translating from the orbits DataFrame to the elements we want in this array.
        if sso is None:
            sso = self.orbits
        # Do we have a single item (Series) or multiples (Dataframe)?
        if isinstance(sso, pd.Series):
            # Passed a single SSO in Series.
            nSso = 1
        elif isinstance(sso, pd.DataFrame):
            # Multiple SSO in dataframe.
            nSso = len(sso)
        else:
            if len(sso.shape) == 0:
                # Single SSO, in a numpy array.
                nSso = 1
            else:
                # Multiple SSSO in numpy array (or something else?).
                nSso = len(sso)
        if nSso == 1:
            orbids = 0
            elem_type = 2
            epoch_type = 3
        else:
            orbids = np.arange(0, nSso, 1)
            elem_type = np.zeros(nSso) + 2
            epoch_type = np.zeros(nSso) + 3
        # Convert to format for pyoorb, INCLUDING converting inclination, node, argperi to RADIANS
        oorbElem = np.column_stack((orbids, sso['q'], sso['e'], np.radians(sso['inc']),
                                     np.radians(sso['node']), np.radians(sso['argPeri']),
                                     sso['tPeri'], elem_type, sso['epoch'], epoch_type,
                                     sso['H'], sso['g']))
        return oorbElem

    def setupOorb(self):
        """
        Initialize oorb. (call once)
        """
        oo.pyoorb.oorb_init(ephemeris_fname="")

    def _generateOorbEphs(self, oorbArray, ephTimes=None, obscode=807):
        """
        Generate ephemerides.
        """
        if ephTimes is None:
            ephTimes = self.ephTimes
        oorbephems, err = oo.pyoorb.oorb_ephemeris(in_orbits = oorbArray, in_obscode=obscode, in_date_ephems=ephTimes)
        if err != 0:
            print 'Oorb returned error %s' %(err)
        return oorbephems

    def _unpackOorbEphs(self, oorbephems, byObject=True):
        """
        Given oorb ephemeris array (shape = object / times / eph@time),
        Return a numpy array aranged with
         columns = ['delta', 'ra', 'dec', 'mag', 'time', 'timescale', 'dradt', 'ddecdt', 'phase', 'solarelon']
         as the second
        grouped either by object (i.e. length of ra array == length of times) (default)
        or grouped by time (i.e. length of ra array == number of objects) (if byObject not true).
        """
        # Python oorb ephem array format:
        #   [objid][time][ephemeris information @ that time]
        # 0 = distance (geocentric distance)
        # 1 = ra (deg)
        # 2 = dec (deg)
        # 3 = mag
        # 4 = ephem mjd
        # 5 = ephem mjd timescale
        # 6 = dra/dt (deg/day) sky motion
        # 7 = ddec/dt (deg/day) sky motion
        # 8 = phase angle (deg)
        # 9 = solar elongation angle (deg)
        # So usually we want to swap the axes at least, so that instead of all the ephemeris information @ a particular time
        # being the accessible bit of information, we have all the RA values over time for a single object ('byObject')
        # Alternatively, we may want all the RA values for all objects at one time.
        #     This is also an option, by setting 'byObject' to False.
        ephs = np.swapaxes(oorbephems, 2, 0)
        # oorbcols=['delta', 'ra', 'dec', 'magV', 'time', 'timescale', 'dradt', 'ddecdt', 'phase', 'solarelon']
        velocity = np.sqrt(ephs[6]**2 + ephs[7]**2)
        if byObject:
            ephs = np.swapaxes(ephs, 2, 1)
            velocity = np.swapaxes(velocity, 1, 0)
        # Create a numpy recarray. We're deviating from the DataFrame here probably mostly due to history.
        ephs = np.rec.fromarrays([ephs[0], ephs[1], ephs[2], ephs[3], ephs[4],
                                  ephs[6], ephs[7], ephs[8], ephs[9], velocity],
                                  names=['delta', 'ra', 'dec', 'magV', 'time', 'dradt',
                                         'ddecdt', 'phase', 'solarelon','velocity'])
        return ephs

    def generateEphs(self, sso=None, ephTimes=None):
        """
        Combines several private methods to pack and unpack oorb-format orbital elements and ephemerides,
        into a single easy-access point to generate numpy arrays of ephemerides.

        Use the individual private methods if you want to unpack the ephemerides in a manner other than 'byObject'.
        """
        if sso is None:
            sso = self.orbits
        oorbelems = self._packOorbElem(sso=sso)
        oorbephs = self._generateOorbEphs(oorbelems, ephTimes=ephTimes)
        ephs = self._unpackOorbEphs(oorbephs, byObject=True)
        return ephs

    # Linear interpolation
    def interpolateEphs(self, ephs, i=0):
        """
        Generate linear interpolations between the quantities in ephs over time.
        """
        interpfuncs = {}
        for n in ephs.dtype.names:
            if n == 'time':
                continue
            interpfuncs[n] = interpolate.interp1d(ephs['time'][i], ephs[n][i], kind='linear',
                                                  assume_sorted=True, copy=False)
        return interpfuncs

    def _setupCamera(self):
        """
        Initialize camera mapper, etc.
        """
        self.mapper = LsstSimMapper()
        self.camera = self.mapper.camera
        self.epoch = 2000.0
        self.cameraFov=np.radians(2.1)

    def ssoInFov(self, interpfuncs, simdata, rFov=np.radians(1.75),
                 useCamera=True,
                 simdataRaCol = 'fieldRA', simdataDecCol='fieldDec'):
        """
        Return the indexes of the simdata observations where the object was inside the fov.
        """
        # See if the object is within 'rFov' of the center of the boresight.
        raSso = np.radians(interpfuncs['ra'](simdata['expMJD']))
        decSso = np.radians(interpfuncs['dec'](simdata['expMJD']))
        sep = haversine(raSso, decSso, simdata[simdataRaCol], simdata[simdataDecCol])
        if not useCamera:
            idxObsRough = np.where(sep<rFov)[0]
            return idxObsRough
        # Or go on and use the camera footprint.
        try:
            self.camera
        except AttributeError:
            self._setupCamera()
        idxObs = []
        idxObsRough = np.where(sep<self.cameraFov)[0]
        for idx in idxObsRough:
            mjd = simdata[idx]['expMJD']
            obs_metadata = ObservationMetaData(unrefractedRA=np.degrees(simdata[idx][simdataRaCol]),
                                               unrefractedDec=np.degrees(simdata[idx][simdataDecCol]),
                                               rotSkyPos=np.degrees(simdata[idx]['rotSkyPos']),
                                               mjd=simdata[idx]['expMJD'])
            raObj = np.radians(np.array([interpfuncs['ra'](simdata[idx]['expMJD'])]))
            decObj = np.radians(np.array([interpfuncs['dec'](simdata[idx]['expMJD'])]))
            raObj, decObj = observedFromICRS(raObj, decObj, obs_metadata=obs_metadata, epoch=self.epoch)
            chipNames = findChipName(ra=raObj,dec=decObj, epoch=self.epoch, camera=self.camera, obs_metadata=obs_metadata)
            if chipNames != [None]:
                idxObs.append(idx)
        idxObs = np.array(idxObs)
        return idxObs


    def _calcColors(self, sedname='C.dat'):
        """
        Calculate the colors for a moving object with sed 'sedname'.
        """
        # Do we need to read in the LSST bandpasses?
        try:
            self.lsst
        except AttributeError:
            filterdir = os.getenv('LSST_THROUGHPUTS_BASELINE')
            filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
            self.lsst ={}
            for f in filterlist:
                self.lsst[f] = Bandpass()
                self.lsst[f].readThroughput(os.path.join(filterdir, 'total_'+f+'.dat'))
            self.vband = Bandpass()
            self.vband.readThroughput('harris_V.dat')
            self.colors = {}
        # See if the sed's colors are in memory already.
        if sedname not in self.colors:
            moSed = Sed()
            moSed.readSED_flambda(sedname)
            vmag = moSed.calcMag(self.vband)
            self.colors[sedname] = {}
            for f in filterlist:
                self.colors[sedname][f] = moSed.calcMag(self.lsst[f]) - vmag
        return self.colors[sedname]


    def _calcMagLosses(self, velocity, seeing, texp=30.):
        """
        Calculate the magnitude losses due to trailing and not matching the point-source detection filter.
        """
        a_trail = 0.76
        b_trail = 1.16
        a_det = 0.42
        b_det = 0.00
        x = velocity * texp / seeing / 24.0
        dmagTrail = 1.25 * np.log10(1 + a_trail*x**2/(1+b_trail*x))
        dmagDetect = 1.25 * np.log10(1 + a_det*x**2 / (1+b_det*x))
        return dmagTrail, dmagDetect

    def _openOutput(self, outfileName):
        self.outfile = open(outfileName, 'w')
        self.wroteHeader = False

    def writeObs(self, objId, interpfuncs, simdata, idxObs, outfileName='out.txt',
                 sedname='C.dat', tol=1e-8,
                 seeingCol='finSeeing', expTimeCol='visitExpTime'):
        """
        Call for each object; write out the observations of each object.
        """
        # Return if there's nothing to write out.
        if len(idxObs) == 0:
            return
        # Open file if needed.
        try:
            self.outfile
        except AttributeError:
            self._openOutput(outfileName)
        # Calculate the ephemerides for the object, using the interpfuncs, for the times in simdata[idxObs].
        tvis = simdata['expMJD'][idxObs]
        ephs = np.recarray([len(tvis)], dtype=([('delta', '<f8'), ('ra', '<f8'), ('dec', '<f8'),
                                                ('magV', '<f8'), ('time', '<f8'), ('dradt', '<f8'), ('ddecdt', '<f8'),
                                                ('phase', '<f8'), ('solarelon', '<f8'), ('velocity', '<f8')]))
        for n in interpfuncs:
            ephs[n] = interpfuncs[n](tvis)
        ephs['time'] = tvis
        # Calculate the extra columns we want to write out (dmag due to color, trailing loss, and detection loss)
        # First calculate and match the color dmag term.
        dmagColor = np.zeros(len(idxObs), float)
        dmagColorDict = self._calcColors(sedname)
        filterlist = np.unique(simdata[idxObs]['filter'])
        for f in filterlist:
            if f not in dmagColorDict:
                raise UserWarning('Could not find filter %s in calculated colors!' %(f))
            match = np.where(simdata[idxObs]['filter'] == f)[0]
            dmagColor[match] = dmagColorDict[f]
        magFilter = ephs['magV'] + dmagColor
        # Calculate trailing and detection loses.
        dmagTrail, dmagDetect = self._calcMagLosses(ephs['velocity'], simdata[seeingCol][idxObs],
                                                    simdata[expTimeCol][idxObs])
        # Turn into a recarray so it's easier below.
        dmags = np.rec.fromarrays([magFilter, dmagColor, dmagTrail, dmagDetect],
                                  names=['magFilter', 'dmagColor', 'dmagTrail', 'dmagDetect'])

        outCols = ['objId',] + list(ephs.dtype.names) + list(simdata.dtype.names) + list(dmags.dtype.names)

        if not self.wroteHeader:
            writestring = ''
            for col in outCols:
                writestring += '%s ' %(col)
            self.outfile.write('%s\n' %(writestring))
            self.wroteHeader = True

        # Write results.
        for eph, simdat, dm in zip(ephs, simdata[idxObs], dmags):
            writestring = '%s ' %(objId)
            for col in ephs.dtype.names:
                writestring += '%s ' %(eph[col])
            for col in simdat.dtype.names:
                writestring += '%s ' %(simdat[col])
            for col in dm.dtype.names:
                writestring += '%s ' %(dm[col])
            self.outfile.write('%s\n' %(writestring))
        self.outfile.flush()


## Function to link the above class methods to generate an output file with moving object observations.

def runMoObs(orbitfile, outfileName, opsimfile,
            dbcols=None, tstep=2./24., nyears=None,
            rFov=np.radians(1.75), useCamera=True):

    from lsst.sims.maf.db import OpsimDatabase

    # Read orbits.
    moogen = MoObs()
    moogen.readOrbits(orbitfile)
    print "Read orbit information from %s" %(orbitfile)

    # Check rfov/camera choices.
    if useCamera:
        print "Using camera footprint"
    else:
        print "Not using camera footprint; using circular fov with %f degrees radius" %(np.degrees(rFov))

    # Read opsim database.
    opsdb = OpsimDatabase(opsimfile)
    if dbcols is None:
        dbcols = ['expMJD', 'night', 'fieldRA', 'fieldDec', 'rotSkyPos', 'filter',
                  'finSeeing', 'fiveSigmaDepth', 'visitExpTime', 'solarElong']
    else:
        reqcols = ['expMJD', 'night', 'fieldRA', 'fieldDec', 'rotSkyPos', 'filter',
                   'visitExpTime', 'finSeeing', 'fiveSigmaDepth']
        for col in reqcols:
            if col not in dbcols:
                dbcols.append(col)
    if nyears is not None:
        ndays = nyears * 365
        sqlconstraint = 'night<%d' %(ndays)
    else:
        ndays = 3650.0
        sqlconstraint = ''
    simdata = opsdb.fetchMetricData(dbcols, sqlconstraint=sqlconstraint)
    print "Queried data from opsim %s, fetched %f days worth of visits." %(opsimfile, ndays/365.)

    moogen.setTimes(timestep=tstep, ndays=ndays, timestart=simdata['expMJD'].min())
    print "Will generate ephemerides on grid of %f day timesteps, then extrapolate to opsim times." %(tstep)

    moogen.setupOorb()
    for i, sso in moogen.orbits.iterrows():
        ephs = moogen.generateEphs(sso)
        interpfuncs = moogen.interpolateEphs(ephs)
        idxObs = moogen.ssoInFov(interpfuncs, simdata, rFov=rFov, useCamera=useCamera)
        moogen.writeObs(sso['objId'], interpfuncs, simdata, idxObs, outfileName=outfileName)
    print "Wrote output observations to file %s" %(outfileName)

# Test example:
if __name__ == '__main__':
    runMoObs('pha20141031.des', 'test_allObs.txt', 'enigma_1189_sqlite.db', nyears=1)
