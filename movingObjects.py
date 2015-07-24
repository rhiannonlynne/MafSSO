"""

This class is intended to replace movingObject.py and
movingObjectList.py. (having two classes seemed too bulky for most
uses, and movingObject.py by itself wasn't actually used except for a
not-strictly-necessary data structure).

To simplify things, this class can hold either a single moving object
or many ... most moving object actions are more efficient when applied
to a large group of moving objects (such as calculating the
ephemerides of all the objects at a particular time, or the magnitudes
of many objects with very little variety in colors). The orbital
elements are not stored as records of a single object, but as
dictionaries of numpy arrays

To support uses beyond LSST's catalogs_generation, additional
functionality has been added to this class (such as reading orbits
from a file or database directly, as well as being able to be supplied
via a list, or supporting translation between common orbital element
formats). 

Data stored:
'orbits' - dictionary of numpy arrays for orbits, with keys appropriate for orbits:
   objid, q [perihelion distance AU], e [eccentricity], inc [inclination RADIANS],
   node [longitude of node RADIANS], argperi [argument of perihelion RADIANS],
   timeperi [time of perihelion MJD, _ORBITTIMEUNIT],
   H [absolute magnitude V], epoch [time of orbit MJD, _ORBITTIMEUNIT]
     but can also be given in Keplerian format (a/e/inc/Omega/omega,M,H,epoch). 
'ephems' = (aka 'positions') - dictionary of dictionaries (of numpy arrays) for
ephemeris positions, keyed to 'MJD',
     then for each time:
   ra, dec [DEGREES], magV, distance [AU]
   dra/dt (sky velocity in RA direction), ddec/dt (sky velocity in Dec) [DEGREES/DAY], 

Methods:
    def set_timeunits(self, orbittimeunit='TT', ephemtimeunit='TAI'):
    def set_orbits_Kep(self, objid, a, e, inc, node, argperi, M, magHv, epoch,
                        phaseV=None, objtype=None, sedname=None):
    def set_orbits_COM(self, objid, q, e, inc, node, argperi, timeperi, magHv, epoch,
                      phaseV=None, objtype=None, sedname=None):
    def get_orbitsDB(self, uiObj, objtype=None, nobjs=None, objid=None, sqlconstraint=None):
    def read_orbits(self, filename, fileformat='COM'):
    def write_orbits(self, outfilename, fileformat='COM'):
    def orbits_KeptoCOM(self):
    def orbits_COMtoKep(self):
    def pack_oorbArray(self):
    def unpack_oorbArray(self):
    def setup_oorb(self):
    def propagate_orbits(self, new_epoch, nbody=True):
    def make_mjdStr(self, mjd):
    def generate_ephemerides(self, mjdlist, obscode=807):
    def set_ephemeris(self, mjd, distance, ra, dec, magV, dradt, ddecdt, phase_angle=None):
    def calcDist_vincenty(self, RA1, Dec1, RA2, Dec2):
    def isInFieldOfView(self, mjd, ra_center, dec_center, radius_fov=1.8, make_new=False):
    def setup_calcMags(self, filterlist, rootSEDdir, rootFILTERdir=None):
    def calcMags(self, mjd, filt, m5=None, withErrors=True):
    def convert_radec2ec(self, mjd):
    def calc_phasefunction(self, mjdlist, gval=0.15):
    def calc_distanceE(self):

"""

"""
import time
def dtime(time_prev):
   return (time.time() - time_prev, time.time())
"""

import copy
import numpy
import pyoorb as oo
import inputOutput as io

## constants that are useful.
_deg2rad = numpy.pi/180.0
_rad2deg = 180.0/numpy.pi
# apparent r band magnitude of the sun. 
# see http://www.ucolick.org/~cnaw/sun.html for apparent magnitudes in other bands.
_mag_sun = -27.1
_km_per_au = 1.496e8

## keys for the class dictionaries.
# Dictionary keys for orbits 
_COM_FORMAT = ['objid', 'q', 'e', 'inc', 'node', 'argperi', 'timeperi', 'magHv', 'epoch',
               'phaseV', 'objtype', 'sedname']
_COM_FORMAT_REQ = ['objid', 'q', 'e', 'inc', 'node', 'argperi', 'timeperi', 'magHv', 'epoch', 'phaseV']
_COM_FORMAT_OTHER = ['objtype', 'sedname']
_COM_FILEFORMAT = ['objid', 'format', 'q', 'e', 'inc', 'node', 'argperi', 'timeperi', 'magHv', 'epoch']
_KEP_FILEFORMAT = ['objid', 'a', 'e', 'inc', 'node', 'argperi', 'M', 'magHv', 'epoch']
# Omega == longitude of node (node)
# omega = argument of perihelion (argperi)
# ALL ELEMENT ANGLES ARE IN DEGREES.


## for OpenOrb.
# UTC/UT1/TT/TAI = 1/2/3/4
_TIMEUNITS = {'UTC':1, 'UT1':2, 'TT':3, 'TAI':4}
# orbittimeunit is the unit of time used for orbital parameters.
# ephemtimeunit is the unit of time for the ephemerides.

class MovingObjects:
    """ Class for holding orbit and ephemeris information on many or a single moving object. """
    def __init__(self, orbits=None, orbittimeunit='TT', ephemtimeunit='TAI'):
        """Initialize movingObjects - which can contain many moving objects orbits' and ephemerides or just one."""
        # But, the initial orbits can be passed in, in a dictionary format
        # This is the best way to get orbital information in, if querying from a database.
        #  (use useful_input.py and query db externally, then pass resulting dictionary back in).
        self.orbits=None
        self.oorbArray=None
        if orbits != None:
            # Test that required COM format elements are present in dictionary.
            for key in _COM_FORMAT_REQ:
                if key not in orbits.keys():
                    raise Exception("KeyError: Necessary orbit key %s not present in input dictionary" %(key))
            self.orbits = copy.deepcopy(orbits)
        # The initial ephemerides are always an empty dictionary.
        self.ephems = {}
        # Set time units.
        self.set_timeunits(orbittimeunit=orbittimeunit, ephemtimeunit=ephemtimeunit)
        return

    def set_timeunits(self, orbittimeunit='TT', ephemtimeunit='TAI'):
        """Set the units of time for the object; orbit time units (default TT) and ephem time units (default TAI)."""
        self.orbittimeunit = orbittimeunit
        self.ephemtimeunit = ephemtimeunit
        return

    ## Some methods to get data into the orbits dictionary.
    def get_orbitsDB(self, uiObj, objtype=None, nobjs=None, objid=None, sqlconstraint=None):
        """Get [nobjs] orbits of [objtype] from standard cometary format orbit database.
        Can also specify objid - which can be a list or numpy array."""
        # This method is convenient, but for more complex orbit DB queries just pass in the final dictionary to init.
        # Unset oorb array.
        self.oorbArray = None
        # Go get the data .. on wookie, is mysql db (db = ssm, table=ssm_orbits).
        #  and on tusken it is a postgres db (db=mops_y5, but table=ssm_orbits
        #  and on fatboy, it is sqlserver & orbits_[x]
        #    where x = 50813, 508443, 50873, 509093, 50933, 50963, 50993, 51023, 51053, 51083, 5113, 51143, 51173
        # Build sql query to database, using the keys above.
        if uiObj.engine.name == 'mssql':
            mssql_keys = ['objid', 'q' ,'e', 'inclination', 'node', 'argperi', 't_peri',
                          'h_g', 'epoch', 'phasegv', 'otype', 'sed_filename']
            # If we are just looking for a number or range of objects ...
            if objid == None:
                query = 'select '
                if nobjs != None:
                    query = query + ' top %d' %(nobjs)
                for k in mssql_keys:
                    query = query + ' %s,' %(k)
                # Remove the last comma, for valid sql.
                query = query[:-1]
                # Add the rest of the query.
                query = query + ' from orbits_50813'
                if objtype!=None:
                    query = query + ' where otype=%d' %(objtype)
                if sqlconstraint != None:
                    if objtype!=None:
                        query = query + ' and %s' %(sqlconstraint)
                    else:
                        query = query + ' where %s' %(sqlconstraint)
            # Otherwise, we were looking for the orbit for a particular objid (or objids)
            else:
                query = 'select '
                for k in mssql_keys:
                    query = query + ' %s,' %(k)
                # Remove the last comma, for valid sql.
                query = query[:-1]
                # Add the rest of the query.
                query = query + ' from orbits_50813'
                # And put in the objid restriction.
                if isinstance(objid, int) | isinstance(objid, float):
                    query = query + ' where objid=%d' %(objid)
                else:
                    if isinstance(objid, numpy.ndarray):
                        query = query + ' where objid in %s' %(objid.tolist())
                    else:
                        query = query + ' where objid in %s' %(objid)
        else:
            keys = ['objid', 'q', 'e', 'i', 'node', 'argperi', 't_peri', 'h_v', 'epoch',
                    'g_value', 'objtype', 'sed_filename']
            # If we are just looking for a number or range of objects ...
            if objid==None:
                query = 'select '
                for k in keys:
                    query = query + ' %s,'%(k)
                query = query[:-1]
                query = query + ' from ssm_orbits'
                if objtype != None:
                    query = query + ' where objtype="%d"' %(objtype)
                if sqlconstraint != None:
                    if objtype !=None:
                        query = query + ' and %s' %(sqlconstraint)
                    else:
                        query = query + ' where %s ' %(sqlconstraint)
            if nobjs != None:
                query = query + ' limit %d' %(nobjs)
            # Otherwise we're looking for a specific objid or objids.
            else:
                query = 'select '
                for k in keys:
                    query = query + ' %s,' %(k)
                # Remove the last comma, for valid sql.
                query = query[:-1]
                # Add the rest of the query.
                query = query + ' from ssm_orbits'
                # And put in the objid restriction
                if isinstance(objid, int) | isinstance(objid, float):
                    query = query + ' where objid=%d' %(objid)
                else:
                    if isinstance(objid, numpy.ndarray):
                        query = query + ' where objid in %s' %(objid.tolist())
                    else:
                        query = query + ' where objid in %s' %(objid)
        #print query
        results = uiObj.sqlQuery(query)
        orbitkeys = _COM_FORMAT
        self.orbits = uiObj.assignResults(results, keys=orbitkeys)
        return

    def set_orbits_Kep(self, objid, a, e, inc, node, argperi, M, magHv, epoch,
                      phaseV=None, objtype=None, sedname=None):
        """Set up the orbits dictionary, using input arrays for Keplerian orbital parameters.

        Translates to COM format"""
        # Check that the inputs are all the same length.
        numobjs = len(objid)
        if ((len(a)!=numobjs) | (len(e)!=numobjs) | (len(inc)!=numobjs) | (len(node)!=numobjs)
            | (len(argperi)!=numobjs) | (len(M)!=numobjs) | (len(magHv)!=numobjs) | (len(epoch)!=numobjs)):
            raise Exception('ValueError: Input values must have the same length.')
        # Unset oorb array.
        self.oorbArray = None
        # Set other values.
        self.orbits = {}
        self.orbits['objid'] = numpy.copy(objid)
        self.orbits['inc'] = numpy.copy(inc)
        self.orbits['node'] = numpy.copy(node)
        self.orbits['argperi'] = numpy.copy(argperi)
        self.orbits['magHv'] = numpy.copy(magHv)
        self.orbits['epoch'] = numpy.copy(epoch)
        self.orbits['e'] = numpy.copy(e)
        # Save these 'keplerian' elements, even though we will translate them.
        self.orbits['a'] = numpy.copy(a)
        self.orbits['M'] = numpy.copy(M)
        self.orbits_KeptoCOM()
        if phaseV == None:
            self.orbits['phaseV'] = numpy.zeros(len(self.orbits['q']), 'float') + 0.15
        else:
            if (len(phaseV)==numobjs):
                self.orbits['phaseV'] = numpy.copy(phaseV)
            elif (isinstance(phaseV, 'float')):
                self.orbits['phaseV'] = numpy.ones(len(self.orbits['q'], 'float')) * phaseV
            else:
                raise Exception('ValueError: PhaseV must have same length as other input elements.')
        if objtype != None:
            self.orbits['objtype'] = numpy.copy(objtype)
        if sedname != None:
            self.orbits['sedname'] = numpy.copy(sedname)
        return

    def set_orbits_COM(self, objid, q, e, inc, node, argperi, timeperi, magHv, epoch,
                       phaseV=None, objtype=None, sedname=None):
        """Set up the orbits dictionary, using input arrays for COM format orbital parameters."""
        # Check that the inputs are all the same length.
        numobjs = len(objid)
        if ((len(q)!=numobjs) | (len(e)!=numobjs) | (len(inc)!=numobjs) | (len(node)!=numobjs) |
            (len(argperi)!=numobjs) | (len(timeperi)!=numobjs) | (len(magHv)!=numobjs) | (len(epoch)!=numobjs)):
            raise Exception('ValueError: Input values must have the same length.')
        # Unset oorb array.
        self.oorbArray = None
        # Set other values.
        self.orbits = {}
        self.orbits['objid'] = numpy.copy(objid)
        self.orbits['q'] = numpy.copy(q)
        self.orbits['e'] = numpy.copy(e)
        self.orbits['inc'] = numpy.copy(inc)
        self.orbits['node'] = numpy.copy(node)
        self.orbits['argperi'] = numpy.copy(argperi)
        self.orbits['timeperi'] = numpy.copy(timeperi)
        self.orbits['magHv'] = numpy.copy(magHv)
        self.orbits['epoch'] = numpy.copy(epoch)
        if phaseV == None:
            self.orbits['phaseV'] = numpy.zeros(len(self.orbits['q']), 'float') + 0.15
        else:
            if (len(phaseV)==numobjs):
                self.orbits['phaseV'] = numpy.copy(phaseV)
            elif (isinstance(phaseV, float)):
                self.orbits['phaseV'] = numpy.ones(len(self.orbits['q'], 'float')) * phaseV
            else:
                raise Exception('ValueError: PhaseV must have same length as other input elements.')
        if objtype != None:
            self.orbits['objtype'] = numpy.copy(objtype)
        if sedname != None:
            self.orbits['sedname'] = numpy.copy(sedname)
        return

    def read_orbits(self, filename, fileformat='COM'):
        """Read the orbital information from a file, setting up the orbits dictionary."""
        # Unset oorbArray.
        self.oorbArray = None
        ui = io.InputOutput()
        # If file on disk is in COM format (aka PS SSM files)
        if fileformat=='COM':
            print 'Reading COM format file from %s' %(filename)
            # This format line is intended to match PS SSM files, so there is a 'format' column
            #   that we will later dump.
            # _COM_FILEFORMAT = ['objid', 'format', 'q', 'e', 'inc', 'node', 'argperi', 'timeperi', 'magHv', 'epoch']
            keys = _COM_FILEFORMAT
            # Specify expected data types in each column.
            keytypes = ['string', 'string', 'double', 'double', 'double', 'double',
                        'double', 'double', 'float', 'float']
            # Read the data.
            self.orbits = ui.readDatafile(filename, keys, keytypes)
            # Delete the unnecessary dictionary key.
            del self.orbits['format']
            # And add the necessary phaseV information. Set to default G value in solar system.
            self.orbits['phaseV'] = numpy.zeros(len(self.orbits['q']), 'float') + 0.15
        # Otherwise, the file on disk is in Keplerian format, and we'll assume it does not have this 'format' line.
        elif ((fileformat=='KEP') | (fileformat=='Kep')):
            print 'Reading KEP format file from %s' %(filename)
            # _KEP_FILEFORMAT = ['objid', 'a', 'e', 'inc', 'node', 'argperi', 'M', 'magHv', 'epoch']
            keys = _KEP_FILEFORMAT
            keytypes = ['string', 'double', 'double', 'double', 'double', 'double', 'double', 'float', 'float']
            # Read the data.
            self.orbits = ui.readDatafile(filename, keys, keytypes)
            # Translate to COM format.
            self.orbits_KeptoCOM()
            # And add the necessary phaseV information. Set to default G value.
            self.orbits['phaseV'] = numpy.zeros(len(self.orbits['q']), 'float') + 0.15
        else:
            print 'Did not understand fileformat %s' %(fileformat)
        return

    def write_orbits(self, outfilename, fileformat='COM'):
        """Write out a COM format file containing these orbital elements."""
        # Open output file.
        f = open(outfilename, 'w')
        # Write in Keplerian orbit format. (Might have to translate first).
        if ((fileformat=='Kep') | (fileformat=='KEP')):
            try:
                self.orbits['a']
            except AttributeError:
                self.orbits_COMtoKep()
            print >>f, "!!ObjID FORMAT a e i Omega/node omega/argperi t_p magHv t_0 "
            for i in range(len(self.orbits['a'])):
                print >>f, "%s %s %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e" \
                       %( i, "KEP", self.orbits['a'][i], self.orbits['e'][i],
                          self.orbits['inc'][i], self.orbits['node'][i],
                          self.orbits['argperi'][i], self.orbits['timeperi'][i],
                          self.orbits['magHv'][i], self.orbits['epoch'][i])
        # Else, write orbits in COM orbit format.
        else:
            print >>f, "!!ObjID FORMAT q e i Omega/node omega/argperi t_p magHv t_0 INDEX N_PAR MOID COMPCODE"
            for i in range(len(self.orbits['q'])):
                print >>f, "%s %s %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %d %d %d %s" \
                      %(i, "COM", self.orbits['q'][i], self.orbits['e'][i], self.orbits['inc'][i], \
                        self.orbits['node'][i], self.orbits['argperi'][i], self.orbits['timeperi'][i], \
                        self.orbits['magHv'][i], self.orbits['epoch'][i], 1, 6, 10, 'PYOORB')
        f.close()
        return

    ## Translate orbital elements between COM and KEP formats.
    def orbits_KeptoCOM(self):
        """Translate orbital elements from a/e/M/epoch to COM format elements,

        returning q/e/TimePeri/epoch."""
        # q = a(1-e)
        self.orbits['q'] = self.orbits['a']*(1-self.orbits['e'])
        # mean anomaly = 0 at perihelion
        # M / 2pi = (t - T) / P
        #  where M = mean anomaly at time t (epoch), and T = time of perihelion
        #  P = period = sqrt(a**3) [years .. so must convert to days]
        # So, T (time of perihelion) = t - (M/2pi)*P
        P = numpy.sqrt(self.orbits['a']**3)
        self.orbits['timeperi'] = self.orbits['epoch'] - (self.orbits['M']/360.0)*(P*365.25)
        return

    def orbits_COMtoKep(self):
        """Translate orbital elements from q/e/T_peri/epoch to Keplerian format elements.

        Returns a/e/M/epoch."""
        self.orbits['a'] = self.orbits['q'] / (1-self.orbits['e'])
        P = numpy.sqrt(self.orbits['a']**3)
        self.orbits['M'] = 360.0*(self.orbits['epoch']-self.orbits['timeperi'])/(P*365.25)
        return

    ## Functions for propagating orbits or predicting positions.
    def pack_oorbArray(self):
        """Translate orbital element dictionary (easy for humans) into pyoorb-suitable input orbit array."""
        # Translate orbital elements into array that pyoorb will like.
        # PyOrb wants ::
        # 0: orbitId
        # 1 - 6: orbital elements, using radians for angles
        # 7: element type code, where 2 = cometary - means timescale is TT, too
        # 8: epoch
        # 9: timescale for the epoch; 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
        # 10: magHv
        # 11: G
        elem_type = numpy.zeros(len(self.orbits['q'])) + 2
        epoch_type = numpy.zeros(len(self.orbits['q'])) + _TIMEUNITS[self.orbittimeunit]
        # Also, the orbitID has to be a float, rather than a string, so substitute if needed.
        if ((isinstance(self.orbits['objid'][0], float) == True) |
            (isinstance(self.orbits['objid'][0], int) == True)):
            orbids = self.orbits['objid']
        else:
            orbids = numpy.arange(0, len(self.orbits['objid']), 1)
        # Convert to format for pyoorb, INCLUDING converting inclination, node, argperi to RADIANS
        # Note that the first column *should* be objid ... but pyoorb does not use it, and if this is
        #  a string, it will actually cause poorb to fail. SO, substitute something else.
        self.oorbArray = numpy.column_stack((orbids, self.orbits['q'], self.orbits['e'], self.orbits['inc']*_deg2rad,
                                              self.orbits['node']*_deg2rad, self.orbits['argperi']*_deg2rad,
                                              self.orbits['timeperi'],
                                              elem_type, self.orbits['epoch'], epoch_type,
                                              self.orbits['magHv'], self.orbits['phaseV']))
        return

    def unpack_oorbArray(self):
        """Translate pyoorb-style orbit array back into orbital element dictionary (easier for humans)."""
        t = numpy.swapaxes(self.oorbArray, 0, 1)
        self.orbits['q'] = t[1]
        self.orbits['e'] = t[2]
        self.orbits['inc'] = t[3] * _rad2deg
        self.orbits['node'] = t[4] * _rad2deg
        self.orbits['argperi'] = t[5] * _rad2deg
        self.orbits['timeperi'] = t[6]
        self.orbits['epoch'] = t[8]
        self.orbits['magHv'] = t[10]
        self.orbits['phaseV'] = t[11]
        return

    def setup_oorb(self):
        oo.pyoorb.oorb_init(ephemeris_fname="")
        return

    def propagate_orbits(self, new_epoch, nbody=True):
        """Propagate (all) orbits from previous epoch (orbits['epoch']) to new_epoch."""
        # Check to see if oorbArray is ready.
        if (self.oorbArray == None):
            print "# Setting up oorbArray."
            self.pack_oorbArray()
        # Set up new orbit epoch information (including timescale).
        orb_epoch = numpy.zeros([1,2], dtype='double', order='F')
        orb_epoch[0][:] = [new_epoch, _TIMEUNITS[self.orbittimeunit]]
        # Now propagate orbits.
        #oo.pyoorb.oorb_init()
        if nbody:
            newoorbArray, err = oo.pyoorb.oorb_propagation_nb(in_orbits = self.oorbArray,
                                                              in_epoch = orb_epoch)
        else:
            newoorbArray, err = oo.pyoorb.oorb_propagation_2b(in_orbits=self.oorbArray,
                                                              in_epoch = orb_epoch)
        self.oorbArray = newoorbArray
        self.unpack_oorbArray()
        return

    def make_mjdStr(self, mjd):
        """ Convert float/number mjdTai to string so can use for dictionary lookup.

        mjdTai should be a single floating point number, which is then returned as fixed-format string."""
        # A utility function, for more reliably finding dictionary keys for ephemeris dates.
        mjdstr = "%.8f" %(mjd)
        return mjdstr

    def generate_ephemerides(self, mjdlist, obscode=807):
        """Generate ephemerides(positions) for all objects at times in mjdTaiList, for observatory 'obscode'.

        Resulting ephemerides are stored in dictionaries, keyed to 'time' and then values ('ra', 'dec', etc.).
        Currently, all dates are ASSUMED to be MJD_TAI  (need to change 'timescale' if not)
        """
        # Comment on speed of this method::
        #   Some tests have shown that generating ephemerides with pyoorb (OO v1.0.1) scales linearly with
        #   both the number of objects which need ephemerides and the number of mjd's
        #   (within a fixed MJD time period) that are requested.
        #   Increasing the time period between the first/last ephemeride is also a linear scaling.
        #   Also, within this function, at least up to the point where memory might become swapped,
        #   the VAST majority of time is taken doing the actual ephemeride generation.
        #   For example: for 400 objects with 101 mjd's within 5 days, the time required
        #    to set up for oorb was 0.001s, the time required to generate the ephemerides was 14s,
        #    and the time to swap the axes around/set up ephem dictionaries was 0.001s (on my macbook pro).
        #   So : ~3-3.5e-4 s / ephemeride.
        #   For ~1000 objects with 100 mjd's (within one night) ~ 30s  (10 years -> 30 hrs)
        #       or 100 objects with 1000 mjd's (within one night) ~30s  (10 years -> 30 hours)
        #
        # Convert float mjdTai's into strings for dictionary lookup
        mjdlistStr = []
        if ((isinstance(mjdlist, float) == True) | (isinstance(mjdlist, int) ==True)):
            mjdlist = [mjdlist]
        for mjd in mjdlist:
            mjdlistStr.append(self.make_mjdStr(mjd))
        # Set up array to hold ephemeris date information for all moving objects (the dates you want positions predicted).
        ephem_dates = numpy.zeros([len(mjdlist),2], dtype='double', order='F')
        for i in range(len(mjdlist)):
            ephem_dates[i][:] = [mjdlist[i], _TIMEUNITS[self.ephemtimeunit]]
        # Make sure that oorbArray is ready.
        if (self.oorbArray is None):
            self.pack_oorbArray()
        # Now do ephemeris generation for all objects on all dates
        #oo.pyoorb.oorb_init(ephemeris_fname="")
        newephems, err = oo.pyoorb.oorb_ephemeris(in_orbits = self.oorbArray,
                                                  in_obscode = obscode,
                                                  in_date_ephems = ephem_dates)
        if (err != 0):
            raise Exception("pyoorb.oorb_ephemeris encountered an error")
        # Returned ephems contain a 3-D Fortran array of ephemerides, the axes are:
        #   [objid][time][ephemeris information element]
        # the ephemeris information elements are (in order):
        # distance, ra, dec, mag, ephem mjd, ephem mjd timescale, dradt(sky), ddecdt(sky)
        # per object, per date, 8 elements (array shape is OBJ(s)/DATE(s)/VALUES)
        # Note that ra/dec, dradt, etc. are all in DEGREES.
        # First: (to arrange ephems for easier later use)
        # Swap the order of the axes: DATE / Objs / values
        e = numpy.swapaxes(newephems, 0, 1)
        i = 0
        for mjdstr in mjdlistStr:
            # Swap again so we can get easy arrays of 'ra', etc.
            eph = numpy.swapaxes(e[i], 0, 1)
            i = i + 1
            self.ephems[mjdstr] = {}
            self.ephems[mjdstr]['objid'] = self.orbits['objid']
            self.ephems[mjdstr]['dist'] = eph[0]
            self.ephems[mjdstr]['ra'] = eph[1]
            self.ephems[mjdstr]['dec'] = eph[2]
            self.ephems[mjdstr]['magV'] = eph[3]
            self.ephems[mjdstr]['dradt'] = eph[6]
            self.ephems[mjdstr]['ddecdt'] = eph[7]
            self.ephems[mjdstr]['phase_angle'] = eph[8]
            self.ephems[mjdstr]['solar_elon'] = eph[9]
        # done calculating ephemerides. ephemerides stored in dictionary with movingObjects
        return

    def set_ephemeris(self, mjd, distance, ra, dec, magV, dradt, ddecdt, phase_angle=None, solar_elon=None):
        """Set values for an ephemeris."""
        mjdstr = self.make_mjdStr(mjd)
        self.ephems[mjdstr] = {}
        self.ephems[mjdstr]['dist'] = numpy.copy(distance)
        self.ephems[mjdstr]['ra'] = numpy.copy(ra)
        self.ephems[mjdstr]['dec'] = numpy.copy(dec)
        self.ephems[mjdstr]['magV'] = numpy.copy(magV)
        self.ephems[mjdstr]['dradt'] = numpy.copy(dradt)
        self.ephems[mjdstr]['ddecdt'] = numpy.copy(ddecdt)
        if phase_angle!=None:
            self.ephems[mjdstr]['phase_angle'] = numpy.copy(phase_angle)
        if solar_elon!=None:
            self.ephems[mjdstr]['solar_elon'] = numpy.copy(solar_elon)
        return

    def calcDist_vincenty(self, RA1, Dec1, RA2, Dec2):
        """Calculates distance on a sphere using the Vincenty formula.
        Give this function RA/Dec values in radians. Returns angular distance(s), in radians.
        Note that since this is all numpy, you could input arrays of RA/Decs."""
        D1 = (numpy.cos(Dec2)*numpy.sin(RA2-RA1))**2 + \
            (numpy.cos(Dec1)*numpy.sin(Dec2) - \
             numpy.sin(Dec1)*numpy.cos(Dec2)*numpy.cos(RA2-RA1))**2
        D1 = numpy.sqrt(D1)
        D2 = (numpy.sin(Dec1)*numpy.sin(Dec2) + \
              numpy.cos(Dec1)*numpy.cos(Dec2)*numpy.cos(RA2-RA1))
        D = numpy.arctan2(D1,D2)
        return D

    def isInFieldOfView(self, mjd, ra_center, dec_center, radius_fov=1.8, make_new=False):
        """Calculate distance between each ephemeris and the center of the fov (all in DEG)

        Sets true/false array in self.ephems[mjd]['fov'] that can be used as a mask.
        If make_new=True, then will return a new movingObjects object which only
         contains objects within fov."""
        mjdstr = self.make_mjdStr(mjd)
        try:
            self.ephems[mjdstr]
        except AttributeError:
            warning.warn('Warning: movingObjects does not have ephemeris on this date; will compute.')
            self.generate_ephemerides(mjd)
        # Calculate the distance between all objects and the center of the field of view.
        dist = self.calcDist_vincenty(self.ephems[mjdstr]['ra']*_deg2rad,
                                      self.ephems[mjdstr]['dec']*_deg2rad,
                                      ra_center*_deg2rad, dec_center*_deg2rad)
        dist = dist * _rad2deg
        # Set True/False values for those within fov.
        self.ephems[mjdstr]['fov'] = numpy.where(dist<radius_fov, True, False)
        # If want to return entire new movingObjects instance with the objects in fov only:
        if make_new:
            newMO = movingObjects()
            condition = self.ephems[mjdstr]['fov']
            # Set up orbital information in new movingObjects.
            newMO.setOrbits_COM(self.orbits['objid'][condition], self.orbits['q'][condition],
                                self.orbits['e'][condition],
                                self.orbits['inc'][condition],
                                self.orbits['node'][condition],
                                self.orbits['argperi'][condition],
                                self.orbits['timeperi'][condition],
                                self.orbits['magV'][condition],
                                self.orbits['epoch'][condition],
                                self.orbits['phaseV'][condition],
                                self.orbits['objtype'][condition],
                                self.orbits['sedname'][condition])
            # Set up ephemeris information (for this mjd) in new movingObjects.
            newMO.setEphemeris(mjd, self.ephems[mjdstr]['dist'], self.ephems['ra'],
                               self.ephems['dec'], self.ephems['magV'],
                               self.ephems['dradt'], self.ephems['ddecdt'],
                               self.ephems['phase_angle'])
            return newMO
        # Otherwise, just return as normal.
        return

    def setup_calcMags(self, filterlist, rootSEDdir, rootFILTERdir=None):
        """ Calculate the magnitude of all objects in the movingObjects at time mjd, in a particular filter. 

        This uses the magV at a particular time to generate the color in a particular filter,
         given the sedname of each object."""
        # Set up data needed to calculate magnitudes for each moving object.
        # First find location of lsst filter throughput curves, if not specified.
        import os
        if rootFILTERdir == None:
            # (assumes throughputs is setup)
            # Just use the default directory of the throughputs as this is probably correct.
            rootFILTERdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
        # Now rootFILTERdir should be defined.
        if rootFILTERdir == None:
            raise Exception("Ack: rootFILTERdir is undefined and it seems that 'throughputs' is not setup.")
        # Import Sed and Bandpass for calculating magnitudes.
        #  (assumes catalogs_measures is setup)
        import lsst.sims.photUtils.Sed as Sed
        import lsst.sims.photUtils.Bandpass as Bandpass
        # read in and set up filter files
        bandpass = {}
        # Add V and imsim bandpasses to filterlist
        filterlist.insert(0, 'V')
        filterlist.append('imsim')
        for f in filterlist:
            if f == 'V':
                filename = os.path.join(rootSEDdir,'harris_V.dat')
                # Instantiate and read the throughput file.
                bandpass[f] = Bandpass()
                bandpass[f].readThroughput(filename)
            elif f=='imsim':
                filename = 'imsim'
                bandpass[f] = Bandpass()
                bandpass[f].imsimBandpass()
            else:
                filename = os.path.join(rootFILTERdir,'total_' + f + '.dat')
                # Instantiate and read the throughput file.
                bandpass[f] = Bandpass()
                bandpass[f].readThroughput(filename)
        # Read in and set up sed files. Assumes asteroid SEDS end in .dat
        possiblesedtypes = os.listdir(rootSEDdir)
        sedtypes = []
        for p in possiblesedtypes:
            if p.endswith('.dat'):
                sedtypes.append(p)
        sed={}
        sedmag = {}
        sedcol = {}
        for sedfile in sedtypes:
            # Read sed files.
            filename = os.path.join(rootSEDdir, sedfile)
            sed[sedfile] = Sed()
            sed[sedfile].readSED_flambda(filename)
            # Set up magnitudes for color calculation.
            sedmag[sedfile] = {}
            sedcol[sedfile] = {}
            for f in filterlist:
               sedmag[sedfile][f] = sed[sedfile].calcMag(bandpass[f])
               sedcol[sedfile][f] = sedmag[sedfile][f] - sedmag[sedfile]['V']
        # Save the info needed later.
        self.sedcol = sedcol
        return

    def calcMags(self, mjd, filt, m5=None, withErrors=True):
        try:
            self.sedcol
        except AttributeError:
            raise Exception("Error: Run setup_calcMags first, to generate color information for every SED type.")
        # Do a few checks on self.orbits['sedname'] .. does it exist (just set to S.dat if not).
        try:
            self.orbits['sedname']
        except AttributeError:
            warning.warn("Warning: Did not find sednames for this movingObjects, so setting to 'S.dat'")
            self.orbits['sedname'] = []
            for i in range(len(self.orbits['objid'])):
                self.orbits['sedname'].append('S.dat')
            self.orbits['sedname'] = numpy.array(self.orbits['sedname'], (string, 6))
        # If it does exist, check that self.orbit['sedname'] only contains SEDS matching our sed files.
        used_sedlist = numpy.unique(self.orbits['sedname'])
        for s in used_sedlist:
            if s not in self.sedcol.keys():
                raise Exception("The sedname, %s, is not present in the SEDs I have colors for (%s)" %(s, self.sedcol.keys()))
        # Set up mjdstr for access to ephemeris dictionaries.
        mjdstr = self.make_mjdStr(mjd)
        # Check that ephems exist at this time.
        try:
            self.ephems[mjdstr]
        except AttributeError:
            warning.warn('Warning: computing ephemeris for this date.')
            self.generate_ephemerides(mjd)
        # Calculate the actual magnitudes.
        self.ephems[mjdstr]['filtmag'] = numpy.empty([len(self.orbits['objid']),], float)
        self.ephems[mjdstr]['imsimmag'] = numpy.empty([len(self.orbits['objid']),], float)
        for u in used_sedlist:
            condition = (self.orbits['sedname'] == u)
            self.ephems[mjdstr]['filtmag'][condition] = self.sedcol[u][filt] + \
                                                        self.ephems[mjdstr]['magV'][condition]
            self.ephems[mjdstr]['imsimmag'][condition] = self.sedcol[u]['imsim'] + \
                                                         self.ephems[mjdstr]['magV'][condition] 
        # Add SNR measurement if given 5-sigma limiting mag for exposure.
        if m5 != None:
            flux_ratio = numpy.power(10, 0.5*(m5 - self.ephems[mjdstr]['filtmag']))
            snr = 5 * (flux_ratio)
            self.ephems[mjdstr]['snr'] = snr
            # Calculate approx errors in ra/dec/mag from magnitude/m5.
            if withErrors:
                # calculate error in ra/dec
                rgamma = 0.039
                # average seeing is 0.7" (or 700 mas)
                error_rand = numpy.sqrt((0.04-rgamma)*flux_ratio +
                                        rgamma*flux_ratio*flux_rtio)
                ast_error_rand = 700.0 * error_rand
                ast_error_sys = 10.0
                astrom_error = numpy.sqrt(ast_error_sys**2 + ast_error_rand**2)
                # convert from mas to deg
                astrom_error = astrom_error / 100.0 / 60.0/ 60.0
                self.ephems[mjdstr]['astrom_err'] = astrom_err
                mag_error_sys = 0.005
                mag_error = numpy.sqrt(error_rand**2 + mag_error_sys**2)
                self.ephems[mjdstr]['mag_err'] = mag_err
        return


    ## These are some possibly interesting functions,
    ##   more like utility functions for alternative calculations.
    ## If not desired, just comment out this whole section.

    def convert_radec2ec(self, mjd):
        # Calculate ecliptic lat/lon for ephemerides on mjd
        #   (which can be a list or array or single value).
        mjdlistStr = []
        if ((isinstance(mjdlist, float) == True) | (isinstance(mjdlist, int) ==True)):
            mjdlist = [mjdlist]
        for mjd in mjdlist:
            mjdlistStr.append(self.make_mjdStr(mjd))
        # Use RA and DEC on [mjd] to calculate eclat and eclon,
        #   adding into the ephems[mjd] dictionary. 
        ecnode = 0.0
        ecinc = 23.439291*_deg2rad
        for mj in mjdlistStr:
            x = numpy.cos(self.ephems[mjd]['ra']*_deg2rad) * \
                numpy.cos(self.ephems[mjd]['dec']*_deg2rad)
            y = numpy.sin(self.ephems[mjd]['ra']*_deg2rad) * \
                numpy.cos(self.ephems[mjd]['dec']*_deg2rad)
            z = numpy.sin(self.ephems[mjd]['dec']*_deg2rad)
            xp = x
            yp = numpy.cos(ecinc)*y + numpy.sin(ecinc)*z
            zp = -numpy.sin(ecinc)*y + numpy.cos(ecinc)*z
            self.ephems[mjd]['eclat'] = numpy.arcsin(zp)*_rad2deg
            self.ephems[mjd]['eclon'] = numpy.arctan2(yp, xp)*_rad2deg
            self.ephems[mjd]['eclon'] = self.ephems[mjd]['eclon'] % 360
        return

    def calc_phasefunction(self, mjdlist, gval=0.15):
        # using the phase angle (in DEGREES), compute how the brightness of
        #  the object varies
        # albedo*phase = how much sunlight is reflected
        mjdlistStr = []
        if ((isinstance(mjdlist, float) == True) | (isinstance(mjdlist, int) ==True)):
            mjdlist = [mjdlist]
        for mjd in mjdlist:
            mjdlistStr.append(self.make_mjdStr(mjd))
        for mj in mjdlistStr:
            phione = 10**(-3.33 * numpy.power(numpy.tan(self.ephems[mjd]['phase_angle']*
                                                        _deg2rad/2), 0.63))
            phitwo = 10**(-1.87 * numpy.power(numpy.tan(self.ephems[mjd]['phase_angle']*
                                                        _deg2rad/2), 1.22))
            phaseval = (1-gval)*phione + gval*phitwo
        return phaseval

    # The next functions calculate some *approximate* information about the
    # location of the objects  at a particular time, from their orbital elements only.
    def calc_distanceE(self):
        # Calculate heliocentric distance (at epoch),
        #  for all objects using only Keplerian orbit information.
        import scipy
        # calculate eccentric anomaly from mean anomaly
        #  M = E - e*sin(E)
        def func(E, e, M):
            return (E - e*numpy.sin(E) - M)
        E = numpy.zeros(len(self.orbits['e']), dtype='float')
        for i in range(len(self.orbits['e'])):
            E[i] = scipy.optimize.newton(func, x0=M[i], args=(self.orbits['e'][i],
                                                              self.orbits['M'][i]*_deg2rad),
                                         tol=1e-10)
        # calculate distance from eccentric anomaly
        r = self.orbits['a'] * (1 -self.orbits['e'] * numpy.cos(E))
        E = E * _rad2deg
        # r = distance from the Sun
        return r, E

"""
    def calc_sundist(self, earth_dist, solarelongation):
        # given the solar elongation (single value, in RADIANS) and the distance from Earth to the object,
        # calculate and return the distance to the sun for the same object
        # uses the law of cosines
        sun_dist = numpy.sqrt(earth_dist**2 + 1 - 2*earth_dist*numpy.cos(solarelongation))
        return sun_dist

    def calc_earthdist(self, sun_dist, solarelongation):
        # given the solar elongation (single value, radians) and the distance from the Sun to the object,
        # return the distance between Earth and the object
        phaseangle = calc_phaseangle(solarelongation, sun_dist)
        sangle = numpy.pi - solarelongation - phaseangle
        earth_dist = numpy.sqrt(1 + sun_dist**2 - 2*sun_dist*numpy.cos(sangle))
        earth_dist = numpy.where(numpy.isnan(earth_dist), 0, earth_dist)
        return earth_dist

    def calc_phaseangle(self, solarelongation, sun_dist):
        # given the solar elongation (in RADIANS) and the distance from Sun to object
        # calculate and return the phase angle for the object (could be array)
        # using law of sines sin(angle)/opposite side = sin(angle2) / opposite side 2
        phaseangle = numpy.arcsin(numpy.sin(solarelongation) / sun_dist)
        return phaseangle

    def calc_Hmag(self, obsmag, earth_dist,  sun_dist, solarelongation, gval=0.15):
        # H magnitude is 'absolute' magnitude, = mag object would have at 1 AU and 0 phase angle
        # H = m - 2.5 log(D**2 d**2 / phasefunc)
        phaseangle = calc_phaseangle(solarelongation, sun_dist)
        phasefunc = calc_phasefunction(phaseangle, gval)
        Hval = obsmag - 2.5 * numpy.log10(sun_dist**2 * earth_dist**2 / phasefunc)
        return Hval

    def calc_diam(self, Hval, albedo):
        # calculate diameter of an object given the H magnitude and albedo
        # H = -2.5 log(Sunflux * albedo * radius**2 / radius(sphere)**2)
        # so, diam = sqrt(2/pi*10^[(H-Msun+2.5log(albedo))/2.5])
        tmpval = mag_sun - Hval - 2.5*numpy.log10(albedo)
        diam = numpy.sqrt(10**(tmpval/2.5)) * 2.0
        # convert units of diam from AU to km
        diam = diam * km_per_au
        return diam
"""
