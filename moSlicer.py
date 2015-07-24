import os
import numpy as np
import pandas as pd

from moObs import MoOrbits

__all__ = ['MoSlicer']

class MoSlicer(MoOrbits):

    def __init__(self, orbitfile, Hrange=None):
        """
        Instantiate the MoSlicer object.

        orbitfile = the file with the orbit information on the objects.

        If Hrange is not None (instead, set to a numpy array), then
        each orbit will be cloned to the H values specified by Hrange.

        Iteration over the MoSlicer will go as:
          - iterate over each orbit;
            - if Hrange is not None, for each orbit, iterate over Hrange.
        """
        # Read orbits (inherited from MoOrbits).
        self.readOrbits(orbitfile)
        # See if we're cloning orbits, and set slicer shape (metric shape) accordingly.
        self.Hrange = Hrange
        if self.Hrange is not None:
            self.slicerShape = [self.nSso, len(Hrange)]
        else:
            self.slicerShape = [self.nSso, 1]
        # Set observations to None.
        self.ssoObs = None


    def readObs(self, obsfile):
        """
        Read observations created by moObs.
        """
        # For now, just read all the observations (should be able to chunk this though).
        self.allObs = pd.read_table(obsfile, delim_whitespace=True)
        # We may have to rename the first column from '#objId' to 'objId'.
        if self.allObs.columns.values[0].startswith('#'):
            newcols = self.allObs.columns.values
            newcols[0] = newcols[0].replace('#', '')
            self.allObs.columns = newcols
        if 'magFilter' not in self.allObs.columns.values:
            self.allObs['magFilter'] = self.allObs['magV'] + self.allObs['dmagColor']
        self.subsetObs()

    def subsetObs(self, pandasConstraint=None):
        """
        Choose a subset of all the observations, such as those in a particular time period.
        """
        if pandasConstraint is None:
            self.obs = self.allObs
        else:
            self.obs = self.allObs.query(pandasConstraint)

    def _sliceObs(self, idx):
        """
        Return the observations of ssoId.
        For now this works for any ssoId; in the future, this might only work as ssoId is
         progressively iterated through the series of ssoIds (so we can 'chunk' the reading).
        """
        # Find the matching orbit.
        orb = self.orbits.iloc[idx]
        # Find the matching observations.
        obs = self.obs.query('objId == %d' %(orb['objId']))
        # Return the values for H to consider for metric.
        if self.Hrange is not None:
            Hvals = self.Hrange
        else:
            Hvals = orb['H']
        return {'obs': obs.to_records(),
                'orb': orb,
                'Hvals': Hvals}

    def __iter__(self):
        """
        Iterate through each of the ssoIds.
        """
        self.idx = 0
        return self

    def next(self):
        """
        Returns result of self._getObs when iterating over moSlicer.
        """
        if self.idx >= 10: #self.nSso:
            raise StopIteration
        idx = self.idx
        self.idx += 1
        return self._sliceObs(idx)
