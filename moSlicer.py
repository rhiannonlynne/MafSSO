import os
import numpy as np
import pandas as pd

from .moObs import MoOrbits

class MoSlicer(MoOrbits):

    def __init__(self):
        self.orbits = None
        self.ssoIds = None
        self.ssoObs = None

    # def readOrbits(self, orbitfile)

    def readObs(self, obsfile):
        """
        Read observations created by moObs.
        """
        # For now, just read all the observations (should be able to chunk this though).
        self.ssoObs = pd.read_table(obsfile, sep='\s*', engine='python')
        # We may have to rename the first column from '#objId' to 'objId'.
        if self.ssoObs.columns.values[0].startswith('#'):
            newcols = self.ssoObs.columns.values
            newcols[0] = newcols[0].replace('#', '')
            self.ssoObs.columns = newcols
        # self.ssoObs = self.ssoObs.to_records()

    def _getObs(self, ssoId_idx):
        """
        Return the observations of ssoId.
        For now this works for any ssoId; in the future, this might only work as ssoId is
         progressively iterated through the series of ssoIds (so we can 'chunk' the reading).
        """
        ssoId = self.ssoIds[ssoId_idx]
        obs = self.ssoObs.query('objId' == ssoId)
        return obs.to_records()

    def __iter__(self):
        """
        Iterate through each of the ssoIds.
        """
        self.idx = 0
        return

    def next(self):
        """
        Returns result of self._getObs when iterating over moSlicer.
        """
        if self.idx >= self.nSso:
            raise StopIteration
        idx = self.idx
        self.idx += 1
        return self._getObs(idx)
