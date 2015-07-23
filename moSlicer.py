import os
import numpy as np
import pandas as pd

from .moo import MoOrbits

class MoSlicer(MoOrbits):

    def __init__(self):
        self.orbits = None
        self.ssoObs = None

    def setupSlicer(self, orbitfile, obsfile):
        """
        Read orbits, then iterate through the orbits
        and read the observation data from the obsfile.
        """
        # Read the orbits.
        self.readOrbits(orbitfile)
