import os
import numpy as np
import matplotlib.pyplt as plt

class MoMetricBundle(object):
    def __init__(self, metric, slicer, pandasConstraint=''):
        self.metric = metric
        self.slicer = slicer
        self.pandasConstraint = pandasConstraint
        self.metricValues = None

    def write(self, filename):
        self.slicer.write(filename, self)


####

class MoMetricBundleGroup(object):
    def __init__(self, bundleDict, outDir='.'):
        self.bundleDict = bundleDict
        self.outDir = outDir
        if not os.path.isdir(self.outDir):
            os.makedirs(self.outDir)
        self.pandasConstraints = list(set([b.pandasConstraint for b in bundleDict.values()]))

    def _setCurrent(self, pandasConstraint):
        """
        Private utility to set the currentBundleDict (i.e. a set of metricBundles with the same constraint).
        """
        self.currentBundleDict = {}
        for k, b in self.bundleDict.iteritems():
            if b.pandasConstraint == pandasConstraint:
                self.currentBundleDict[k] = b

    def _checkCompatible(self, b1, b2):
        """
        Check if two bundles are compatible.
        This basically means they have the same slicer and the same pandasconstraint.
        """
        if b1.slicer == b2.slicer:
            if b1.pandasConstraint == b2.pandasConstraint:
                return True
        else:
            return False

    def runAll(self):

