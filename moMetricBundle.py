import os
import numpy as np
import numpy.ma as ma


class MoMetricBundle(object):
    def __init__(self, metric, slicer, pandasConstraint='',
                 runName='opsim', metadata='',
                 fileRoot='metricOutput',
                 plotDict=None, plotFuncs=None,
                 secondaryMetrics=None):

        self.metric = metric
        self.slicer = slicer
        self.pandasConstraint = pandasConstraint
        self.metricValues = None
        self.fileRoot = fileRoot
        self.runName = runName
        self.metadata = metadata
        if plotDict is not None:
            self.plotDict = plotDict
        else:
            self.plotDict = {}
        if self.plotFuncs is not None:
            self.plotFuncs = plotFuncs
        else:
            self.plotFuncs = [pFunc() for pFunc in self.slicer.plotFuncs]
        if secondaryMetrics is not None:
            self.secondaryMetrics = secondaryMetrics

    def _setupMetricValues(self):
        """
        Set up the numpy masked array to store the metric value data.
        """
        dtype = self.metric.metricDtype
        # Can't store some mask values in an int array.
        if dtype == 'int':
            dtype = 'float'
        self.metricValues = ma.MaskedArray(data = np.empty(self.slicer.slicerShape, dtype),
                                            mask = np.zeros(self.slicer.slicerShape, 'bool'),
                                            fill_value= self.slicer.badval)

    def runSecondaryMetric(self):
        """
        Calculate secondary-type metrics, such as completeness or integrate over H distribution.
        These return a new metric bundle.
        """
        
        newBundle = MoMetricBundle(Metric
        
    def write(self):
        # This doesn't really do the full job yet.
        self.slicer.write(self.fileRoot, self)


####

class MoMetricBundleGroup(object):
    def __init__(self, bundleDict, outDir='.'):
        self.bundleDict = bundleDict
        self.outDir = outDir
        if not os.path.isdir(self.outDir):
            os.makedirs(self.outDir)
        self.slicer = self.bundleDict.itervalues().next().slicer
        for b in self.bundleDict.itervalues():
            if b.slicer != self.slicer:
                raise ValueError('Currently, the slicers for the MoMetricBundleGroup must be equal - using the same observations and Hvals.')
        self.pandasConstraints = list(set([b.pandasConstraint for b in bundleDict.values()]))

    def _setCurrent(self, pandasConstraint):
        """
        Private utility to set the currentBundleDict (i.e. a set of metricBundles with the same constraint).
        """
        self.currentBundleDict = {}
        for k, b in self.bundleDict.iteritems():
            if b.pandasConstraint == pandasConstraint:
                self.currentBundleDict[k] = b

    def runCurrent(self, pandasConstraint):
        """
        Calculate the metric values for a set of compatible slicers.
        """
        for b in self.currentBundleDict.itervalues():
            b._setupMetricValues()
        for i, slicePoint in enumerate(self.slicer):
            ssoObs = slicePoint['obs']
            for j, Hval in enumerate(slicePoint['Hvals']):
                for b in self.currentBundleDict.itervalues():
                    if len(ssoObs) == 0:
                        b.metricValues.mask[i][j] = True
                    else:
                        b.metricValues.data[i][j] = b.metric.run(ssoObs, slicePoint['orbit'], Hval)

    def runAll(self):
        """
        Run all pandasConstraints and metrics for these moMetricBundles.
        """
        for pandasConstraint in self.pandasConstraints:
            self._setCurrent(pandasConstraint)
            self.runCurrent(pandasConstraint)
        print 'Calculated all metrics.'

    def plotAll(self):
        """
        Generate 
