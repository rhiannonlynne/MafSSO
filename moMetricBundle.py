import os
import numpy as np
import numpy.ma as ma

import moPlots as moPlots
import lsst.sims.maf.utils as utils
import lsst.sims.maf.plots as plots

class MoMetricBundle(object):
    def __init__(self, metric, slicer, constraint=None,
                 runName='opsim', metadata=None,
                 fileRoot=None,
                 plotDict=None, plotFuncs=None,
                 secondaryMetrics=None):
        """
        Instantiate moving object metric bundle, save metric/slicer/constraint, etc.
        """
        self.metric = metric
        self.slicer = slicer
        if constraint == '':
            constraint = None
        self.constraint = constraint
        # For compatibility with plotHandler/etc until we reconcile this better.
        self.sqlconstraint = constraint
        if self.sqlconstraint is None:
            self.sqlconstraint = ''
        self.runName = runName
        self._buildMetadata(metadata)
        # Set output file root name.
        self._buildFileRoot(fileRoot)
        self.plotDict = {'units':'@H'}
        self.setPlotDict(plotDict)
        self.setPlotFuncs(plotFuncs)
        if secondaryMetrics is not None:
            self.secondaryMetrics = secondaryMetrics
        # Set up metric value storage.
        self.metricValues = None

    def _buildFileRoot(self, fileRoot=None):
        """
        Build an auto-generated output filename root (i.e. minus the plot type or .npz ending).
        """
        if fileRoot is None:
            # Build basic version.
            self.fileRoot = '_'.join([self.runName, self.metric.name, self.metadata])
        else:
            self.fileRoot = fileRoot
        # Sanitize output name if needed.
        self.fileRoot = utils.nameSanitize(self.fileRoot)

    def _buildMetadata(self, metadata):
        """
        Combine any provided metadata and constraint.
        """
        # Use obsfile name for metadata if none provided.
        if metadata is not None:
            self.metadata = metadata
        else:
            self.metadata = self.slicer.obsfile.replace('.txt', '').replace('_allObs', '').replace('.dat', '')
        # And modify by constraint.
        if self.constraint is not None:
            self.metadata += ' ' + self.constraint

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

    def setPlotDict(self, plotDict):
        """
        Set or update any property of plotDict.
        """
        # Don't auto-generate anything here - the plotHandler does it.
        if plotDict is not None:
            self.plotDict.update(plotDict)

    def setPlotFuncs(self, plotFuncs=None):
        """
        Set or reset the plotting functions.
        Default is to use all the plotFuncs associated with a slicer.
        """
        if plotFuncs is not None:
            if plotFuncs is isinstance(plotFuncs, plots.BasePlotter):
                self.plotFuncs = [plotFuncs]
            else:
                self.plotFuncs = []
                for pFunc in plotFuncs:
                    if not isinstance(pFunc, plots.BasePlotter):
                        raise ValueError('plotFuncs should contain instantiated lsst.sims.maf.plotter objects.')
                    self.plotFuncs.append(pFunc)
        else:
            # Moving object slicers keep instantiated plotters in the self.slicer.plotFuncs.
            self.plotFuncs = [pFunc for pFunc in self.slicer.plotFuncs]

    def runSecondaryMetric(self, secondaryMetric):
        """
        Calculate secondary-type metrics, such as completeness or integrate over H distribution.
        These return a new metric bundle.
        """
        newBundle = MoMetricBundle(metric=secondaryMetric, slicer=self.slicer,
                                   constraint=self.constraint,
                                   runName=self.runName, metadata=self.metadata,
                                   plotDict=self.plotDict, plotFuncs=self.plotFuncs)
        if secondaryMetric.name == 'Completeness':
            newBundle.plotFuncs = [moPlots.MetricVsH()]
            newBundle.setPlotDict({'ylabel':'Completeness @H'})
        if secondaryMetric.name == 'IntegrateOverH':
            newBundle.metric.name = self.metric.name + ' (cumulative H)'
            newBundle.setPlotDict({'units':'<=H'})
            if 'ylabel' in newBundle.plotDict:
                newBundle.plotDict['ylabel'] = newBundle.plotDict['ylabel'].replace('@H', '<=H')
        # Calculate new metric values.
        newBundle.metricValues, Hrange = secondaryMetric.run(self.metricValues, self.slicer.slicePoints['H'])
        if len(Hrange) != len(self.slicer.slicePoints['H']):
            newBundle.slicer.slicePoints['H'] = Hrange
            if secondaryMetric.name == 'Completeness':
                newBundle.slicer.slicerShape = [len(Hrange), len(Hrange)]
        return newBundle

    def plot(self, plotHandler=None, plotFunc=None, outfileSuffix=None, savefig=False):
        """
        Create all plots available from the slicer. plotHandler holds the output directory info, etc.
        """
        # Generate a plotHandler if none was set.
        if plotHandler is None:
            plotHandler = plots.PlotHandler(savefig=savefig)
        # Make plots.
        if plotFunc is not None:
            if isinstance(plotFunc, plots.BasePlotter):
                plotFunc = plotFunc
            else:
                plotFunc = plotFunc()

        plotHandler.setMetricBundles([self])
        plotHandler.setPlotDicts(plotDicts=[self.plotDict], reset=True)
        madePlots = {}
        if plotFunc is not None:
            # We haven't updated plotHandler to know about these kinds of plots yet.
            # and we want to automatically set some values for the ylabel for metricVsH.
            tmpDict = {}
            if plotFunc.plotType == 'MetricVsH':
                if 'ylabel' not in self.plotDict:
                    tmpDict['ylabel'] = self.metric.name
            fignum = plotHandler.plot(plotFunc, plotDicts=tmpDict, outfileSuffix=outfileSuffix)
            madePlots[plotFunc.plotType] = fignum
        else:
            for plotFunc in self.plotFuncs:
                # We haven't updated plotHandler to know about these kinds of plots yet.
                # and we want to automatically set some values for the ylabel for metricVsH.
                tmpDict = {}
                if plotFunc.plotType == 'MetricVsH':
                    if 'ylabel' not in self.plotDict:
                        tmpDict['ylabel'] = self.metric.name
                fignum = plotHandler.plot(plotFunc, plotDicts=tmpDict, outfileSuffix=outfileSuffix)
                madePlots[plotFunc.plotType] = fignum
        return madePlots

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
        self.constraints = list(set([b.constraint for b in bundleDict.values()]))

    def _setCurrent(self, constraint):
        """
        Private utility to set the currentBundleDict (i.e. a set of metricBundles with the same constraint).
        """
        self.currentBundleDict = {}
        for k, b in self.bundleDict.iteritems():
            if b.constraint == constraint:
                self.currentBundleDict[k] = b

    def runCurrent(self, constraint):
        """
        Calculate the metric values for set of bundles using the same constraint and slicer.
        """
        self._setCurrent(constraint)
        self.slicer.subsetObs(constraint)
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
        Run all constraints and metrics for these moMetricBundles.
        """
        for constraint in self.constraints:
            self.runCurrent(constraint)
        print 'Calculated all metrics.'

    def plotCurrent(self, constraint):
        self._setCurrent(constraint)
        for b in self.currentBundleDict.itervalues():
            b.plot()

    def plotAll(self):
        """
        Make a few generically desired plots. This needs more flexibility in the future.
        """
        for constraint in self.constraints:
            self.plotCurrent(constraint)
        print 'Plotted all metrics.'

