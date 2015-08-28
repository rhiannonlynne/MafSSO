import os
from copy import deepcopy
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

from moSlicer import MoSlicer
from moMetrics import BaseMoMetric
import moPlots as moPlots
import moStackers as moStackers
import lsst.sims.maf.utils as utils
from lsst.sims.maf.plots import PlotHandler, BasePlotter


class MoMetricBundle(object):
    def __init__(self, metric, slicer, constraint=None,
                 runName='opsim', metadata=None,
                 fileRoot=None,
                 plotDict=None, plotFuncs=None,
                 summaryMetrics=None):
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
        self.setSummaryMetrics(summaryMetrics)
        # Set up metric value storage.
        self.metricValues = None
        self.summaryValues = None

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

    def setSummaryMetrics(self, summaryMetrics):
        """
        Set (or reset) the summary metrics for the metricbundle.
        """
        if summaryMetrics is not None:
            if isinstance(summaryMetrics, BaseMoMetric):
                self.summaryMetrics = [summaryMetrics]
            else:
                self.summaryMetrics = []
                for s in summaryMetrics:
                    if not isinstance(s, BaseMoMetric):
                        raise ValueError('SummaryStats must only contain instantiated moving object metric objects')
                    self.summaryMetrics.append(s)
        else:
            self.summaryMetrics = []

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
            if plotFuncs is isinstance(plotFuncs, BasePlotter):
                self.plotFuncs = [plotFuncs]
            else:
                self.plotFuncs = []
                for pFunc in plotFuncs:
                    if not isinstance(pFunc, BasePlotter):
                        raise ValueError('plotFuncs should contain instantiated lsst.sims.maf.plotter objects.')
                    self.plotFuncs.append(pFunc)
        else:
            # Moving object slicers keep instantiated plotters in the self.slicer.plotFuncs.
            self.plotFuncs = [pFunc for pFunc in self.slicer.plotFuncs]

    def reduceMetric(self, reduceFunc, reducePlotDict=None, reduceDisplayDict=None):
        """
        Run reduce methods on the metric bundle, such as completeness or integrate over H distribution.
        These return a new metric bundle.
        """
        rName = reduceFunc.__name__.replace('reduce', '')
        reduceName = self.metric.name + '_' + rName
        newmetric = deepcopy(self.metric)
        newmetric.name = reduceName
        newBundle = MoMetricBundle(metric=newmetric, slicer=self.slicer,
                                   constraint=self.constraint,
                                   runName=self.runName, metadata=self.metadata,
                                   plotDict=self.plotDict, plotFuncs=self.plotFuncs,
                                   summaryMetrics=self.summaryMetrics)
        if rName == 'Completeness':
            newBundle.metric.name = 'Completeness'
            newBundle.plotFuncs = [moPlots.MetricVsH()]
            newBundle.setPlotDict({'ylabel':'Completeness @H'})
        if rName == 'CumulativeH':
            newBundle.metric.name = self.metric.name + ' (cumulative H)'
            newBundle.setPlotDict({'units':'<=H'})
            if 'ylabel' in newBundle.plotDict:
                newBundle.plotDict['ylabel'] = newBundle.plotDict['ylabel'].replace('@H', '<=H')
        newBundle._buildFileRoot()
        # Calculate new metric values.
        newBundle.metricValues, Hvals = reduceFunc(self.metricValues, self.slicer.slicePoints['H'])
        if newBundle.metricValues.shape[1] != self.slicer.slicerShape[1]:
            # Then something happened with reduceFunction --
            #  .. usually, this would be with 'completeness'. It can reshape the metric values so that
            #  (a) if we cloned H, we now go from multiple metricvalues per H value to a single value per H [(nSso, nHrange) -> (1, nHrange)]
            #  (b) if we did not clone H, we now have a binned metricvalues - one per H bin, instead of one per nsso [(nSso, 1) -> (nHrange, nHrange)]
            # and we don't really care (in general) because we won't be re-slicing the observations. But we do need to update the slicePoints['H'],
            #  and for that we need a new slicer. Rereading the orbit file is a bit of a pain, but not too bad.
            newslicer = MoSlicer(orbitfile=self.slicer.orbitfile, Hrange=Hvals)
            newBundle.slicer = newslicer
        return newBundle

    def computeSummaryStats(self, resultsDb=None):
        """
        Compute summary statistics on metricValues, using summaryMetrics (metricbundle list).
        So far, the only summary metric that is possible to use would be 'ValueAtHMetric'.
        """
        if self.summaryValues is None:
            self.summaryValues = {}
        if self.summaryMetrics is not None:
            # Build array of metric values, to use for (most) summary statistics.
            for m in self.summaryMetrics:
                summaryName = m.name
                summaryVal = m.run(self.metricValues, self.slicer.slicePoints['H'])
                self.summaryValues[summaryName] = summaryVal
                # Add summary metric info to results database, if applicable.
                if resultsDb:
                    metricId = resultsDb.updateMetric(self.metric.name, self.slicer.slicerName,
                                                      self.runName, self.constraint, self.metadata, None)
                    resultsDb.updateSummaryStat(metricId, summaryName=summaryName, summaryValue=summaryVal)

    def plot(self, plotHandler=None, plotFunc=None, outfileSuffix=None, savefig=False):
        """
        Create all plots available from the slicer. plotHandler holds the output directory info, etc.
        """
        # Generate a plotHandler if none was set.
        if plotHandler is None:
            plotHandler = PlotHandler(savefig=savefig)
        # Make plots.
        if plotFunc is not None:
            if isinstance(plotFunc, BasePlotter):
                plotFunc = plotFunc
            else:
                plotFunc = plotFunc()

        plotHandler.setMetricBundles([self])
        # The plotDict will be automatically accessed when the plotHandler calls the plotting method.
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
    def __init__(self, bundleDict, outDir='.', resultsDb=None, verbose=True):
        # Not really handling resultsDb yet.
        self.verbose = verbose
        self.bundleDict = bundleDict
        self.outDir = outDir
        if not os.path.isdir(self.outDir):
            os.makedirs(self.outDir)
        self.resultsDb = resultsDb

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
        allStackers = moStackers.AllStackers()
        self.slicer.subsetObs(constraint)
        for b in self.currentBundleDict.itervalues():
            b._setupMetricValues()
        for i, slicePoint in enumerate(self.slicer):
            ssoObs = slicePoint['obs']
            for j, Hval in enumerate(slicePoint['Hvals']):
                # Run stackers to add extra columns (that depend on H)
                ssoObs = allStackers.run(ssoObs, slicePoint['orbit']['H'], Hval)
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
            self._setCurrent(constraint)
            self.runCurrent(constraint)
        if self.verbose:
            print 'Calculated all metrics.'

    def plotCurrent(self, savefig=True, outfileSuffix=None, figformat='pdf', dpi=600, thumbnail=True,
                closefigs=True):
        plotHandler = PlotHandler(outDir=self.outDir, resultsDb=self.resultsDb,
                                  savefig=savefig, figformat=figformat, dpi=dpi, thumbnail=thumbnail)
        for b in self.currentBundleDict.itervalues():
            b.plot(plotHandler=plotHandler, outfileSuffix=outfileSuffix, savefig=savefig)
            if closefigs:
                plt.close('all')
        if self.verbose:
            print 'Plotting complete.'

    def plotAll(self, savefig=True, outfileSuffix=None, figformat='pdf', dpi=600, thumbnail=True,
                closefigs=True):
        """
        Make a few generically desired plots. This needs more flexibility in the future.
        """
        for constraint in self.constraints:
            self._setCurrent(constraint)
            self.plotCurrent(savefig=savefig, outfileSuffix=outfileSuffix, figformat=figformat, dpi=dpi,
                             thumbnail=thumbnail, closefigs=closefigs)
        if self.verbose:
            print 'Plotted all metrics.'

    def reduceCurrent(self, updateSummaries=True):
        """
        Run all reduce functions for the metricbundle in the currentBundleDict.
        """
        # Create a temporary dictionary to hold the reduced metricbundles.
        reduceBundleDict = {}
        for b in self.currentBundleDict.itervalues():
            # If there are no reduce functions associated with the metric, skip this metricBundle.
            if len(b.metric.reduceFuncs) > 0:
                # Apply reduce functions, creating a new metricBundle in the process (new metric values).
                for reduceFunc in b.metric.reduceFuncs.itervalues():
                    newmetricbundle = b.reduceMetric(reduceFunc)
                    # Add the new metricBundle to our metricBundleGroup dictionary.
                    name = newmetricbundle.metric.name
                    if name in self.bundleDict:
                        name = newmetricbundle.fileRoot
                    reduceBundleDict[name] = newmetricbundle
                # Remove summaryMetrics from top level metricbundle if desired.
                if updateSummaries:
                    b.summaryMetrics = []
        # Add the new metricBundles to the MetricBundleGroup dictionary.
        self.bundleDict.update(reduceBundleDict)
        # And add to to the currentBundleDict too, so we run as part of 'summaryCurrent'.
        self.currentBundleDict.update(reduceBundleDict)

    def reduceAll(self, updateSummaries=True):
        """
        Run the reduce methods for all metrics in bundleDict.
        This assumes that 'clearMemory' was false.
        """
        for constraint in self.constraints:
            self._setCurrent(constraint)
            self.reduceCurrent(updateSummaries=updateSummaries)

    def summaryCurrent(self):
        """
        Run summary statistics on all the metricBundles in currentBundleDict.
        """
        for b in self.currentBundleDict.itervalues():
            b.computeSummaryStats(self.resultsDb)

    def summaryAll(self):
        """
        Run the summary statistics for all metrics in bundleDict.
        This assumes that 'clearMemory' was false.
        """
        for constraint in self.constraints:
            self._setCurrent(constraint)
            self.summaryCurrent()
