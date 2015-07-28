import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection

from lsst.sims.maf.plots import BasePlotter


class MetricVsH(BasePlotter):
    """
    Plot metric values versus H.
    Marginalize over metric values in each H bin using 'npReduce'.
    """
    def __init__(self):
        self.plotType = 'MetricVsH'
        self.objectPlotter = False
        self.defaultPlotDict = {'title':None, 'xlabel':'H (mag)', 'ylabel':None, 'label':None,
                                'color':'b', 'linestyle':'-', 'npReduce':np.mean, 'nbins':30}

    def __call__(self, metricValue, slicer, userPlotDict, fignum=None):
        fig = plt.figure(fignum)
        plotDict = {}
        plotDict.update(self.defaultPlotDict)
        plotDict.update(userPlotDict)
        Hvals = slicer.slicePoints['H']
        if len(Hvals) == slicer.slicerShape[1]:
            # Using cloned H distribution.
            # Apply 'npReduce' method directly to metric values, and plot at matching H values.
            mVals = plotDict['npReduce'](metricValue, axis=0)
        else:
            # Probably each object has its own H value.
            stepsize = (Hvals.max() - Hvals.min())  / float(plotDict['nbins'])
            bins = np.arange(Hvals.min(), Hvals.max() + stepsize/2.0, stepsize)
            # In each bin of H, calculate the 'npReduce' value of the corresponding metricValues.
            inds = np.digitize(Hvals, bins)
            inds = inds-1
            mVals = np.zeros(len(bins), float)
            for i in range(len(bins)):
                mVals[i] = plotDict['npReduce'](metricValue[inds == i].filled())
            Hvals = bins
        plt.plot(Hvals, mVals, color=plotDict['color'], linestyle=plotDict['linestyle'],
                label=plotDict['label'])
        plt.xlabel(plotDict['xlabel'])
        plt.ylabel(plotDict['ylabel'])
        if 'xMin' in plotDict:
            plt.xlim(xmin = plotDict['xMin'])
        if 'xMax' in plotDict:
            plt.xlim(xmax = plotDict['xMax'])
        return fig.number


class MetricVsOrbit(BasePlotter):
    """
    Plot metric values (at a particular H value) vs. orbital parameters.
    Marginalize over metric values in each orbital bin using 'npReduce'.
    """
    def __init__(self):
        self.plotType = 'MetricVsOrbit'
        self.objectPlotter = False
        self.defaultPlotDict = {'title':None, 'xlabel':None, 'ylabel':None,
                                'label':None, 'cmap':cm.cubehelix,
                                'xaxis':'q', 'yaxis':'e', 'npReduce':np.mean,
                                'nxbins':100, 'nybins':100,
                                'Hval':None, 'Hwidth':1.0}

    def __call__(self, metricValue, slicer, userPlotDict, fignum=None):
        fig = plt.figure(fignum)
        plotDict = {}
        plotDict.update(self.defaultPlotDict)
        plotDict.update(userPlotDict)
        xvals = slicer.slicePoints['orbits'][plotDict['xaxis']]
        yvals = slicer.slicePoints['orbits'][plotDict['yaxis']]
        if plotDict['xlabel'] is None:
            plotDict['xlabel'] = plotDict['xaxis']
        if plotDict['ylabel'] is None:
            plotDict['ylabel'] = plotDict['yaxis']
        # Set x/y bins.
        if 'xbins' in plotDict:
            xbins = plotDict['xbins']
        else:
            xbinsize = (xvals.max() - xvals.min())/float(plotDict['nxbins'])
            xbins = np.arange(xvals.min(), xvals.max() + xbinsize/2.0, xbinsize)
        if 'ybins' in plotDict:
            ybins = plotDict['ybins']
        else:
            ybinsize = (yvals.max() - yvals.min())/float(plotDict['nybins'])
            ybins = np.arange(yvals.min(), yvals.max() + ybinsize/2.0, ybinsize)
        nxbins = len(xbins)
        nybins = len(ybins)
        # Identify the relevant metricValues for the Hvalue we want to plot.
        Hvals = slicer.slicePoints['H']
        if plotDict['Hval'] is None:
            if len(Hvals) == slicer.slicerShape[1]:
                Hidx = len(Hvals) / 2
                Hval = Hvals[Hidx]
            else:
                Hval = np.median(Hvals)
                Hidx = np.where(np.abs(Hvals - Hval) <= plotDict['Hwidth']/2.0)[0]
        if len(Hvals) == slicer.slicerShape[1]:
            mVals = np.swapaxes(metricValue, 1, 0)[Hidx].filled()
        else:
            mVals = metricValue[Hidx].filled()
        # Calculate the mean metric values at each x/y bin.
        mSum = np.zeros((nybins, nxbins), dtype='float')
        mNums = np.zeros((nybins, nxbins), dtype='int')
        xidxs = np.digitize(xvals, xbins) - 1
        yidxs = np.digitize(yvals, ybins) - 1
        for iy in range(nybins):
            ymatch = np.where(yidxs == iy)[0]
            for ix in range(nxbins):
                xmatch = np.where(xidxs[ymatch] == ix)[0]
                val = plotDict['npReduce'](mVals[ymatch][xmatch])
        xi, yi = np.meshgrid(xbins, ybins)
        if 'colorMin' in plotDict:
            vMin = plotDict['colorMin']
        else:
            vMin = val.min()
        if 'colorMax' in plotDict:
            vMax = plotDict['colorMax']
        else:
            vMax = val.max()
        levels = np.arange(vMin, vMax, (vMax-vMin)/200.0)
        plt.contourf(xi, yi, val, levels, extend='max', zorder=0)
        cbar = plt.colorbar()
        plt.xlabel(plotDict['xlabel'])
        plt.ylabel(plotDict['ylabel'])
        return fig.number

class MetricVsOrbitPoints(BasePlotter):
    """
    Plot metric values (at a particular H value) as function of orbital parameters,
    using points for each metric value.
    """
    def __init__(self):
        self.plotType = 'MetricVsOrbit'
        self.objectPlotter = False
        self.defaultPlotDict = {'title':None, 'xlabel':None, 'ylabel':None,
                                'label':None, 'cmap':cm.cubehelix,
                                'xaxis':'q', 'yaxis':'e',
                                'Hval':None, 'Hwidth':1.0,
                                'foregroundPoints':True, 'backgroundPoints':False}

    def __call__(self, metricValue, slicer, userPlotDict, fignum=None):
        fig = plt.figure(fignum)
        plotDict = {}
        plotDict.update(self.defaultPlotDict)
        plotDict.update(userPlotDict)
        xvals = slicer.slicePoints['orbits'][plotDict['xaxis']]
        yvals = slicer.slicePoints['orbits'][plotDict['yaxis']]
        if plotDict['xlabel'] is None:
            plotDict['xlabel'] = plotDict['xaxis']
        if plotDict['ylabel'] is None:
            plotDict['ylabel'] = plotDict['yaxis']
        # Identify the relevant metricValues for the Hvalue we want to plot.
        Hvals = slicer.slicePoints['H']
        if plotDict['Hval'] is None:
            if len(Hvals) == slicer.slicerShape[1]:
                Hidx = len(Hvals) / 2
                Hval = Hvals[Hidx]
            else:
                Hval = np.median(Hvals)
                Hidx = np.where(np.abs(Hvals - Hval) <= plotDict['Hwidth']/2.0)[0]
        if len(Hvals) == slicer.slicerShape[1]:
            mVals = np.swapaxes(metricValue, 1, 0)[Hidx]
        else:
            mVals = metricValue[Hidx]
        if 'colorMin' in plotDict:
            vMin = plotDict['colorMin']
        else:
            vMin = mVals.min()
        if 'colorMax' in plotDict:
            vMax = plotDict['colorMax']
        else:
            vMax = mVals.max()
        if plotDict['backgroundPoints']:
            # This isn't quite right for the condition .. but will do for now.
            condition = np.where(mVals == 0)
            plt.plot(xvals[condition], yvals[condition], 'r.', markersize=4, alpha=0.5, zorder=3)
        if plotDict['foregroundPoints']:
            plt.scatter(xvals, yvals, c=mvals, vmin=vMin, vmax=vMax, s=15, alpha=0.8, zorder=0)
            cb = plt.colorbar()
        plt.xlabel(plotDict['xlabel'])
        plt.ylabel(plotDict['ylabel'])
        return fig.number
