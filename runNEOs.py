import os
import numpy as np
import matplotlib.pyplot as plt

import moObs as moObs

import lsst.sims.maf.plots as plots
from moSlicer import MoSlicer
import moMetrics as moMetrics
from moSummaryMetrics import ValueAtHMetric
import moPlots as moPlots
import moMetricBundle as mmb



def makeObs(orbitfile, obsfile, opsimdb, rFov, useCamera):
    """
    Generates and writes observations of objects in orbitfile into obsfile, using observations from opsimdb.
    This is assuming that the population is NEO-ish.
    """
    # rFov = fov in radians
    # useCamera = True/False (use camera footprint)
    moObs.runMoObs(orbitfile, obsfile, opsimdb, rFov=rFov, useCamera=useCamera, tstep=2./24.)

def calcCompleteness(orbitfile, obsfile, outDir):
    """
    Calculate a basic set of metrics on the NEOs, including completeness.
    Saves the plots to disk.
    """
    # Read data back from disk.
    mos = MoSlicer(orbitfile, Hrange=np.arange(13, 26, 0.5))
    mos.readObs(obsfile)

    # Nobs
    metric = MoMetrics.NObsMetric()
    slicer = mos
    pandasConstraint = None
    plotDict = {'nxbins':100, 'nybins':100}
    nobs = mmb.MoMetricBundle(metric, slicer, pandasConstraint,
                          runName=runName, metadata=metadata, plotDict=plotDict)

    # Calculate completeness. First we must calculate "DiscoveryChances".
    # Set up an example metric bundle.
    metric = MoMetrics.DiscoveryChancesMetric()
    slicer = mos
    pandasConstraint = None
    discovery = mmb.MoMetricBundle(metric, slicer, pandasConstraint,
                                    runName=runName, metadata=metadata, plotDict=plotDict)


    bdict = {'nobs':nobs, 'discovery':discovery}
    bg = mmb.MoMetricBundleGroup(bdict, outDir=outDir)
    bg.runAll()
    bg.plotAll()

    ph = plots.PlotHandler(outDir=outDir)
    # Then calculate 'completeness' as function of H, as a secondary metric.
    completeness = discovery.reduceMetric(discovery.metric.reduceFuncs['Completeness'])

    # And we can make an 'integrated over H distribution' version.
    completenessInt = completeness.reduceMetric(completeness.metric.reduceFuncs['CumulativeH'])

    Hmark = 21.0
    for c in [completeness, completenessInt]:
        summaryMetric = ValueAtHMetric(Hmark=Hmark)
        c.setSummaryMetrics(summaryMetric)
        c.computeSummaryStats()
        label = "Completeness at H=%.1f: %.2f" %(Hmark, c.summaryValues['Value At H=%.1f' %Hmark])
        c.setPlotDict({'label':label})
        c.plot(plotFunc = moPlots.MetricVsH())
        plt.axvline(Hmark, color='r', linestyle=':')
        plt.axhline(c.summaryValues['Value At H=%.1f' %(Hmark)], color='r', linestyle='-')
        plt.legend(loc=(0.9, 0.2))
        plt.savefig(os.path.join(outDir, c.fileRoot + '_vsH' + '.png'), format='png')


if __name__ == '__main__':

    outDir = 'out'
    orbitfile = 'pha20141031.des'
    obsfile = os.path.join(outDir, 'pha_obs.txt')
    opsimdb = 'enigma_1189_sqlite.db'
    rFov = np.radians(1.75)
    useCamera = True

    makeObs(orbitfile, obsfile, opsimdb, rFov, useCamera)
    calcCompleteness(orbitfile, obsfile, outDir)
    print "All Done!"
