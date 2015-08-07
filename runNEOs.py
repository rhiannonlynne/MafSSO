import os, argparse
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

def calcCompleteness(orbitfile, obsfile, outDir, runName, metadata=None):
    """
    Calculate a basic set of metrics on the NEOs, including completeness.
    Saves the plots to disk.
    """
    # Read data back from disk.
    mos = MoSlicer(orbitfile, Hrange=np.arange(15, 27, 0.25))
    mos.readObs(obsfile)

    # Nobs
    metric = moMetrics.NObsMetric()
    slicer = mos
    pandasConstraint = None
    plotDict = {'nxbins':100, 'nybins':100}
    nobs = mmb.MoMetricBundle(metric, slicer, pandasConstraint,
                          runName=runName, metadata=metadata, plotDict=plotDict)

    # Calculate completeness. First we must calculate "DiscoveryChances".
    # Set up an example metric bundle.
    discovery = {}
    nObsList = [2, 3, 4]
    for nObs in nObsList:
        md = metadata + ' %d visits/night' %(nObs)
        metric = moMetrics.DiscoveryChancesMetric(nObsPerNight=nObs, Hindex=0.3)
        slicer = mos
        pandasConstraint = None
        discovery[nObs] = mmb.MoMetricBundle(metric, slicer, pandasConstraint,
                                        runName=runName, metadata=md, plotDict=plotDict)

    bdict = {'discovery_%s' %(nObs):discovery[nObs] for nObs in nObsList }
    bdict['nobs'] = nobs

    bg = mmb.MoMetricBundleGroup(bdict, outDir=outDir)
    bg.runAll()
    bg.plotAll()

    ph = plots.PlotHandler(outDir=outDir)
    # Then calculate 'completeness' as function of H, as a secondary metric.
    completeness = {}
    completenessInt = {}
    for nObs in nObsList:
        completeness[nObs] = discovery[nObs].reduceMetric(discovery[nObs].metric.reduceFuncs['Completeness'])
        # And we can make an 'integrated over H distribution' version.
        completenessInt[nObs] = completeness[nObs].reduceMetric(completeness[nObs].metric.reduceFuncs['CumulativeH'])

    Hmark = 22.0
    for nObs in nObsList:
        for c in [completeness[nObs], completenessInt[nObs]]:
            summaryMetric = ValueAtHMetric(Hmark=Hmark)
            c.setSummaryMetrics(summaryMetric)
            c.computeSummaryStats()
            label = "Completeness at H=%.1f: %.2f" %(Hmark, c.summaryValues['Value At H=%.1f' %Hmark])
            c.setPlotDict({'label':label})
            c.plot(plotFunc = moPlots.MetricVsH())
            plt.axvline(Hmark, color='r', linestyle=':')
            plt.axhline(c.summaryValues['Value At H=%.1f' %(Hmark)], color='r', linestyle='-')
            plt.legend(loc=(0.1, 0.2))
            plt.savefig(os.path.join(outDir, c.fileRoot + '_vsH' + '.png'), format='png')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Python script to run MAF NEO completeness analysis')
    parser.add_argument('dbFile', type=str, default=None,help="full file path to the opsim sqlite file")
    parser.add_argument('--orbitfile', type=str, default='pha20141031.des', help="full file path to the orbit file")
    parser.add_argument("--outDir",type=str, default=None, help='Output directory for MAF outputs. Default opsimdb_orbitroot')
    parser.add_argument('--plotOnly', dest='plotOnly', action='store_true',
                        default=False, help="Reload the metric values from disk and re-plot them.")
    parser.set_defaults()
    args, extras = parser.parse_known_args()

    orbitfile = args.orbitfile
    orbitroot = os.path.split(orbitfile)[-1].replace('.txt', '').replace('.des', '').replace('.dat', '')
    obsfile = orbitroot + '_allObs.txt'
    opsimdb = args.dbFile
    runName = os.path.split(opsimdb)[-1].replace('_sqlite.db', '')
    metadata = '%s' %(orbitroot)
    rFov = np.radians(1.75)
    useCamera = True

    outDir = args.outDir
    if outDir is None:
        outDir = runName + '_' + orbitroot
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    if not args.plotOnly:
        makeObs(orbitfile, os.path.join(outDir, obsfile), opsimdb, rFov, useCamera)
        print 'Made observations'
    calcCompleteness(orbitfile, os.path.join(outDir, obsfile), outDir, runName, metadata)
    print "All Done!"
