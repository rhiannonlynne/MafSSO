import os, argparse
import numpy as np
import matplotlib.pyplot as plt

import moObs as moObs

import lsst.sims.maf.plots as plots
from moSlicer import MoSlicer
import moMetrics as moMetrics
import moSummaryMetrics as moSummaryMetrics
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
    plotDict = {'nxbins':200, 'nybins':200}
    nobs = mmb.MoMetricBundle(metric, slicer, pandasConstraint,
                          runName=runName, metadata=metadata, plotDict=plotDict)

    # Calculate completeness. First we must calculate "DiscoveryChances".
    # Set up an example metric bundle.
    discovery = {}
    nObsList = [2, 3, 4]
    for nObs in nObsList:
        md = metadata + ' %d visits/night' %(nObs)
        metric = moMetrics.DiscoveryMetric(nObsPerNight=nObs)
        slicer = mos
        pandasConstraint = None
        discovery[nObs] = mmb.MoMetricBundle(metric, slicer, pandasConstraint,
                                             runName=runName, metadata=md,
                                             plotDict=plotDict, plotFuncs=[moPlots.MetricVsH()])

    discovery['4nights'] = mmb.MoMetricBundle(metric=moMetrics.DiscoveryMetric(nObsPerNight=2, nNightsPerWindow=4),
                                              slicer=mos, constraint=None, runName=runName,
                                              metadata = metadata+'4 nights/track',
                                              plotDict=plotDict, plotFuncs=[moPlots.MetricVsH()])

    bdict = {'discovery_%s' %(nObs):discovery[nObs] for nObs in nObsList }
    bdict['4nights'] = discovery['4nights']
    bdict['nobs'] = nobs

    bg = mmb.MoMetricBundleGroup(bdict, outDir=outDir)
    bg.runAll()
    bg.plotAll()

    Hmark = 22

    completeness = {}
    completenessInt = {}
    hVals = {}
    fig1 = plt.figure()
    plt.title('Completeness %s' %(runName))
    fig2 = plt.figure()
    plt.title('Cumulative completeness %s' %(runName))
    for nObs in nObsList:
        discChances = discovery[nObs].childBundles['N_Chances']
        discChances.setSummaryMetrics([moSummaryMetrics.CompletenessMetric(), moSummaryMetrics.CumulativeCompletenessMetric()])
        discChances.computeSummaryStats()
        completeness[nObs] = discChances.summaryValues['Completeness'][0]
        hVals[nObs] = discChances.summaryValues['Completeness'][1]
        completenessInt[nObs] = discChances.summaryValues['CumulativeCompleteness'][0]
        print discChances.summaryValues['Completeness']
        print nObs
        print hVals[nObs]
        print completeness[nObs]
        print completenessInt[nObs]
        Hidx = np.where(hVals[nObs] == Hmark)[0]
        plt.figure(fig1.number)
        cval = completeness[nObs][Hidx]
        plt.plot(hVals[nObs], completeness[nObs], label='%s per tracklet, %.2f @H=%.1f' %(nObs, cval, Hmark))
        plt.figure(fig2.number)
        cval = completenessInt[nObs][Hidx]
        plt.plot(hVals[nObs], completenessInt[nObs], label='%s per tracklet, %.2f @H=%.1f' %(nObs, cval, Hmark))
    for nObs in ['4nights']:
        discChances = discovery[nObs].childBundles['N_Chances']
        discChances.setSummaryMetrics([moSummaryMetrics.CompletenessMetric(), moSummaryMetrics.CumulativeCompletenessMetric()])
        discChances.computeSummaryStats()
        completeness[nObs] = discChances.summaryValues['Completeness'][0]
        hVals[nObs] = discChances.summaryValues['Completeness'][1]
        completenessInt[nObs] = discChances.summaryValues['CumulativeCompleteness'][0]
        Hidx = np.where(hVals[nObs] == Hmark)[0]
        plt.figure(fig1.number)
        cval = completeness[nObs][Hidx]
        plt.plot(hVals[nObs], completeness[nObs], label='%s pairs per track, %.2f @H=%.1f' %(nObs, cval, Hmark))
        plt.figure(fig2.number)
        cval = completenessInt[nObs][Hidx]
        plt.plot(hVals[nObs], completenessInt[nObs], label='%s pairs per track, %.2f @H=%.1f' %(nObs, cval, Hmark))
    plt.figure(fig1.number)
    plt.axvline(Hmark, color='r', linestyle=':')
    plt.legend(loc='upper right', fancybox=True, numpoints=1, fontsize='smallest')
    plt.savefig(os.path.join(outDir, 'completeness.pdf'), format='pdf')
    plt.figure(fig2.number)
    plt.axvline(Hmark, color='r', linestyle=':')
    plt.legend(loc='upper right', fancybox=True, numpoints=1, fontsize='smallest')
    plt.savefig(os.path.join(outDir, 'completenessInt.pdf'), format='pdf')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Python script to run MAF NEO completeness analysis')
    parser.add_argument('--orbitFile', type=str, default='pha20141031.des', help="full file path to the orbit file")
    parser.add_argument("--outDir",type=str, default=None, help='Output directory for MAF outputs. Default opsimdb_orbitroot')
    parser.add_argument('--plotOnly', dest='plotOnly', action='store_true',
                        default=False, help="Reload the metric values from disk and re-plot them.")
    parser.add_argument('dbFile', type=str, help="full file path to the opsim sqlite file")
    parser.set_defaults()
    args, extras = parser.parse_known_args()

    orbitfile = args.orbitFile
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
