import os, argparse
import numpy as np
import matplotlib.pyplot as plt

from moObs import MoObs

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
    #moObs.runMoObs(orbitfile, obsfile, opsimdb, rFov=rFov, useCamera=useCamera, tstep=2./24.)
    # Read orbits.
    moogen = MoObs()
    moogen.readOrbits(orbitfile)
    print "Read orbit information from %s" %(orbitfile)

    # Check rfov/camera choices.
    if useCamera:
        print "Using camera footprint"
    else:
        print "Not using camera footprint; using circular fov with %f degrees radius" %(np.degrees(rFov))

    # Read data for visits from numpy file on disk (generated from CombineObs)
    import pandas as pd
    simdata = pd.read_csv(opsimdb)
    simdata = simdata.to_records()
    print "Read data from opsim %s, fetched %d visits." %(opsimdb, len(simdata['expMJD']))

    tstep=2./24.
    moogen.setTimesRange(timeStep=tstep, timeStart=simdata['expMJD'].min(), timeEnd=simdata['expMJD'].max())
    print "Will generate ephemerides on grid of %f day timesteps, then extrapolate to opsim times." %(tstep)

    moogen.setupOorb()
    for i, sso in moogen.orbits.iterrows():
        ephs = moogen.generateEphs(sso)
        interpfuncs = moogen.interpolateEphs(ephs)
        idxObs = moogen.ssoInFov(interpfuncs, simdata, rFov=rFov, useCamera=useCamera)
        moogen.writeObs(sso['objId'], interpfuncs, simdata, idxObs,  sedname=sso['sed_filename'],  outfileName=obsfile)
    print "Wrote output observations to file %s" %(obsfile)


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
    nyears = [3, 5, 10, 12, 15]
    for yr in nyears:
        discovery[yr] = {}
        for nObs in nObsList:
            md = metadata + ' %d visits/night after %d years' %(nObs, yr)
            #metric = moMetrics.DiscoveryMetric(tMin=0, tMax=90./60/24., nObsPerNight=nObs, nNightsPerWindow=3, tWindow=30)
            metric = moMetrics.DiscoveryChancesMetric(tNight=90./60./24., nObsPerNight=nObs, nNightsPerWindow=3, tWindow=15)
            slicer = mos
            pandasConstraint = 'night <= %d' %(yr*365)
            discovery[yr][nObs] = mmb.MoMetricBundle(metric, slicer, pandasConstraint,
                                                    runName=runName, metadata=md,
                                                    plotDict=plotDict, plotFuncs=[moPlots.MetricVsH()])

        #metric = moMetrics.DiscoveryMetric(tMin=0, tMax=90./60/24., nObsPerNight=2, nNightsPerWindow=4, tWindow=30)
        metric = moMetrics.DiscoveryChancesMetric(tNight=90./60./24., nObsPerNight=2, nNightsPerWindow=4, tWindow=30)
        discovery[yr]['4nights'] = mmb.MoMetricBundle(metric=metric,
                                                slicer=mos, constraint=pandasConstraint, runName=runName,
                                                metadata = metadata+'4 nights/track after %d years' %(yr),
                                                plotDict=plotDict, plotFuncs=[moPlots.MetricVsH()])

    bdict = {}
    for yr in nyears:
        for nObs in nObsList:
            bdict['discovery_%s_%d' %(nObs, yr)] = discovery[yr][nObs]
        bdict['4nights_%d' %(yr)] = discovery[yr]['4nights']
    bdict['nobs'] = nobs

    bg = mmb.MoMetricBundleGroup(bdict, outDir=outDir)
    bg.runAll()
    bg.plotAll()

    Hmark = 22

    completeness = {}
    completenessInt = {}
    hVals = mos.Hrange
    for yr in nyears:
        completeness[yr] = {}
        completenessInt[yr] = {}
        for nObs in nObsList:
        #discChances = discovery[nObs].childBundles['N_Chances']
            discChances = discovery[yr][nObs]
            discChances.setSummaryMetrics([moSummaryMetrics.CompletenessMetric(), moSummaryMetrics.CumulativeCompletenessMetric()])
            discChances.computeSummaryStats()
            completeness[yr][nObs] = discChances.summaryValues['Completeness'][0]
            completenessInt[yr][nObs] = discChances.summaryValues['CumulativeCompleteness'][0]
        for nObs in ['4nights']:
        #discChances = discovery[nObs].childBundles['N_Chances']
            discChances = discovery[yr][nObs]
            discChances.setSummaryMetrics([moSummaryMetrics.CompletenessMetric(), moSummaryMetrics.CumulativeCompletenessMetric()])
            discChances.computeSummaryStats()
            completeness[yr][nObs] = discChances.summaryValues['Completeness'][0]
            completenessInt[yr][nObs] = discChances.summaryValues['CumulativeCompleteness'][0]

    # Make a figure of completeness with the different discovery patterns, for one year.
    yr = 12
    Hidx = np.where(hVals == Hmark)[0]
    plt.figure()
    plt.title('Completeness %s' %(runName))
    for nObs in nObsList:
        cval = completeness[yr][nObs][Hidx]
        plt.plot(hVals, completeness[yr][nObs], label='%d obs/tracklet, %.0f%s H @ %.1f' %(nObs, cval*100, '%', Hmark))
    for nObs in ['4nights']:
        cval = completeness[yr][nObs][Hidx]
        plt.plot(hVals, completeness[yr][nObs], label='4 pairs/track, %.0f%s H @ %.1f' %(cval*100, '%', Hmark))
    plt.axvline(Hmark, color='r', linestyle=':')
    plt.xlabel('H')
    plt.ylabel('Completeness @ H')
    plt.legend(loc='upper right', fancybox=True, numpoints=1, fontsize='smaller')
    plt.savefig(os.path.join(outDir, 'completeness.pdf'), format='pdf')

    plt.figure()
    plt.title('Cumulative completeness %s' %(runName))
    for nObs in nObsList:
        cval = completenessInt[yr][nObs][Hidx]
        plt.plot(hVals, completenessInt[yr][nObs], label='%d obs/tracklet, %.0f%s H <= %.1f' %(nObs, cval*100, '%', Hmark))
    for nObs in ['4nights']:
        cval = completenessInt[yr][nObs][Hidx]
        plt.plot(hVals, completenessInt[yr][nObs], label='4 pairs/track, %.0f%s H <= %.1f' %(cval*100, '%', Hmark))
    plt.axvline(Hmark, color='r', linestyle=':')
    plt.xlabel('H')
    plt.ylabel('Completeness <= H')
    plt.legend(loc='upper right', fancybox=True, numpoints=1, fontsize='smaller')
    plt.savefig(os.path.join(outDir, 'completenessInt.pdf'), format='pdf')

    # And make a figure for one set of discovery criterion, over time.
    nObs = 2
    plt.figure()
    plt.title('Completeness over time %s' %(runName))
    for yr in nyears:
        cval = completeness[yr][nObs][Hidx]
        plt.plot(hVals, completeness[yr][nObs], label='Year %d: %.0f @ H=%.1f' %(yr, cval*100, Hmark))
    plt.axvline(Hmark, color='r', linestyle=':')
    plt.xlabel('H')
    plt.ylabel('Completeness @ H')
    plt.legend(loc='upper right', fancybox=True, numpoints=1, fontsize='smaller')
    plt.savefig(os.path.join(outDir, 'completenessTime.pdf'), format='pdf')

    nObs = 2
    plt.figure()
    plt.title('Completeness over time %s' %(runName))
    for yr in nyears:
        cval = completenessInt[yr][nObs][Hidx]
        plt.plot(hVals, completenessInt[yr][nObs], label='Year %d: %.0f <= H=%.1f' %(yr, cval*100, Hmark))
    plt.axvline(Hmark, color='r', linestyle=':')
    plt.xlabel('H')
    plt.ylabel('Completeness <= H')
    plt.legend(loc='upper right', fancybox=True, numpoints=1, fontsize='smaller')
    plt.savefig(os.path.join(outDir, 'completenessIntTime.pdf'), format='pdf')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Python script to run MAF NEO completeness analysis')
    parser.add_argument('--orbitFile', type=str, default='pha20141031.des', help="full file path to the orbit file")
    parser.add_argument('--obsFile', type=str, default=None, help='full file path to the observation file')
    parser.add_argument("--outDir",type=str, default=None, help='Output directory for MAF outputs. Default opsimdb_orbitroot')
    parser.add_argument('--plotOnly', dest='plotOnly', action='store_true',
                        default=False, help="Reload the metric values from disk and re-plot them.")
    parser.add_argument('dbFile', type=str, help="full file path to the opsim sqlite file")
    parser.set_defaults()
    args, extras = parser.parse_known_args()

    orbitfile = args.orbitFile
    orbitroot = os.path.split(orbitfile)[-1].replace('.txt', '').replace('.des', '').replace('.dat', '')
    if args.obsFile is None:
        obsfile = orbitroot + '_allObs.txt'
    else:
        obsfile = args.obsFile
    opsimdb = args.dbFile
    runName = os.path.split(opsimdb)[-1].replace('_summary.csv', '')
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
