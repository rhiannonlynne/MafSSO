import inspect
import numpy as np
import numpy.ma as ma
from lsst.sims.maf.metrics import BaseMetric

__all__ = ['BaseMoMetric', 'NObsMetric', 'NObsNoSinglesMetric',
           'NNightsMetric', 'ObsArcMetric',
           'DiscoveryMetric', 'Discovery_N_ChancesMetric', 'Discovery_N_ObsMetric',
           'Discovery_TimeMetric', 'Discovery_RADecMetric', 'Discovery_EcLonLatMetric',
           'ActivityOverTimeMetric', 'ActivityOverPeriodMetric',
           'DiscoveryChancesMetric', 'MagicDiscoveryMetric',
           'HighVelocityMetric', 'HighVelocityNightsMetric']


class BaseMoMetric(BaseMetric):
    """Base class for the moving object metrics."""

    def __init__(self, cols=None, metricName=None, units='#', badval=0,
                 comment=None, childMetrics=None,
                 appMagCol='appMag', m5Col='magLimit',
                 nightCol='night', expMJDCol='expMJD',
                 snrCol='SNR',  visCol='vis',
                 raCol='ra', decCol='dec', seeingCol='finSeeing',
                 expTimeCol='visitExpTime'):
        # Set metric name.
        self.name = metricName
        if self.name is None:
            self.name = self.__class__.__name__.replace('Metric', '', 1)
        # Set badval and units, leave space for 'comment' (tied to displayDict).
        self.badval = badval
        self.units = units
        self.comment = comment
        # Set some commonly used column names.
        self.m5Col = m5Col
        self.appMagCol = appMagCol
        self.nightCol = nightCol
        self.expMJDCol = expMJDCol
        self.snrCol = snrCol
        self.visCol = visCol
        self.raCol = raCol
        self.decCol = decCol
        self.seeingCol = seeingCol
        self.expTimeCol = expTimeCol
        self.colsReq = [self.appMagCol, self.m5Col,
                        self.nightCol, self.expMJDCol,
                        self.snrCol, self.visCol]
        if cols is not None:
            for col in cols:
                self.colsReq.append(col)

        if childMetrics is None:
            try:
                if not isinstance(self.childMetrics, dict):
                    raise ValueError('self.childMetrics must be a dictionary (possibly empty)')
            except AttributeError:
                self.childMetrics = {}
                self.metricDtype = 'float'
        else:
            if not isinstance(childMetrics, dict):
                raise ValueError('childmetrics must be provided as a dictionary.')
            self.childMetrics = childMetrics
            self.metricDtype = 'object'

    def run(self, ssoObs, orb, Hval):
        raise NotImplementedError


class NObsMetric(BaseMoMetric):
    """
    Count the total number of observations where an object was 'visible'.
    """
    def __init__(self, snrLimit=None, **kwargs):
        """
        @ snrLimit .. if snrLimit is None, this uses the _calcVis method/completeness
                      if snrLimit is not None, this uses that value as a cutoff instead.
        """
        super(NObsMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit

    def run(self, ssoObs, orb, Hval):
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
            return vis.size
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
            return vis.size


class NObsNoSinglesMetric(BaseMoMetric):
    """
    Count the number of observations for an object, but don't
    include any observations where it was a single observation on a night.
    """
    def __init__(self, snrLimit=None, **kwargs):
        super(NObsNoSinglesMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit

    def run(self, ssoObs, orb, Hval):
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return 0
        nights = ssoObs[self.nightCol][vis]
        nights = nights.astype('int')
        ncounts = np.bincount(nights)
        nobs = ncounts[np.where(ncounts > 1)].sum()
        return nobs


class NNightsMetric(BaseMoMetric):
    """
    Count the number of distinct nights an object is observed.
    """
    def __init__(self, snrLimit=None, **kwargs):
        """
        @ snrLimit : if SNRlimit is None, this uses _calcVis method/completeness
                     else if snrLimit is not None, it uses that value as a cutoff.
        """
        super(NNightsMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit

    def run(self, ssoObs, orb, Hval):
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return 0
        nights = len(np.unique(ssoObs[self.nightCol][vis]))
        return nights

class ObsArcMetric(BaseMoMetric):
    """
    Calculate the difference between the first and last observation of an object.
    """
    def __init__(self, snrLimit=None, **kwargs):
        super(ObsArcMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit

    def run(self, ssoObs, orb, Hval):
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return 0
        arc = ssoObs[self.expMJDCol][vis].max() - ssoObs[self.expMJDCol][vis].min()
        return arc

class DiscoveryMetric(BaseMoMetric):
    """Identify the discovery opportunities for an object."""
    def __init__(self, nObsPerNight=2,
                 tMin=5./60.0/24.0, tMax=90./60./24.0,
                 nNightsPerWindow=3, tWindow=15,
                 snrLimit=None, **kwargs):
        """
        @ nObsPerNight = number of observations per night required for tracklet
        @ tMin = min time start/finish for the tracklet (days)
        @ tMax = max time start/finish for the tracklet (days)
        @ nNightsPerWindow = number of nights with observations required for track
        @ tWindow = max number of nights in track (days)
        @ snrLimit .. if snrLimit is None then uses 'completeness' calculation in 'vis' column.
                   .. if snrLimit is not None, then uses this SNR value as a cutoff.
        """
        self.childMetrics = {'N_Chances':Discovery_N_ChancesMetric(self),
                             'N_Obs':Discovery_N_ObsMetric(self),
                             'Time':Discovery_TimeMetric(self),
                             'RADec':Discovery_RADecMetric(self),
                             'EcLonLat':Discovery_EcLonLatMetric(self)}
        super(DiscoveryMetric, self).__init__(childMetrics = self.childMetrics, **kwargs)
        self.snrLimit = snrLimit
        self.nObsPerNight = nObsPerNight
        self.tMin = tMin
        self.tMax = tMax
        self.nNightsPerWindow = nNightsPerWindow
        self.tWindow = tWindow
        self.badval = None

    def run(self, ssoObs, orb, Hval):
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return self.badval
        # Identify discovery opportunities.
        #  Identify visits where the 'night' changes.
        visSort = np.argsort(ssoObs[self.expMJDCol][vis])
        nights = ssoObs[self.nightCol][vis][visSort]
        #print 'all nights', nights
        n = np.unique(nights)
        # Identify all the indexes where the night changes in value.
        nIdx = np.searchsorted(nights, n)
        #print 'nightchanges', nights[nIdx]
        # Count the number of observations per night (except last night)
        obsPerNight = (nIdx - np.roll(nIdx, 1))[1:]
        # Add the number of observations on the last night.
        obsLastNight = np.array([len(nights) - nIdx[-1]])
        obsPerNight = np.concatenate((obsPerNight, obsLastNight))
        # Find the nights with more than nObsPerNight.
        nWithXObs = n[np.where(obsPerNight >= self.nObsPerNight)]
        nIdxMany = np.searchsorted(nights, nWithXObs)
        nIdxManyEnd = np.searchsorted(nights, nWithXObs, side='right') - 1
        # Check that nObsPerNight observations are within tMin/tMax
        timesStart = ssoObs[self.expMJDCol][visSort][nIdxMany]
        timesEnd = ssoObs[self.expMJDCol][visSort][nIdxManyEnd]
        # Identify the nights with 'clearly good' observations.
        good = np.where((timesEnd - timesStart >= self.tMin) & (timesEnd - timesStart <= self.tMax), 1, 0)
        # Identify the nights where we need more investigation (a subset of the visits may be within the interval).
        check = np.where((good==0) & (nIdxManyEnd + 1 - nIdxMany > self.nObsPerNight) & (timesEnd-timesStart > self.tMax))[0]
        for i, j, c in zip(visSort[nIdxMany][check], visSort[nIdxManyEnd][check], check):
            t = ssoObs[self.expMJDCol][i:j+1]
            dtimes = (np.roll(t, 1- self.nObsPerNight) - t)[:-1]
            tidx = np.where((dtimes >= self.tMin) & (dtimes <= self.tMax))[0]
            if len(tidx) > 0:
                good[c] = 1
        # 'good' provides mask for observations which could count as 'good to make tracklets' against ssoObs[visSort][nIdxMany]
        # Now identify tracklets which can make tracks.
        goodIdx = visSort[nIdxMany][good == 1]
        goodIdxEnds = visSort[nIdxManyEnd][good == 1]
        #print 'good tracklets', nights[goodIdx]
        if len(goodIdx) < self.nNightsPerWindow:
            return self.badval
        deltaNights = np.roll(ssoObs[self.nightCol][vis][goodIdx], 1 - self.nNightsPerWindow) - ssoObs[self.nightCol][vis][goodIdx]
        # Identify the index in ssoObs[vis][goodIdx] (sorted by expMJD) where the discovery opportunity starts.
        startIdxs = np.where((deltaNights >= 0) & (deltaNights <= self.tWindow))[0]
        # Identify the index where the discovery opportunity ends.
        endIdxs = np.zeros(len(startIdxs), dtype='int')
        for i, sIdx in enumerate(startIdxs):
            inWindow = np.where(ssoObs[self.nightCol][vis][goodIdx] - ssoObs[self.nightCol][vis][goodIdx][sIdx] <= self.tWindow)[0]
            endIdxs[i] = np.array([inWindow.max()])
        # Convert back to index based on ssoObs[vis] (sorted by expMJD).
        startIdxs = goodIdx[startIdxs]
        endIdxs = goodIdxEnds[endIdxs]
        #print 'start', startIdxs,  nights[startIdxs]#, orb['objId'], Hval
        #print 'end', endIdxs, nights[endIdxs]#, orb['objId'], Hval
        return {'start':startIdxs, 'end':endIdxs, 'trackletNights':ssoObs[self.nightCol][goodIdx]}


class Discovery_N_ChancesMetric(BaseMoMetric):
    """
    Child metric to be used with DiscoveryMetric.
    Calculates total number of discovery opportunities.
    """
    def __init__(self, parentDiscoveryMetric, **kwargs):
        super(Discovery_N_ChancesMetric, self).__init__(**kwargs)
        self.parentMetric = parentDiscoveryMetric

    def run(self, ssoObs, orb, Hval, metricValues):
        """
        Return the number of different discovery chances we had for each object/H combination.
        """
        nchances = metricValues['start'].size
        return nchances

class Discovery_N_ObsMetric(BaseMoMetric):
    """
    Calculates the number of observations in the i-th discovery track.
    """
    def __init__(self, parentDiscoveryMetric, i=0, **kwargs):
        super(Discovery_N_ObsMetric, self).__init__(**kwargs)
        self.parentMetric = parentDiscoveryMetric
        # The number of the discovery chance to use.
        self.i = i

    def run(self, ssoObs, orb, Hval, metricValues):
        """
        Return the number of observations in the i-th discovery opportunity.
        """
        if self.i >= len(metricValues['start']):
            return 0
        startIdx = metricValues['start'][self.i]
        endIdx = metricValues['end'][self.i]
        nobs = endIdx - startIdx
        return nobs

class Discovery_TimeMetric(BaseMoMetric):
    """
    Returns the time of the i-th discovery opportunity.
    """
    def __init__(self, parentDiscoveryMetric, i=0, **kwargs):
        super(Discovery_TimeMetric, self).__init__(**kwargs)
        self.parentMetric = parentDiscoveryMetric
        self.i = i

    def run(self, ssoObs, orb, Hval, metricValues):
        """
        Return the time of the i-th discovery opportunity.
        """
        if self.i>=len(metricValues['start']):
            return 0
        startIdx = metricValues['start'][self.i]
        return ssoObs[self.expMJDCol][startIdx]

class Discovery_RADecMetric(BaseMoMetric):
    """
    Returns the RA/Dec of the i-th discovery opportunity.
    """
    def __init__(self, parentDiscoveryMetric, i=0, **kwargs):
        super(Discovery_RADecMetric, self).__init__(**kwargs)
        self.parentMetric = parentDiscoveryMetric
        self.metricDtype = 'object'
        self.i = i

    def run(self, ssoObs, orb, Hval, metricValues):
        """
        Return the RA/Dec of the i-th discovery opportunity.
        """
        if self.i>=len(metricValues['start']):
            return (-999, -999)
        startIdx = metricValues['start'][self.i]
        return (ssoObs[self.raCol][startIdx], ssoObs[self.decCol][startIdx])

class Discovery_EcLonLatMetric(BaseMoMetric):
    """
    Returns the ecliptic lon/lat (in degrees) of the i-th discovery opportunity.
    """
    def __init__(self, parentDiscoveryMetric, i=0, **kwargs):
        super(Discovery_EcLonLatMetric, self).__init__(**kwargs)
        self.parentMetric = parentDiscoveryMetric
        self.i = i
        self.metricDtype = 'object'

    def run(self, ssoObs, orb, Hval, metricValues):
        if self.i>=len(metricValues['start']):
            return (-999, -999)
        startIdx = metricValues['start'][self.i]
        return (ssoObs['ecLon'][startIdx], ssoObs['ecLat'][startIdx])


class ActivityOverTimeMetric(BaseMoMetric):
    """
    Count the time periods where we would have a chance to detect activity on
    a moving object.
    Splits observations into time periods set by 'window', then looks for observations within each window,
    and reports what fraction of the total windows receive 'nObs' visits.
    """
    def __init__(self, window, snrLimit=5, surveyYears=10.0, metricName=None, **kwargs):
        if metricName is None:
            metricName = 'Chance of detecting activity lasting %.0f days' %(window)
        super(ActivityOverTimeMetric, self).__init__(metricName=metricName, **kwargs)
        self.snrLimit = snrLimit
        self.window = window
        self.surveyYears = surveyYears
        self.windowBins = np.arange(0, self.surveyYears*365 + self.window/2.0, self.window)
        self.nWindows = len(self.windowBins)
        self.units = '%f Day Windows' %(self.window)

    def run(self, ssoObs, orb,  Hval):
        # For cometary activity, expect activity at the same point in its orbit at the same time, mostly
        # For collisions, expect activity at random times
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return self.badval
        n, b = np.histogram(ssoObs[vis][self.nightCol], bins=self.windowBins)
        activityWindows = np.where(n>0)[0].size
        return activityWindows / float(self.nWindows)


class ActivityOverPeriodMetric(BaseMoMetric):
    """
    Count the fraction of the orbit (when split into nBins) that receive
    observations, in order to have a chance to detect activity.
    """
    def __init__(self, binsize, snrLimit=5,
                 aCol='a', tPeriCol='tPeri', metricName=None, **kwargs):
        """
        @ binsize : size of orbit slice, in degrees.
        """
        if metricName is None:
            metricName = 'Chance of detecting activity in %.1f of the orbit' %(window)
        super(ActivityOverPeriodMetric, self).__init__(metricName=metricName, **kwargs)
        self.aCol = aCol
        self.tPeriCol = tPeriCol
        self.snrLimit = snrLimit
        self.binsize = np.radians(binsize)
        self.anomalyBins = np.arange(0, 2*np.pi + self.binsize/2.0, self.binsize)
        self.nBins = len(self.anomalyBins)
        self.units = '%.1f deg' %(np.degrees(self.binsize))

    def run(self, ssoObs, orb, Hval):
        # For cometary activity, expect activity at the same point in its orbit at the same time, mostly
        # For collisions, expect activity at random times
        period = np.power(orb[self.aCol], 3./2.) * 365.25
        anomaly = ((ssoObs[self.expMJDCol] - orb[self.tPeriCol]) / period) % (2*np.pi)
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return self.badval
        n, b = np.histogram(anomaly[vis], bins=self.anomalyBins)
        activityWindows = np.where(n>0)[0].size
        return activityWindows / float(self.nBins)


class DiscoveryChancesMetric(BaseMoMetric):
    """
    Count the number of discovery opportunities for an object.
    """
    def __init__(self, nObsPerNight=2, tNight=90./60./24.,
                 nNightsPerWindow=3, tWindow=15, snrLimit=None,
                 **kwargs):
        """
        @ nObsPerNight = number of observations per night required for tracklet
        @ tNight = max time start/finish for the tracklet (days)
        @ nNightsPerWindow = number of nights with observations required for track
        @ tWindow = max number of nights in track (days)
        @ snrLimit .. if snrLimit is None then uses 'completeness' calculation,
                   .. if snrLimit is not None, then uses this value as a cutoff.
        """
        super(DiscoveryChancesMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit
        self.nObsPerNight = nObsPerNight
        self.tNight = tNight
        self.nNightsPerWindow = nNightsPerWindow
        self.tWindow = tWindow
        self.gamma = 0.038
        self.sigma = 0.12
        self.badval = 0

    def run(self, ssoObs, orb, Hval):
        """SsoObs = Dataframe, orb=Dataframe, Hval=single number."""
        # Calculate visibility for this orbit at this H.
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return self.badval
        else:
            # Now to identify where observations meet the timing requirements.
            #  Identify visits where the 'night' changes.
            visSort = np.argsort(ssoObs[self.nightCol])[vis]
            n = np.unique(ssoObs[self.nightCol][visSort])
            # Identify all the indexes where the night changes (swap from one night to next)
            nIdx = np.searchsorted(ssoObs[self.nightCol][visSort], n)
            # Count the number of observations per night (except last night)
            obsPerNight = (nIdx - np.roll(nIdx, 1))[1:]
            # Add the number of observations on the last night.
            obsLastNight = np.array([len(visSort) - nIdx[-1]])
            obsPerNight = np.concatenate((obsPerNight, obsLastNight))
            # Find the nights with at least nObsPerNight visits.
            nWithXObs = n[np.where(obsPerNight >= self.nObsPerNight)]
            nIdxMany = np.searchsorted(ssoObs[self.nightCol][visSort], nWithXObs)
            nIdxManyEnd = np.searchsorted(ssoObs[self.nightCol][visSort], nWithXObs, side='right') - 1
            # Check that nObsPerNight observations are within tNight
            timesStart = ssoObs[self.expMJDCol][visSort][nIdxMany]
            timesEnd = ssoObs[self.expMJDCol][visSort][nIdxManyEnd]
            # Identify the nights with 'clearly good' observations.
            good = np.where(timesEnd - timesStart <= self.tNight, 1, 0)
            # Identify the nights where we need more investigation (a subset of the visits may be within the interval).
            check = np.where((good==0) & (nIdxManyEnd + 1 - nIdxMany > self.nObsPerNight) & (timesEnd-timesStart > self.tNight))[0]
            for i, j, c in zip(visSort[nIdxMany][check], visSort[nIdxManyEnd][check], check):
                t = ssoObs[self.expMJDCol][i:j+1]
                dtimes = (np.roll(t, 1- self.nObsPerNight) - t)[:-1]
                if np.any(dtimes <= self.tNight):
                    good[c] = 1
            # 'good' provides mask for observations which could count as 'good to make tracklets' against ssoObs[visSort][nIdxMany]
            # Now identify tracklets which can make tracks.
            goodIdx = visSort[nIdxMany][good == 1]
            #print 'good tracklet nights', ssoObs[self.nightCol][goodIdx]
            # Now (with indexes of start of 'good' nights with nObsPerNight within tNight),
            # look at the intervals between 'good' nights (for tracks)
            if len(goodIdx) < self.nNightsPerWindow:
                discoveryChances = 0
            else:
                dnights = (np.roll(ssoObs[self.nightCol][goodIdx], 1-self.nNightsPerWindow) - ssoObs[self.nightCol][goodIdx])
                discoveryChances = len(np.where((dnights >= 0) & (dnights <= self.tWindow))[0])
        return discoveryChances


class MagicDiscoveryMetric(BaseMoMetric):
    """
    Count the number of discovery opportunities with very good software.
    """
    def __init__(self, nObs=6, tWindow=30, snrLimit=None, **kwargs):
        """
        @ nObs = the total number of observations required for 'discovery'
        @ tWindow = the timespan of the discovery window.
        @ snrLimit .. if snrLimit is None then uses 'completeness' calculation,
                   .. if snrLimit is not None, then uses this value as a cutoff.
        """
        super(MagicDiscoveryMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit
        self.nObs = nObs
        self.tWindow = tWindow
        self.badval = 0

    def run(self, ssoObs, orb, Hval):
        """SsoObs = Dataframe, orb=Dataframe, Hval=single number."""
        # Calculate visibility for this orbit at this H.
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return self.badval
        tNights = np.sort(ssoObs[self.nightCol][vis])
        deltaNights = np.roll(tNights, 1-self.nObs) - tNights
        nDisc = np.where((deltaNights < self.tWindow) & (deltaNights >= 0))[0].size
        return nDisc

class HighVelocityMetric(BaseMoMetric):
    """
    Count the number of times an asteroid is observed with a velocity high enough to make it appear
    trailed by a factor of (psfFactor)*PSF - i.e. velocity >= psfFactor * seeing / visitExpTime.
    Simply counts the total number of observations with high velocity.
    """
    def __init__(self, psfFactor=2.0,  snrLimit=None, **kwargs):
        """
        @ psfFactor = factor to multiply seeing/visitExpTime by (velocity(deg/day) >= 24*psfFactor*seeing(")/visitExptime(s))
        """
        super(HighVelocityMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit
        self.psfFactor = psfFactor
        self.badval = 0

    def run(self, ssoObs, orb, Hval):
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return self.badval
        highVelocityObs = np.where(ssoObs['velocity'][vis] >=
                                   (24.*  self.psfFactor * ssoObs[self.seeingCol][vis] / ssoObs[self.expTimeCol][vis]))[0]
        return highVelocityObs.size

class HighVelocityNightsMetric(BaseMoMetric):
    """
    Count the number of times an asteroid is observed with a velocity high enough to make it appear
    trailed by a factor of (psfFactor)*PSF - i.e. velocity >= psfFactor * seeing / visitExpTime,
    where we require nObsPerNight observations within a given night.
    Counts the total number of nights with enough high-velocity observations.
    """
    def __init__(self, psfFactor=2.0, nObsPerNight=2, snrLimit=None, **kwargs):
        """
        @ psfFactor = factor to multiply seeing/visitExpTime by (velocity(deg/day) >= 24*psfFactor*seeing(")/visitExptime(s))
        @ nObsPerNight = number of observations required per night
        """
        super(HighVelocityNightsMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit
        self.psfFactor = psfFactor
        self.nObsPerNight = nObsPerNight
        self.badval = 0

    def run(self, ssoObs, orb, Hval):
        if self.snrLimit is not None:
            vis = np.where(ssoObs[self.snrCol] >= self.snrLimit)[0]
        else:
            vis = np.where(ssoObs[self.visCol] > 0)[0]
        if len(vis) == 0:
            return self.badval
        highVelocityObs = np.where(ssoObs['velocity'][vis] >=
                                   (24.*  self.psfFactor * ssoObs[self.seeingCol][vis] / ssoObs[self.expTimeCol][vis]))[0]
        if len(highVelocityObs) == 0:
            return self.badval
        nights = ssoObs[self.nightCol][vis][highVelocityObs]
        n = np.unique(nights)
        nIdx = np.searchsorted(nights, n)
        # Count the number of observations per night (except last night)
        obsPerNight = (nIdx - np.roll(nIdx, 1))[1:]
        # Add the number of observations on the last night.
        obsLastNight = np.array([len(nights) - nIdx[-1]])
        obsPerNight = np.concatenate((obsPerNight, obsLastNight))
        # Find the nights with at least nObsPerNight visits (this is already looking at only high velocity observations).
        nWithXObs = n[np.where(obsPerNight >= self.nObsPerNight)]
        return nWithXObs.size
