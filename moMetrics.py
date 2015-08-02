import inspect
import numpy as np
import numpy.ma as ma
from lsst.sims.maf.metrics import MetricRegistry

__all__ = ['BaseMoMetric', 'NObsMetric', 'DiscoveryChancesMetric',
           'NNightsMetric', 'ObsArcMetric',
           'ActivityOverTimeMetric', 'ActivityOverPeriodMetric']


class BaseMoMetric(object):
    """Base class for the moving object metrics."""
    __metaclass__ = MetricRegistry

    def __init__(self, metricName=None, units='#', badval=0,
                 m5Col='fiveSigmaDepth', lossCol='dmagDetect',
                 magFilterCol='magFilter',
                 nightCol='night', expMJDCol='expMJD',
                 Hindex = 0.3):
        self.metricDtype = float
        self.name = metricName
        if self.name is None:
            self.name = self.__class__.__name__.replace('Metric', '', 1)
        self.badval = badval
        self.units = units
        self.comment = None
        # Set some commonly used column names.
        self.m5Col = m5Col
        self.lossCol = lossCol
        self.magFilterCol = magFilterCol
        self.nightCol = nightCol
        self.expMJDCol = expMJDCol
        self.snrLimit = None
        self.colsReq = [self.m5Col, self.lossCol, self.magFilterCol,
                        self.nightCol, self.expMJDCol]
        # Set parameters used for reduce methods.
        self.Hindex = Hindex

        # Set up dictionary of reduce functions (may be empty).
        self.reduceFuncs = {}
        self.reduceOrder = {}
        self.reduceUnits = {}
        for i, r in enumerate(inspect.getmembers(self, predicate=inspect.ismethod)):
            if r[0].startswith('reduce'):
                reducename = r[0].replace('reduce', '', 1)
                self.reduceFuncs[reducename] = r[1]
                self.reduceOrder[reducename] = i
                try:
                    self.reduceUnits[reducename] = r[1].units
                except AttributeError:
                    self.reduceUnits[reducename] = '@H'
        self.reduceOrder['CumulativeH'] = len(self.reduceFuncs)


    def _calcAppMag(self, ssoObs, Hval, Href):
        """
        Adjust the apparent magnitude of the object in this filter for any changes to H
         (in case of cloning the objects in this orbit).
        """
        return ssoObs[self.magFilterCol] + Hval - Href

    def _calcMagLimit(self, ssoObs):
        """
        Calculate the effective magnitude limit, accounting for the velocity of the moving object.
        The mag limit must be adjusted by 'dmagDetect' if detection is done on a 'stationary
        object likelihood image' (this is the maximum loss).
        The mag limit should be adjusted by 'dmagTrailing' if only accounting for SNR loss
        due to increased number of sky pixels.
        """
        return ssoObs[self.m5Col] - ssoObs[self.lossCol]

    def _calcSNR(self, appMag, magLimit, gamma=0.038):
        """
        Calculate the SNR of a source with appMag in an image with 'magLimit'.
        """
        xval = np.power(10, 0.5*(appMag - magLimit))
        snr = 1.0 / np.sqrt((0.04 - gamma)*xval + gamma*xval*xval)
        return snr

    def _calcVis(self, appMag, magLimit, sigma=0.12):
        """
        Calculate whether an object is visible according to
        Fermi-Dirac completeness formula (see SDSS, eqn 24, Stripe82 analysis:
         http://iopscience.iop.org/0004-637X/794/2/120/pdf/apj_794_2_120.pdf).
        Calculate estimated completeness/probability of detection,
        then evaluates if this object could be visible.
        """
        completeness = 1.0 / (1 + np.exp((appMag - magLimit)/sigma))
        probability = np.random.random_sample(len(appMag))
        vis = np.where(probability <= completeness)[0]
        return vis

    def _prep(self, ssoObs, orb, Hval):
        """
        We will almost always need to convert the incoming observation + Hval
        into apparent magnitude / see what's visible, etc.
        This is a convenience function for that.
        """
        if len(ssoObs) == 0:
            return ValueError('No data here')
        Href = orb['H']
        if Hval is None:
            Hval = Href
        appMag = self._calcAppMag(ssoObs, Hval, Href)
        magLimit = self._calcMagLimit(ssoObs)
        if self.snrLimit is None:
            snr = None
            vis = self._calcVis(appMag, magLimit)
        else:
            snr = self._calcSNR(appMag, magLimit)
            vis = np.where(snr >= self.snrLimit)[0]
        return appMag, magLimit, vis, snr

    def run(self, ssoObs, orb, Hval):
        raise NotImplementedError


    def reduceCumulativeH(self, metricVals, Hvals):
        """
        Take a calculated (differential H, "H @ X") metric value and
        integrate over a size distribution to return an "H<=X" value.
        Currently supports only a single power law distribution.
        """
        self.units = '<= H'
        # Set expected H distribution.
        # dndh = differential size distribution (number in this bin)
        dndh = np.power(10., self.Hindex*(Hvals-Hvals.min()))
        # dn = cumulative size distribution (number in this bin and brighter)
        intVals = np.cumsum(metricVals*dndh, axis=1)/np.cumsum(dndh)
        return intVals, Hvals


class NObsMetric(BaseMoMetric):
    """
    Count the number of observations for an object.
    """
    def __init__(self, snrLimit=None, **kwargs):
        """
        @ snrLimit .. if snrLimit is None, this uses the _calcVis method/completeness
                      if snrLimit is not None, this uses that value as a cutoff instead.
        """
        super(NObsMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit

    def run(self, ssoObs, orb, Hval):
        try:
            appMag, magLimit, vis, snr = self._prep(ssoObs, orb, Hval)
            return vis.size
        except ValueError:
            return 0


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
        try:
            appMag, magLimit, vis, snr = self._prep(ssoObs, orb, Hval)
        except ValueError:
            return 0
        if len(vis) > 0:
            nights = len(np.unique(ssoObs[self.nightCol][vis]))
        else:
            nights = 0
        return nights

class ObsArcMetric(BaseMoMetric):
    """
    Calculate the difference between the first and last observation of an object.
    """
    def __init__(self, snrLimit=None, **kwargs):
        super(ObsArcMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit

    def run(self, ssoObs, orb, Hval):
        try:
            appMag, magLimit, vis, snr = self._prep(ssoObs, orb, Hval)
        except ValueError:
            return 0
        if len(vis) > 0:
            arc = ssoObs[self.expMJDCol][vis].max() - ssoObs[self.expMJDCol][vis].min()
        else:
            arc = 0
        return arc


class DiscoveryChancesMetric(BaseMoMetric):
    """
    Count the number of discovery opportunities for an object.
    """
    def __init__(self, nObsPerNight=2, tNight=90.*60.,
                 nNightsPerWindow=3, tWindow=15, snrLimit=None,
                 requiredChances=1, **kwargs):
        """
        @ nObsPerNight = number of observations per night required for tracklet
        @ tNight = max time start/finish for the tracklet (seconds)
        @ nNightsPerWindow = number of nights with observations required for track
        @ tWindow = max number of nights in track (days)
        @ snrLimit .. if snrLimit is None then uses 'completeness' calculation,
                   .. if snrLimit is not None, then uses this value as a cutoff.

        Parameters for reduce method (Completeness)
        @ requiredChances = number of possible discovery chances required to count an object as 'found'
        @ nBins = number of bins to split "H" into, if not using cloned H distribution.
        """
        super(DiscoveryChancesMetric, self).__init__(**kwargs)
        self.snrLimit = snrLimit
        self.nObsPerNight = nObsPerNight
        self.tNight = tNight
        self.nNightsPerWindow = nNightsPerWindow
        self.tWindow = tWindow
        self.requiredChances = requiredChances
        # If H is not a cloned distribution, then we need to specify how to bin these values.
        self.nbins = 20
        self.minHrange = 1.0

    def run(self, ssoObs, orb, Hval):
        """SsoObs = Dataframe, orb=Dataframe, Hval=single number."""
        # Calculate visibility for this orbit at this H.
        try:
            appMag, magLimit, vis, snr = self._prep(ssoObs, orb, Hval)
        except ValueError:
            return 0
        # Calculate number of discovery chances.
        if len(vis) == 0:
            discoveryChances = 0
        else:
            # Now to identify where observations meet the timing requirements.
            #  Identify visits where the 'night' changes.
            visSort = np.argsort(ssoObs[self.nightCol])[vis]
            n = np.unique(ssoObs[self.nightCol][visSort])
            # Identify all the indexes where the night changes (swap from one night to next)
            nIdx = np.searchsorted(ssoObs[self.nightCol][visSort], n)
            # Add index pointing to last observation.
            nIdx = np.concatenate([nIdx, np.array([len(visSort)-1])])
            # Find the nights & indexes where there were more than nObsPerNight observations.
            obsPerNight = (nIdx - np.roll(nIdx, 1))[1:]
            nWithXObs = n[np.where(obsPerNight >= self.nObsPerNight)]
            nIdxMany = np.searchsorted(ssoObs[self.nightCol][visSort], nWithXObs)
            nIdxManyEnd = np.searchsorted(ssoObs[self.nightCol][visSort], nWithXObs, side='right') - 1
            # Check that nObsPerNight observations are within tNight
            timesStart = ssoObs[self.expMJDCol][visSort][nIdxMany]
            timesEnd = ssoObs[self.expMJDCol][visSort][nIdxManyEnd]
            # Identify the nights where the total time interval may exceed tNight
            # (but still have a subset of nObsPerNight which are within tNight)
            check = np.where((timesEnd - timesStart > self.tNight) & (nIdxManyEnd + 1 - nIdxMany > self.nObsPerNight))[0]
            bad = []
            for i, j, c in zip(visSort[nIdxMany][check], visSort[nIdxManyEnd][check], check):
                t = ssoObs[self.expMJDCol][i:j+1]
                dtimes = (np.roll(t, 1-nObsPerNight) - t)[:-1]
                if np.all(dtimes > self.tNight+eps):
                    bad.append(c)
            goodIdx = np.delete(visSort[nIdxMany], bad)
            # Now (with indexes of start of 'good' nights with nObsPerNight within tNight),
            # look at the intervals between 'good' nights (for tracks)
            if len(goodIdx) < self.nNightsPerWindow:
                discoveryChances = 0
            else:
                dnights = (np.roll(ssoObs[self.nightCol][goodIdx], 1-self.nNightsPerWindow) - ssoObs[self.nightCol][goodIdx])
                discoveryChances = len(np.where((dnights >= 0) & (dnights <= self.tWindow))[0])
        return discoveryChances

    def reduceCompleteness(self, discoveryChances, Hvals):
        """
        Take the discoveryChances metric results and turn it into
        completeness estimate (relative to the entire population).
        Require at least 'requiredChances' to count an object as "found".
        """
        nSsos = discoveryChances.shape[0]
        nHval = len(Hvals)
        discoveriesH = discoveryChances.swapaxes(0, 1)
        if nHval == discoveryChances.shape[1]:
            # Hvals array is probably the same as the cloned H array.
            completeness = ma.MaskedArray(data = np.zeros([1, nHval], float),
                                          mask = np.zeros([1, nHval], 'bool'),
                                          fill_value = 0.0)
            for i, H in enumerate(Hvals):
                completeness.data[0][i] = np.where(discoveriesH[i].filled(0) >= self.requiredChances)[0].size
            completeness = completeness / float(nSsos)
        else:
            # The Hvals are spread more randomly among the objects (we probably used one per object).
            hrange = Hvals.max() - Hvals.min()
            minH = Hvals.min()
            if hrange < self.minHrange:
                hrange = self.minHrange
                minH = Hvals.min() - hrange/2.0
            stepsize = hrange / float(self.nbins)
            bins = np.arange(minH, minH + hrange + stepsize/2.0, stepsize)
            Hvals = bins[:-1]
            n_all, b = np.histogram(discoveriesH[0], bins)
            condition = np.where(discoveriesH[0] >= self.requiredChances)[0]
            n_found, b = np.histogram(discoveriesH[0][condition], bins)
            completeness = ma.MaskedArray(data = np.zeros([1, len(Hvals)], float),
                                          mask = np.zeros([1, len(Hvals)], bool),
                                          fill_value = 0.0)
            completeness.data[0] = n_found.astype(float) / n_all.astype(float)
            completeness.mask[0] = np.where(n_all==0, True, False)
        return completeness, Hvals


class ActivityOverTimeMetric(BaseMoMetric):
    """
    Count the time periods where we would have a chance to detect activity on
    a moving object.
    Splits observations into time periods set by 'window', then looks for observations within each window,
    and reports what fraction of the total windows receive 'nObs' visits.
    """
    def __init__(self, window, snrLimit=5, surveyYears=10.0, **kwargs):
        metricName = 'ActivityOver%.1fDays' %(window)
        super(ActivityOverTimeMetric, self).__init__(metricName=metricName, **kwargs)
        self.snrLimit = snrLimit
        self.window = window
        self.surveyYears = surveyYears
        self.units = '%f Day Windows' %(self.window)

    def run(self, ssoObs, orb,  Hval):
        # For cometary activity, expect activity at the same point in its orbit at the same time, mostly
        # For collisions, expect activity at random times
        windowBins = np.arange(0, self.surveyYears*365 + self.window/2.0, self.window)
        nWindows = len(windowBins)
        try:
            appMag, magLimit, vis, snr = self._prep(ssoObs, orb, Hval)
        except ValueError:
            return 0
        if len(vis) == 0:
            return 0
        else:
            n, b = np.histogram(ssoObs[vis][self.nightCol], bins=windowBins)
            activityWindows = np.where(n>0)[0].size
        return activityWindows / float(nWindows)


class ActivityOverPeriodMetric(BaseMoMetric):
    """
    Count the fraction of the orbit (when split into nBins) that receive
    observations, in order to have a chance to detect activity.
    """
    def __init__(self, window, snrLimit=5, nBins=10,
                 aCol='a', tPeriCol='tPeri', **kwargs):
        super(ActivityOverPeriodMetric, self).__init(**kwargs)
        self.aCol = aCol
        self.tPeriCol = tPeriCol
        self.snrLimit = snrLimit
        self.nBins = nBins
        self.units = '%d bins' %(self.nBins)

    def run(self, ssoObs, orb, Hval):
        # For cometary activity, expect activity at the same point in its orbit at the same time, mostly
        # For collisions, expect activity at random times
        period = np.power(orb[self.aCol], 3./2.) * 365.25
        anomaly = ((ssoObs[self.expMJDCol] - orb[self.tPeriCol]) / period) % (2*np.pi)
        binsize = 2*np.pi / float(self.nBins)
        anomalyBins = np.arange(0, 2*np.pi + binsize/2.0, binsize)
        try:
            appMag, magLimit, vis, snr = self._prep(ssoObs, orb, Hval)
        except ValueError:
            return 0
        if len(vis) == 0:
            return 0
        else:
            n, b = np.histogram(anomaly[vis], bins=anomalyBins)
            activityWindows = np.where(n>0)[0].size
        return activityWindows / float(nBins)
