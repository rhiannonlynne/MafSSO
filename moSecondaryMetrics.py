import numpy as np
import numpy.ma as ma

from moMetrics import BaseMoMetric

__all__ = ['MoCompletenessMetric', 'IntegrateOverHMetric']

class MoCompletenessMetric(BaseMoMetric):
    """
    Take the discoveryChances metric results and turn it into
    completeness estimate (relative to the entire population).
    Require at least 'requiredChances' to count an object as "found".
    """
    def __init__(self, requiredChances=1, nbins=30, metricName='Completeness', **kwargs):
        super(MoCompletenessMetric, self).__init__(metricName=metricName, **kwargs)
        self.requiredChances = requiredChances
        self.nbins = nbins
        self.minHrange = 1.0

    def run(self, discoveryChances, Hvals):
        nSsos = discoveryChances.shape[0]
        discoveriesH = discoveryChances.swapaxes(0, 1)
        if len(Hvals) == len(discoveriesH):
            # Then we specified Hvals as an array, so use this for bins.
            completeness = ma.MaskedArray(data = np.zeros(len(Hvals), float),
                                          mask = np.zeros(len(Hvals), 'bool'),
                                          fill_value = 0.0)
            for i, H in enumerate(Hvals):
                completeness.data[i] = np.where(discoveriesH[i].filled(0) >= self.requiredChances)[0].size
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
            data = n_found.astype(float) / n_all.astype(float)
            mask = np.where(n_all==0, True, False)
            completeness = ma.MaskedArray(data = data,
                                          mask = mask,
                                          fill_value = 0.0)
        return completeness, Hvals


class IntegrateOverHMetric(BaseMoMetric):
    """
    Take a calculated (differential H, "H @ X") metric value and
    integrate over a size distribution to return an "H<=X" value.
    """
    def __init__(self, Hindex=0.3, **kwargs):
        """
        Currently only supports single power law H distribution.
        """
        super(IntegrateOverHMetric, self).__init__(*kwargs)
        self.Hindex = Hindex
        self.units = '<= H'

    def run(self, metricVals, Hvals):
        # Set expected H distribution.
        # dndh = differential size distribution (number in this bin)
        dndh = np.power(10., self.Hindex*(Hvals-Hvals.min()))
        # dn = cumulative size distribution (number in this bin and brighter)
        dn = np.cumsum(dndh)
        intVals = np.cumsum(metricVals*dndh)/dn
        return intVals, Hvals
