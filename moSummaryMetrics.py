import numpy as np

from moMetrics import BaseMoMetric

__all__ = ['MoCompletenessMetric', 'IntegrateOverHMetric']

class MoCompletenessMetric(BaseMoMetric):
    """
    Take the discoveryChances metric results and turn it into
    completeness estimate (relative to the entire population).
    Require at least 'requiredChances' to count an object as "found".
    """
    def __init__(self, requiredChances=1, nbins=30, **kwargs):
        super(MoCompletenessMetric, self).__init__(**kwargs)
        self.requiredChances = requiredChances
        self.nbins = nbins

    def run(self, discoveryChances, Hvals):
        nSsos = discoveryChances.shape[0]
        discoveries = discoveryChances.swapaxes(0, 1)
        if len(Hvals) == discoveryChances.shape[1]:
            # Then we specified Hvals as an array, so use this for bins.
            completeness = np.zeros(len(Hvals), float)
            for i, H in enumerate(Hvals):
                completeness[i] = np.where(discoveries[i] >= self.requiredChances)[0].size
        else:
            # The Hvals are spread more randomly among the objects (we probably used one per object).
            stepsize = (Hvals.max() - Hvals.min()) / self.nbins
            bins = np.arange(Hvals.min(), Hvals.max() + stepsize, stepsize)
            n_all, b = np.histogram(discoveries[0], bins)
            condition = np.where(discoveries[0] >= self.requiredChances)[0]
            n_found, b = np.histogram(discoveries[0][condition], bins)
            completeness = float(n_found) / float(n_all)
        return completeness


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

    def run(metricVals, Hrange):
        # Set expected H distribution.
        # dndh = differential size distribution (number in this bin)
        dndh = np.power(10., Hindex*(Hrange-Hrange.min()))
        # dn = cumulative size distribution (number in this bin and brighter)
        dn = np.cumsum(dndh)
        intVals = np.cumsum(metricVals*dndh)/dn
        return intVals
