import numpy as np
import numpy.ma as ma
import warnings

from moMetrics import BaseMoMetric

__all__ = ['ValueAtHMetric']

class ValueAtHMetric(BaseMoMetric):
    """
    Return the value of a metric at a given H.
    """
    def __init__(self, Hmark=22, **kwargs):
        metricName = 'Value At H=%.1f' %(Hmark)
        super(ValueAtHMetric, self).__init__(metricName=metricName, **kwargs)
        self.Hmark = Hmark

    def run(self, metricVals, Hvals):
        # Check if desired H value is within range of H values.
        if (self.Hmark < Hvals.min()) or (self.Hmark > Hvals.max()):
            warnings.warn('Desired H value of metric outside range of provided H values.')
            return None
        nHvals = len(Hvals)
        nHMetricVals = metricVals.shape[1]
        if nHvals == nHMetricVals:
            # Hvals matched the points where the metric values were calculated (clone H distribution).
            eps = 1.0e-6
            # Hvals is an array used for each metric value,
            #  we have to pick out the particular metricValues to use.
            diffs = np.abs(self.Hmark - Hvals)
            Hidx = np.where(diffs == diffs.min())[0]
            result = metricVals.swapaxes(0,1)[Hidx]
            Hmark = Hvals[Hidx]
            self.name = 'Value At H=%.1f' %(Hmark)
        else:
            # We have a range of metric values, one per Hval.
            result = np.interpolate([self.Hmark], Hvals, metricVals.swapaxes(0, 1))
        return result
