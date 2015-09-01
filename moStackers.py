import numpy as np
from lsst.sims.maf.stackers import BaseStacker

class AppMagStacker(object):
    """
    Add the apparent magnitude of an object with a given Hval to the observations dataframe.
    """
    def __init__(self, magFilterCol='magFilter'):
        self.magFilterCol = magFilterCol
        self.colsReq = [self.magFilterCol]
        self.colsAdded = ['appMag']

    def run(self, ssoObs, Href, Hval):
        ssoObs['appMag'] = ssoObs[self.magFilterCol] + Hval - Href
        return ssoObs


class MagLimitStacker(object):
    """
    Add the apparent magnitude limit with trailing or detection losses to the observations dataframe.
    """
    def __init__(self, m5Col='fiveSigmaDepth', lossCol='dmagDetect'):
        self.m5Col = m5Col
        self.lossCol = lossCol
        self.colsReq = [self.m5Col, self.lossCol]
        self.colsAdded = ['magLimit']

    def run(self, ssoObs, Href=None, Hval=None):
        ssoObs['magLimit'] = ssoObs[self.m5Col] - ssoObs[self.lossCol]
        return ssoObs

class SNRStacker(object):
    """
    Add the SNR to the observations dataframe.
    """
    def __init__(self, magLimitCol='magLimit', appMagCol='appMag',gamma=0.038):
        self.appMagCol = appMagCol
        self.magLimitCol = magLimitCol
        self.colsReq = [self.appMagCol, self.magLimitCol]
        self.colsAdded = ['SNR']
        self.gamma = gamma

    def run(self, ssoObs, Href=None, Hval=None):
        xval = np.power(10, 0.5*(ssoObs[self.appMagCol] - ssoObs[self.magLimitCol]))
        ssoObs['SNR'] = 1.0 / np.sqrt((0.04 - self.gamma)*xval + self.gamma*xval*xval)
        return ssoObs

class VisStacker(object):
    """
    Calculate whether an object is visible according to
    Fermi-Dirac completeness formula (see SDSS, eqn 24, Stripe82 analysis:
    http://iopscience.iop.org/0004-637X/794/2/120/pdf/apj_794_2_120.pdf).
    Calculate estimated completeness/probability of detection,
    then evaluates if this object could be visible.
    """
    def __init__(self, magLimitCol='magLimit',
                appMagCol='appMag', sigma=0.12):
        self.magLimitCol = magLimitCol
        self.appMagCol = appMagCol
        self.sigma = sigma
        self.colsReq = [self.magLimitCol, self.appMagCol]
        self.colsAdded = ['vis']

    def run(self, ssoObs, Href=None, Hval=None):
        completeness = 1.0 / (1 + np.exp((ssoObs[self.appMagCol] - ssoObs[self.magLimitCol])/self.sigma))
        probability = np.random.random_sample(len(ssoObs[self.appMagCol]))
        ssoObs['vis'] = np.where(probability <= completeness, 1, 0)
        return ssoObs


class AllStackers(object):
    """
    Since for moving objects we usually want to run all of these at once/in order,
    provide a convenient way to do it.
    """
    def __init__(self):
        self.appMag = AppMagStacker()
        self.magLimit = MagLimitStacker()
        self.snr = SNRStacker()
        self.vis = VisStacker()

    def run(self, ssoObs, Href, Hval=None):
        if Hval is None:
            Hval = Href
        ssoObs = self.appMag.run(ssoObs, Href, Hval)
        ssoObs = self.magLimit.run(ssoObs, Href, Hval)
        ssoObs = self.snr.run(ssoObs, Href, Hval)
        ssoObs = self.vis.run(ssoObs, Href, Hval)
        return ssoObs
