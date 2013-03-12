import numpy as np
import tractor
import tractor.sdss_galaxy as sdss
from astrometry.util.util import Tan


def star2dev():



def star2exp(tracobj):
    '''
    Expand point sources into an exponential form galaxy if it works out
    better for the model.

    '''
    cat = tracobj.getCatalog()

    for src in cat:
        '''Initial log probability.'''
        lnp0 = tr.getLogProb()
        '''Switch source i to a galaxy.'''
        oldsrc = cat[i]
        '''
        args are r_e (in arcsec), b/a ratio,
        position angle phi (in deg)
        '''
        shape = sdss.GalaxyShape(1., 10., 0.)
        gal = sdss.ExpGalaxy(oldsrc.getPosition(), oldsrc.getBrightness(),
                             shape)
        cat[i] = gal

        # try to get brightness thawed
        cat[i].freezeAllBut('brightness')


        print 'Thawed parameters:'
        for nm in tr.getNamedParams():
            print ' ', nm

        tracobj.optimize()

        lnp1 = tracobj.getLogProb()
        penalty = 3    # MAGIC_NUM - let's see about this
        if lnp1 - penalty > lnp0:
            '''Keep galaxy'''
            pass
        else:
            '''lave her be!'''
            cat[i] = oldsrc






def exp2composite():




