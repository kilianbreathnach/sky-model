import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import os
import logging
import sys
import numpy as np
import random as rd
from glob import glob
import pyfits
import tractor
import tractor.sdss_galaxy as sdss
from astrometry.util.util import Tan
import pic_compare as bpl


# base branches to the directory for web publishing

baseweb = "/home/kilian/public_html/try_galaxy-cm/by_hand/optplot"

# extra bits for this attempt

extrabits = "try4/psf5"


# The default:
#lvl = logging.INFO
lvl = logging.DEBUG
logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

# Bit o' logging n' bookkeeping

class Logger(object):
    def __init__(self, filename, err_type):
        self.terminal = err_type       # sys.stdout or sys.stderr
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

sys.stdout = Logger("{0}/{1}/log.txt".format(baseweb, extrabits),
        sys.stdout)
sys.stderr = Logger("{0}/{1}/error.txt".format(baseweb, extrabits),
        sys.stderr)


def optloop(tr):
    while True:
        dlnp, X, alpha = tr.optimize()
        if dlnp < 1.00e-3:
            break

def optplot(image, tr, pub_name, comp, step, param):
            residual = abs(image.getImage() - tr.getModelImage(image))
            bpl.wristband(image.getImage(), tr.getModelImage(image),
                    residual, tr.getChiImage(),
                    "{0}/{1}/{2}-{3}_step{4}_{5}.png".format(baseweb,
                        extrabits, pub_name, comp, step, param))


MAGIC_NO = 1.0484
'''For converting MAD to standard deviation.'''

def MAD(image):
    '''
    Get Median Absolute Deviation of image.
    image should be an np array.
    '''
    row, col = image.shape
    sample_list = []
    for i in range(50):
        x, y = rd.randrange(0, row - 6), rd.randrange(0, col - 5)
        sample_list += [ abs(np.float64(image[x, y]) - np.float64(image[x + 5,
            y + 5])) ]

    return np.median(np.array(sample_list))

def Initial_PSF(FWHM,double=False):

    # NB. FWHM of PSF is given in pixels.

    if not double:
        # Single Gaussian default:
        w = np.array([1.0])
        mu = np.array([[0.0,0.0]])               # centroid position in pixels
        var = (FWHM/2.35)**2.0
        cov = np.array([[[var,0.0],[0.0,var]]])  # pixels^2, covariance matrix

    else:
        # Double Gaussian alternative:
        w = np.array([0.75,0.25])
        mu = np.array([[0.0,0.0],[0.0,0.0]])
        var = (FWHM/2.35)**2.0
        cov = np.array([[[var,0.0],[0.0,var]],[[4*var,0.0],[0.0,4*var]]])

    return tractor.GaussianMixturePSF(w,mu,cov)



wcsfile = "/data2/nova/BACKUP-jobs/00000046/wcs.fits"

folder = os.path.dirname(wcsfile)
pub_name = folder.split("/")[-1]    # to save on the interwebs
rdlsfiles = glob("./rdls*fits")
log = folder + "/job.log"

words = open(log).readlines()[3].split()

imagef_ind = words.index("--image") + 1    # wrong filename for the image
badfile = words[imagef_ind]
getpath = badfile.split("/")[-4:]
goodfile = "/home/nova/BACKUP-data/{0}/{1}/{2}/{3}".format(*tuple(getpath))
'''Now we have the right filename for the original.'''


# Now finally the real bacon

try:
    fits = pyfits.open(goodfile)
    picture = fits[0].data
    if type(picture) != np.ndarray:
        picture = fits[1].data
except (IOError):
    picture = mpimg.imread(goodfile)[::-1]


if len(picture.shape) == 3:
    dim = picture.shape[2]
else:
    dim = 1

for comp in range(dim):
    '''this loop works on each component of the RGB image'''
    if dim == 1:
        pic = picture
    else:
        pic = picture[:,:,comp]

    pic = pic.astype('float64')
    lenX, lenY = pic.shape
    med = np.median(pic)
    std_dev = MAGIC_NO * MAD(pic)
    if std_dev <= 0.000001:
        print pub_name
        break
    invvar = (1. / std_dev**2) * np.ones((lenX, lenY))
    '''Making the inverse variance matrix guess.'''


    # Make a wcs object

    ext = 0
    wcs = tractor.FitsWcs(Tan(wcsfile, ext))


    # PSF

    psf = Initial_PSF(5., double=True)


    # Get a background guess for the sky noise

    sky = tractor.ConstantSky(med)


    # Get a crude photocalibration

    photocal = tractor.MagsPhotoCal("r",12.)


    # Now to generate the image...

    image = tractor.Image(data=pic, invvar=invvar, psf=psf, wcs=wcs,
                          sky=sky, photocal=photocal, name="pic")

    # and the Catalogue...

    catalog =[]

    for rdfn in rdlsfiles:
        rdls = pyfits.open(rdfn)
        posdat = rdls[1].data

        for j in range(len(posdat)):
            ra, dec = posdat[j]
            t_radec = tractor.RaDecPos(ra, dec)
            src_xy = wcs.positionToPixel(t_radec)
            cent = np.rint(src_xy).astype("int")
            if cent[0] < 0 or cent[0] >= lenY or cent[1] < 0 or cent[1] >= lenX:
                continue
            else:
                mag = tractor.Mags(r=6.)
                radec = tractor.RaDecPos(posdat[j][0], posdat[j][1])
                ps = tractor.PointSource(radec, mag)
                catalog.append(ps)

    print "sources found:"
    print catalog


    # Belt on th'auld Tractor:

    tr = tractor.Tractor(images=[image], catalog=catalog)


    # Now to optimise

    print "STEP 0:"
    step = 0

    tr.thawParam('catalog')
    tr.freezeParamsRecursive('images')

    for src in tr.getCatalog():
        src.thawAllRecursive()
        src.freezeAllBut('brightness')

    print "And the thawed params for brightness optimisation are:"
    for nm in tr.getParamNames():
        print " ", nm

#    image.freezeAllBut('sky')
    optloop(tr)
    param = "brightness" # "skynbright"
    optplot(image, tr, pub_name, comp, step, param)
    step += 1


    print "STEP 1"

    cat = tr.getCatalog()
    '''Initial log probability.'''
#    lnp0 = tr.getLogProb()
    '''Switch source i to a galaxy.'''
    problist = []
    for i in range(10): # len(cat)):

        oldsrc = cat[i]

        initprob = tr.getLogProb()
        print "Source ", i, " has initial probability: ", initprob

        print "initially, galaxy is ", oldsrc
        '''
        args are r_e (in arcsec), b/a ratio,
        position angle phi (in deg)
        '''
        shape = sdss.GalaxyShape(32., 1., 0.)
        bri = oldsrc.getBrightness()
        print "brightness is ", bri
        for j, p in enumerate(bri.getParams()):
            bri.setParam(j, p - 2)
        print "brightness is ", bri
        gal = sdss.ExpGalaxy(oldsrc.getPosition(), bri, shape)
        cat[i] = gal

        print "after we change, galaxy is ", cat[i]

        tr.freezeParamsRecursive('images')
        tr.thawParam('catalog')

        for src in tr.getCatalog():
            src.freezeAllBut('brightness')

        cat[i].thawAllRecursive()
        cat[i].freezeParam('pos')

        print tr.getParamNames()
        tr.optimize()


        print "and after optimising, galaxy is: ", cat[i]

        finprob = tr.getLogProb()
        print "Source ", i, " has final probability: ", finprob
        problist.append(finprob - initprob)
        print "Probability difference is: ", problist[i]

        if problist[i] < 12:
            cat[i] = oldsrc
        else:
            pass


    print "FINISHED SWITCHING TO GALAXIES"

    tr.freezeParamsRecursive('images')
    tr.thawParam('catalog')

    for src in tr.getCatalog():
        src.thawAllRecursive()
        src.freezeParam('pos')

    optloop(tr)
    param = "allsrcbutpos"
    optplot(image, tr, pub_name, comp, step, param)
    step += 1


    print "STEP 2:"

    '''sky + flux + gal shape'''
    tr.freezeParamsRecursive('catalog')
    tr.freezeParamsRecursive('images')

    for src in tr.getCatalog():
        src.thawAllRecursive()
        src.freezeParam('pos')

    image.freezeAllBut('sky')
    optloop(tr)
    param = "skynbrightngal"
    optplot(image, tr, pub_name, comp, step, param)
    step += 1


    '''psf'''
    tr.freezeParamsRecursive('catalog')
    tr.thawParams('images')
    image.freezeAllBut('psf')

    print "Starting psf optimisation"
    print tr.getNamedParams()
    optloop(tr)
    print "End of psf optimisation"
    print tr.getNamedParams()
    param = "psf"
    optplot(image, tr, pub_name, comp, step, param)

