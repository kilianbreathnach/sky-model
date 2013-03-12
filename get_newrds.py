from glob import glob
import os
from astrometry.util.util import Tan


wcsfn = "/data2/nova/BACKUP-jobs/00000046/wcs.fits"

axyfn = wcsfn.replace('wcs.fits', 'job.axy')
dirname = os.path.basename(os.path.dirname(wcsfn))


wcs = Tan(wcsfn)
r1,d1 = wcs.radec_center()
R1 = wcs.radius()

print "RA is ", r1
print "Dec is ", d1
print "radius is ", R1

inds = glob("/data1/INDEXES/4000/index-4001*")

for j, ind in enumerate(inds):

    cmdline = 'query-starkd -v -o rdls4001-{0}.fits -r {1} -d {2} -R {3} {4}'.format(j, r1, d1, R1, ind)

    os.system(cmdline)
