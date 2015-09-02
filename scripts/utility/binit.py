import pyfits
import sys
import numpy as np

f = pyfits.open(sys.argv[1])

a_new = []
for idx, i in enumerate(f[0].data):
    if idx % 2 == 1:
        l2 = np.array(i)
        a_new.append((l1+l2))
    else:
        l1 = np.array(i)

a_new = np.swapaxes(a_new,0,1)
a_new2 = []
for idx, i in enumerate(a_new):
    if idx % 2 == 1:
        l2 = np.array(i)
        a_new2.append((l1+l2))
    else:
        l1 = np.array(i)

a_new2 = np.swapaxes(a_new2,0,1)

# change header keys
a_new2_header = f[0].header

a_new2_header['CCDXBIN'] = 2
a_new2_header['CCDYBIN'] = 2

hdu = pyfits.PrimaryHDU(np.asarray(a_new2), header=f[0].header)
hdu.writeto(sys.argv[2])
f.close()
