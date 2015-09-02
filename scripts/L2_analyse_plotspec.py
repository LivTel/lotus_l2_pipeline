'''
  - Produces a plot of the final spectrum.
'''

import matplotlib
from matplotlib import pyplot as plt
import pyfits
import sys
import numpy as np
from optparse import OptionParser

PLOT_PADDING = 20

def execute(f_in, hdu, f_out, plt_title, color, legend=False, leg_title="", save=True, hold=False):
  hdulist = pyfits.open(f_in)
  data = hdulist[hdu].data
  hdrs = hdulist[hdu].header
  
  CRVAL1        = hdrs['CRVAL1']
  CDELT1        = hdrs['CDELT1']
  NAXIS1        = hdrs['NAXIS1']
  
  x = []
  for i in range(0, NAXIS1, 1):
      x.append(CRVAL1 + i*CDELT1)
  
  y = data.reshape(NAXIS1)
        
  # plot
  y_min = 0 # min(y) - PLOT_PADDING
  y_max = max(y)*2 # max(y) + PLOT_PADDING
                  
  plt.plot(x, y, label=leg_title, color=color)
  plt.ylim([y_min, y_max])
  plt.title(plt_title)
  plt.xlabel("Wavelength ($\AA$)")
  plt.ylabel("Intensity (counts)")
  if not isinstance(data[0], np.ndarray):    # this occurs if the hdu contains no data..   
      ax = plt.gca()
      ax.spines['bottom'].set_color('red')
      ax.spines['top'].set_color('red') 
      ax.spines['right'].set_color('red')
      ax.spines['left'].set_color('red')
      plt.annotate("NO DATA", xy=(0.5,0.5), xycoords="axes fraction", fontsize=42, horizontalalignment='center', verticalalignment='center')

  if save:
      plt.savefig(f_out)
  if not hold:
      plt.clf() 
  if legend:
      leg = plt.legend(loc="upper left", fontsize=10) 
      leg.get_frame().set_alpha(0.5)
  
  return 0

if __name__ == "__main__":
  parser = OptionParser()
  parser.add_option('--f', dest='f_in', action='store', default='out.fits')
  parser.add_option('--hdu', dest='hdu', action='store')
  parser.add_option('--o', dest='f_out', action='store', default='plot.png')
  parser.add_option('--ot', dest='plt_title', action='store', default='plot')
  (options, args) = parser.parse_args()
  
  f_in		= str(options.f_in)
  hdu 		= str(options.hdu)
  f_out		= str(options.f_out)
  plt_title	= str(options.plt_title)
  
  execute(f_in, hdu, f_out, plt_title)
  

                  