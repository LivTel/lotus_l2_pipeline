'''
  - Produces a plot of the detected spectrum locations and corresponding fit.
  - Checks the curvature of the spectrum is within a defined limit.
'''

from matplotlib import pyplot as plt
import pyfits
import sys
import numpy as np
from optparse import OptionParser

PLOT_PADDING = 20

def execute(f_in, dat_peaks, dat_traces, f_out, plt_title, max_curvature, save=True, hold=False):
  hdulist = pyfits.open(f_in)
  data = hdulist[0].data
  
  plt.imshow(data, aspect='auto', vmin=np.median(data), vmax=np.percentile(data, 99.5))
  
  x = []
  y = []
  with open(dat_peaks) as f:
    for line in f:
      if not line.startswith('#') and line.strip() != "" and not line.startswith('-1'):
        x.append(float(line.split('\t')[0]))
        y.append(float(line.split('\t')[1]))
        
        coeffs = []
        with open(dat_traces) as f:
          for line in f:
            if line.startswith('# Polynomial Order'):
              ncoeffs = int(line.split('\t')[1]) + 1  
            elif not line.startswith('#') and line.strip() != "" and not line.startswith('-1'):
              for i in range(ncoeffs):
                coeffs.append(float(line.split('\t')[i]))
                  
  x_fitted = range(0, len(data[0]))
  y_fitted = np.polyval(coeffs[::-1], x_fitted)
        
  # plot
  y_min = min(y_fitted) - PLOT_PADDING
  y_max = max(y_fitted) + PLOT_PADDING
                  
  plt.plot(x, y, 'ko')
  plt.plot(x_fitted, y_fitted, 'r-')
  plt.colorbar()
  plt.ylim([y_min, y_max])
  plt.title(plt_title)
  plt.xlabel("x")
  plt.ylabel("y")
  
  if save:
      plt.savefig(f_out)
  if not hold:
      plt.clf()  
  
  # check deviation and set return code
  deviation = np.max(y_fitted) - np.min(y_fitted)
  if deviation > float(max_curvature):
    return 1
  else:
    return 0            

if __name__ == "__main__":
  parser = OptionParser()
  parser.add_option('--f', dest='f_in', action='store', default='out.fits')
  parser.add_option('--p', dest='dat_peaks', action='store', default='lofind_peaks.dat')
  parser.add_option('--t', dest='dat_traces', action='store', default='lotrace_traces.dat')
  parser.add_option('--o', dest='f_out', action='store', default='plot.png')
  parser.add_option('--ot', dest='plt_title', action='store', default='plot')
  parser.add_option('--c', dest='max_curvature', action='store', default=1.0)
  (options, args) = parser.parse_args()
  
  f_in		= str(options.f_in)
  dat_peaks  	= str(options.dat_peaks)
  dat_traces 	= str(options.dat_traces)
  f_out		= str(options.f_out)
  plt_title	= str(options.plt_title)
  max_curvature = float(options.max_curvature)
  
  execute(f_in, dat_peaks, dat_traces, f_out, plt_title, max_curvature)
  
  
