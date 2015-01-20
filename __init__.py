# Licensed under a 3-clause BSD style license - see LICENSE.rst 

__version__ = "0.1"
__author__ = "John Antoniadis"

from .genutils import polyfit2d, polyval2d, bestfit, sfile, doppler
from .fitdaspectrum import DAModel, Fit, fitwdmodel
from .findvelocity import findvel, minchi2
from .plotting import plot_da_fit
