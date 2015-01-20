'''Last Change 24/11/2013'''

from __future__ import division 

import numpy as np
from numpy.polynomial.polynomial import polyfit 
from numpy.polynomial.polynomial import polyval

from scipy.interpolate import interp1d
import scipy.ndimage.filters as flt

from matplotlib import pyplot as plt

from astropy.table import Table, Column

from .findvelocity import *
from .genutils import *


import os


class DAModel(object):
      def __init__(self,teff,logg,seeing=0.75,slit_width=1.5,resolution=0.15,dispersion=0.47):
                 self._teff = teff
                 self._logg = logg     
                 self._slit = slit_width
                 self._seeing = seeing
                 self._res = resolution
                 self._disp = dispersion
                 self.getModelSpectrum()
      
      @property
      def teff(self):
                 return self._teff

      @teff.setter
      def teff(self,value):
                 if ((value < 6000.) | (value > 20000.)):
                       raise ValueError("Requested Value Outside of Possible Range")
                 else:
                       self._teff = value
                       self.getModelSpectrum()

      @property
      def logg(self):
                 return self._logg

      @logg.setter
      def logg(self,value):
                 if ((value < 5.0) | (value > 9.5)):
                       raise ValueError("Requested Value Outside of Supported Range")
                 else:
                       self._logg = value
                       self.getModelSpectrum()

      @property
      def slit(self):
                 return self._slit

      @slit.setter
      def slit(self,value):
                 if ((value < 0.0)):
                       raise ValueError("Value has to be Possitive")
                 else:
                       self._slit = value
                       self.getModelSpectrum()

      @property
      def seeing(self):
                 return self._seeing

      @seeing.setter
      def seeing(self,value):
                 if ((value < 0.0)):
                       raise ValueError("Value has to be Possitive")
                 else:
                       self._seeing = value
                       self.getModelSpectrum()
     
      @property
      def res(self):
                 return self._res

      @res.setter
      def res(self):
                 if ((value < 0.0)):
                       raise ValueError("Value has to be Possitive")
                 else:
                       self._res = value
                       self.getModelSpectrum()

      @property
      def disp(self):
                 return self._disp

      @disp.setter
      def disp(self):
                 if ((value < 0.0)):
                       raise ValueError("Value has to be Possitive")
                 else:
                       self._disp = value
                       self.getModelSpectrum()

                 
      def getModelSpectrum(self):
                 if (np.mod(self.teff,250)==0 and np.mod(self.logg,0.25)==0):
                   if (self.teff < 10000.0 and self.teff >= 6000.):
                    self.modelname = 'da0' + str(np.int(self.teff)) + '_' + str(np.int(100.*self.logg)) + '.dk'
                    self.openFile()
                   else:
                    self.modelname = 'da' + str(np.int(self.teff)) + '_' + str(np.int(100.*self.logg)) + '.dk'
                    self.openFile()
                 else:
                    print "Check Teff and logg and try again"


      def openFile(self):
                 try:
                   name = os.path.join(os.path.dirname(__file__), 'data', self.modelname)
                   self.spectrum = np.genfromtxt(open(name),skip_header=35)
                   self.base = np.arange(self.spectrum[0,0],self.spectrum[-1,0],0.23)
                 except IOError:
                   print "File not found in the target folder. Please try again"



      def __call__(self,v):
            # First we need to find the slit width in units of the model-spectrum resolution
            # Slit width in pixels is: slit_in_pix = slit_in_arcsec/instrument_resolution
            # Slit width in dispersion units is: slit_in_Ang = slit_in_pix*dispersion

                 f = np.interp(self.base,self.spectrum[:,0]*(1+v/299792.458),self.spectrum[:,1])
                 dspectrum = np.array([self.base,f]).T
                 gres = np.min(np.abs(np.diff(self.base)))
                 slitinsres = (self.slit/self.res)*self.disp/gres
                 self.fwhm_slit = np.round(slitinsres/2.)

            # Same for fwhm seeing
                 self.seeing = (self.seeing/self.res)*self.disp # in A

            # Corresponding sigma in dispersion units: 
                 self.sigma = (self.seeing/(8.*np.log(2.))**0.5)/gres

            # Now we create a gaussian convolution kernel truncated at the fwhm slit width:      

                 self.kernel = np.arange(-self.fwhm_slit/2,self.fwhm_slit/2+1.0)
                 self.con_kern = np.exp(-((self.kernel)**2)/(2*self.sigma**2))
                 self.con_kern /= self.con_kern.sum()

            # ...and convolve our model spectrum with it

                 convolved = (np.array([dspectrum[:,0],flt.convolve(dspectrum[:,1],self.con_kern)])).T

                 return convolved



def fitwdmodel(s,tgrid,ggrid,vgrid=None,degree=3,chi2range=1e99,plot=True,*args, **kwargs):
    
    t,g = np.meshgrid(tgrid,ggrid)
    t,g = t.T,g.T

    if vgrid == None:
          vgrid = np.array([0.])

#table to hold the data

    fittable = Table([t,g,np.zeros_like(t),np.zeros_like(t),np.zeros_like(t)],names=['teff','logg','chi2','chi2red','prob'])

    if plot:
        from matplotlib.pylab import subplots,close
        fig,ax = subplots(1,1)
        ax.hold(True)
        plt.ion()
        plt.show()

#loop

    for i in np.arange(len(tgrid)):

        for j in np.arange(len(ggrid)):

            m = DAModel(t[i,j],g[i,j],*args,**kwargs)
            fit,model = findvel(s,m,vgrid,plot=plot,pold=degree,sigmarange=0.5)
            fittable['chi2'][i,j] = fit.meta['bestchi2']

            if plot:

                          plt.cla()
                          ax.errorbar(s[:,0],s[:,1],yerr=s[:,2])
                          ax.plot(s[:,0],model)
                          plt.draw()


    fittable.meta['dof'] = len(s[:,0]) -2.
    fittable['chi2red'] = fittable['chi2']/float(fittable.meta['dof'])

    i=np.where(fittable['chi2'] < fittable['chi2'].min() + chi2range)
    x = fittable['teff'][i].flatten()
    tt = fittable['teff']
    y = fittable['logg'][i].flatten()
    gg = fittable['logg']
    z = ((fittable['chi2'][i]-fittable['chi2'][i].min()).flatten())/(fittable.meta['dof']*fittable['chi2'][i].min())
    m = polyfit2d(x,y,z)

    m = m*fittable.meta['dof']*fittable['chi2'][i].min()
    fittable.meta['pol2dcoeffs'] = m

    xbest = -(m[1]*m[4]-2.*m[3]*m[2])/(m[4]**2. -4.*m[6]*m[2])
    ybest = -(m[3]*m[4]-2*m[1]*m[6])/(m[4]**2. -4.*m[6]*m[2])
    sigmax = np.sqrt((-4*(m[2])/(m[4]**2. -4*m[6]*m[2])))
    sigmay = np.sqrt((-4*(m[6])/(m[4]**2. -4*m[6]*m[2])))

    fittable.meta['BestFit'] = np.array([xbest,ybest,sigmax,sigmay,2*sigmax,2*sigmay,3*sigmax,3*sigmay])

    fittable['chi2anal'] = m[0] + m[1]*gg + m[2]*gg*gg + m[3]*tt + m[4]*tt*gg + m[6]*tt*tt
    fittable['probanal'] = 0.5*np.exp(-(fittable['chi2anal'])/2.)

    return fittable



def bayes(fittable):

    fittable['prob'] = np.exp(-(fittable['chi2']-np.min(fittable['chi2']))/2.)

    fittable.meta['ncon']=fittable['prob'].sum()
    fittable['prob'] = fittable['prob']/float(fittable.meta['ncon'])

    fittable.meta['TeffProb'] = fittable['prob'].sum(1)
    fittable.meta['loggProb'] = fittable['prob'].sum(0)

    fittable.meta['TeffCDF'] = np.cumsum(fittable.meta['TeffProb'])
    fittable.meta['loggCDF'] = np.cumsum(fittable.meta['loggProb'])

#Get precentiles from inverse CDF

    tcdfi=interp1d(fittable.meta['TeffCDF'],tgrid)
    gcdfi=interp1d(fittable.meta['loggCDF'],ggrid)

    fittable.meta['TeffPercentiles'] = np.array([np.max(tcdfi(0.5)),np.max(tcdfi(0.15832)),np.min(tcdfi(1-0.15832)),np.max(tcdfi(0.02252)),np.min(tcdfi(1-0.02252)),np.max(tcdfi(0.00137)),np.min(tcdfi(1-0.00137))])

    fittable.meta['loggPercentiles'] = np.array([np.max(gcdfi(0.5)),np.max(gcdfi(0.15832)),np.min(gcdfi(1-0.15832)),np.max(gcdfi(0.02252)),np.min(gcdfi(1-0.02252)),np.max(gcdfi(0.00137)),np.min(gcdfi(1-0.00137))])
    

    return fittable
