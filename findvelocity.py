from astropy.table import Table
from numpy.polynomial.polynomial import polyfit 
from numpy.polynomial.polynomial import polyval
import numpy as np
from matplotlib import pyplot as plt


class ErrorContinue(Exception):
       pass

def findvel(s,m,vgrid,pold=3,sigma=0.5,plot=False):
#todo. cover case len(vgrid)=1=v
       chi2 = Table([vgrid,np.zeros_like(vgrid),np.zeros_like(vgrid)],names=['v','chi2','chi2red'])
       chi2['chi2'] = np.array([Fit(s,m(v),degree=pold)()[0] for v in vgrid])
       chi2.meta['ndof'] = len(s[:,1]) - 2. - pold
       chi2['chi2red'] = chi2['chi2']/chi2.meta['ndof']



       if len(vgrid) < 3:
            ibest = chi2['chi2'].argmin()
            vbest, bestchi2,bestchi2red = chi2[ibest]
            chi2.meta['vbest'] = vbest
            chi2.meta['verr'] = 0.
            chi2.meta['bestchi2'] = bestchi2
            sol = 0
       else:
            chi2.meta['vbest'], chi2.meta['verr'], chi2.meta['bestchi2'], sol = minchi2(chi2,sigmarange=sigma)


            if plot:

                plt.plot(chi2['v'],polyval(chi2['v'],sol),color='r')
                plt.scatter(chi2['v'], chi2['chi2'])
                plt.show()


       return chi2, Fit(s,m(chi2.meta['vbest']))()[1], sol



def minchi2(chi2,sigmarange=None):
       if sigmarange is None:
           sigmarange = 0.5
       ibest = chi2['chi2'].argmin()
       i = np.where((chi2['chi2'] < chi2['chi2'][ibest] + (sigmarange)))
       sol = polyfit(chi2['v'][i],chi2['chi2'][i],2)

       vbest = -sol[1]/2./sol[2]
       verr = np.sqrt(polyval(vbest,sol)/chi2.meta['ndof']/sol[2])
       bestchi2 = polyval(vbest,sol)

       return vbest, verr, bestchi2, sol


class Fit(object):
         def __init__(self,spectrum,model,degree=3):
                 self.s = spectrum
                 self.m = model
                 self.degree = degree
                 self.fde = self.s[:,1]/self.s[:,2]
                 self.Vp = np.polynomial.polynomial.polyvander(self.s[:,0]/self.s[:,0].mean()-1., self.degree)
                 self.rcond = len(self.s[:,1])*np.finfo(self.s[:,1].dtype).eps
                 self.sol = None
                 self.interpolateToObservedWavelengths()

         def __call__(self):
                 V = self.Vp*(self.mf/self.s[:,2])[:,np.newaxis] 
                 sq =np.sqrt((V*V).sum(0))
                 self.sol, res, rank, s = np.linalg.lstsq(V/sq, self.fde, self.rcond)
                 self.sol = (self.sol.T/sq).T
                 fit = np.dot(V, self.sol)*self.s[:,2]
                 chi2 = np.sum(((self.s[:,1]-fit)/self.s[:,2])**2.)
                 return chi2,fit

         def interpolateToObservedWavelengths(self):
                 try:
                   self.mf = np.interp(self.s[:,0],self.m[:,0],self.m[:,1]) 
                 except ValueError:
                   print "Something is wrong with the observed wavelengths"


