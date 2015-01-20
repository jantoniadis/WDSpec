import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from numpy import interp as interp
from astropy.table import Table
from scipy.interpolate import griddata

from .genutils import getnearest

class MCmyDAFit(object):
    def __init__(self,chi2,full=False):
        self.chi2 = chi2
        self.full = full
    

    def __call__(self,n):
        if n == 1:
            print "Can't MC with one point"
            return 0.
        else: 
            if self.full:
                mctable=self.mc2d(n)
                return mctable
            else:
                mctable=self.mc1d(n)
                return mctable


    def mc1d(self,n,tstep=25.,gstep=0.01):
        tlow,tmax,gmin,gmax=self.chi2['teff'][:,0][0],self.chi2['teff'][:,0][-1],self.chi2['logg'][0,:][0], self.chi2['logg'][0,:][-1]

        #Resample the marginalized PDFs 
        t = np.arange(tlow,tmax+10e-09,tstep)
        g = np.arange(gmin,gmax+10e-09,gstep)
        tpdf = interp(t,self.chi2['teff'][:,0],self.chi2.meta['TeffProb'])
        gpdf = interp(g,self.chi2['logg'][0,:],self.chi2.meta['loggProb'])
        
        # Calculate CDFs and re-normalize 
        tcdf = np.cumsum(tpdf)
        tcdf = tcdf/tcdf.max()
        gcdf = np.cumsum(gpdf)
        gcdf = gcdf/gcdf.max()

        #Calculate inverse cdfs
        itcdf = interp1d(tcdf,t)
        igcdf = interp1d(gcdf,g)
        mctable=Table([itcdf(np.random.random(n)),igcdf(np.random.random(n))],names=['teff','logg'])
        return mctable

#Experimental. Doesn't work well for sparcely sampled data.


    def mc2d(self,n,tstep=10.,gstep=0.01):
        tlow,tmax,gmin,gmax=self.chi2['teff'][:,0][0],self.chi2['teff'][:,0][-1],self.chi2['logg'][0,:][0], self.chi2['logg'][0,:][-1]


        seed = np.random.random(n)

        #Resample the marginalized PDF 
        t = np.arange(tlow,tmax+1e-09,tstep)
        g = np.arange(gmin,gmax+1e-09,gstep)
        tt,gg = np.meshgrid(t,g)
        tt,gg = tt.T,gg.T

        m = self.chi2.meta['pol2dcoeffs']
        chi = m[0] + m[1]*gg + m[2]*gg*gg + m[3]*tt + m[4]*tt*gg + m[6]*tt*tt
        chi = chi - chi.min()

        grid = np.exp(-chi/2.)
        grid = grid/grid.sum()
        tpdf = grid.sum(1)

        tcdf = np.cumsum(tpdf)
        itcdf = interp1d(tcdf,t)      

        tprob = itcdf(seed)

        # For each element of the array get nearest coordinate on the original grid
        tpos = (getnearest(tprob,tstep) - tlow)/tstep
        
        mctable=Table([tprob,np.zeros_like(tprob)],names=['teff','logg'])

        self.chi2.meta['resampled_tgrid'] = tt
        self.chi2.meta['resampled_ggrid'] = gg
        self.chi2.meta['resampled_chi2'] = chi
        self.chi2.meta['resampled_prob'] = np.exp(-chi/2.)/np.exp(-chi/2.).sum()

        gcdf = np.cumsum(grid,axis=1)

        gcdf = (gcdf.T/gcdf.max(1)).T

        #Now loop (this seems to be the fastest way of doing it)
        mctable['logg']=np.array([interp(np.random.random(),gcdf[i,:],g) for i in tpos])

        return mctable

        

        

        
