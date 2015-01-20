import numpy as np
import itertools
from numpy import interp as interp

def sfile(file):
      spectrum = np.genfromtxt(open(file))
      return spectrum


def bestfit(x,y,sigmarange=None):
       if sigmarange is None:
           sigmarange = 1e+10
       ibest = y.argmin()
       i = np.where((y < y[ibest]*(1.+sigmarange)))
       sol = polyfit(x[i],y[i],2)
       ybest = -sol[1]/2./sol[2]
       yerr = np.sqrt(polyval(ybest,sol)/(len(y)-2.)/sol[2])
       return ybest, yerr


def doppler(v):
    return 1. + v/299792.458


def polyfit2d(x, y, z):
    V=np.polynomial.polynomial.polyvander2d(x,y,[2,2])
    V[:,5] = 0
    V[:,7] = 0
    V[:,8] = 0
    m,_,_,_ = np.linalg.lstsq(V, z)
    return m


def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


def getnearest(a, MinClip):
      return np.round(a / MinClip) * MinClip


def interpll(x,xi,yi,base=2.718281828459045):
      lx,lxi,lyi = np.log(x)/np.log(base),np.log(xi)/np.log(base),np.log(yi)/np.log(base)
      return base**interp(lx,lxi,lyi)


def polygrid(a,bins):
      # 
      # http://stackoverflow.com/questions/13542402/how-to-create-multidimensional-array-with-numpy-mgrid
      D = len(a)
      amin = a-a
      amax = 3*a + a/bins
      n_bins =  bins*np.ones(D)
      k = np.int(np.prod(n_bins))
      bounds = np.array([amin,amax]).T
      agrid = (np.mgrid[[slice(row[0], row[1], n*1j) for row, n in zip(bounds, n_bins)]]).reshape(D,k)
      return agrid, k



      
      
