import numpy as np
from astropy.table import Table

###Functions for making publication-quality plots of fitted spectrum ####
def plot_da_fit(chi2,mc=None,cl=True):
    from matplotlib import pyplot as plt
    import matplotlib.colors
    from matplotlib import cm
    from matplotlib.ticker import NullFormatter

    nullfmt = NullFormatter()

    left,width = 0.1,0.65
    bottom,height = 0.1,0.65
    bottom_h = left_h = left+width+0.02
    left = 0.03
    left_h = left+0.02+0.2


    rect_scatter = [left_h,bottom,width,height]
    rect_histx = [left_h,bottom_h,width,0.2]
    rect_histy = [left,bottom, 0.2,height]


    plt.figure(1, figsize=(8,8))

    chi = plt.axes(rect_scatter)
    teff = plt.axes(rect_histx)
    logg = plt.axes(rect_histy)

    #Format Labels and Ticks 
    teff.xaxis.set_major_formatter(nullfmt)
    teff.yaxis.set_major_formatter(nullfmt)
    logg.xaxis.set_major_formatter(nullfmt)
    logg.yaxis.set_major_formatter(nullfmt)
    chi.yaxis.set_ticks_position("right")
    chi.yaxis.set_label_position("right")

    
    #Range for the plot ~5 times the 1sigma error 
    tmin,tmax = chi2['teff'][:,0][0],chi2['teff'][:,0][-1]
    gmin,gmax = chi2['logg'][0,:][0],chi2['logg'][0,:][-1]
#    tmin,tmax = chi2.meta['TeffPercentiles'][5]-200., chi2.meta['TeffPercentiles'][6]+200.
#    gmin,gmax = chi2.meta['loggPercentiles'][5]-0.3,chi2.meta['loggPercentiles'][6]+0.3


    #Plot Projected Probabilities
    teff.plot(chi2['teff'][:,0],chi2.meta['TeffProb'])
    teff.scatter(chi2['teff'][:,0],chi2.meta['TeffProb'])
    teff.axis([tmin,tmax,0.0,1.0])

    logg.plot(chi2.meta['loggProb'],chi2['logg'][0,:])
    logg.scatter(chi2.meta['loggProb'],chi2['logg'][0,:])
    logg.axis([0.0,1.0,gmin,gmax])

    #Plot x2 map

    chi.contour(chi2['teff'],chi2['logg'],(chi2['chi2']-np.min(chi2['chi2'])),[2.3,6.2,11.8])
    if cl:
        if mc == None:
            chi.contourf(chi2['teff'],chi2['logg'],chi2['chi2'],20,cmap='winter')
        else: 
        #    chi.contourf(chi2['teff'],chi2['logg'],(chi2['chi2']-np.min(chi2['chi2'])),20,cmap='winter')
            chi.hist2d(mc['teff'],mc['logg'],20)
    else:
        if mc != None:
            chi.scatter(mc['teff'],mc['logg'],s=0.001)
    chi.axis([tmin,tmax,gmin,gmax])


    # Labels
    chi.set_xlabel('Temperature (K)', family='serif',size='x-large')
    chi.set_ylabel('Surface gravity ($\log g$)', family='serif',size='x-large')


    plt.show()
