import numpy as np
import numexpr as ne 

def ml17_to_3d_teff(teff,logg):
    a1 = -1.0461690E-03
    a2 = -2.6846737E-01
    a3 = 3.0654611E-01
    a4 = 1.8025848E+00
    a5 = 1.5006909E-01
    a6 = 1.0125295E-01
    a7 = -5.2933335E-02
    a8 = -1.3414353E-01

    teff0 = ne.evaluate("(teff-10000.0)/1000.00")
    logg0 = ne.evaluate("(logg-8.00000)")
    shift = ne.evaluate("a1+(a2+a3*teff0+(a6+a7*teff0+a8*logg0)*logg0)*exp(-a3*((teff0-a4)**2))")

    return ne.evaluate("teff + shift*1000.00")

def ml17_to_3d_logg(teff,logg):
    a0  = 1.1922481E-03
    a1  = -2.7230889E-01
    a2  = -6.7437328E-02
    a3  = -8.7753624E-01
    a4  = 1.4936511E-01
    a5  = -1.9749393E-01
    a6  = 4.1687626E-01
    a7  = 3.8195432E-01
    a8  = -1.4141054E-01
    a9  = -2.9439950E-02
    a10 = 1.1908339E-01

    teff0 = ne.evaluate("(teff-10000.0)/1000.00")
    logg0 = ne.evaluate("logg-8.0")

    shift= ne.evaluate("a0 +a1*exp(-(a2 +(a4+a6*exp(-a7*((teff-a8)**2.)))*teff0+(a5+a9*teff0+a10*logg0)*logg0)**2.*((teff0-a3)**2))")
    
    return ne.evaluate("logg+shift")

def ml18_to_3d_teff(teff,logg):
    a1 = 1.0947335E-03
    a2 = -1.8716231E-01
    a3 = 1.9350009E-02
    a4 = 6.4821613E-01
    a5 = -2.2863187E-01
    a6 = 5.8699232E-01
    a7 = -1.0729871E-01
    a8 = 1.1009070E-01
    teff0 = ne.evaluate("(teff-10000.0)/1000.0")
    logg0 = ne.evaluate("logg-8.0")

    shift = ne.evaluate("a1 +(a2+a7*teff0+a8*logg0)*exp(-a3+a5*teff0+a6*logg0)**2. *((teff0-a4)**2.)")

    return ne.evaluate("teff + shift*1000.")

def ml18_to_3d_logg(teff,logg):
    a1  = 7.5209868E-04
    a2  = -9.2086619E-01
    a3  = 3.1253746E-01
    a4  = -1.0348176E+01
    a5  = 6.5854716E-01
    a6  = 4.2849862E-01
    a7  = -8.8982873E-02
    a8  = 1.0199718E+01
    a9  = 4.9277883E-02
    a10 = -8.6543477E-01
    a11 = 3.6232756E-03
    a12 = -5.8729354E-02

    teff0 = ne.evaluate("(teff-10000.0)/1000.0")
    logg0 = ne.evaluate("logg-8.0")

    shift=ne.evaluate("(a1+a5*exp(-a6*((teff0-a7)**2)))+a2*exp(-a3*((teff0-(a4+a8*exp(-(a9+a11*teff0+a12*logg0)**2. *((teff0-a10)**2))))**2))")

    return ne.evaluaate("logg+shift")
