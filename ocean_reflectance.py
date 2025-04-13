import numpy as np

def ocean_dielectric(sst, sss, f):
    """
    ocean_dieletric

    written by Chris Ruf, CLaSP Dept., University of Michigan, cruf@umich.edu
    
    translated into python by Aronne Merrelli (merrelli@umich.edu)
    
    ********************************************************************************
     OCEAN_DIELECTRIC
       Determines the complex relative dielectric constant of the ocean surface 
         as a function of its temperature and salinity and the measurement frequency
       The outputs are zeps = eps_real - j*eps_imag, unitless
    ********************************************************************************
    
    USAGE: zeps = ocean_dielectric(sst, sss, f)
    INPUTS:
        sst - sea surface temperature [Celsius]
        sss - sea surface salinity [ppt]
        f - Frequency being observed [GHz]
        the inputs sst, sss, f can be scalars or numpy arrays.
    
    OUTPUTS:
        zeps - complex relative dielectric constant [unitless]

        output will be scalar or the broadcasted shape if numpy arrays
        were used as inputs.
    ********************************************************************************

    """

    #fresh water permittivity:
    eps = (-82.168 * sst + 37088.6) / (sst + 421.854)
    #fresh water 1st debye relaxation time [ns]:
    tau1pi = (0.7246 * sst + 255.04) / ((49.25 + sst) * (45.0 + sst))

    #conductivity:
    if (sss > 0):
        p = (0.014409 * sss + 5.45216) * sss + 37.5109
        q = (sss + 182.283) * sss + 1004.75
        r15 = sss * p / q
        sig35 = (((4.3047E-09 * sst - 2.991E-06) * sst + 4.738817E-04) * sst + 0.08607) * sst + 2.903602
        p = (-0.099486 * sss + 3.2841) * sss + 6.9431
        q = (sss + 69.024) * sss + 84.85
        rt = 1+(sst - 15)*p/(((0.00198 * sss - 0.2276) * sss + 49.843 + sst)*q)
        #conductivity at actual sst&sss
        sigma = sig35 * rt * r15
        sig = 17.97510357 * sigma / f
        a = (0.02012 + 0.002781 * sss) * (34.16 + sst)/((5.074 + sss)*(26.2+sst))
        #correct fresh water permittivity for salinity:
        eps = eps * (1 - sss * a)
        b = (0.04034 + 0.002703 * sss) / (9.001 + sss) - (0.00181 * sst + 0.00562) * sst / ((sst + 7.35) * sst + 127.0)
        #correct fresh water relaxation time for salinity [ns]
        tau1pi = tau1pi * (1.0 - sss * b)
    else:
        #No conductivity if no salinity:
        sig = 0

    #Change units:
    eps1 = 0.0787 * eps

    #High freq (optical) limit to dielectric constant of water:
    epsinf = 0.0186 * sst + 4.05;

    #2nd debye relaxation time [ns] (assumes no salinity dependence):
    tau2pi = 0.00628;

    x1 = f * tau1pi
    x2 = f * tau2pi
    term1 = (eps - eps1) / (x1 * x1 + 1)
    term2 = (eps1 - epsinf) / (x2 * x2 + 1)

    #real part of dielectric constant
    epsreal = term1 + term2 + epsinf

    #absolute value of imaginary part of dielectric constant:
    epsimag = term1 * x1 + term2 * x2 + sig;

    #complex relative dielectric constant is epsreal - j*epsimag
    zeps = epsreal - 1j*epsimag;

    return zeps


def ocean_R(sst, sss, f, theta):
    """

    by Chris Ruf, CLaSP Dept., University of Michigan, cruf@umich.edu

    translated into python by Aronne Merrelli (merrelli@umich.edu)

    *******************************************************************************
    OCEAN_R
	Determines the power reflection coefficient of the ocean surface as a function of 
    its temperature and salinity
    and as a function of measurement frequency, incidence angle and polarization
    The outputs are v- and h-pol reflectivity, unitless
    *********************************************************************************/
    USAGE: [Rv, Rh] = ocean_R(sst,sss,f,theta,pol)

    INPUTS:
        sst - sea surface temperature [Celsius]
        sss - sea surface salinity [ppt]
        f - Frequency being observed [GHz]
        theta - incidence angle being observed [deg]
    OUTPUTS:
        Rv - vertical polarization ocean surface reflectivity [unitless, 0 <= emis <= 1]
        Rh - horizontal polarization ocean surface reflectivity [unitless, 0 <= emis <= 1]

        output will be scalar or the broadcasted shape if numpy arrays
        were used as inputs.
    
    *******************************************************************************

    """

    #get complex relative dielectric constant from function call
    zeps = ocean_dielectric(sst, sss, f)

    #Compute electric field Fresnel reflection coefficients, reh and rev:
    zsin = np.sin(np.deg2rad(theta))
    zcos = np.cos(np.deg2rad(theta))
    ztmp = np.sqrt(zeps - zsin**2)
    reh = (zcos - ztmp) / (zcos + ztmp)
    rev = (zeps*zcos - ztmp) / (zeps*zcos + ztmp)

    #H-pol reflectivity:
    Rh = (reh * reh.conj()).real

    #V-pol reflectivity:
    Rv = (rev * rev.conj()).real

    return Rv, Rh
