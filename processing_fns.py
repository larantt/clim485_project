"""
processing_fns.py

Script containing all the functions needed to for the Climate 485 project.

QUESTIONS:
1. how to get uncertainty due to the clear sky transmission and could we
fold this into the brightness temperature calculation

2. how do we get delta_LWP and delta_TClear

3. how specific do we need to be for the orbit - what assumptions can we
make to get the NESR assuming we want to derive the orbit from the ideal 
frequency not the other way around.

4. How do we do the math for an eliptical orbit, where do we need to start
and how is this different from a polar orbit?
"""

## IMPORTS ##
import numpy as np
import pandas as pd
from ocean_reflectance import ocean_R

## GLOBALS ##
GAS_ABS = ['freq_GHz', 'dry_air_tau', 'wv_abs']
LIQ_CLD = ['freq_GHz', 'mass_abs']

GAS_ABS_FP = "/Users/laratobias-tarsh/Documents/wn25/clim485/project/mw_gas_absorption.txt"
LIQ_CLD_FP = "/Users/laratobias-tarsh/Documents/wn25/clim485/project/liquid_cloud_ka_0degC.txt" # assumption cloud is 0c

def read_text_file(filepath,col_names):
    """
    Reads text file for absorbption coefficient or optical depths
    
    Parameters
    ----------
    filepath : str
        path to the file to read in
    col_names : list(str)
        names of the columns in the file

    Returns
    -------
    df : pd.DataFrame
        pandas dataframe containing text file data
    """
    df = pd.read_csv(
        filepath,
        delim_whitespace=True,
        comment='#',
        header=None,
        names=col_names
    )
    return df

def total_atmospheric_trans(freq, lwp, df_liq_water, df_wv_gases, cwv = 14, theta = 0):
    """
    Function to determine the total atmospheric transmission by using 
    the given mass extinction coefficients at a given frequency. 
    Incorporates optical depth due to gases (WV and others) and liquid water at 0C.

    Parameters
    ----------
    freq : float
        Frequency being sensed in GHz

    lwp : float
        liquid water path in kg/m^2
    
    df_liq_water : pd.DataFrame
        table of mass extinction coefficients by frequency for liquid water at 0C
    
    df_wv_gases : pd.Dataframe
        table of mass extinction coefficients by frequency for water vapor and other gases 
    
    cwv : float
        the column water vapor in kg/m^2. Default is 14 for a U.S. Standard Atmosphere
    
    theta : float
        the angle off nadir that the satellite is detecting radiance from, in degrees. Defaults to 0 (nadir)

    Returns
    -------
    t_star : float
        total atmospheric transmission
    
    """
    # Cloud mass extinction coefficient
    k_liq_water = np.interp(freq, df_liq_water['freq_GHz'], df_liq_water['mass_abs'])

    # Water vapor mass extinction coefficient
    k_wv = np.interp(freq, df_wv_gases['freq_GHz'], df_wv_gases['wv_abs'])

    # Dry air optical depth
    tau_dry_air = np.interp(freq, df_wv_gases['freq_GHz'], df_wv_gases['dry_air_tau'])

    # convert degrees to radians
    theta = np.radians(theta)
    mu = np.cos(theta)

    # Calculate total atmospheric transmission (t_star)
    t_star = np.exp(-1 * (k_liq_water*lwp + k_wv*cwv + tau_dry_air)/mu)

    return t_star


def tb_polarization(lst, freq, t_star, cloud_temp = 273.15, lss = 0, theta = 0):
    '''Calculate the brightness temperatures of the lake surface for different polarizations
   
    PARAMETERS
    ----------
    sst: float
        lake surface temperature (Celsius)
    cloud_temp: float
        temperature of the cloud deck (isothermal assumption, in Kelvin)
    freq: float or numpy array
        frequency observed (GHz)
    t_star : float or numpy array
        total atmospheric transmission
    lss: float
        lake surface salinity (ppt)
    theta: float
        incidence angle observed (degrees)

    For the purposes of this project, lss and theta will always be zero.
    Returns the brightness temperature of the lake surface for each polarization'''

    # Calculate the reflectances for the polarizations
    try:
        Rv, Rh = ocean_R(lst, lss, np.array(freq), theta)
    except:
        Rv, Rh = ocean_R(lst, lss, freq, theta)

    tb_Rv = (1-Rv)*(lst+273.15)*t_star + (1-t_star)*cloud_temp
    tb_Rh = (1-Rh)*(lst+273.15)*t_star + (1-t_star)*cloud_temp

    return tb_Rv, tb_Rh

def tb_sensitivity_lst(lst, delta_lst, freq, t_star):
    '''Calculate the sensitivity of brightness temperature to salinity for each polarization
    
    PARAMETERS
    ----------
    sst: float
        lake surface temperature (Celsius)
    delta_lst : float
        small change in lst used to compute the central difference derivative
    freq: float or numpy array
        frequency observed (GHz)
    t_star : float or numpy array
        total atmospheric transmission
    
    Returns the partial derivatives of brightness temperature with respect to LST for each polarization'''

    # Calculate brightness temperatures at each frequency for positive delta_lst and negative delta_lst
    tb_Rv_posdlst, tb_Rh_posdlst = tb_polarization(lst+delta_lst, np.array(freq), t_star)
    tb_Rv_negdlst, tb_Rh_negdlst = tb_polarization(lst-delta_lst, np.array(freq), t_star)

    # Setup finite difference approach to estimate the sensitivity of the brightness temperature to lst
    dtb_dlst_Rv =  (tb_Rv_posdlst - tb_Rv_negdlst)/(2*delta_lst)
    dtb_dlst_Rh =  (tb_Rh_posdlst - tb_Rh_negdlst)/(2*delta_lst)

    return dtb_dlst_Rv, dtb_dlst_Rh

def tb_sensitivity_lwp(lst, lwp, delta_lwp, freq, df_liq_water, df_wv_gases, cwv = 14, theta = 0):
    '''Calculate the sensitivity of brightness temperature to salinity for each polarization
    
    PARAMETERS
    ----------
    sst: float
        lake surface temperature (Celsius)
    delta_lwp : float
        small change in lwp used to compute the central difference derivative
    freq: float or numpy array
        frequency observed (GHz)
    df_liq_water : pd.DataFrame
        table of mass extinction coefficients by frequency for liquid water at 0C
    df_wv_gases : pd.Dataframe
        table of mass extinction coefficients by frequency for water vapor and other gases 
    cwv : float
        the column water vapor in kg/m^2. Default is 14 for a U.S. Standard Atmosphere
    theta : float
        the angle off nadir that the satellite is detecting radiance from, in degrees. Defaults to 0 (nadir)
    
    Returns the partial derivatives of brightness temperature with respect to LST for each polarization'''

    # Calculate transmission for positive and negative delta_lwp
    t_star_posdlwp = total_atmospheric_trans(freq, lwp+delta_lwp, df_liq_water, df_wv_gases, cwv = 14, theta = 0)
    t_star_negdlwp = total_atmospheric_trans(freq, lwp-delta_lwp, df_liq_water, df_wv_gases, cwv = 14, theta = 0)

    # Calculate brightness temperatures at each frequency for positive and negative delta_lwp
    tb_Rv_posdlwp, tb_Rh_posdlwp = tb_polarization(lst, np.array(freq), np.array(t_star_posdlwp))
    tb_Rv_negdlwp, tb_Rh_negdlwp = tb_polarization(lst, np.array(freq), np.array(t_star_negdlwp))

    # Setup finite difference approach to estimate the sensitivity of the brightness temperature to lwp
    dtb_dlwp_Rv =  (tb_Rv_posdlwp - tb_Rv_negdlwp)/(2*delta_lwp)
    dtb_dlwp_Rh =  (tb_Rh_posdlwp - tb_Rh_negdlwp)/(2*delta_lwp)

    return dtb_dlwp_Rv, dtb_dlwp_Rh


def total_uncertainty(NEDT,dTB_dLST,dTB_dLWP,delta_LWP):
    """
    Function calculates the total additive uncertainty for the 
    brightness temperature

    Parameters
    ----------
    NEDT : float
        NEDT for the sensor
    dTB_dLST : float
        sensitivity of brightness temperature to surface temperature
    dTB_dLWP : float
        sensitivity of brightness temperature to liquid water path
    delta_LWP : float

    Returns
    -------
    sigma_LST : float
        total uncertainty on the lake surface temperature measurement
    """
    term_1 = (NEDT**2) * (dTB_dLST**-2)
    term_2 = (delta_LWP**2) * (dTB_dLWP**2) * (dTB_dLST**-2)
    sigma_LST = np.sqrt(term_1 + term_2)

    return sigma_LST


def fwd_model_rad(L_TOA,)