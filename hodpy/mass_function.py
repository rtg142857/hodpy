#! /usr/bin/env python
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit
import h5py
import yaml

from hodpy.power_spectrum import PowerSpectrum
from hodpy.cosmology import CosmologyMXXL, CosmologyAbacus, CosmologyFlamingo
from hodpy import lookup


class MassFunction(object):
    """
    Class for fitting a Sheth-Tormen-like mass function to measurements from a simulation snapshot

    Args:
        cosmology: hodpy.Cosmology object, in the cosmology of the simulation
        redshift: redshift of the simulation snapshot
        [fit_params]: if provided, sets the best fit parameters to these values. 
                      fit_params is an array of [dc, A, a, p]
        [measured_mass_function]: mass function measurements to fit to, if provided.
                      measured_mass_function is an array of [log mass bins, mass function]
    """
    def __init__(self, cosmology, redshift, fit_params=None, measured_mass_function=None):
        
        self.cosmology = cosmology
        self.power_spectrum = PowerSpectrum(self.cosmology)
        self.redshift = redshift
        
        if not fit_params is None:
            self.dc, self.A, self.a, self.p = fit_params
            
        if not measured_mass_function is None:
            self.__mass_bins = measured_mass_function[0]
            self.__mass_func = measured_mass_function[1]
            
     
    def __func(self, sigma, dc, A, a, p):
        # Sheth-Tormen mass function
        mf = A * np.sqrt(2*a/np.pi)
        mf *= 1 + (sigma**2 / (a * dc**2))**p
        mf *= dc / sigma
        mf *= np.exp(-a * dc**2 / (2*sigma**2))

        return np.log10(mf)
        
    
    def get_fit(self):
        """
        Fits the Sheth-Tormen mass function to the measured mass function, returning the
        best fit parameters
        
        Returns:
            an array of [dc, A, a, p]
        """
        sigma = self.power_spectrum.sigma(10**self.__mass_bins, self.redshift)
        mf = self.__mass_func / self.power_spectrum.cosmo.mean_density(0) * 10**self.__mass_bins
        
        popt, pcov = curve_fit(self.__func, sigma, np.log10(mf), p0=[1,0.1,1.5,-0.5])
        
        self.update_params(popt)
        
        return popt

    
    def update_params(self, fit_params):
        '''
        Update the values of the best fit params
        
        Args:
            fit_params: an array of [dc, A, a, p]
        '''     
        self.dc, self.A, self.a, self.p = fit_params
        
    
    def mass_function(self, log_mass, redshift=None):
        '''
        Returns the halo mass function as a function of mass and redshift
        (where f is defined as Eq. 4 of Jenkins 2000)

        Args:
            log_mass: array of log_10 halo mass, where halo mass is in units Msun/h
        Returns:
            array of halo mass function
        '''        
        
        sigma = self.power_spectrum.sigma(10**log_mass, self.redshift)

        return 10**self.__func(sigma, self.dc, self.A, self.a, self.p)

    
    def number_density(self, log_mass, redshift=None):
        '''
        Returns the number density of haloes as a function of mass and redshift

        Args:
            log_mass: array of log_10 halo mass, where halo mass is in units Msun/h
        Returns:
            array of halo number density in units (Mpc/h)^-3
        '''  
        mf = self.mass_function(log_mass)

        return mf * self.power_spectrum.cosmo.mean_density(0) / 10**log_mass
        

        

class MassFunctionMXXL(object):
    """
    Class containing the fits to the MXXL halo mass function

    Args:
        mf_fits_file: Tabulated file of the best fit mass function parameters
    """
    def __init__(self, mf_fits_file=lookup.mxxl_mass_function):
        
        #self.power_spectrum = power_spectrum
        self.cosmology = CosmologyMXXL()
        self.power_spectrum = PowerSpectrum(self.cosmology)
        
        # read in MXXL mass function fit parameters
        snap, redshift, A, a, p = \
                   np.loadtxt(mf_fits_file, skiprows=1, unpack=True)
        
        # interpolate parameters
        self._A = RegularGridInterpolator((redshift,), A, bounds_error=False, 
                                          fill_value=None)

        self._a = RegularGridInterpolator((redshift,), a, bounds_error=False, 
                                          fill_value=None)

        self._p = RegularGridInterpolator((redshift,), p, bounds_error=False, 
                                          fill_value=None)

    def A(self, redshift):
        return self._A(redshift)

    def a(self, redshift):
        return self._a(redshift)

    def p(self, redshift):
        return self._p(redshift)

    def mass_function(self, log_mass, redshift):
        '''
        Returns the halo mass function as a function of mass and redshift
        (where f is defined as Eq. 4 of Jenkins 2000)

        Args:
            log_mass: array of log_10 halo mass, where halo mass is in units Msun/h
            redshift: array of redshift
        Returns:
            array of halo mass function
        '''        
        
        sigma = self.power_spectrum.sigma(10**log_mass, redshift)
        
        # sigma(z) evolution was incorrect when doing these fits
        # apply correction to this
        sigma = sigma * (self.power_spectrum.delta_c(redshift)/self.power_spectrum.delta_c(0))**2

        dc=1
        A = self.A(redshift)
        a = self.a(redshift)
        p = self.p(redshift)
        
        mf = A * np.sqrt(2*a/np.pi)
        mf *= 1 + (sigma**2 / (a * dc**2))**p
        mf *= dc / sigma
        mf *= np.exp(-a * dc**2 / (2*sigma**2))

        return mf

    def number_density(self, log_mass, redshift):
        '''
        Returns the number density of haloes as a function of mass and redshift

        Args:
            log_mass: array of log_10 halo mass, where halo mass is in units 
                      Msun/h
            redshift: array of redshift
        Returns:
            array of halo number density in units (Mpc/h)^-3
        '''  
        mf = self.mass_function(log_mass, redshift)

        return mf * self.power_spectrum.cosmo.mean_density(0) / 10**log_mass
        


class MassFunctionAbacus(object):
    """
    Class containing the fits to the AbacusSummit halo mass function

    Args:
        cosmo:   AbacusSummit cosmology number
        [phase]: AbacusSummit phase. Default is 0.
        [mf_file]: File containing the AbacusSummit mass function measurements.
                   If none is provided, will read the file specified in 
                   lookup.abacus_mass_function. Default value is None
    """
    def __init__(self, cosmo, phase=0, mf_file=None):
        
        self.cosmology = CosmologyAbacus(cosmo)
        
        if mf_file is None:
            # use default file specified in lookup.py
            mf_file = lookup.abacus_mass_function.format(cosmo,phase)
        
        self.__logM, self.__z, self.__logn = self.__read_mf_file(mf_file, num_snap=15)
        
        self.__mass_function_interpolator = RegularGridInterpolator((self.__logM,self.__z),
                            self.__logn, bounds_error=False, fill_value=None)
        
        
    def __read_mf_file(self, mf_file, num_snap=15):
        # read the file of mass functions measured from the simulation
        f = h5py.File(mf_file, 'r')
        logM = f['0/log_mass'][...]
        zsnap = np.zeros(num_snap)
        logn = np.zeros((len(logM),num_snap))
        
        for i in range(num_snap):
            zsnap[i] = f['%i/z'%i][...][0]
            logn[:,i] = f['%i/log_n'%i][...]
        f.close()
        
        return logM, zsnap, logn
        

    def mass_function(self, log_mass, redshift):
        '''
        Returns the halo mass function as a function of mass and redshift
        (where f is defined as Eq. 4 of Jenkins 2000)

        Args:
            log_mass: array of log_10 halo mass, where halo mass is in units Msun/h
            redshift: array of redshift
        Returns:
            array of halo mass function
        '''     
        n = self.number_density(log_mass, redshift)
        
        return mf * 10**log_mass / self.cosmology.mean_density(0)

    
    def number_density(self, log_mass, redshift):
        '''
        Returns the number density of haloes as a function of mass and redshift

        Args:
            log_mass: array of log_10 halo mass, where halo mass is in units 
                      Msun/h
            redshift: array of redshift
        Returns:
            array of halo number density in units (Mpc/h)^-3
        '''  
        
        return 10**self.__mass_function_interpolator((log_mass,redshift))
    
class MassFunctionFlamingo(object):
    """
    Class containing the fits to the Flamingo halo mass function

    Args:
        path_config_fiename:   yml file containing config stuff
    Old args:
        cosmo:
        [phase]: AbacusSummit phase. Default is 0.
        [mf_file]: File containing the AbacusSummit mass function measurements.
                   If none is provided, will read the file specified in 
                   lookup.abacus_mass_function. Default value is None
    """
    def __init__(self, path_config_filename, mf_file=None):
        with open(path_config_filename, "r") as file:
            path_config = yaml.safe_load(file)
        label = path_config["Labels"]["sim_label"]
        
        self.cosmology = CosmologyFlamingo(path_config_filename)
        
        if mf_file is None:
            # use default file specified in lookup.py
            mf_file = lookup.flamingo_mass_function.format(label)
        
        num_redshifts = path_config["Params"]["num_redshifts"]
        self.__logM, self.__z, self.__logn = self.__read_mf_file(mf_file, num_snap=num_redshifts)
        
        self.__mass_function_interpolator = RegularGridInterpolator((self.__logM,self.__z),
                            self.__logn, bounds_error=False, fill_value=None)
        
        
    def __read_mf_file(self, mf_file, num_snap=15):
        # read the file of mass functions measured from the simulation
        f = h5py.File(mf_file, 'r')
        logM = f['0/log_mass'][...]
        zsnap = np.zeros(num_snap)
        logn = np.zeros((len(logM),num_snap))
        
        for i in range(num_snap):
            zsnap[i] = f['%i/z'%i][...][0]
            logn[:,i] = f['%i/log_n'%i][...]
        f.close()
        
        return logM, zsnap, logn
        

    def mass_function(self, log_mass, redshift):
        '''
        Returns the halo mass function as a function of mass and redshift
        (where f is defined as Eq. 4 of Jenkins 2000)

        Args:
            log_mass: array of log_10 halo mass, where halo mass is in units Msun/h
            redshift: array of redshift
        Returns:
            array of halo mass function
        '''     
        n = self.number_density(log_mass, redshift)
        
        return mf * 10**log_mass / self.cosmology.mean_density(0)

    
    def number_density(self, log_mass, redshift):
        '''
        Returns the number density of haloes as a function of mass and redshift

        Args:
            log_mass: array of log_10 halo mass, where halo mass is in units 
                      Msun/h
            redshift: array of redshift
        Returns:
            array of halo number density in units (Mpc/h)^-3
        '''  
        
        return 10**self.__mass_function_interpolator((log_mass,redshift))

