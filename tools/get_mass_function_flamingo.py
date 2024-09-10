#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.interpolate import interp1d
import h5py
import yaml
import os

import sys
sys.path.append('..')
from hodpy.mass_function import MassFunction
from hodpy.cosmology import CosmologyFlamingo
from hodpy.power_spectrum import PowerSpectrum

#from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
#from abacusnbody.data.read_abacus import read_asdf


def measure_mass_function_box(path_config_filename, snapshot, bin_size=0.01):
    """
    Measure the halo mass function of the cubic box
    Args:
        path_config_filename: yml file containing paths to simulation
        snapshot: soap file for a given redshift
        #snapshot_folder: directory name within the soap_meta_path containing the soap files for a given redshift
        [bin_size]: 0.01 by default
    Returns:
        Array of halo mass bin centres, in log10(mass)
        Array of log(n)
    """
    with open(path_config_filename, "r") as file:
        path_config = yaml.safe_load(file)
    soap_meta_path = path_config["Paths"]["soap_meta_path"]
    #soap_path = soap_meta_path + snapshot_folder
    soap_path = soap_meta_path + snapshot
    L = path_config["Params"]["L"]

    # initialize Flamingo cosmology
    #cosmology = CosmologyFlamingo(path_config_filename)
    
    # input path of the simulation
    #abacus_path = "/global/cfs/cdirs/desi/cosmosim/Abacus/"
    #mock = "AbacusSummit_%s_c%03d_ph%03d"%(simulation, cosmo, phase)
    #input_file = abacus_path+mock+"/halos/z%.3f/halo_info/halo_info_%03d.asdf"

    # location of the snapshots
    #soap_files_list = os.listdir(soap_path)
    #log_mass = [None]*len(soap_files_list)
    # for file_name in soap_files_list:
    #     file_number = file_name.split(".")[1]
    #     print("Reading files ({}%)".format(int(100*file_number/len(soap_files_list))), end='\r')
    #     input_file = soap_path + file_name

    #     halo_cat = h5py.File(input_file, "r")
    #     halo_mass = np.array(halo_cat["SO"]["200_mean"]["DarkMatterMass"])
    #     log_mass[file_number] = np.log10(halo_mass[halo_mass>0])
        
        #print(file_number, len(log_mass[file_number]))

    print("Reading halo file...")
    halo_cat = h5py.File(soap_path, "r")
    halo_mass = np.array(halo_cat["SO"]["200_mean"]["DarkMatterMass"])*1e10
    log_mass = np.log10(halo_mass[halo_mass>0])

    # log_mass = [None]*Nfiles
    
    # for file_number in range(Nfiles):
    #     print("Reading files ({}%)".format(int(100*file_number/Nfiles)), end='\r')
        
    #     halo_cat = CompaSOHaloCatalog(input_file%(redshift, file_number), 
    #                                 cleaned=True, fields=['N'])
    #     m_par = halo_cat.header["ParticleMassHMsun"]
    #     halo_mass = np.array(halo_cat.halos["N"])*m_par
    #     log_mass[file_number] = np.log10(halo_mass[halo_mass>0])

    print("Reading files (100%)")
    print("Getting mass function")
    
    #log_mass = np.concatenate(log_mass)

    # get number densities in mass bins  
    mass_bins = np.arange(10,16,bin_size)
    mass_binc = mass_bins[:-1]+bin_size/2.
    hist, bins = np.histogram(log_mass, bins=mass_bins)
    n_halo = hist/bin_size/L**3

    # remove bins with zero haloes
    keep = n_halo > 0
    return mass_binc[keep], n_halo[keep]


def get_mass_functions(path_config_filename, mass_function_file, snapshot_redshifts=None):
    """
    Measure the mass function at each snapshot, and smooth to remove noise
    """
    # if no redshifts provided, use all the snapshots with z < 1.0
    #if snapshot_redshifts is None:
    #    snapshot_redshifts = 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, \
    #                         0.450, 0.500, 0.575, 0.650, 0.725, 0.800, 0.875, 0.950
    with open(path_config_filename, "r") as file:
        path_config = yaml.safe_load(file)
    soap_meta_path = path_config["Paths"]["soap_meta_path"]
    output_list_path = path_config["Paths"]["output_list_path"] # file with list of redshifts

    snapshot_list = os.listdir(soap_meta_path)
    snapshot_list = [filename for filename in snapshot_list if filename.endswith(".hdf5")]
    snapshot_redshifts = np.genfromtxt(output_list_path, comments="#")

    if len(snapshot_list) != len(snapshot_redshifts):
        raise Exception("Snapshot folder has a different number of redshifts ("+str(len(snapshot_list))+") than output_list.txt suggests ("+str(len(snapshot_redshifts))+")")
    #for snapshot_redshift in redshifts:
    for i in range(len(snapshot_list)):
        snapshot_redshift = snapshot_redshifts[i]
        if snapshot_redshift > 7.0:
            continue

        print("z = %.3f"%snapshot_redshift)
        #logM, n = measure_mass_function_box(path_config_filename, snapshot_folder=snapshot_list[i], bin_size=0.01)
        logM, n = measure_mass_function_box(path_config_filename, snapshot=snapshot_list[i], bin_size=0.01)

        cosmology = CosmologyFlamingo(path_config_filename)

        mf = MassFunction(cosmology=cosmology, redshift=snapshot_redshift, 
                          measured_mass_function=[logM, n])

        # get fit to mass function
        fit_params = mf.get_fit()
        print("Fit parameters: "+str(fit_params))

        # take ratio of measured mass function to the fit, and smooth
        # smooth the ratio since it covers a smaller dynamic range
        ratio = np.log10(n/mf.number_density(logM))

        kernel = norm.pdf(np.arange(-25,26), loc=0, scale=8)
        kernel /= np.sum(kernel)
        ratio_convolved = np.convolve(ratio, kernel, mode='same')

        start, stop = 20, -20
        logM = logM[start:stop]
        logn = np.log10(10**ratio_convolved[start:stop] * mf.number_density(logM))

        logM_bins = np.arange(10,16,0.01)
        f_interp = interp1d(logM[30:-1], logn[30:-1], bounds_error=False, fill_value='extrapolate', kind='linear')
        logn_interp = f_interp(logM_bins)

        f = h5py.File(mass_function_file,'a')
        f.create_dataset('%i/z'%i, data=np.array([snapshot_redshift,]))
        f.create_dataset('%i/log_mass'%i, data=logM_bins, compression='gzip')
        f.create_dataset('%i/log_n'%i, data=logn_interp, compression='gzip')
        f.close()
        
        
if __name__ == '__main__':
    path_config_filename = sys.argv[1]
    with open(path_config_filename, "r") as file:
        path_config = yaml.safe_load(file)
    label = path_config["Labels"]["sim_label"]
    
    # Output file to save the mass functions
    mass_function_file = '../lookup/flamingo/flamingo_mass_functions_%s.hdf5'%(label)

    # Measure the mass functions. By default for all snapshots z<1.0
    get_mass_functions(path_config_filename, mass_function_file)
