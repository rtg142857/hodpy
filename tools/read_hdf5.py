import numpy as np
import h5py
import swiftsimio as sw
import yaml
import os

from colossus.cosmology import cosmology as colossusCosmology
colossusCosmology.setCosmology('planck18')
from colossus.halo import mass_defs

def read_soap_log_mass(input_file, UnitMass_in_cgs, h, redshift, cosmology):
    """
    Get an array representing the base-10 logarithm of the 200-crit dark matter field halo masses from a SOAP file, in units of log solar mass/h.
    Args:
        input_file: Path to the SOAP file.
        UnitMass_in_cgs: The mass unit used in the SOAP calc (same as the snapshot units) in grams.
        h: Reduced Hubble constant.
    """
    halo_cat = h5py.File(input_file, "r")
    is_not_subhalo = np.array(halo_cat["InputHalos"]["HBTplus"]["Depth"]) == 0
    # Some (field) halos have 0 mass and 0 radius; these are probably the orphan halos that were ejected from their host halo; we discard these
    is_not_0mass = np.array(halo_cat["SO"]["200_crit"]["DarkMatterMass"]) != 0
    relevant_field_halos = np.logical_and(is_not_0mass, is_not_subhalo)

    UnitMass_in_Msol_h = UnitMass_in_cgs * h / 1.98841e33
    M200c = np.array(halo_cat["SO"]["200_crit"]["DarkMatterMass"])[relevant_field_halos] * UnitMass_in_Msol_h
    rvmax = np.array(halo_cat["BoundSubhalo"]["MaximumDarkMatterCircularVelocityRadius"])[relevant_field_halos] * h
    rho = cosmology.critical_density(redshift)
    r200c = (3./(800*np.pi) * M200c / rho)**(1./3) * (1.+redshift)
    conc = 2.16 * r200c / rvmax
    M200m, r200m, c200m = mass_defs.changeMassDefinition(M200c, conc, redshift, "200c", "200m", profile="nfw")

    log_mass = np.log10(M200m)
    return log_mass

def read_hbt_log_mass(input_file, UnitMass_in_cgs, h):
    """
    Get an array representing the base-10 logarithm of the 200-crit dark matter field halo masses from an HBT file, in units of log solar mass/h.
    Args:
        input_file: Path to the HBT file.
        UnitMass_in_cgs: The mass unit used in the HBT calc (same as the snapshot units) in grams.
        h: Reduced Hubble constant.
    """
    halo_cat = h5py.File(input_file, "r")
    is_not_subhalo = np.array(halo_cat["Subhalos"]["Rank"]) == 0
    # Some (field) halos have 0 mass and 0 radius; these are probably the orphan halos that were ejected from their host halo; we discard these
    is_not_0mass = np.array(halo_cat["Subhalos"]["BoundM200Crit"]) != 0
    relevant_field_halos = np.logical_and(is_not_0mass, is_not_subhalo)

    UnitMass_in_Msol_h = UnitMass_in_cgs * h / 1.98841e33
    log_mass = np.log10(np.array(halo_cat["Subhalos"]["BoundM200Crit"])[relevant_field_halos] * UnitMass_in_Msol_h)
    return log_mass

def find_field_particles_snapshot_file(input_file, group_id_default, particle_rate):
    """
    Identify the field particles (DM particles not in FOF halos) in a Flamingo snapshot file.
    Return a boolean array representing whether or not each particle in the file is a field particle.
    Args:
        input_file: Path to the Flamingo snapshot file.
        group_id_default: Halo ID given to field particles (particles not in halos).
        particle_rate: If particle_rate is N, only consider every Nth particle. Done to avoid out of memory errors.
    """
    data = sw.load(input_file)
    DM_IDs = data.dark_matter.fofgroup_ids[::particle_rate]
    field_boolean = (DM_IDs == group_id_default)
    del data
    del DM_IDs
    return field_boolean

def get_average_dm_particle_mass(path_config_filename):
    """ 
    Returns the average mass of a dark matter particle in a hdf5 file, in whatever mass is used internally
    Args: path_config_filename
    """
    with open(path_config_filename, "r") as file:
        path_config = yaml.safe_load(file)
    snapshot_path = path_config["Paths"]["snapshot_path"]

    if ".hdf5" in snapshot_path:
        file_path = snapshot_path
    else:
        snapshot_files_list = os.listdir(snapshot_path)
        snapshot_files_list = [file for file in snapshot_files_list if file.count(".") == 2]
        file_path = snapshot_path + snapshot_files_list[0]

    f = h5py.File(file_path, "r")
    mass_array = f["DMParticles"]["Masses"]
    return np.mean(mass_array)

def get_log_min_halo_mass(path_config_filename):
    """
    Gets the minimum mass of the halos found by a group finder, in log Msol/h.
    """
    with open(path_config_filename, "r") as file:
        path_config = yaml.safe_load(file)
    soap_path = path_config["Paths"]["soap_path"]
    try:
        halo_type = path_config["Misc"]["halo_type"]
    except:
        halo_type = "soap"
    used_parameters_path = path_config["Paths"]["params_path"]
    with open(used_parameters_path, "r") as file:
        params = yaml.safe_load(file)

    unit_mass = params["Snapshots"]["UnitMass_in_cgs"] / 1.98841e33
    particle_mass_Msol_h = get_average_dm_particle_mass(path_config_filename) * unit_mass * params["Cosmology"]["h"]
    
    if halo_type == "soap":
        halo_cat = h5py.File(soap_path, "r")
        min_halo_mass_Msol = particle_mass_Msol_h * halo_cat["BoundSubhalo"]["MaximumDarkMatterCircularVelocityRadius"].attrs["Mask Threshold"]
        return np.log10(min_halo_mass_Msol)
    elif halo_type == "peregrinus":
        # if ".hdf5" in soap_path:
        #     halo_cat = h5py.File(soap_path, "r")
        # else:
        #     soap_files_list = os.listdir(soap_path)
        #     soap_files_list = [file for file in soap_files_list if "Catalogue" in file]
        #     halo_cat = h5py.File(soap_path + soap_files_list, "r")
        print("Assuming 20 particles in smallest halos (double-check in hbt_params.txt, \"MinNumPartOfSub\")")
        MinNumPartOfSub = 20
        return np.log10(particle_mass_Msol_h * MinNumPartOfSub)

#def count_field_particles_snapshot_directory(input_directory, group_id_default):
    """
    Count the number of field particles (particles not in halos) across all Flamingo snapshot files in a directory,
    not including the virtual file containing every particle.
    Use this if you encounter out-of-memory errors.
    Args:
        input_file: Path to the Flamingo snapshot file.
        group_id_default: Halo ID given to field particles (particles not in halos).
    """