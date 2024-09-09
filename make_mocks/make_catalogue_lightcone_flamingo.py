#! /usr/bin/env python
from __future__ import print_function
import numpy as np
from astropy.table import Table, vstack
import yaml
import os

import sys
sys.path.append("..")
from hodpy.halo_catalogue import FlamingoCatalogue
from hodpy.galaxy_catalogue import BGSGalaxyCatalogueFlamingo
from hodpy.cosmology import CosmologyFlamingo
#from hodpy.mass_function import MassFunctionMXXL WRONG
from hodpy.hod_bgs_flamingo import HOD_BGS
from hodpy.k_correction import GAMA_KCorrection
from hodpy.colour import ColourDESI
from hodpy import lookup


def main(input_file, output_file, path_config_filename, photsys, mag_faint=20.0):

    import warnings
    warnings.filterwarnings("ignore")

    # to be safe, re-make the magnitude lookup tables on first loop iteration
    file_number = int(input_file[-8:-5])
    replace_lookup = file_number==0

    # create halo catalogue
    halo_cat = FlamingoCatalogue(input_file)

    # empty galaxy catalogue
    gal_cat  = BGSGalaxyCatalogueFlamingo(halo_cat, path_config_filename)

    # use hods to populate galaxy catalogue
    hod = HOD_BGS(path_config_filename, redshift_evolution=True,
                  replace_central_lookup=replace_lookup, replace_satellite_lookup=replace_lookup)
    gal_cat.add_galaxies(hod)

    # position galaxies around their haloes
    gal_cat.position_galaxies()

    # add g-r colours
    col = ColourDESI(path_config_filename=path_config_filename, photsys=photsys, hod=hod)
    gal_cat.add_colours(col)

    # use colour-dependent k-correction to get apparent magnitude
    kcorr = GAMA_KCorrection(CosmologyFlamingo(path_config_filename))
    gal_cat.add_apparent_magnitude(kcorr)

    # cut to galaxies brighter than apparent magnitude threshold
    gal_cat.cut(gal_cat.get("app_mag") <= mag_faint)

    # save catalogue to file
    gal_cat.save_to_file(output_file, format="fits_BGS")
    

    
if __name__ == "__main__":
    path_config_filename = sys.argv[1] # Config file path

    with open(path_config_filename, "r") as file:
        path_config = yaml.safe_load(file)
    
    input_file = path_config["Paths"]["lightcone_path"]
    #input_file = "input/halo_catalogue_small.hdf5"
    output_file = "lightcone_catalogue.hdf5"
    mag_faint = path_config["Params"]["mag_faint"] # faintest apparent magnitude
    photsys = path_config["Params"]["photsys"]
    
    main(input_file, output_file, path_config_filename, photsys=photsys, mag_faint = mag_faint)
