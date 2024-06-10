#! /usr/bin/env python
from __future__ import print_function
import numpy as np
from astropy.table import Table, vstack
import yaml
import os

import sys
sys.path.append('..')
from hodpy.halo_catalogue import FlamingoSnapshot
from hodpy.galaxy_catalogue_snapshot import BGSGalaxyCatalogueSnapshotFlamingo
from hodpy.hod_bgs_abacus import HOD_BGS
from hodpy.colour import ColourDESI
from hodpy import lookup

def main(input_file, output_file, path_config_filename, photsys, snapshot_redshift=0.0, mag_faint=-18):
    '''
    Create a cubic box BGS mock
    '''
    import warnings
    warnings.filterwarnings("ignore")
    
    # to be safe, re-make the magnitude lookup tables on first loop iteration
    file_number = int(input_file[-8:-5])
    replace_lookup = file_number==0
    
    # create halo catalogue
    halo_cat = FlamingoSnapshot(input_file, path_config_filename=path_config_filename, snapshot_redshift=snapshot_redshift)

    # empty galaxy catalogue
    gal_cat  = BGSGalaxyCatalogueSnapshotFlamingo(halo_cat, path_config_filename)

    # use hods to populate galaxy catalogue
    hod = HOD_BGS(path_config_filename, photsys=photsys, mag_faint_type='absolute', mag_faint=mag_faint, redshift_evolution=True,
                  replace_central_lookup=replace_lookup, replace_satellite_lookup=replace_lookup)
    gal_cat.add_galaxies(hod)

    # position galaxies around their haloes
    gal_cat.position_galaxies()

    # add g-r colours
    col = ColourDESI(photsys=photsys, hod=hod)
    gal_cat.add_colours(col)

    # cut to galaxies brighter than absolute magnitude threshold
    gal_cat.cut(gal_cat.get("abs_mag") <= mag_faint)
    
    # save catalogue to file
    gal_cat.save_to_file(output_file, format="fits_BGS")
    
    
def join_files(path, photsys):
    '''
    Combine the cubic box outputs into a single file
    '''
    output_files_list = os.listdir(path)
    for input_file in output_files_list:
        print(input_file)
        input_file_number = input_file.split(".")[1]

        #table_i = Table.read(path+'BGS_box_%s_%03d.fits'%(photsys,file_number))
        table_i = Table.read(input_file)
        table_i['FILE_NUM'][:] = input_file_number

        #table_i = table_i['R_MAG_ABS', 'G_R_REST', 'HALO_MASS', 'cen', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'FILE_NUM', 'HALO_ID']

        if input_file_number==0:
            table = table_i
        else:
            table = vstack([table, table_i])

        del table_i
    
    # write the new table
    table.write(output_path+'BGS_box_%s.fits'%(photsys), format="fits")
    
        
if __name__ == "__main__":
    path_config_filename = sys.argv[1] # Config file path

    with open(path_config_filename, "r") as file:
        path_config = yaml.safe_load(file)

    soap_path = path_config["Paths"]["soap_path"]
    photsys = path_config["Params"]["photsys"]
    mag_faint = path_config["Params"]["mag_faint"]
    snapshot_redshift = path_config["Params"]["redshift"]

    # location of the snapshots
    soap_files_list = os.listdir(soap_path)

    output_path = "./output_files/"

    # each snapshot is split into 34 files
    for input_file in soap_files_list:
        print("Populating galaxies for snapshot "+input_file+"...")
        input_file_number = input_file.split(".")[1]
        output_file = output_path+'BGS_box_%s.%03d.fits'%(photsys, input_file_number)

        # populate the haloes with galaxies
        main(soap_path+input_file, output_file, path_config_filename=path_config_filename, photsys=photsys, snapshot_redshift=snapshot_redshift, mag_faint=mag_faint)

    print("Joining output files into a single file...")
    # join the 34 outputs into a single file
    join_files(output_path, photsys)

