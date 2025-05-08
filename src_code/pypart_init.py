######################################
#pypart_init.py
'''
CREATE INPUT DICTIONARY FOR 
PARTICLE RELEASE THAT INCLUDES
ALL RELELVANT PARAMETERS INITIALIZED
IN pypart_params.py

AS WELL AS INITIAL PARTICLE
LOCATIONS THAT ARE SET USING
SPECIFICED SEEDING FUNCTION FROM
SEEDING LIBRARY
'''
#DANIEL DAUHAJRE UCLA OCTOBER 2016
######################################

############################################
#IMPORT MODULES
import os
import scipy.io as sio
import numpy as np
from netCDF4 import Dataset as netcdf
import pypart_funcs as PF
from utils import *
from part  import *
import ROMS_solutions_paths as ROMS_out_paths
import pickle as pickle
###########################################

####################################
# INITIALIZE PARAMETERS
#####################################
code_path = './src_code/'
execfile(code_path + 'pypart_params.py')


###########################################
# OBTIAN ROMS OUTPUT PATHS 
###########################################
#ROMS_ID = 'L4PV_nowec_P3_1min'
print 'Setting ROMS solutions paths for: ' + ROMS_ID
ROMS_obj = ROMS_out_paths.ROMS_run(ROMS_ID)
ROMS_obj.set_paths()
#################add by zsm
print 'hiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii*******************************'
print ROMS_ID
print ROMS_obj.path_grid + ROMS_obj.grid_name
print ROMS_obj.path_output + ROMS_obj.out_base_name
########################end
#####################################
#SET OUTPUT AND GRID FILE PATHS/NAMES
####################################
grd_name = ROMS_obj.path_grid + ROMS_obj.grid_name
dat_name = ROMS_obj.path_output + ROMS_obj.out_base_name


########################################
#INITIALIZE INPUT PARTICLE DICTIONARY
########################################
init_dict = {}
param_keys = ['ROMS_ID','run_name','fnum', 'frpf', 'nfr', 'fr', 'dfr','sub','veloc_adv', 'n_releases', 'tim_opt', 'omega_online', 't_interp','file_type']
param_key_types = [ROMS_ID, run_name,fnum, frpf, nfr, fr, dfr, sub, veloc_adv, n_releases,tim_opt, omega_online,t_interp,file_type]
init_dict = PF.add_keys_dict(param_keys, param_key_types, var_dict = init_dict) 
#################################add by zsm
print 'fuck it fuck it *******************************************************'
print param_keys
print ROMS_obj.path_output + ROMS_obj.out_base_name + '%04d' %fnum + '.nc','r'
####################################end

######################################################
# INITIALIZE PARTICLE LOCATIONS BASED ON SEEDING
# STRATEGY CHOSEN BY USER
#####################################################
nc_out = netcdf(ROMS_obj.path_output + ROMS_obj.out_base_name + '%04d' %fnum + '.nc','r')
nc_grd = netcdf(grd_name, 'r')

######################################
# SET HORIZONTAL LOCATIONS FIRST
####################################
if seed_choice == 'preload_CC_masks_and_part':
   dict_load = PF.load_pickle_file(path_mask + fname_mask)
   dict_part = PF.load_pickle_file(path_mask + fname_init)
   init_dict['dest_masks'] = dict_load['mask_sites']
   px_temp = dict_part['px']
   py_temp = dict_part['py']
if seed_choice == 'preload_CC_masks':
   dict_load = PF.load_pickle_file(path_mask + fname_mask)
   init_dict['dest_masks'] = dict_load['mask_sites']
   px_temp,py_temp = PF.fill_px_py_by_masks(init_dict['dest_masks'],nc_grd.variables['mask_rho'][:,:],nq_sites)
if seed_choice == 'box':
   px_temp,py_temp = PF.box_seed(nq,sx_0,sx_1,sy_0,sy_1)
   print sx_0
   print sx_1
   print sy_0
   print sy_1
if seed_choice == 'isobath_fill':
   px_temp,py_temp = PF.isobath_fill_seed(nc_grd,h1,h2,nq_isos)
if seed_choice == 'isobath_fill_sectors':
   px_temp,py_temp = PF.isobath_fill_sectors(nc_grd,h1,h2,nq_isos,i0,i1,n_sectors)
if seed_choice == 'radial_feature':
   #ACCESS NETCDF OUTPUT AND GRID FILES TO USE FOR SEEDING VISUALIZATION
   px_temp,py_temp = PF.radial_feature_seed(nc_grd, nc_out, fr,field_look=field_seed)
if seed_choice == 'multi_site_circle':
   px_temp, py_temp = PF.multi_site_circles(nc_grd, lat_sites, lon_sites, nq_sites, rad_sites, fig_check_seed = check_seed) 

if seed_choice == 'multi_site_along_iso':
   px_temp,py_temp = PF.multi_site_iso_contour_circles(nc_grd,h_seed,n_sites,nq_sites,rad_sites,fig_check_seed=check_seed)

if seed_choice == 'fill_multi_isos':
   px_temp, py_temp,dest_masks = PF.seed_between_isos(nc_grd, iso_list, n_sites, nq_sites, rad_sites, fig_check_seed = check_seed,check_dest_map=dest_check)
if seed_choice == 'mask_zsm':
    datap=sio.loadmat(pathmat + fname_mask)
    px_temp=datap['px']
    py_temp=datap['py']
    # px_temp
##########################################
# CHECK SEEDING
############################################
if check_seed:
   PF.check_fig_seed(nc_grd,px_temp,py_temp)


####################################
# VERTICAL LOCATIONS
####################################
if veloc_adv == '3D':
   if pz_seed == 'k':
      pz_temp = PF.fill_zs(px_temp,py_temp,nc_out,nc_grd,pz_choice=pz_seed,k_seed=k_seed,t_seed=fr)
   if pz_seed == 'z':
      pz_temp = PF.fill_zs(px_temp,py_temp,nc_out,nc_grd,pz_choice=pz_seed,z_seed=z_seed,t_seed=fr)

#################################
# FINALIZE init_dict 
# arrays with multiple releases
# if necessary
#################################
if n_releases >1:
   if veloc_adv == '3D':
      init_dict['px'], init_dict['py'], init_dict['pz'] = PF.multi_seed_3D(px_temp,py_temp,pz_temp,n_releases)
   else:
      init_dict['px'], init_dict['py'] = PF.multi_seed_2D(px_temp,py_temp,n_releases)
else:
   if veloc_adv == '3D':
      init_dict['px'] = px_temp 
      init_dict['py'] = py_temp
      init_dict['pz'] = pz_temp 
   else:
      init_dict['px'] = px_temp
      init_dict['py'] = py_temp

print px_temp.shape
print py_temp.shape
###################################################
#	SAVE INPUT DICTIONARY AS PICKLE FILE
##################################################
PF.change_dir_py('part_in')
PF.save_to_pickle(init_dict,run_name + '_in')
os.chdir('../')






