###################################
#pypart_params.py
'''
SET PARAMETERS FOR PARTICLE RELEASE
SIMULATION
'''
###################################

run_name = 'SBC_SITES_L4_V5'


##########################################
# SET ROMS ID AND FILES TO INTEGRATE OVER
"""
ROMS_ID corresponds to cases 
in ROMS_solutions_paths.py class, add on to this
when using your own solution in the class
"""
#########################################
ROMS_ID = 'L4_MidCal_SBC'


################################
'''
fnum ---> first file index
frpf --> file index spacing (change this by multiples to change temporal resolution)
nfr --> number of frames to integrate over
fr --> starting frame in velocity file (zero index based)
dfr --> time-step increment per netcdf file (default is 1)
sub --> number of time-steps between frames
'''
##################################
fnum    = 144 
frpf    = 24
fr      = 0 
dfr     = 1  
nfr     = 48*40
sub     = 10 

##########################
# timing option
'''
True --> time processes
False --> no timing
'''
###########################
tim_opt = True



#####################################
# TYPE OF VELOCITY TO ADVECT
#####################################
veloc_adv = '3D'
file_type = 'roms' #roms or uv depending on file output
omega_online = True
z_online = True
t_interp = 'linear'
###############################
# set output file name
#################################
#out_part = run_name + '.nc'


############################################
# PARTICLE SEEDING PARAMETERS
'''
SEED WITH DESTINATION MASKS BASED ON h_bounds

but save a modified version of destination
masks that extends to shoreline for use
in connectivity analysis

in this manner, particles are sourced more offshore
than they can be trapped at a certain destination
'''
###########################################
seed_choice = 'preload_CC_masks_and_part'
path_mask = '/home/dauhajre/MidCal/py_particles/particle_runs/SBC_SITES_L4_V5/'
fname_mask = 'L4_CC_masks.p'
fname_init = 'L4_SBC_init_part_oct31.p'
nq_sites = 100
n_releases = 40*3
check_seed = True
pz_seed = 'z'
z_seed = -2 


