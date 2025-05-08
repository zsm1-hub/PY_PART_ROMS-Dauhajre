###################################
#pypart_params.py
'''
SET PARAMETERS FOR PARTICLE RELEASE
SIMULATION
'''
###################################

run_name = 'NOTIDE'


##########################################
# SET ROMS ID AND FILES TO INTEGRATE OVER
"""
ROMS_ID corresponds to cases 
in ROMS_solutions_paths.py class, add on to this
when using your own solution in the class
"""
#########################################
ROMS_ID = 'notide'


################################
'''
fnum ---> first file index
frpf --> file index spacing (change this by multiples to change temporal resolution) how many tindex in a netcdf file
nfr --> number of frames to integrate over
fr --> starting frame in velocity file (zero index based)
dfr --> time-step increment in per netcdf file (default is 1)
sub --> number of time-steps between frames, time intervals keyi tiaobu, or we can say, dt=delt/sub is real time-step in model
'''
##################################
fnum    = 16 
frpf    = 24
fr      = 0 
dfr     = 1  
nfr     = 16*24 
sub     = 1 

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
seed_choice = 'box'
nq=50
# # same to netcdf  read way y equals j
sy_0=300
sy_1=350
sx_0=150
sx_1=190
n_releases =1
check_seed = True
pz_seed = 'z'
z_seed = -2 

# seed_choice = 'isobath_fill'


# nq_isos=50
# h1=20
# h2=30
# n_releases =1
# check_seed = True

# pz_seed = 'z'
# z_seed = -2 

# seed_choice = 'mask_zsm'
# pathmat='/home/zsm/PY_PART_ROMS-Dauhajre/part_in/'
# fname_mask='boxpart.mat'
# n_releases =1
# check_seed = True

# pz_seed = 'z'
# z_seed = -2 