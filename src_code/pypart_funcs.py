#######################################
#pypart_funcs.py
'''
LIBRARY OF FUNCTIONS USED IN PYTHON
PARTICLE CODE
'''
######################################

#########################################
#IMPORT MODULES
import sys
sys.path.append('./src_code/f2py_modules/')
import pickle as pickle
import numpy as np
import os
from pylab import *
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf
import time as clock
import omega_fort as FO
import ROMS_depths as RD
######################################



##################################################
# FUNCTION TO ACCESS ROMS VELOCITIES
##################################################
def get_uvW_nc(nc_obj,tind,timing=False,z_r=[],z_w=[], pm=[],pn=[], nc_omega=[],omega_online=True,t_interp='spline'):
    '''
    LOAD u,v,W and respective time-derivatives dt_u, dt_v, dt_W
    at time-specified
    '''
    if timing:
       start_time = clock.time()
    u_temp = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables['u'][tind,:,:,:]),-1,0))
    v_temp = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables['v'][tind,:,:,:]),-1,0))
    if t_interp =='spline':
       dt_u_temp  = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables['dt_u'][tind,:,:,:]),-1,0))
       dt_v_temp  = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables['dt_v'][tind,:,:,:]),-1,0))
    if timing:
       print '  u,v reading time: ' + str(clock.time() - start_time)
       start_time = clock.time()
    
    #EXTRAPOLATE u,v to boundaries at surface and bottom for 
    # index-space advection within free-surface and bottom (defined by w-levels)
    [nx_u,ny,nz] = u_temp.shape
    [nx,ny_v,nz] = v_temp.shape
    u_out = np.zeros([nx_u,ny,nz+2])
    v_out = np.zeros([nx,ny_v,nz+2])
    u_out[:,:,1:-1] = np.copy(u_temp)
    v_out[:,:,1:-1] = np.copy(v_temp)
    #HOMOGENOUS NEUMAN B.C.
    u_out[:,:,-1]   = u_temp[:,:,-1]
    u_out[:,:,0]    = u_temp[:,:,0]
    v_out[:,:,-1]   = v_temp[:,:,-1]
    v_out[:,:,0]    = v_temp[:,:,0]

    u_out[u_out > 1e30] = 0.
    v_out[v_out > 1e30] = 0.

    if t_interp =='spline':
       dt_u = np.zeros([nx_u,ny,nz+2])
       dt_v = np.zeros([nx,ny_v,nz+2])
       dt_u[:,:,1:-1] = np.copy(dt_u_temp)
       dt_v[:,:,1:-1] = np.copy(dt_v_temp)
       dt_u[:,:,-1]   = dt_u_temp[:,:,-1]
       dt_u[:,:,0]    = dt_u_temp[:,:,0]
       dt_v[:,:,-1]   = dt_v_temp[:,:,-1]
       dt_v[:,:,0]    = dt_v_temp[:,:,0]

       dt_u[dt_u > 1e30] =  0.
       dt_v[dt_v > 1e30] = 0.
 
    #return u_out,v_out,0,0,0,0
    #PRE-CALCULATED OMEGA FIELDS
    if omega_online==False:
       W_out = np.swapaxes(np.squeeze(nc_omega.variables['omega'][tind,:,:,:]),-1,0)
       if t_interp =='spline':
          dt_W = np.swapaxes(np.squeeze(nc_omega.variables['dt_omega'][tind,:,:,:]),-1,0)
       if timing:
          print '  W reading time: ' + str(clock.time() - start_time)

    # ONLINE OMEGA CALCULATION
    if omega_online:
       z_w[z_w>1e30] = 0.
       z_r[z_r>1e30] = 0.
       if timing:
          start_time = clock.time()
       W_out = FO.get_omega(u_out[:,:,1:-1],v_out[:,:,1:-1],np.swapaxes(z_r,-1,0),np.swapaxes(z_w,-1,0),np.swapaxes(pm,-1,0),np.swapaxes(pn,-1,0))
       if t_interp == 'spline':
          dt_W = FO.get_omega(dt_u[:,:,1:-1],dt_v[:,:,1:-1],np.swapaxes(z_r,-1,0),np.swapaxes(z_w,-1,0),np.swapaxes(pm,-1,0),np.swapaxes(pn,-1,0))
       if timing:
          print '   Online W, dt_W calculation: ' + str(clock.time() - start_time)

    # ZERO OUT MASKED VALUES
    W_out[W_out > 1e30] = 0.
    W_out[np.isnan(W_out)] = 0.
    if t_interp == 'spline':
       dt_W[dt_W>1e30] = 0.
       dt_W[np.isnan(dt_W)] = 0.
       return u_out,v_out, W_out,dt_u,dt_v,dt_W
    else:
       return u_out,v_out,W_out
    ######################################################################################

def get_uv_2D(nc_obj,tind,k_load=-1,timing=False,t_interp='linear',keys=['u','v']):
    '''
    LOAD 2D u,v and respective time-derivatives dt_u, dt_v
    at time-specified

    t_interp --> 'linear' (linear interpolation, only load u,v)
                 'spline' (spline interplation, load u,v,dt_u,dt_v-->pre-calculated)
    '''
    dt_keys = ['dt_' + n for n in keys]

    ##########################################
    # LOAD FROM k-level
    ##########################################
    if len(nc_obj.variables[keys[0]].shape)>3:
       if timing:
          start_time = clock.time()
       u_out = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables[keys[0]][tind,k_load,:,:]),-1,0))
       v_out = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables[keys[1]][tind,k_load,:,:]),-1,0))
       if t_interp=='spline':
          dt_u  = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables[dt_keys[0]][tind,k_load,:,:]),-1,0))
          dt_v  = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables[dt_keys[1]][tind,k_load,:,:]),-1,0))
       if timing:
          print '  u,v reading time: ' + str(clock.time() - start_time)
       u_out[u_out > 1e30] =  np.nan
       v_out[v_out > 1e30] = np.nan
       if t_interp == 'spline': 
          dt_u[dt_u > 1e30] =  np.nan
          dt_v[dt_v > 1e30] = np.nan
          return u_out, v_out, dt_u, dt_v
 
    ##########################################
    # LOAD 2D u,v: either bar or z-sliced
    #####################################   
    else:
       if timing:
          start_time = clock.time()
       u_out = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables['u'][tind,:,:]),-1,0))
       v_out = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables['v'][tind,:,:]),-1,0))
       if t_interp=='spline':
          dt_u  = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables['dt_u'][tind,:,:]),-1,0))
          dt_v  = np.asarray(np.swapaxes(np.squeeze(nc_obj.variables['dt_v'][tind,:,:]),-1,0))
       if timing:
          print '  u,v reading time: ' + str(clock.time() - start_time)
       u_out[u_out > 1e30] =  np.nan
       v_out[v_out > 1e30] = np.nan
       if t_interp == 'spline': 
          dt_u[dt_u > 1e30] =  np.nan
          dt_v[dt_v > 1e30] = np.nan
          return u_out, v_out, dt_u, dt_v
    
    return u_out, v_out
    ######################################################################################

def get_zs_nc(nc_obj,tind,timing=False):
    '''
    LOAD z_r, z_w t time-specified
    '''
    if timing:
       start_time = clock.time()

    z_r = nc_obj.variables['z_r'][tind,:,:,:]
    z_w = nc_obj.variables['z_w'][tind,:,:,:]
    print '  z reading time: ' + str(clock.time() - start_time)
    z_r[z_r>1e30] = np.nan
    z_w[z_w>1e30] = np.nan
    return z_r, z_w 
    ######################################################################################


	####################################################################
	#	          SET UP NETCDF OUTPUT SEQUENCE OF FILES 
        ######################################################################
def setup_out_nc(base_name, tstart,tend,tfile, nq,veloc_adv='3D'):
    '''
    CREATE NETCDF FILES THAT WILL STORE OUTPUT
    TO MATCH UP WITH ROMS-HISTORY FILE SEQUENCING (e.g.,
    his.0000.nc will match up with part_out.0000.nc)
    '''
    #print 'Moving to directory ' + path_nc
    #os.chdir(path_nc)
    nc_name_list = []
    for t in range(tstart,tend,tfile):
        nc_name = base_name + '.' +  str(t).zfill(4) + '.nc'
        nc_name_list.append(nc_name)
        print 'Creating output file: ' + nc_name
        if veloc_adv == '3D':
           create_nc_part_3d(nc_name,nq)
        else:
           create_nc_part_2d(nc_name,nq)
    return nc_name_list
    #################################################

def create_nc_part_2d(nc_name, nq):
    '''
    CREATE OUTPUT FILE 
    '''
    nc = netcdf(nc_name, 'w')
    nc.createDimension('time',0)
    nc.createDimension('nq',nq)
    
    nct  = nc.createVariable('time','d',('time',))
    ncf  = nc.createVariable('frame','i',('time',))
    ncpx = nc.createVariable('px','d',('time','nq'))
    ncpy = nc.createVariable('py','d',('time','nq'))
    nc.close()
    #################################################

def create_nc_part_3d(nc_name, nq):
    '''
    CREATE OUTPUT FILE 
    '''
    nc = netcdf(nc_name, 'w')
    nc.createDimension('time',0)
    nc.createDimension('nq',nq)
    
    nct  = nc.createVariable('time','d',('time',))
    ncf  = nc.createVariable('frame','i',('time',))
    ncpx = nc.createVariable('px','d',('time','nq'))
    ncpy = nc.createVariable('py','d',('time','nq'))
    ncpz = nc.createVariable('pz','d',('time','nq'))
    nc.close()
    #################################################



def write_data_nc_2d(nc_obj, px, py,time, frame, idx):
    nc_obj.variables['time'][idx]   = time
    nc_obj.variables['px'][idx,:]   = px
    nc_obj.variables['py'][idx,:]   = py
    nc_obj.variables['frame'][idx]  = frame
    ############################################



def write_data_nc_3d(nc_obj, px, py, pz,time, frame, idx):
    nc_obj.variables['time'][idx]   = time
    nc_obj.variables['px'][idx,:]   = px
    nc_obj.variables['py'][idx,:]   = py
    nc_obj.variables['pz'][idx,:]   = pz
    nc_obj.variables['frame'][idx]  = frame
    ############################################



	####################################################################
	#	   DICTIONARY SETUP I/O FUNCTIONS (FOR INPUT DICT)
        #####################################################################

def add_keys_dict(keys_add, keys_type, var_dict={}): 
    """
    var_dict ---> dictionary to add keys to (default is empty {} or
                  can send a pre-existing dictionary)
   
    keys_add --> list strings of keys to add to dictionary

    keys_type --> list of data types corresponding to each key (i.e, empty list, numpy array, constant)
    """
    for ke in range(len(keys_add)):
	var_dict[keys_add[ke]] = keys_type[ke]

    return var_dict
    ####################################################


def save_to_pickle(var_dict, out_name):
    """
    SAVE DICTIONARY WITH PICKLE
    var_dict ---> dictionary
    out_name ---> file name to save as
    """
    pickle.dump(var_dict,open(out_name + ".p", "wb"))
    print 'Saved dictionary as: ' + out_name + '.p'
    ###################################################

def load_pickle_file(path_full):
    """
    LOAD PICKLE DICTIONARY
    path_full --> full path  including file name and  ".p" suffix
    """
    print ''
    print 'LOADING : ' + path_full
    print ''
    return pickle.load(open(path_full, "rb"))
    #############################################



def change_dir_py(dir_name):
    if not os.path.exists(dir_name):
       print 'Creating directory: ' + dir_name
       os.makedirs(dir_name)
    print 'Moving into directory: ' + dir_name
    os.chdir(dir_name)
    ##############################################


	####################################################################
	#	   PARTICLE SEEDING FUNCTIONS
        #####################################################################
def z_to_k(z_in,px,py,z):
    '''
    Convert depth to k-level index
    to seed particles at specific depth
    '''
    pz = np.zeros(len(px))
    for ip in range(len(px)):
        i = int(np.floor(px[ip]+0.5))
        j = int(np.floor(py[ip]+0.5))
        z_arr = z[i,j,:]
        k_arr = np.arange(len(z_arr))
        #simple interpolation of z to k 
        pz[ip] =  np.interp(z_in,z_arr,k_arr)
    return pz
    #####################################

def fill_zs(px,py,ncr,ncg,pz_choice='k',k_seed=0,z_seed=0,t_seed=0):
    '''
    Fill pz array with either k-index (k_seed)
    or k-index defined by a depth (z_seed<0)
    '''
    pz = np.zeros(len(px))
    if pz_choice == 'k':
       pz[:] = k_seed
    if pz_choice == 'z':
       #GET DEPTHS
       [Ly,Lx] = ncg.variables['pm'].shape
       #z_r = ncz.variables['z_r'][t_seed,:,:,:]
       #z_w = ncz.variables['z_w'][t_seed,:,:,:]
       z_r,z_w = RD.get_zr_zw_tind(ncr,ncg,t_seed,[0,Ly,0,Lx])
       pz[:] = z_to_k(z_seed,px,py,np.swapaxes(z_w,-1,0)) 
    return pz
    ########################################################

def multi_seed_3D_new(px_in, py_in, pz_in, n_releases):
    '''
    CREATE  px, py that has initial particle positions
    spaced every n_releases

    e.g., 
    if n_releases = 4
    px = [px0_0, px1_0, px_2_0, px_0_1, px1_1...]
         where px0_0 is position of particle 0 for release=0
	       px0_1 is position of particle 0 for release=1


    px_in, py_in --> initial particle positions made assuming a single release
    '''
    nq = len(px_in)
    px = np.zeros(nq*n_releases, order='F')
    py = np.zeros(nq *n_releases, order='F')
    pz = np.zeros(nq*n_releases,order='F')
    
    ind=0
    for n in range(n_releases):
        px[ind:ind+nq] = px_in
        py[ind:ind+nq] = py_in
        pz[ind:ind+nq] = pz_in
        ind+=nq

    return px, py, pz
    ##################################################################


def multi_seed_3D(px_in, py_in, pz_in, n_releases):
    '''
    CREATE  px, py that has initial particle positions
    spaced every n_releases

    e.g., 
    if n_releases = 4
    px = [px0_0, px0_1, px0_2, px0_3, px1_0, ...]
         where px0_0 is position of particle 0 for release=0
	       px0_1 is position of particle 0 for release=1


    px_in, py_in --> initial particle positions made assuming a single release
    '''
    nq = len(px_in)
    px = np.zeros(nq*n_releases, order='F')
    py = np.zeros(nq *n_releases, order='F')
    pz = np.zeros(nq*n_releases,order='F')

    for n in range(nq):
	px[n*n_releases:(n*n_releases) + n_releases] = px_in[n]
	py[n*n_releases:(n*n_releases) + n_releases] = py_in[n]
	pz[n*n_releases:(n*n_releases) + n_releases] = pz_in[n]


    return px, py, pz
    ##################################################################


def multi_seed_2D(px_in, py_in, n_releases):
    '''
    CREATE  px, py that has initial particle positions
    spaced every n_releases

    e.g., 
    if n_releases = 4
    px = [px0_0, px0_1, px0_2, px0_3, px1_0, ...]
         where px0_0 is position of particle 0 for release=0
	       px0_1 is position of particle 0 for release=1


    px_in, py_in --> initial particle positions made assuming a single release
    '''
    nq = len(px_in)
    px = np.zeros(nq*n_releases, order='F')
    py = np.zeros(nq *n_releases, order='F')
    for n in range(nq):
	px[n*n_releases:(n*n_releases) + n_releases] = px_in[n]
	py[n*n_releases:(n*n_releases) + n_releases] = py_in[n]
    return px, py
    ##################################################################


def seed_between_isos(nc_grd,iso_list,ns,nqs,rs,fig_check_seed=True,check_dest_map=True):
    '''
    SEED IN CIRCLES AT MUTLIPLE ISOBATHS
    '''
    px_temp = []
    py_temp = []

    nisos = len(iso_list)
    nparticles = nisos * ns * nqs
    px_temp = []
    py_temp = []
    #px_arr = np.zeros(nparticles)
    #py_arr = np.zeros(nparticles)
    pind = 0
    #dp = ns * nqs
    x_sites = []
    y_sites = []
    for i in range(nisos):
        #print str(pind + dp)
        #px_arr[pind:pind+dp], py_arr[pind:pind+dp] = multi_site_iso_contour_circles(nc_grd,iso_list[i],ns,nqs, rs, fig_check_seed=False)
        #pind+=dp           
        px, py,xcent,ycent = multi_site_iso_contour_circles(nc_grd,iso_list[i],ns,nqs, rs, fig_check_seed=False)
        px_temp.extend(px.tolist())
        py_temp.extend(py.tolist())

        x_sites.extend(xcent)
        y_sites.extend(ycent)

    
    if check_dest_map:
       nsites = len(x_sites)
       [Ly,Lx] = nc_grd.variables['h'][:,:].shape
       dest_masks = [np.zeros([Ly,Lx]) for n in range(nsites)]
       r_grd = (rs *1E3) * (np.mean((nc_grd.variables['pm'][:,:] + nc_grd.variables['pn'][:,:])   /2))
       print '     CREATING MASKED ARRAY FOR DESTINATIONS'
       for s in range(nsites):
           print '     site #' + str(s+1)
           #CENTER-POINT OF SITES
           x_s = x_sites[s]
           y_s = y_sites[s]
           for j in range(Ly):
               for i in range(Lx):
                   dist = np.sqrt( (x_s - i)**2 + (y_s - j)**2)
                   if dist <= r_grd:
                      dest_masks[s][j,i] = 1.

    else:
       dest_masks = []
    #print 'len = ' + str(len(px_temp))
    px_arr = np.asarray(px_temp)
    py_arr = np.asarray(py_temp)
    if fig_check_seed:
    	mask = nc_grd.variables['mask_rho'][:,:]
	h = nc_grd.variables['h'][:,:].T
	mask_rho = np.copy(mask)
        mask_rho[mask_rho==0] = np.nan
        plt.ion() 
        fig_check = plt.figure(figsize=[24,8])
        #h_con = plt.contour(h.T,colors='k',linewidths=2.5)
        plt.imshow(h.T*mask_rho,origin='lower',cmap=cm.jet)
        cbar_plot = plt.colorbar()
        cbar_plot.set_label(r'$h(m)$',fontsize=22)
        plt.plot(px_arr,py_arr, 'o', color='k',markersize=2.5)
        plt.title('INITIAL PARTICLE LOCATIONS')
	ax = plt.gca()
        ax.set_xlim([0,mask_rho.shape[1]])
        ax.set_ylim([0,mask_rho.shape[0]])
        move_on = raw_input('Press enter to continue...')
	plt.ioff()
	plt.close('all')

        

    return px_arr, py_arr, dest_masks
    #################################################


def multi_site_iso_contour_circles(nc_grd,h_seed,ns,nqs, rs, fig_check_seed=True):
    '''
    SEED CIRCLES OF PARTICLES AT MULTIPLE SITES
    DEFINED AS POINTS ALONG AN ISOBATH
    
    h_seed --> isobath to take as contour line to seed particles along

    ns  --> number of sites along contour line
    
    nqs --> number of particles per each site

    rs --> radius in km of particles for each site
    '''
    [Ly,Lx] = nc_grd.variables['h'][:,:].shape
    # GET COORDINATES OF ISOBATH CONTOUR
    CD = plt.contour(nc_grd.variables['h'][:,:], [h_seed])
    len_cons = [len(CD.collections[0].get_paths()[i]) for i in range(len(CD.collections[0].get_paths()))]
    ind_max_len = np.asarray(len_cons).argmax()
    p = CD.collections[0].get_paths()[ind_max_len]
    x_iso = p.vertices[:,0]
    y_iso = p.vertices[:,1]
   
    #GET COORDINATES OF EACH SITE EVENLY DISTRIBUTED ALONG THE ISOBATH
    dpt = int(len(x_iso)/ns)
    #len_temp = len(x_iso[0::dpt])
    #diff_len = abs(len_temp - ns)
    #x_site = x_iso[::dpt][0:len_temp-diff_len]
    #y_site = y_iso[::dpt][0:len_temp-diff_len]

    #slicing below avoids asymmetry at lateral boundaries
    x_site = x_iso[dpt:-1:dpt]
    y_site = y_iso[dpt:-1:dpt]
     
    nsites = len(x_site)
    #print nsites
    nparticles = len(x_site)*nqs
    #print nparticles
    px_arr = np.zeros(nparticles,order='F')
    py_arr = np.zeros(nparticles, order='F')
    p_ind = 0 
    for s in range(nsites):
        ##################################
	# FIND SITE CENTERS IN GRIDPOINTS
	####################################
 	r_grd = (rs *1E3) * (np.mean((nc_grd.variables['pm'][:,:] + nc_grd.variables['pn'][:,:])   /2))
       

        ####################################
	# CREATE CIRCLE OF PARTICLES
	####################################
	px_arr[p_ind:p_ind+nqs], py_arr[p_ind:p_ind+nqs] = part_circle(x_site[s], y_site[s], r_grd, nqs)
        p_ind+=nqs

    buff = 3
    px_arr[px_arr>=Lx-1] = Lx- buff
    py_arr[py_arr>=Ly-1] = Ly - buff

    #REMOVE PARTICLES ON MASK
    px_arr,py_arr = nan_parts_on_mask(px_arr,py_arr,nc_grd)    
    ####################################
    # SHOW USER FLOAT LOCATIONS ON GRID
    ####################################
    if fig_check_seed:
	mask = nc_grd.variables['mask_rho'][:,:]
	h = nc_grd.variables['h'][:,:].T
	mask_rho = np.copy(mask)
        mask_rho[mask_rho==0] = np.nan
        plt.ion() 
        fig_check = plt.figure(figsize=[24,8])
        #h_con = plt.contour(h.T,colors='k',linewidths=2.5)
        plt.imshow(h.T*mask_rho,origin='lower',cmap=cm.jet)
        cbar_plot = plt.colorbar()
        cbar_plot.set_label(r'$h(m)$',fontsize=22)
        plt.plot(px_arr,py_arr, 'o', color='k',markersize=2.5)
        plt.title('INITIAL PARTICLE LOCATIONS')
	ax = plt.gca()
        ax.set_xlim([0,mask_rho.shape[1]])
        ax.set_ylim([0,mask_rho.shape[0]])
        move_on = raw_input('Press enter to continue...')
	plt.ioff()
	plt.close('all')


    return px_arr, py_arr, x_site, y_site
    #####################################



def nan_parts_on_mask(px,py,nc_grd):
    npart = len(px)
    px_new = []
    py_new = []
    mask_rho = nc_grd.variables['mask_rho'][:,:]
    for ip in range(npart):
        i = int(np.floor(px[ip]+0.5))
        j = int(np.floor(py[ip]+0.5))
        if mask_rho[j,i] == 1.:
           px_new.append(px[ip])
           py_new.append(py[ip])
    return np.asarray(px_new),np.asarray(py_new)

def multi_site_circles(nc_grd, lat_sites, lon_sites, nq_sites, rad_sites, fig_check_seed=True):
    '''
    SEED CIRCLES OF PARTICLES AT MULTIPLE SITES
    SPECIFIED BY USER IN LAT / LON

    lat_sites --> list of latitude points of sites
    lon_sites --> list of longitude points of sites

    nq_site --> (list) number of particles per each site

    rad_site --> (list) radius in km of particles for each site
    '''
    
    ###############################################
    # FIND GRID POINTS CORRESPONDING TO LAT / LON
    ###############################################
    lon_nc = nc_grd.variables['lon_rho'][:,:]
    lat_nc = nc_grd.variables['lat_rho'][:,:]
   
    nsites = len(lat_sites)
    nparticles = sum(nq_sites)

    px_arr = np.zeros(nparticles,order='F')
    py_arr = np.zeros(nparticles, order='F')

    p_ind = 0
    for s in range(nsites):
        ##################################
	# FIND SITE CENTERS IN GRIDPOINTS
	####################################
	min_1D = np.abs( (lat_nc - lat_sites[s])**2 + (lon_nc - lon_sites[s])**2)
        y_site, x_site = np.unravel_index(min_1D.argmin(), min_1D.shape)
	r_grd = (rad_sites[s] *1E3) * np.mean((nc_grd.variables['pm'][:,:] + nc_grd.variables['pn'][:,:])   /2)
 
        ####################################
	# CREATE CIRCLE OF PARTICLES
	####################################
	px_arr[p_ind:p_ind+nq_sites[s]], py_arr[p_ind:p_ind+nq_sites[s]] = part_circle(x_site, y_site, r_grd, nq_sites[s])
        p_ind+=nq_sites[s]


    ####################################
    # SHOW USER FLOAT LOCATIONS ON GRID
    ####################################
    if fig_check_seed:
	mask = nc_grd.variables['mask_rho'][:,:]
	h = nc_grd.variables['h'][:,:].T
	mask_rho = np.copy(mask)
        mask_rho[mask_rho==0] = np.nan
        plt.ion() 
        fig_check = plt.figure(figsize=[24,8])
        #h_con = plt.contour(h.T,colors='k',linewidths=2.5)
        plt.imshow(h.T*mask_rho,origin='lower',cmap=cm.gist_earth_r)
        cbar_plot = plt.colorbar()
        cbar_plot.set_label(r'$h(m)$',fontsize=22)
        plt.plot(px_arr-1,py_arr-1, 'o', color='k',markersize=2.5)
        plt.title('INITIAL PARTICLE LOCATIONS')
	ax = plt.gca()
        ax.set_xlim([0,mask_rho.shape[1]])
        ax.set_ylim([0,mask_rho.shape[0]])
        move_on = raw_input('Press enter to continue...')
	plt.ioff()
	plt.close('all')




    return px_arr, py_arr
    #####################################


def mask_sites_circles(cent_sites,r_sites,mask_land):
    '''
    Create destination masked sites
    based on circles of particular radius
    that also do not intersect 2D array mask_land

    mask_land does not have to be mask_rho, in sense
    that particles can end up on actual land mask, so mask_land
    can portrude a bit into actual land mask
    '''
    nsites = len(cent_sites)
    [Ly,Lx] = mask_land.shape
    masks_temp = [np.zeros([Ly,Lx]) for s in range(nsites)]
    print 'Creating masked sites...'
    for i in range(Lx):
        if (i+1)%10 == 0: print round(100.*(i+1)/Lx-1), '%'
        for j in range(Ly):
            for s in range(nsites):
                x_s = cent_sites[s][1]
                y_s = cent_sites[s][0]
                dist = np.sqrt( (x_s - i)**2 + (y_s-j)**2)
                if dist <=r_sites:
                   masks_temp[s][j,i] = 1.

    #NOW MULTIPLY dest_masks times mask_land
    masks_out = [masks_temp[s] * mask_land for s in range(nsites)]
    return masks_out
    ###########################################################

def part_circle(x_cent, y_cent, r_grdpt, nq):
    '''
    MAKE A CIRCLE OF PARTICLES FOR GIVEN RADIUS
    AROUND CENTER COORDINATES

    x_cent --> x-gridpoint
    y_cent --> y-gridpoint
    r_grdpt   ---> radius of circle in grid points
    nq     ---> number of particles within circle
    '''
    #print 'r = ' + str(r_grdpt)
    rho = np.sqrt(np.random.uniform(0,r_grdpt**2,nq))
    phi = np.random.uniform(0,2*np.pi,nq)
    px_arr = x_cent + (rho*np.cos(phi))
    py_arr = y_cent + (rho*np.sin(phi))

    
    '''
    px_list = []
    py_list = []
    for n in range(nq):
        q = np.random.random() * 2 * np.pi
	r = np.sqrt(np.random.random())
        px_list.append(x_cent + (r_grdpt * r) * np.cos(q))
	py_list.append(y_cent + (r_grdpt * r) * np.sin(q))
    
    px_arr =  np.asarray(px_list,dtype='float64') 
    py_arr = np.asarray(py_list,dtype='float64') 
    '''
    return px_arr, py_arr
    ###################################################




def box_seed(nq,i0,i1,j0,j1):
     """
     SEED PARTICLES IN A BOX 

     nq ---> number of particles
     i0 --> x-index lower-bound of box
     j0 --> y-index lower-bound of box
     i1 --> x-index upper-bound of box
     j1 --> y-index upper-bound of box

     n_ind --> index multiplier

     """
     from random import randint
     px = np.zeros(nq,order='F')
     py = np.zeros(nq,order='F')
     px[:] = np.nan
     py[:] = np.nan
     
     # RANDOMLY DISTRIBUTE PARTICLES WITHIN BOUNDS OF BOX
     for n in range(nq):
         px[n] = randint(i0,i1)
         py[n] = randint(j0,j1)

     print '###################################################'
     print 'DEFAULT SEEDING:'
     print ' nq = ' + str(nq)  + ' total particles'
     print ''
     print '####################################################'
     return px,py
     ###################################################

def create_source_dest_masks(h,h_bounds,n_sites,i0,i1):
    '''
    CREATE A LIST OF MASKED ARRAYS
    THAT DEFINE SOURCE/DESTINATION SITES
    FOR PARTICLES THAT LIE BETWEEN 2 isobaths (h1, and h2)

    h --> [Ly,Lx] array of depth
    h_bounds --> list of isobath cutoffs [ [h1,h2], [h1,h2], etc]
              where  h1 --> isobath 1 cutoff (m)
                     h2 --> isobath 2 cutoff (m)
               and      h1 < h2
    n_sites --> number of sites to split up along xi-axis
    i0,i1 --> (int), xi-axis lateral boundaries for start/stop of sites
    '''
    [Ly,Lx] = h.shape
    out_masks = []
    for c in range(len(h_bounds)):
        [h1,h2] = h_bounds[c]
        bool_h = (h>=h1) & (h<=h2)
        mask_H = np.zeros([Ly,Lx])
        mask_H[bool_h] = 1.

        #DESTINATION MASK GENERATION
        dest_masks = [np.zeros([Ly,Lx]) for n in range(n_sites)]
        di = (i1-i0) / n_sites
        bounds = range(i0,i1+di,di)
        for b in range(len(bounds)-1):
            dest_masks[b][:,bounds[b]:bounds[b+1]] = mask_H[:,bounds[b]:bounds[b+1]]
        out_masks.extend(dest_masks)
    
    return out_masks
    ######################################



def fill_px_py_by_masks(dest_masks,mask_rho,nq_site):
    '''
    FILL PARTICLE ARRAYS WITH
    RANDOM DISTRIBUTION OF PARTICLES
    WITHIN A SECTOR DEFINED BY dest_masks

    Corraborate particle placement with mask_rho
    if dest_masks overlaps with land_mask so particles
    are not initialized on land points
    '''
    px_l = []
    py_l = []
    n_sites = len(dest_masks)
    for n in range(n_sites):
        y,x = np.where(dest_masks[n]*mask_rho==1)
        inds = np.random.randint(0,len(x),size=nq_site)
        px_l.extend([x[inds[i]] for i in range(nq_site)])
        py_l.extend([y[inds[i]] for i in range(nq_site)])
    return np.asarray(px_l), np.asarray(py_l)
    ###################################################


def get_iso_vertices(h,iso1,iso2,diso=0.1):
    '''
    GET A LIST OF x,y POINTS
    THAT LIE BETWEEN 2 ISOBATHS FOR
    PARTICLE SEEDING

    WILL GENERATE A LARGE AMOUNT OF (x,y)
    POINTS THAT ARE THEN RANDOMLY INDEXED
    FOR PARTICLE SEEDING
    '''
    iso_hd = np.arange(iso1,iso2,diso)
    CD = plt.contour(h,iso_hd)
    x = []
    y = []
    for c in range(len(CD.collections)):
        len_cons = [len(CD.collections[c].get_paths()[i]) for i in range(len(CD.collections[c].get_paths()))]
        ind_max_len = np.asarray(len_cons).argmax()
        p = CD.collections[c].get_paths()[ind_max_len]
        x.extend(p.vertices[:,0].tolist())
        y.extend(p.vertices[:,1].tolist())
    return x,y

def isobath_fill_seed(nc_grd,h1,h2,nq):
    """
    FILL AREA BETWEEN 2 ISOBATHS WITH PARTICLES
    
    nc_grd ---> netcdf path of ROMS grid

    h1 ---> shallowest isobath boundary
    h2 ---> deepest isobath boundary

    nq ---> number of particles (estimate)
    """
    #DECLARE LISTS, # OF PARTICLES IS NOT SPECIFIED, SO APPEND POSTIONS TO LIST FIRST
    px_l = []
    py_l = []
    h = nc_grd.variables['h'][:,:]
    x_verts,y_verts = get_iso_vertices(h,h1,h2)
    inds = np.random.randint(0,len(x_verts),size=nq)
    px_l = [x_verts[inds[i]] for i in range(nq)]
    py_l = [y_verts[inds[i]] for i in range(nq)]
    return np.asarray(px_l), np.asarray(py_l)
    ###########################################
    


def isobath_fill_sectors(nc_grd,h1,h2,nq,i0,i1,n_secs):
    '''
    Fill area between 2 isobaths with particles

    same as isobath_fill_seed, but 
    divide up isobath area laterally in X-coordinate
    and fill each sector with (nq) particles
    '''
    px_l = []
    py_l = []
    h = nc_grd.variables['h'][:,:]
    [Ly,Lx] = h.shape
    x_verts,y_verts = get_iso_vertices(h,h1,h2)
    x_arr = np.asarray(x_verts)
    y_arr = np.asarray(y_verts)
    di = (i1-i0) / n_secs
    bounds = range(i0,i1+di,di)
    for b in range(len(bounds)-1):
        # COORDINATES ALLOWED IN SECTOR
        bool_x = (x_arr >=bounds[b]) & (x_arr<bounds[b+1])
        x_in = x_arr[bool_x].tolist()
        y_in = y_arr[bool_x].tolist()
        inds = np.random.randint(0,len(x_in),size=nq)
        px_l.extend( [x_in[inds[i]] for i in range(nq)])
        py_l.extend([y_in[inds[i]] for i in range(nq)])
    
    return np.asarray(px_l), np.asarray(py_l)
    ###########################################
 



def check_fig_seed(nc_grd,px,py,c_contour=[]):
    ####################################
    # SHOW USER FLOAT LOCATIONS ON GRID
    ####################################
    h = nc_grd.variables['h'][:,:]
    xi_rho=nc_grd.variables['xi_rho'][:,:]
    eta_rho=nc_grd.variables['eta_rho'][:,:]
    mask = nc_grd.variables['mask_rho'][:,:]
    mask_rho = np.copy(mask)
    mask_rho[mask_rho==0] = np.nan
    plt.ion() 
    fig_check = plt.figure(figsize=[24,8])
    if len(c_contour)==0:
       h_con = plt.contour(xi_rho,eta_rho,h,colors='k',linewidths=2.5)
    else:
       h_con = plt.contour(xi_rho,eta_rho,h,c_contour,colors='k',linewidth=2.5)
    plt.contourf(xi_rho,eta_rho,h*mask_rho,origin='lower',cmap=cm.gist_earth_r)
    cbar_plot = plt.colorbar()
    cbar_plot.set_label(r'$h(m)$',fontsize=22)
    plt.plot(px,py, 'o', color='slategrey',markersize=2.5)
    plt.title('INITIAL PARTICLE LOCATIONS')
    ax = plt.gca()
    # ax.set_xlim([0,mask_rho.shape[1]])
    # ax.set_ylim([0,mask_rho.shape[0]])
    move_on = raw_input('Press enter to continue...')
    plt.ioff()
    plt.close('all')
    #################################################################################################


def radial_feature_seed(nc_grd,nc_out,nc_tind,field_look='div',radial_inputs=[]):
    """
    SEED PARTICLES RADIALLY AROUND A 
    FEATURE (i.e., front or filament)

    USES GUI-INTERACTIVE PLOTTING WHERE USER
    CLICKS ON FEATURE AND INPUTS RADIUS OF PARTICLE
    CIRCLE AND NUMBER OF PARTICLES TO SEED WITHIN
    THAT RADIAL DISTRIBUTION
    
   
    nc_grd --> NETCDF object of grid-file

    nc_out --> NETCDF object of file to use to look at
	       2D map of fields to place particles
 
    nc_tind --> time-step w/in netcdf file to look at fields (ZERO-BASED INDEXING)

    field_look --> 'temp', 'salt', 'div', 'vort'
                    field to look at to place particles
		    (default is surface divergence), all other
		    options are for surface fields



    radial_inputs---> list containing 2 components
                      [r_km, nparticles]
		        radius, number of particles

			defualt is [] (empty) and in default
			mode, the user will be prompted to enter
			in these values
			
    """

    #############################################
    # 		LOAD GRID DATA
    ##############################################
    mask = nc_grd.variables['mask_rho'][:,:]
    mask_rho = np.copy(mask)
    mask_rho[mask_rho==0] = np.nan
    pm       = nc_grd.variables['pm'][:,:]
    pn       = nc_grd.variables['pn'][:,:]
    h        = nc_grd.variables['h'][:,:]
    f        = nc_grd.variables['f'][:,:]

    #############################################
    # LOAD OUTPUT DATA TO VISUALIZE FOR SEEDING
    #############################################
    if field_look == 'div':
       try:
          div = calc_div(nc_out.variables['u'][nc_tind,-1,:,:],nc_out.variables['v'][nc_tind,-1,:,:],pm,pn,f)
       except:
	  div = calc_div(nc_out.variables['ubar'][nc_tind,:,:],nc_out.variables['vbar'][nc_tind,:,:],pm,pn,f)
       
       field_plot = (div * mask_rho) / f
       cmap_plot = cm.seismic
       cbar_label = r'$\delta / f$'
       min_max = [-5,5]

    if field_look == 'vort':
      try:
         vort = calc_vort(nc_out.variables['u'][nc_tind,-1,:,:],nc_out.variables['v'][nc_tind,-1,:,:],pm,pn,f)
      except:
          vort = calc_vort(nc_out.variables['ubar'][nc_tind,:,:],nc_out.variables['vbar'][nc_tind,:,:],pm,pn,f)
       
      field_plot = (psi2rho(vort) * mask_rho) / f
      cmap_plot = cm.seismic
      cbar_label = r'$\zeta / f$'
      min_max = [-10,10]      

    if field_look == 'temp':
       field_plot = nc_out.variables['temp'][nc_tind,-1,:,:] * mask_rho 
       cmap_plot = cm.rainbow
       cbar_label = r'$SST \left(^{\circ} C\right)'
       min_max = [np.min(field_plot), np.max(field_plot)]


    if field_look == 'salt':
       field_plot = nc_out.variables['salt'][nc_tind,-1,:,:] *mask_rho
       cmap_plot = cm.rainbow
       cbar_label = r'$SSS'
       min_max = [np.min(field_plot), np.max(field_plot)]

   
    ##############################################
    # PLOT AND PROMPT USER TO CLICK ON PLOT 
    
    '''
    ALL OF THIS HAPPENS IN WHILE LOOP
    UNTIL USER IS SATISFIED WITH PARTICLE DISTRIBUTION
    SO USER CAN RE-SELECT FEATURE, RADIUS, AND # OF PARTICLES
    UNTIL THEY HAVE DISTRIBUTION THEY WANT
    '''
    ##############################################
    choose_again = 'Y'
    while choose_again == 'Y':
         plt.ion()
         plt.figure(figsize=[14,7])
         h_con = plt.contour(h,colors='k',linewidths=2.5)
         plt.imshow(field_plot,origin='lower',cmap=cmap_plot,vmin=min_max[0],vmax=min_max[1])
         cbar_plot = plt.colorbar()
         cbar_plot.set_label(cbar_label,fontsize=22)
         plt.title('INSTRUCTIONS FOR CLICKING: ' + '\n' + 'Select point: left click' + '\n' + 'Cancel last input: right click' + '\n' + 'Stop click recording: middle click')
         coords = plt.ginput(n=-1,show_clicks=True,mouse_add=1,mouse_pop=3,mouse_stop=2)

         #######################################################
         # ASK USER HOW MANY PARTICLES AND WHAT RADIUS THEY WANT
         ########################################################
	 if radial_inputs==[]:
	    print '################################################'
	    print ''
            print ''
            print ''
            r_km                = input('Radius of circle (in km)?  ')
            print ''
	    nparticles = input('How many total particles?  ')

	    print '##################################################'
            #CONVERT RADIUS TO GRID-POINT UNITS
         
         else:
             r_km = radial_inputs[0]
	     nparticles = radial_inputs[1]


	 r_grdpt = (r_km *1E3) * np.mean((pm+pn)/2)
    
         ###############################################
         # DISTRIBUTE PARTICLES INSIDE CRICLE
         ################################################
         x_center = coords[0][0]
         y_center = coords[0][1]

         px_list = []
         py_list = []
         for n in range(nparticles):
	     q = np.random.random() * 2 * np.pi
	     r = np.sqrt(np.random.random())
	     px_list.append(x_center + (r_grdpt * r) * np.cos(q))
	     py_list.append(y_center + (r_grdpt * r) * np.sin(q))

     
         px_arr =  np.asarray(px_list,dtype='float64') 
         py_arr = np.asarray(py_list,dtype='float64') 

         ####################################
         # SHOW USER FLOAT LOCATIONS ON GRID
         ####################################
         fig_check = plt.figure(figsize=[14,7])
         h_con = plt.contour(h,colors='k',linewidths=2.5)
         plt.imshow(field_plot,origin='lower',cmap=cmap_plot,vmin=min_max[0],vmax=min_max[1])
         cbar_plot = plt.colorbar()
         cbar_plot.set_label(cbar_label,fontsize=22)
         plt.plot(px_arr,py_arr, 'o', color='k',markersize=2.5)
         ax = plt.gca()
         ax.set_xlim([0,mask_rho.shape[1]])
         ax.set_ylim([0,mask_rho.shape[0]])


         ##########################################
         # PROMPT USER TO ENTER INPUTS AGAIN
	 ##########################################
         user_in_bad = True    
	 while user_in_bad:
	     temp_in = raw_input('Choose feature / particle distribution again (Y/N)? ' )
	     if temp_in == 'Y' or temp_in == 'N':
		choose_again = temp_in
		user_in_bad = False
         plt.close('all')
	 ##########################################################################################


    return px_arr-1,py_arr-1
    ############################################################################



	####################################################################
	#	   MISC ROMS CALCULATIONS
        #####################################################################

def calc_div(ubar,vbar,pm,pn,f):
    """
    Calculate 2D divergence
    ubar, vbar --> [neta,nxi] arrays
    pm,pn ---> horizontal spacing [neta,nxi]
    f     ---> [neta,nxi] coriolis parameter
    """
    [neta,nxi] = ubar.shape
    div = np.zeros([neta,nxi])

    A = (pm * pn)
    temp1, dxi_div   = np.gradient( (1/pn) * u2rho(ubar))
    deta_div, temp2  = np.gradient( (1/pm) * v2rho(vbar))
    div_horiz = A * (dxi_div + deta_div)

    return div_horiz
    ############################################################


def calc_vort(ubar, vbar, pm, pn, f):
    """
    Calculate 2D divergence
    ubar, vbar --> [neta,nxi] arrays
    pm,pn ---> horizontal spacing [neta,nxi]
    f     ---> [neta,nxi] coriolis parameter
    """
    [neta,nxi] = pm.shape  
    d1 = np.zeros([neta-1,nxi-1])
    d2 = np.zeros([neta-1,nxi-1])
    d1[:,:] = (rho2v(1/pn)[:,1:]*vbar[:,1:]) - (rho2v(1/pn)[:,:-1] * vbar[:,:-1])
    d2[:,:] = (rho2u(1/pm)[1:,:] *ubar[1:,:]) - (rho2u(1/pm)[:-1,:] * ubar[:-1,:])
    xi = rho2psi(pm*pn) * (d1 -d2)
    return xi
    ##############################################################################

#######################################################
#Transfert a field at psi points to rho points
#######################################################

def psi2rho(var_psi):

    if np.rank(var_psi)<3:
        var_rho = psi2rho_2d(var_psi)
    else:
        var_rho = psi2rho_3d(var_psi)

    return var_rho


##############################

def psi2rho_2d(var_psi):

    [M,L]=var_psi.shape
    Mp=M+1
    Lp=L+1
    Mm=M-1
    Lm=L-1

    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,1:L]=0.25*(var_psi[0:Mm,0:Lm]+var_psi[0:Mm,1:L]+var_psi[1:M,0:Lm]+var_psi[1:M,1:L])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]

    return var_rho

#############################

def psi2rho_3d(var_psi):


    [Nz,Mz,Lz]=var_psi.shape
    var_rho=np.zeros((Nz,Mz+1,Lz+1))

    for iz in range(0, Nz, 1):    
        var_rho[iz,:,:]=psi2rho_2d(var_psi[iz,:,:])


    return var_rho

#######################################################
#Transfert a field at rho points to psi points
#######################################################

def rho2psi(var_rho):

    if np.rank(var_rho)<3:
        var_psi = rho2psi_2d(var_rho)
    else:
        var_psi = rho2psi_3d(var_rho)

    return var_psi


##############################

def rho2psi_2d(var_rho):

    var_psi = 0.25*(var_rho[1:,1:]+var_rho[1:,:-1]+var_rho[:-1,:-1]+var_rho[:-1,1:])

    return var_psi

#############################

def rho2psi_3d(var_rho):

    var_psi = 0.25*(var_rho[:,1:,1:]+var_rho[:,1:,:-1]+var_rho[:,:-1,:-1]+var_rho[:,:-1,1:])

    return var_psi


#######################################################
#Transfert a 2 or 3-D field at rho points to u points
#######################################################

def rho2u(var_rho):

    if np.rank(var_rho)<3:
        var_u = rho2u_2d(var_rho)
    else:
        var_u = rho2u_3d(var_rho)

    return var_u

def rho2u_2d(var_rho):

    [Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,0:L]+var_rho[:,1:Lp])

    return var_u


def rho2u_3d(var_rho):

    [N,Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,:,0:L]+var_rho[:,:,1:Lp])

    return var_u

#######################################################
#Transfert a 3-D field at rho points to v points
#######################################################

def rho2v(var_rho):

    if np.rank(var_rho)<3:
        var_v = rho2v_2d(var_rho)
    else:
        var_v = rho2v_3d(var_rho)

    return var_v

#######################################################

def rho2v_2d(var_rho):

    [Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[0:M,:]+var_rho[1:Mp,:]);

    return var_v

#######################################################

def rho2v_3d(var_rho):

    [N,Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[:,0:M,:]+var_rho[:,1:Mp,:]);

    return var_v


#######################################################
#Transfert a 2-D field at u points to the rho points
#######################################################

def u2rho(var_u):


    if np.rank(var_u)<3:
        var_rho = u2rho_2d(var_u)
    else:
        var_rho = u2rho_3d(var_u)

    return var_rho

#######################################################

def u2rho_2d(var_u):

    [Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[:,1:L]=0.5*(var_u[:,0:Lm]+var_u[:,1:L])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]
    return var_rho

#######################################################

def u2rho_3d(var_u):

    [N,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,:,1:L]=0.5*(var_u[:,:,0:Lm]+var_u[:,:,1:L])
    var_rho[:,:,0]=var_rho[:,:,1]
    var_rho[:,:,Lp-1]=var_rho[:,:,L-1]
    return var_rho


#######################################################
#Transfert a 2 or 2-D field at v points to the rho points
#######################################################

def v2rho(var_v):

    if np.rank(var_v)<3:
        var_rho = v2rho_2d(var_v)
    else:
        var_rho = v2rho_3d(var_v)

    return var_rho

#######################################################

def v2rho_2d(var_v):

    [M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,:]=0.5*(var_v[0:Mm,:]+var_v[1:M,:])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]

    return var_rho

#######################################################

def v2rho_3d(var_v):

    [N,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,1:M,:]=0.5*(var_v[:,0:Mm,:]+var_v[:,1:M,:])
    var_rho[:,0,:]=var_rho[:,1,:]
    var_rho[:,Mp-1,:]=var_rho[:,M-1,:]

    return var_rho




