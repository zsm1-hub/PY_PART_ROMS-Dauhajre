##########################################

#  ROMS_solutions_paths.py 
'''
Streamline setting of solution paths
by calling function in this module
that returns paths / start_days / dx / etc.
for any set of model solutions

Add in a solution set every time with unique ROMS_run_ID


CLASS STRUCTURE



Also has 2 functions that set time conventions
and set UTC and PST hours for a given ocean_time
'''

# DANIEL DAUHAJRE UCLA 2015

##########################################


#######################################

import os
import sys
import numpy as np
import datetime 
import matplotlib.dates as dates
import pytz


########################################
class ROMS_run(object):
    """ Represents a set of solutions """

    def __init__(self, ID_in):
      self.ROMS_run_ID = ID_in
      self.pst_hours   = []
      self.path_output = []
      self.path_grid   = []
      self.path_input  = []
      self.out_base_name = []
      self.grid_name     = []
      self.wind_name     = []
      self.hflux_name    = []
      self.swflux_name   = []
      #self.model_start_day = []
      #self.model_dx        = []
      #self.nc_tstep        = []
      self.utc_hour_start  = []

    def __repr__(self):
      return str(self.__dict__)




    ##################################
    # set_paths()

    '''
    Set paths for a specific ROMS run
    '''
    #################################
    def set_paths(self):
        if self.ROMS_run_ID == 'L4_MidCal_SBC':

           self.path_output   = '/hun/dauhajre/Simulations/MidCal/L4/NO_WEC/Output/'
           self.path_grid = '/hun/dauhajre/Simulations/MidCal/L4/NO_WEC/Input/'
           #self.path_input = '/mnt/inca/akanc/MidCal/L3/02_input/winter2006/'
           self.grid_name     = 'usw4sbc_grd.nc'
           self.out_base_name = 'usw4sbc_avg_3sec.'
           
           self.date_origin = datetime.datetime(1994,1,1,0,0,0,0)
           self.dx        = 36
	   self.nc_tstep        = 24

        if self.ROMS_run_ID == 'L3_MidCal_30min':

           self.path_output   = '/hun/dauhajre/Simulations/MidCal/L3/NO_WEC/Output_30min/'
           self.path_grid = '/hun/dauhajre/Simulations/MidCal/L3/NO_WEC/Input/'
           #self.path_input = '/mnt/inca/akanc/MidCal/L3/02_input/winter2006/'
           self.grid_name     = 'usw3_grd.nc'
           self.out_base_name = 'usw3_avg.'
           
           self.date_origin = datetime.datetime(1994,1,1,0,0,0,0)
           self.dx        = 100
	   self.nc_tstep        = 24

        if self.ROMS_run_ID == 'L2_MidCal_30min':

           self.path_output   = '/mnt/shirma/dauhajre/MidCal/L2_nowec/Output_30min/'
           self.path_grid = '/mnt/shirma/dauhajre/MidCal/L2_nowec/Input/'
           #self.path_input = '/mnt/inca/akanc/MidCal/L3/02_input/winter2006/'
           self.grid_name     = 'usw2_grd_DD.nc'
           self.out_base_name = 'usw2_avg_DPD.'
           
           self.date_origin = datetime.datetime(1994,1,1,0,0,0,0)
           self.dx        = 300
	   self.nc_tstep        = 24


 
        if self.ROMS_run_ID == 'L3_MidCal_Winter2006':

           self.path_output   = '/mnt/inca/akanc/MidCal/L3/03_out_winter2006/'
           self.path_grid = self.path_output
           self.path_input = '/mnt/inca/akanc/MidCal/L3/02_input/winter2006/'
           self.grid_name     = 'usw3_grd.nc'
           self.out_base_name = 'usw3_avg.'
           
           self.date_origin = datetime.datetime(1994,1,1,0,0,0,0)
           #############################
           # netcdf ID number paramters
	   ###########################
	   # Model start day (in model day convetions)
	   self.start_day = 4711
	   
	   # index in first netcdf file of real data
	   # this takes into account zero padding
	   self.nc_step_start = 0

	   #number ID of first netcdf FILE
           self.nc_ID_initial = 0

           ###################################
           # OTHER PARAMTERS
	   ##################################
	   self.dx        = 100
	   self.nc_tstep        = 24
	   self.ROMS_run_parent_set = 'MidCal'
        

        if self.ROMS_run_ID == 'L4_SBC':
           
	   self.path_output   = '/mnt/kamiya/dauhajre/Simulations/MidCal/usw4_sbc/his/'
           self.path_ts_uv = '/mnt/kamiya/dauhajre/Simulations/MidCal/usw4_sbc/ts_uv/'
	   self.path_grid = '/mnt/kamiya/dauhajre/Simulations/MidCal/usw4_sbc/input/'
           self.path_omega = '/data/project3/ddauhajre/MidCal/py_particles/omega_calc_offline/L4_SBC_omega/' 
           self.path_zs    = '/data/project3/ddauhajre/MidCal/usw4_sbc/L4_SBC_zs/'
           self.path_input = self.path_grid 
	   self.out_base_name = 'usw4sbc_his.'
	   self.grid_name     = 'usw4sbc_grd.nc'
           self.wind_name     = 'usw4sbc_wnd.nc'
	   self.hflux_name    = 'usw4sbc_rad.nc'
           self.omega_name    = 'usw4sbc_his_omega.'
           self.zs_name        = 'usw4sbc_his_zs.'
           self.date_origin = datetime.datetime(2004,1,1,0,0,0,0)
           #############################
           # netcdf ID number paramters
	   ###########################
	   # Model start day (in model day convetions)
	   self.start_day = 4711
	   
	   # index in first netcdf file of real data
	   # this takes into account zero padding
	   self.nc_step_start = 0

	   #number ID of first netcdf FILE
           self.nc_ID_initial = 0

           ###################################
           # OTHER PARAMTERS
	   ##################################
	   self.dx        = 36
	   self.nc_tstep        = 48
	   self.ROMS_run_parent_set = 'MidCal'
        
           self.dx        = 100
	   self.nc_tstep        = 24
	   self.ROMS_run_parent_set = 'MidCal'
        #####################
        # change by zsm
        #######################
        if self.ROMS_run_ID == 'child_':
           
	   self.path_output   = '/leader/user/zsm/TWS/res/hourly_child/'
           self.path_ts_uv = '/leader/user/zsm/TWS/res/hourly_child/'
	   self.path_grid = '/leader/user/zsm/TWS/files/'
           self.path_omega = '/leader/user/zsm/TWS/res/hourly_child/' 
           self.path_zs    = '/leader/user/zsm/TWS/res/hourly_child/'
           self.path_input = self.path_grid 
	   self.out_base_name = 'child_his_'
	   self.grid_name     = 'chinav2_ATWS1.nc.1'
           self.wind_name     = 'usw4sbc_wnd.nc'
	   self.hflux_name    = 'usw4sbc_rad.nc'
           self.omega_name    = 'child_his_'
           self.zs_name        = 'child_his_'
           self.date_origin = datetime.datetime(1900,1,1,0,0,0,0)
           #############################
           # netcdf ID number paramters
	   ###########################
	   # Model start day (in model day convetions)
	   self.start_day = 4711
	   
	   # index in first netcdf file of real data
	   # this takes into account zero padding
	   self.nc_step_start = 0

	   #number ID of first netcdf FILE
           self.nc_ID_initial = 0



        ##############################################tide
        if self.ROMS_run_ID == 'tide':
           
	   self.path_output   = '/leader/user/zsm/TWS2_400600/TIDE/U3/'
           self.path_ts_uv = '/leader/user/zsm/TWS2_400600/TIDE/U3/'
	   self.path_grid = '/sugon7/zsm/croco_tools_xmd1204/CROCO_FILES/TWS2/'
           self.path_omega = '/leader/user/zsm/TWS2_400600/TIDE/U3/' 
           self.path_zs    = '/leader/user/zsm/TWS2_400600/TIDE/U3/'
           self.path_input = self.path_grid 
	   self.out_base_name = 'CCL1_his_'
	   self.grid_name     = 'TWS2_rot2.nc'
           self.wind_name     = 'usw4sbc_wnd.nc'
	   self.hflux_name    = 'usw4sbc_rad.nc'
           self.omega_name    = 'CCL1_his_'
           self.zs_name        = 'CCL1_his_'
           self.date_origin = datetime.datetime(1900,1,1,0,0,0,0)
        ############################################notide
        if self.ROMS_run_ID == 'notide':
           
	   self.path_output   = '/leader/user/zsm/TWS2_400600/NOTIDE/'
           self.path_ts_uv = '/leader/user/zsm/TWS2_400600/NOTIDE/'
	   self.path_grid = '/sugon7/zsm/croco_tools_xmd1204/CROCO_FILES/TWS2/'
           self.path_omega = '/leader/user/zsm/TWS2_400600/NOTIDE/' 
           self.path_zs    = '/leader/user/zsm/TWS2_400600/NOTIDE/'
           self.path_input = self.path_grid 
	   self.out_base_name = 'CCL1_his.nc.not_'
	   self.grid_name     = 'TWS2_rot2.nc'
           self.wind_name     = 'usw4sbc_wnd.nc'
	   self.hflux_name    = 'usw4sbc_rad.nc'
           self.omega_name    = 'CCL1_his.nc.not_'
           self.zs_name        = 'CCL1_his.nc.not_'
           self.date_origin = datetime.datetime(1900,1,1,0,0,0,0)
           #############################
           # netcdf ID number paramters
	   ###########################
	   # Model start day (in model day convetions)
	   self.start_day = 4711
	   
	   # index in first netcdf file of real data
	   # this takes into account zero padding
	   self.nc_step_start = 0

	   #number ID of first netcdf FILE
           self.nc_ID_initial = 0

           ###################################
           # OTHER PARAMTERS
	   ##################################
	   self.dx        = 36
	   self.nc_tstep        = 48
	   self.ROMS_run_parent_set = 'MidCal'
       
 
    ###########################################
    # RETURN netcdf file ID for a given model day
    # and time step in the sequence of netcdf files

    '''

    INPUTS
    model_day ---> day in model day units
    t_step ---. starts at 0, ends at nc_tstep -1
    
    
    OUTPUTS
    nc_file_ID --> the number corresponding the netcdf number ID (to acces file like '%04g' %nc_file_ID)
    t_step_in_file ---> the index within that netcdf file to get that specific time point
    '''

    ###########################################
    def get_nc_file_ID(self, model_day, t_step):

	if t_step <0 or t_step >= self.nc_tstep:
           print 'ERROR: t_step out of bounds'
	   return
        
        # get netcdf ID number for the START of that day (so if netcdf files have
	# days that overlap within each file, this will account for that)
	day_look = model_day - self.start_day
	nc_index = (day_look * self.nc_tstep) + (self.nc_ID_initial + self.nc_step_start + t_step)
	
	# now get the specific index of the time step, to know which netcdf file to get (if days overalp within netcdf file)
	nc_file_ID = (nc_index / self.nc_tstep) * self.nc_tstep
        t_step_in_file = nc_index - nc_file_ID
        

	return nc_file_ID, t_step_in_file


    ######################################
    # GET START INDEX FOR WRF forcing index
    # for a parent set of solutions
    # since higher time resolution
    # solutions could use same forcing
    # files sampled at same rate


    '''
    FIND INDEX WHERE sms_time is closest to ocean_time
    and return index in sms_time to access forcing
    that is phase aligned with ocean_time w/out 
    interpolating in time

    ocean_time --> in seconds, actual number, not array at all points
    sms_time ---> array of forcing time
    '''
    #####################################
    def get_forcing_time(self, ocean_time, sms_time):
        return np.abs(sms_time[:] - (ocean_time/86400.)).argmin()




    ################################
    # USING DATE ORIGIN AS REFERENCE
    # RETURN FULL UTC DATE (YEAR,MONTH,DAY,HOUR,MIN,SEC)
    ##################################
    def round_minutes(self,dt,direction,resolution):
        new_minute = (dt.minute//resolution + (1 if direction == 'up' else 0)) * resolution
        return dt + datetime.timedelta(minutes=new_minute-dt.minute)


    def get_UTC_date(self,secs):
	stt = self.date_origin.toordinal()
	tcoef = 86400.
	n1 = dates.num2date(secs/tcoef + stt)
	return n1

    def get_PST_date(self,secs):
	stt = self.date_origin.toordinal()
	tcoef = 86400.
	n1 = dates.num2date(secs/tcoef + stt)
	return n1.astimezone(pytz.timezone('US/Pacific'))
    def get_PST_date_rounded(self,secs):
	stt = self.date_origin.toordinal()
	tcoef = 86400.
	n1 = dates.num2date(secs/tcoef + stt)
	return self.round_minutes(n1.astimezone(pytz.timezone('US/Pacific')),'up',1)




    ######################################
    # RETURN A STRING OF UTC AND PST HOUR
    # for plotting purposes

    # ocean_time --> array of size model_nc_step
    
    '''
    takes time passed since model day
    and converts this time to separate
    UTC and PST hour strings
    '''
    ######################################
    def get_utc_pst_hour(self, ocean_time):
        
	#####################################
	# set initial UTC hour
	#####################################
	if self.ROMS_run_parent_set == 'L4PV':
           utc_initial = 0
        
        

        #######################################################
        #  GET DAY, HOUR, MINUTE FOR OCEAN_TIME OF A FULL DAY
	###################################################

	time_passed_days = np.zeros([len(ocean_time)])
        ocean_time_days = ocean_time / 86400.
	model_day = np.zeros([len(ocean_time)])
        for i in range(len(ocean_time)):
	    model_day[i] = int(ocean_time_days[i])
	    time_passed_days[i] = (ocean_time_days[i] - int(ocean_time_days[i]))
           
	minutes   = np.round(time_passed_days *24 *60, decimals =0)
        
	utc_hours_temp = np.zeros([len(ocean_time)])
	for i in range(len(ocean_time)):
	    utc_hours_temp[i] = int(minutes[i] / 60.)
        
	# CONVERT utc_hours = 0 TO 24

        utc_hours = np.copy(utc_hours_temp)
	for i in range(len(utc_hours_temp)):
	    if utc_hours_temp[i] == 0.0:
	       utc_hours[i] = 24.


        minutes_of_hour = np.zeros([len(minutes)])
	for i in range(len(minutes)):
	    if minutes[i] > 60.:
               minutes_of_hour[i] = (minutes[i] - (60 *utc_hours[i-1])) %60
            else:
	       minutes_of_hour[i] = minutes[i] % 60

        self.minute_arr = minutes_of_hour
	self.utc_hour_arr = utc_hours
        self.model_day_arr = model_day
        


	##### CONVERT HOURS FROM UTC TO PST ####
        ###### UTC 8 HOURS AHEAD OF PST ######

        utc_day = [1, 2, 3, 4, 5, 6, 7, 8 , 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]        
        pst_day = [17, 18, 19, 20, 21, 22, 23, 24, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14,15,16,17]
        
	pst_hour = np.zeros([len(ocean_time)])
	for i in range(len(ocean_time)):
	    hour = self.utc_hour_arr[i]
	    
	    ind = np.where(utc_day == hour)[0][0]

	    pst_hour[i] = pst_day[ind]


        self.pst_hour_arr = pst_hour

        
	###################################
        # CONVERT HOUR/MINUTE TO A STRING 
	# IN A LIST 
	##################################
       
        pst_hour_list = []
        for i in range(len(ocean_time)):
	    pst_hour_i = int(self.pst_hour_arr[i])
	    minute_i   = int(self.minute_arr[i])
	    
	    if pst_hour_i < 10:
	       if minute_i < 10:
		  pst_hour_list.append('0' + str(pst_hour_i) + ':0' + str(minute_i))
    
	       if minute_i >=10:
		  pst_hour_list.append('0' + str(pst_hour_i) + ':' + str(minute_i))

	    if pst_hour_i >= 10:
	       if minute_i <10:
		  pst_hour_list.append(str(pst_hour_i) + ':0' + str(minute_i))
	
	       if minute_i >=10:
		  pst_hour_list.append(str(pst_hour_i) + ':' + str(minute_i))


        #############################################################
        # SET LIST OF STRINGS FOR PST HOURS FOR ROMS ocean_time array
	#############################################################
        self.pst_hour_str_list = pst_hour_list





