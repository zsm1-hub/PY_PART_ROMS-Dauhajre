land_masks = PF.load_pickle_file('L3_land_masks.p')
CI_masks   = PF.load_pickle_file('CI_sites.p')

land_m = land_masks['mask_sites']
CI_m   = CI_masks['mask_sites']

land_cents = land_masks['cents']
CI_cents  = CI_masks['cents']

#mask_all = land_m[::-1] + CI_m

#RE-ARRANGING ORDERS
ind_ords_new = [8,6,7,0,1,2,3,4,5,15,9,10,11,12,18,17,13,14,16]
CI_m_new = [CI_m[n] for n in ind_ords_new]
mask_all = land_m[::-1] + CI_m_new
#mask_all = CI_m_new


import sys
sys.path.append('/home/dauhajre/MidCal/py_ana/py_modules/')
import custom_cmap as cm_cust
import matplotlib.pyplot as plt
h = nc_grd.variables['h'][:,:]
mask_rho = nc_grd.variables['mask_rho'][:,:]
plt.ion()
plt.figure()
mask_nan = np.copy(mask_rho)
mask_nan[mask_nan==0.] = np.nan
#plt.contour(h*mask_nan,[1.5],origin='lower',colors='k',linewidths=2)
plt.imshow(h*mask_nan,origin='lower',cmap=plt.cm.gist_earth)
#plt.imshow(mask_nan,origin='lower')
col_sites = cm_cust.make_cbar_even_space(cm_cust.my_cmaps('my_jet'),len(mask_all))
[plt.contour(mask_all[n],colors=[col_sites[n]],linewidths=3) for n in range(len(mask_all))]
#[plt.contour(mask_all[n],colors='slategrey',linewidths=3) for n in range(len(mask_all))]

plt.contour(h_smth,[15],colors='slategrey',linewidths=5)

#for n in range(len(mask_all)):
#    CS =plt.contour(mask_all[n],colors=[col_sites[n]],linewidths=3)
#    print '     Site index = ' + str(n)
#    move_one = raw_input('press enter to continue')



nq_site = 100
px_l = []
py_l = []
n_sites = len(mask_all)
for n in range(n_sites):
    y,x = np.where(mask_all[n]*mask_rho==1)
    inds = np.random.randint(0,len(x),size=nq_site)
    px_l.extend([x[inds[i]] for i in range(nq_site)])
    py_l.extend([y[inds[i]] for i in range(nq_site)])

plt.figure()
plt.imshow(h*mask_nan,origin='lower',cmap=plt.cm.gist_earth)
col_sites = cm_cust.make_cbar_even_space(cm_cust.my_cmaps('my_jet'),len(mask_all))
[plt.contour(mask_all[n],colors=[col_sites[n]],linewidths=3) for n in range(len(mask_all))]
plt.plot(px_l,py_l,'o',color='slategrey')




'''
#SAVE MASKS
save_dict = {}
save_dict['mask_sites'] = mask_all
PF.save_to_pickle(save_dict, 'L3_CC_masks')
'''
