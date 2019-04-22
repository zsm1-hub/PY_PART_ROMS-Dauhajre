import matplotlib.pyplot as plt
h = nc_grd.variables['h'][:,:]
from scipy import ndimage
#h=ndimage.filters.gaussian_filter(nc_grd.variables['h'][:,:],[5,5])
h_seed = 1.5
ns = 30
plt.ioff()
CD = plt.contour(h, [h_seed])
len_cons = [len(CD.collections[0].get_paths()[i]) for i in range(len(CD.collections[0].get_paths()))]
ind_max_len = np.asarray(len_cons).argmax()

p = CD.collections[0].get_paths()[ind_max_len]
x_iso = p.vertices[:,0]
y_iso = p.vertices[:,1]
dpt = int(len(x_iso)/ns)
x_site = x_iso[dpt:-1:dpt]
y_site = y_iso[dpt:-1:dpt]
cents = [[y_site[n],x_site[n]] for n in range(len(x_site))]
r_sites = 50#gridpoints
mask_rho = nc_grd.variables['mask_rho'][:,:]
#USE 1.5 m isobath as land mask
h_smth =ndimage.filters.gaussian_filter(nc_grd.variables['h'][:,:],[5,5])
h_mask = np.zeros(mask_rho.shape)
h_mask[:,:] = 1.
h_mask[h_smth<=1.8] = 0.
mask_sites = PF.mask_sites_circles(cents,r_sites,h_mask)
plt.ion()
plt.figure()
mask_nan = np.copy(mask_rho)
mask_nan[mask_nan==0.] = np.nan
plt.imshow(h*mask_nan,origin='lower',cmap=plt.cm.jet)
[plt.plot(x_site[n],y_site[n],'o',color='g') for n in range(len(x_site))]
[plt.contour(mask_sites[n],colors='k',linewidths=3) for n in range(len(mask_sites))]

#RECURSIVE RELATION THAT AVOIDS SITE OVERLAPS
import copy
mask_sites_new = copy.deepcopy(mask_sites)
list_all = range(len(mask_sites))
for n in range(len(mask_sites)):
    list_others = copy.copy(list_all)
    list_others.pop(n)
    sum_others = np.copy(mask_sites_new[n])
    for k in list_others:
        sum_others+=mask_sites_new[k]
    mask_sites_new[n][sum_others>1.] = 0.


'''
mask_sites_new = copy.copy(mask_sites) 
for n in range(len(mask_sites)-1):
    sum_site = mask_sites[n] + mask_sites[n+1]
    mask_sites_new[n+1][sum_site==2.] = 0

mask_sites_nn = copy.copy(mask_sites_new) 
for n in range(len(mask_sites)-2):
    sum_site = mask_sites_new[n] + mask_sites_new[n+1] + mask_sites_new[n+2]
    mask_sites_nn[n+2][sum_site==2.] = 0.
'''




[plt.contour(mask_sites_new[n],colors='grey',linewidths=3) for n in range(len(mask_sites))]
#[plt.contour(mask_sites_nn[n],colors='grey',linewidths=3) for n in range(len(mask_sites))]



[Ly,Lx] = mask_nan.shape
mask_sum = np.zeros([Ly,Lx])
mask_sum_nn = np.zeros([Ly,Lx])
for n in range(len(mask_sites)):
    mask_sum = mask_sum + mask_sites_new[n]
    #mask_sum_nn = mask_sum_nn + mask_sites_nn[n]

plt.figure()
#plt.subplot(2,1,1)
plt.imshow(mask_sum,origin='lower')
plt.colorbar()
#plt.subplot(2,1,2)
#plt.imshow(mask_sum_nn,origin='lower')
#plt.colorbar()

