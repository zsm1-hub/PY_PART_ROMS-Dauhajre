import matplotlib.pyplot as plt
h = nc_grd.variables['h'][:,:]
from scipy import ndimage
#h=ndimage.filters.gaussian_filter(nc_grd.variables['h'][:,:],[5,5])

#mask_rho = nc_grd.variables['mask_rho'][:,:]
h_seed = 2
ns =4
#plt.ioff()
CD = plt.contour(h, [h_seed],colors='k')
len_cons = [len(CD.collections[0].get_paths()[i]) for i in range(len(CD.collections[0].get_paths()))]
ind_max_len = np.asarray(len_cons).argmax()
#INDEX 3 = CI 1, 6 big one

#ISLAND 1
ns =3 
p = CD.collections[0].get_paths()[1]
#p = CD.collections[0].get_paths()[ind_max_len]
x_iso = p.vertices[:,0]
y_iso = p.vertices[:,1]
dpt = int(len(x_iso)/ns)
x_site_1 = x_iso[dpt:-1:dpt].tolist()
y_site_1 = y_iso[dpt:-1:dpt].tolist()
x_site_1.append(852.051)
y_site_1.append(137.684)


#ISLAND 2
ns = 6 
p = CD.collections[0].get_paths()[3]
#p = CD.collections[0].get_paths()[ind_max_len]
x_iso = p.vertices[:,0]
y_iso = p.vertices[:,1]
dpt = int(len(x_iso)/ns)
x_site_2 = x_iso[dpt:-1:dpt].tolist()
y_site_2 = y_iso[dpt:-1:dpt].tolist()
x_site_2.pop(-1)
y_site_2.pop(-1)
x_site_2.append(1084.59)
y_site_2.append(174.204)


#ISLAND 3 
ns =8
p = CD.collections[0].get_paths()[6]
#p = CD.collections[0].get_paths()[ind_max_len]
x_iso = p.vertices[:,0]
y_iso = p.vertices[:,1]
#dpt = int(len(x_iso)/ns)
dpt = 148
x_site_3 = x_iso[dpt:-1:dpt].tolist()
y_site_3 = y_iso[dpt:-1:dpt].tolist()



##SMALL LAST ISLAND
#x_site_3 = [1164.07,1410.,     1204.92   ,1347.0780316440523,1399.1332992098774]
#y_site_3 = [407.385,559.01464, 373.41,         442.0,              497.]
#hardcode changes
#x_site_3[3] = 1373.55
#y_site_3[3] = 568.914
x_site_3[3] =1371 
y_site_3[3] = 544

x_site_3[2] = 1320.02
y_site_3[2] = 492.41

x_site_3[6]=1277.78
y_site_3[6]=374.284

x_site_3.pop(4)
y_site_3.pop(4)


#ONE LAST SITE
x_site_3.append(1212.7)
y_site_3.append(368.82)

#SMALL LAST ISLAND
#x_site_3.append(1500.25)
#y_site_3.append(588.075)
x_site_3.append(1512.03)
y_site_3.append(593.219)

x_site_3.append(1406)
y_site_3.append(460)
#x_site_3.append(1396.07)
#y_site_3.append(463.355)

#x_site_3.append(1432.18)
#y_site_3.append(521.713)
x_site_3.append(1421)
y_site_3.append(549)

#x_site = np.concatenate([x_site_1,x_site_2,x_site_3])
#y_site = np.concatenate([y_site_1,y_site_2,y_site_3])
x_site = np.concatenate([x_site_2,x_site_1,x_site_3])
y_site = np.concatenate([y_site_2,y_site_1,y_site_3])







#HARDCODE SITE CENTERS
#x_site = [850.556, 881.042,909.661]
#y_site = [137.06, 191.749, 163.619]
#x_site = [850.556,922.755, 1000.81]
#y_site = [137.06,194.313, 210.824]
#x_site = [891.982 ,1000.81]
#y_site = [173.012, 210.824]



cents = [[y_site[n],x_site[n]] for n in range(len(x_site))]
r_sites = 50#gridpoints
mask_rho = nc_grd.variables['mask_rho'][:,:]
h_smth =ndimage.filters.gaussian_filter(nc_grd.variables['h'][:,:],[5,5])
h_mask = np.zeros(mask_rho.shape)
h_mask[:,:] = 1.
h_mask[h_smth<=2.4] = 0.
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

mask2 = copy.deepcopy(mask_sites)
list_all = range(len(mask_sites))
for n in range(len(mask_sites)):
    list_others = copy.copy(list_all)
    list_others.pop(n)
    sum_others = np.copy(mask2[n]) 
    for k in list_others:
        sum_others+=mask2[k]
    mask2[n][sum_others>1.] = 0.


for n in range(len(mask_sites)-1):
    sum_site = mask_sites[n] + mask_sites[n+1]
    mask_sites_new[n+1][sum_site==2.] = 0

mask_sites_nn = copy.copy(mask_sites_new) 
for n in range(len(mask_sites)-2):
    sum_site = mask_sites_new[n] + mask_sites_new[n+1] + mask_sites_new[n+2]
    mask_sites_nn[n+2][sum_site==2.] = 0.

mask_sites_n3 = copy.copy(mask_sites_nn) 
for n in range(len(mask_sites)-3):
    sum_site = mask_sites_nn[n] + mask_sites_nn[n+1] + mask_sites_nn[n+2] + mask_sites_nn[n+3]
    mask_sites_n3[n+3][sum_site==2.] = 0.

mask_sites_n4 = copy.copy(mask_sites_n3) 
for n in range(len(mask_sites)-4):
    sum_site = mask_sites_n3[n] + mask_sites_n3[n+1] + mask_sites_n3[n+2] + mask_sites_n3[n+3] + mask_sites_n3[n+4]
    mask_sites_n4[n+4][sum_site==2.] = 0.



[plt.contour(mask_sites_new[n],colors='r',linewidths=3) for n in range(len(mask_sites))]
#[plt.contour(mask_sites_nn[n],colors='grey',linewidths=3) for n in range(len(mask_sites))]
#[plt.contour(mask_sites_n3[n],colors='grey',linewidths=3) for n in range(len(mask_sites))]
[plt.contour(mask2[n],colors='grey',linewidths=3) for n in range(len(mask_sites))]





[Ly,Lx] = mask_nan.shape
mask_sum = np.zeros([Ly,Lx])
mask_sum_nn = np.zeros([Ly,Lx])
mask_sum_n3 = np.zeros([Ly,Lx])
mask_sum_n4 = np.zeros([Ly,Lx])
mask2sum   = np.zeros([Ly,Lx])
for n in range(len(mask_sites)):
    mask_sum = mask_sum + mask_sites_new[n]
    mask_sum_nn = mask_sum_nn + mask_sites_nn[n]
    mask_sum_n3 = mask_sum_n3 + mask_sites_n3[n]
    mask_sum_n4 = mask_sum_n4 + mask_sites_n4[n]
    mask2sum    = mask2sum + mask2[n]

plt.figure()
#plt.subplot(3,1,1)
plt.imshow(mask2sum,origin='lower')
plt.colorbar()
#plt.subplot(3,1,2)
#plt.imshow(mask_sum_nn,origin='lower')
#plt.colorbar()
#plt.subplot(3,1,3)
#plt.imshow(mask_sum_nn,origin='lower')
#plt.colorbar()
'''
save_dict = {}
save_dict['mask_sites'] = mask2
save_dict['cents'] = cents
PF.save_to_pickle(save_dict, 'CI_sites')
k''
