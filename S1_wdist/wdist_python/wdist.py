#!/mnt/lfs4/NOS/vdatum/Elena.Tolkova/miniconda3/bin/python3 -u
import netCDF4 as nc
import numpy as np
import sys
import time
import datetime
from datetime import datetime

# call as ./wdist.py TASK_ID gageinfo_filename mesh_filename

print(datetime.now())
tstart=time.time()
# task array ID, also gage # in the list starting with 0
iid=int(sys.argv[1])
#read gage id and node number
fn=sys.argv[2]
pmoe_datum = np.loadtxt(fn, dtype=int)
print('gage info: ', pmoe_datum[iid,:])
# read mesh
fn=sys.argv[3] #input("input file path/name?  ")
with open(fn, 'rt') as fid:
    tit = fid.readline().strip()
    nums = np.loadtxt(fid, dtype=int, max_rows=1)
    numele = nums[0]
    numnod = nums[1]

    BB = np.loadtxt(fid, dtype=float, max_rows=numnod)
    p = BB[:, 1:4]
    B = np.loadtxt(fid, dtype=int, max_rows=numele)
    meshele = B[:, [2, 3, 4]]

x=[row[0] for row in p]
y=[row[1] for row in p]
z=[row[2] for row in p]
x=np.array(x,dtype=np.double)
y=np.array(y,dtype=np.double)
z=np.array(z,dtype=np.double)
    
print('done reading fort.14 ; elapsed time (min): ', (time.time()-tstart)/60)
print('mesh size: ', numnod,numele)
#----- remove dry nodes -----
dry = 0.1
mskdry = (z <= dry)
drynodes = np.arange(1, len(z) + 1, dtype=np.int32)
drynodes = drynodes[mskdry]
i1,i2=np.where(np.isin(meshele, drynodes))
i1=np.unique(i1)
meshele = np.delete(meshele, i1, axis=0)
#----- compute unique segments -------
meshele=np.array(meshele,dtype=np.int32)
A=np.concatenate((meshele[:,0:2],meshele[:,1:3],meshele[:,[2,0]]),axis=0)

Au0, ind = np.unique(np.sort(A, axis=1), axis=0, return_index=True)
nAu = len(Au0)
print(f'Total segments = {nAu}')

#----- compute segments' lengths in km (for x,y being lon,lat) -------
dA0= np.full(nAu, np.nan)
i1=Au0[:,0]-1
i2=Au0[:,1]-1
xx = (x[i1]-x[i2])*np.cos(np.pi*(y[i1]+y[i2])/360)
yy = y[i1]-y[i2]
xx = xx**2
yy = yy**2
dA0=111.111*np.sqrt(xx+yy)

#--------------compute for station iid-------------------------
nn=len(x)  # number of nodes
dis = np.full(nn, 999999.9)
last_node = np.full(nn, -1, dtype=np.int32) # last_node connected to current node
klayer = np.full(nn, -1, dtype=np.int32)
node0 = pmoe_datum[iid,1]
dis[node0-1] = 0
last_node[node0-1] = 0
klayer[node0-1] = 0
t1 = node0
kl = 0
    
while True:
        
    kl += 1
    node1 = t1
        
    # --- find the locations of node1 pairs in segments
    i1,i2=np.where(np.isin(Au0, node1))
    iAu = Au0[i1]  # pairs: node1 - connected node
    iAu0=Au0[i1,i2].copy()
    
    node2 = np.sum(iAu, axis=1) - iAu0  # nodes connected to node0 via node1
    n1 = len(node2) 
    #---- accumulated distance ---
    idis=dA0[i1]+dis[iAu0-1]
    j1 = np.argsort(idis)
    idis=idis[j1]

    node2 = node2[j1]
    node2_1 = iAu0[j1]

    node2u, i1 = np.unique(node2, return_index=True)
    i1=np.sort(i1)
    node2u=node2[i1]
    idis=idis[i1]
    node2_1 = node2_1[i1]
        
    mask = dis[node2u - 1] > idis
    i2 = np.where(mask)[0]
    if len(i2)==0:
        break
            
    dis[node2u - 1] = np.where(mask, idis, dis[node2u - 1])
    t1=node2u[i2]
    last_node[t1-1]=node2_1[i2]
    klayer[t1-1]=kl
    print(kl)
		
print('writing station'+str(pmoe_datum[iid,0])+'; elapsed time (min): ', (time.time()-tstart)/60)

ncid = nc.Dataset('station'+str(pmoe_datum[iid,0])+'.nc', 'w', format='NETCDF4')

mesh_dimid = ncid.createDimension('mesh', numnod)
start_dimid = ncid.createDimension('node0', 1)

dis_var = ncid.createVariable('dis', 'f4', ('mesh',))
lst_var = ncid.createVariable('last_node', 'i4', ('mesh',))
klr_var = ncid.createVariable('klayer', 'i4', ('mesh',))
node0_var = ncid.createVariable('node0', 'i4', ('node0',))

dis_var[:]=dis
lst_var[:]=last_node
klr_var[:]=klayer
node0_var[:]=node0

ncid.close()    

print('done task ',iid, '; elapsed time (min): ', (time.time()-tstart)/60)