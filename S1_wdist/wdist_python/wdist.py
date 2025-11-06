# #!/home/Elena.Tolkova/miniconda3/bin/python3 -u
import netCDF4 as nc
import numpy as np
import sys
import time
import datetime
from datetime import datetime
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

# call as ./wdist_test_vect_mp.py TASK_ID gageinfo_filename mesh_filename

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

#----- compute unique segments -------
meshele=np.array(meshele,dtype=np.int32)
A=np.concatenate((meshele[:,0:2],meshele[:,1:3],meshele[:,[2,0]]),axis=0)
Au0, ind = np.unique(np.sort(A, axis=1), axis=0, return_index=True)
nAu = len(Au0)
print(f'Total segments = {nAu}')
#----- remove dry nodes -----
dry = 0.1
mskdry = (z <= dry)
drynodes = np.arange(1, len(z) + 1, dtype=np.int32)
drynodes = drynodes[mskdry]
i1,i2=np.where(np.isin(Au0, drynodes))
i1=np.unique(i1)
Au0 = np.delete(Au0, i1, axis=0)
nAu = len(Au0)
print(f'Total wet segments = {nAu}')
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
node0 = pmoe_datum[iid,1]
# Create a sparse adjacency matrix
rows = Au0[:, 0] - 1
cols = Au0[:, 1] - 1
graph = csr_matrix((dA0, (rows, cols)), shape=(numnod, numnod)) 
# Use Dijkstra's algorithm from scipy
dis, last_node = dijkstra(csgraph=graph, directed=False, indices=node0-1, return_predecessors=True)
#----------------------------------------------------------------
		
print('writing station'+str(pmoe_datum[iid,0])+'; elapsed time (min): ', (time.time()-tstart)/60)

ncid = nc.Dataset('station'+str(pmoe_datum[iid,0])+'.nc', 'w', format='NETCDF4')

mesh_dimid = ncid.createDimension('mesh', numnod)
start_dimid = ncid.createDimension('node0', 1)

dis_var = ncid.createVariable('dis', 'f4', ('mesh',))
lst_var = ncid.createVariable('last_node', 'i4', ('mesh',))
#klr_var = ncid.createVariable('klayer', 'i4', ('mesh',))
node0_var = ncid.createVariable('node0', 'i4', ('node0',))

dis_var[:]=dis
lst_var[:]=last_node+1
#klr_var[:]=klayer
node0_var[:]=node0

ncid.close()    

print('done task ',iid, '; elapsed time (min): ', (time.time()-tstart)/60)
