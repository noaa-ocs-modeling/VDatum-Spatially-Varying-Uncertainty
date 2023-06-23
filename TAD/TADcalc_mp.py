#!/home/Elena.Tolkova/miniconda3/bin/python3 -u
import netCDF4 as nc
import numpy as np
from numpy import diff, sign, mean
from scipy.signal import butter, filtfilt
import TADfunctions as func
import sys
import time
import datetime
from datetime import datetime

# call as TADcalc_mp.py TASK_ID TASK_COUNT filename

# task array ID
idtask=int(sys.argv[1])
ntasks=int(sys.argv[2])
# start time
tstart=time.time()
# input file
fn=sys.argv[3]      #input("input file path/name?  ")
print("reading  ",fn," task ",idtask," of ",ntasks)
ds=nc.Dataset(fn,"r")
nsta=len(ds['x'][:]) #number of time histories
tm=ds['time'][:]     #time in sec
print('nsta=',nsta)

nsta_loop=int(np.floor(nsta/ntasks))
nrmndr=int(np.rint(nsta-ntasks*nsta_loop))
nsta0=np.empty(ntasks, dtype=np.int32)
nsta0[0:nrmndr]=nsta_loop+1
nsta0[nrmndr:]=nsta_loop
i0=np.sum(nsta0[0:idtask])
nsers=nsta0[idtask]
print('processing ',nsers,' records starting with ',i0)

xx=np.empty((nsers,7),dtype=np.float64)
sta=np.empty(nsers,dtype=np.int32)
hh=np.empty((nsers,480),dtype=np.float64)
ll=np.empty((nsers,480),dtype=np.float64)
xx.fill(-9999.9)
hh.fill(-9999.9)
ll.fill(-9999.9)
    
for i in range(nsers):

    if i%500==0:
        i1=i+i0
        i2=i1+500
        if i2>nsta: i2=nsta
        yy=ds['zeta'][:,i1:i2]
        yy=np.array(yy,dtype=np.float64)
        
    y=yy[:,i%500]
    
    sta[i]=int(i+i0+1)  #station/node number starting with 1
    if min(y)>-99 and max(y)<99 and max(y)-min(y)>0.05:  #Min_Height_Diff=3 cm
            
        datums, LTimes, HTimes, LVals, HVals, nlows, nhighs = func.ComputeDatums(tm,y)

        xx[i,0:7]=datums
        hh[i,0:nhighs]=HTimes
        hh[i,240:(240+nhighs)]=HVals
        ll[i,0:nlows]=LTimes
        ll[i,240:(240+nlows)]=LVals

    #if i%1000==0: 
    print("record ",i+i0,flush=True)    


print("writing xx_hh_ll")
### ------- output netcdf -----
ncid = nc.Dataset("xx_hh_ll_mp" + str(idtask) + ".nc", "w", format="NETCDF3_CLASSIC")

sta_dim = ncid.createDimension("station", nsers)
dtm_dim = ncid.createDimension("datum", 7)
pik_dim = ncid.createDimension("high_low",480)

sta_var=ncid.createVariable("stationN","i4",("station"))
xx_var=ncid.createVariable("datums","f4",("station","datum"))
hh_var=ncid.createVariable("Htime_Hval","f4",("station","high_low"))
ll_var=ncid.createVariable("Ltime_Lval","f4",("station","high_low"))

sta_var[0:nsers]=np.array(sta[:],dtype=np.int32)  #station/node number starting with 1
xx_var[0:nsers,0:7]=np.array(xx[:,:],dtype=np.float32)  #datums[:]
hh_var[0:nsers,0:480]=np.array(hh[:,:],dtype=np.float32)    #HTimes[:],HVals[:]
ll_var[0:nsers,0:480]=np.array(ll[:,:],dtype=np.float32)   #LVals[:]

ncid.close()    

print('done ',idtask, '; elapsed time (min): ', (time.time()-tstart)/60)

