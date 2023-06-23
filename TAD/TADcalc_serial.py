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

# input file
print(datetime.now())
tstart=time.time()
fn=sys.argv[1] #input("input file path/name?  ")
print("reading  ",fn)
ds=nc.Dataset(fn)
nsta=len(ds['x'][:]) #number of time histories
print(nsta, ' records')
tm=ds['time'][:]    #time in sec

ncid = nc.Dataset("xx_hh_ll.nc", "w", format="NETCDF3_CLASSIC")

sta_dim = ncid.createDimension("station", nsta)
dtm_dim = ncid.createDimension("datum", 7)
pik_dim = ncid.createDimension("high_low",480)

sta_var=ncid.createVariable("stationN","i4",("station"))
xx_var=ncid.createVariable("datums","f4",("station","datum"))
hh_var=ncid.createVariable("Htime_Hval","f4",("station","high_low"))
ll_var=ncid.createVariable("Ltime_Lval","f4",("station","high_low"))

for i in range(nsta):
    y=ds['zeta'][:,i]
    y=np.array(y,dtype=np.float64)
    if min(y)>-99 and max(y)<99 and max(y)-min(y)>0.05:  #Min_Height_Diff=3 cm
        datums, LTimes, HTimes, LVals, HVals, nlows, nhighs = func.ComputeDatums(tm,y)
        sta_var[i]=i+1  #station/node number starting with 1
        xx_var[i,0:7]=datums[:]
        hh_var[i,0:nhighs]=HTimes[:]
        hh_var[i,240:(240+nhighs)]=HVals[:]
        ll_var[i,0:nlows]=LTimes[:]
        ll_var[i,240:(240+nlows)]=LVals[:]
    if np.mod(i,1000)==0: print("record ",i) 
    
ncid.close()    
print('done ; elapsed time (min): ', (time.time()-tstart)/60)
