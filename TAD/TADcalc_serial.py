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
iprint=float(sys.argv[2])  # to print or not to print
tcut=-1.0 # optional 3nd argument - time(hr) to ditch extrema at the record end
if (len(sys.argv)==4):
    tcut=float(sys.argv[3])

print("reading  ",fn, "disregard last ",tcut, " hours")
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
    ymn=min(y)
    ymx=max(y)
    print(i,ymn,ymx,np.isnan(y).any())
    if ymn>-99 and ymx<99 and ymx-ymn>0.05:  #Min_Height_Diff=3 cm
        datums, LTimes, HTimes, LVals, HVals, nlows, nhighs = func.ComputeDatums(tm,y,tcut,iprint)
        sta_var[i]=i+1  #station/node number starting with 1
        xx_var[i,0:7]=datums[:]
        hh_var[i,0:nhighs]=HTimes[:]
        hh_var[i,240:(240+nhighs)]=HVals[:]
        ll_var[i,0:nlows]=LTimes[:]
        ll_var[i,240:(240+nlows)]=LVals[:]
    else:
        xx_var[i,0:7]=-9999.9
        
    if iprint>0:
        print("record ",i)
    #if np.mod(i,1000)==0: print("record ",i) 
    
ncid.close()    
print('done ; elapsed time (min): ', (time.time()-tstart)/60)
