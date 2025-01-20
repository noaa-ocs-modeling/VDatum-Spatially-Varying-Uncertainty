import numpy as np
from numpy import diff, sign, mean
from scipy.signal import butter, filtfilt
from scipy.interpolate import Akima1DInterpolator
#####################  definitions  ###############################################
#
# Filter definitions for Datums Calculator Tide picker
#
def butter_lowpass(cutOff, fs, order=5):
    nyq = 0.5 * fs
    normalCutoff = cutOff / nyq
    b, a = butter(order, normalCutoff, btype='low', analog = False)
    return b, a

def butter_lowpass_filter(data, cutOff, fs, order=4):
    b, a = butter_lowpass(cutOff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

####################################################################

def Check_Tide_Order(dt, h, l, iprint):
#   Check that tides are in High-Low order
#   Merge the Highs and Lows into a single time-ordered list
    if iprint > 0:
        print (len(h), 'highs  ', len(l), 'lows   ')
    hi = 0
    li = 0
    tides = []
    tide_types = []
    tide_indexes = []
    while ((hi<len(h)) and (li<len(l))):
        if (dt[h[hi]] <= dt[l[li]]):
            tides.append(h[hi])
            tide_types.append('H')
            tide_indexes.append(hi)
            hi = hi + 1
        else:
            tides.append(l[li])
            tide_types.append('L')
            tide_indexes.append(li)
            li = li + 1
    while hi < len(h):
        tides.append(h[hi])
        tide_types.append('H')
        tide_indexes.append(hi)
        hi = hi + 1
    while li < len(l):
        tides.append(l[li])
        tide_types.append('L')
        tide_indexes.append(li)
        li = li + 1

    ttype = tide_types[0]
    for i in range(1,len(tides)-1):
        if tide_types[i] == ttype:
            print ('Tides are out of order at:', dt[tides[i]])
            return -1
        ttype = tide_types[i]
        i = i+1

    return 1
   
####################################################################   

def Check_Tides(dt, wl, h, l, Units_Factor, iprint):
# Check Tides for Minimum time and height between neighbors
#   Merge the Highs and Lows into a single time-ordered list
#    print (len(h), 'highs  ', len(l), 'lows   ')
    Min_Height_Diff = 0.03 * Units_Factor
    hi = 0
    li = 0
    tides = []
    tide_types = []
    tide_indexes = []
    while ((hi<len(h)) and (li<len(l))):
        if (dt[h[hi]] <= dt[l[li]]):
            tides.append(h[hi])
            tide_types.append('H')
            tide_indexes.append(hi)
            hi = hi + 1
        else:
            tides.append(l[li])
            tide_types.append('L')
            tide_indexes.append(li)
            li = li + 1
    while hi < len(h):
        tides.append(h[hi])
        tide_types.append('H')
        tide_indexes.append(hi)
        hi = hi + 1
    while li < len(l):
        tides.append(l[li])
        tide_types.append('L')
        tide_indexes.append(li)
        li = li + 1

    hi_mask = np.ones(len(h), dtype=bool)
    lo_mask = np.ones(len(l), dtype=bool)

#/* Walk through the tides and mark offending pairs for deletion */

    t1 = 0
    t2 = 1
    aredeletedtides = 0
    while (t2<len(tides)):
#        print t1, t2, dt[tides[t1]], dt[tides[t2]], dt[tides[t2]]-dt[tides[t1]], abs(wl[tides[t2]] - wl[tides[t1]])
#        print 't1 = {0:s} ({1:s})    t2 = {2:s} ({3:s})\n'.format(str(dt[tides[t1]]), tide_types[t1], str(dt[tides[t2]]), tide_types[t2])
        if ((dt[tides[t2]] - dt[tides[t1]]) < 7200 and (tide_types[t1] == tide_types[t2])):
#           Delete second tide 
            if tide_types[t2] == "H":
                if iprint>0:
                    print ('Deleting tide at ', dt[h[tide_indexes[t2]]], ' for min time')
                hi_mask[tide_indexes[t2]] = False
            else:
                if iprint>0:
                    print ('Deleting tide at ', dt[l[tide_indexes[t2]]], ' for min time')
                lo_mask[tide_indexes[t2]] = False
            t2=t2+1
            aredeletedtides = 1;
        elif (((dt[tides[t2]] - dt[tides[t1]]) < 7200) or (abs(wl[tides[t2]]-wl[tides[t1]]) < Min_Height_Diff)) and \
                   (tide_types[t1] != tide_types[t2]):
#           Delete both tides 
            if tide_types[t1] == "H":
                if iprint>0:
                    print (' Deleting 2 tides at {0:s} for min time/range.'.format(str(dt[h[tide_indexes[t1]]])))
                hi_mask[tide_indexes[t1]] = False
            else:
                if iprint>0:
                    print (' Deleting 2 tides at {0:s} for min time/range.'.format(str(dt[l[tide_indexes[t1]]])))
                lo_mask[tide_indexes[t1]] = False
            if tide_types[t2] == "H":
                hi_mask[tide_indexes[t2]] = False
            else:
                lo_mask[tide_indexes[t2]] = False
            t1=t2+1
            t2=t1+1
        else:
            t1=t2
            t2=t1+1
    return hi_mask, lo_mask

####################################################################    

def Highest(h_dts, h_vals, t1, t2):
#   Return the index of the highest value between t1 and t2
    mxindex = -1
    mxval = -99999.99
    for i in range(len(h_dts)):
        if ((h_dts[i] >= t1) and (h_dts[i] <= t2)):
            if (h_vals[i] > mxval):
               mxval = h_vals[i]
               mxindex = i
    return mxindex

####################################################################    

def Lowest(l_dts, l_vals, t1, t2):
#   Return the index of the lowest value between t1 and t2
    mxindex = -1
    minval = 99999.99
    for i in range(len(l_dts)):
        if ((l_dts[i] >= t1) and (l_dts[i] <= t2)):
            if (l_vals[i] < minval):
               mxval = l_vals[i]
               mxindex = i
    return mxindex

####################################################################    

def NextClosest(h_dts, t, nc1):
#   Return the index of the tide nearest to t, following index nc1
    tmin=50.*3600
    for i in range(len(h_dts)):
        if abs(h_dts[i]-t)<tmin:
            nc2=i
            tmin=abs(h_dts[i]-t)
    if nc2==nc1:
        nc2=nc1+1
    return nc2

####################################################################    

def Local_Max(dt, wl, h, t):
# Check wl for max value in +-t around time of h
    win_center = dt[h]
    win_start =  max(win_center - t,dt[1])
    win_end = min(win_center + t,dt[-2])
    max_val = wl[h]
    max_loc = h
    loc = h
    while dt[loc] >= win_start:
        if wl[loc] > max_val:
            max_val = wl[loc]
            max_loc = loc
        loc = loc - 1
    loc = h
    while dt[loc] <= win_end:
        if wl[loc] > max_val:
            max_val = wl[loc]
            max_loc = loc
        loc = loc + 1
    return dt[max_loc], wl[max_loc]

####################################################################   

def Local_Min(dt, wl, h, t):
# Check wl for min value in +-t around time of h
    win_center = dt[h]
    win_start = max(win_center - t,dt[1])
    win_end = min(win_center + t,dt[-2])
    min_val = wl[h]
    min_loc = h
    loc = h
    while dt[loc] >= win_start:
        if wl[loc] < min_val:
            min_val = wl[loc]
            min_loc = loc
        loc = loc - 1
    loc = h
    while dt[loc] <= win_end:
        if wl[loc] < min_val:
            min_val = wl[loc]
            min_loc = loc
        loc = loc + 1
    return dt[min_loc], wl[min_loc]

####################################################################   

def ComputeDatums(x,y,tcut,iprint):
	# Set up Filter parameters.
	fs= 86400 / (x[2]-x[1]) # Interval.seconds 
	order = 6
	cutoff = 5.0  # desired cutoff frequency of the filter, per day
	#cutoff = 8.0
	# SDC_Print(['Sampling Rate: ', fs, ' per day.   Using cutoff frequency of ', cutoff , ' per day'])
	# Filter the data
	filt = butter_lowpass_filter(y, cutoff, fs, order)
	
	# find inflection points (tides) in fitered signal
	highs = (diff(sign(diff(filt))) < 0).nonzero()[0] + 1 # local max
	lows  = (diff(sign(diff(filt))) > 0).nonzero()[0] + 1 # local min
	
	# check potential tides for spacing in time and height
	CFactor=1  # metric units
	highs_mask, lows_mask = Check_Tides(x, y, highs, lows, CFactor, iprint) 
	
	#  Delete bad tides
	highs = highs[highs_mask]
	lows = lows[lows_mask]
	
	MSL=0
	MLLW = 0.0
	MLW = 0.0
	nlows = 0
	nllows = 0
	LVals = []
	LTimes = []		
	MHHW = 0.0
	MHW = 0.0
	nhighs = 0
	nhhighs = 0
	HVals = []
	HTimes = []		
	
	if (len(highs)>3) and (len(lows)>3):
	
		#Check Tide Order
		CTO = Check_Tide_Order(x, highs, lows, iprint)
		#if (CTO < 0):
		#    SDC_Print(["***Warning*** - Tides are out of order"])
		#    OutFile.close
		##    exit(-1)
		
		high_values = []
		high_dts = []
		low_values = []
		low_dts = []
		#   Just pick the highest/lowest point within specified window (+- 30 minutes)
		for i in range(len(highs)):
		    high_dt, high_val = Local_Max(x, y, highs[i], 1800)
		    high_values.append(high_val)
		    high_dts.append(high_dt)
		for i in range(len(lows)):
		    low_dt, low_val = Local_Min(x, y, lows[i], 1800)
		    low_values.append(low_val)
		    low_dts.append(low_dt)

		high_types=analyze_types(high_dts,high_values)
		lows=[]
		for i in range(len(low_values)):
		    lows.append((-1)*low_values[i])
		low_types=analyze_types(low_dts,lows)
		
		# added on Jan 19 2025 - remove extrema within tcut(hr) from the end of the record
		if tcut>0:
		    tstop = x[-1] - 3600 * tcut  # Calculate the cutoff threshold

		    # Find indices where conditions are met
		    ihcut = [i for i, val in enumerate(high_dts) if val > tstop]
		    ilcut = [i for i, val in enumerate(low_dts) if val > tstop]

		    # Print the number of items to be deleted
		    if iprint>0:
		        print("Deleting ", len(ihcut), " highs and ", len(ilcut), " lows in the last ", tcut, " hours")
		else:
		    ihcut = []
		    ilcut = []
    
		#    Calculate Datums by First Reduction
		MSL = np.mean(y, dtype=np.float64)
		#     High Means 
		for i in range(len(highs)-len(ihcut)): 
		    if (high_types[i] == 1):
		        MHHW = MHHW + high_values[i]
		        nhhighs = nhhighs + 1
		        MHW = MHW + high_values[i]
		        nhighs = nhighs + 1
		        HTimes.append(high_dts[i])
		        HVals.append(high_values[i]+100)            
		    if (high_types[i] == 0):
		        MHW = MHW + high_values[i]
		        nhighs = nhighs + 1
		        HTimes.append(high_dts[i])
		        HVals.append(high_values[i])
		        
		MHHW = MHHW / nhhighs
		MHW  = MHW  / nhighs
		
		#     Low Means 
		for i in range(len(lows)-len(ilcut)): 
		    if (low_types[i] == 1):
		        MLLW = MLLW + low_values[i]
		        nllows = nllows + 1
		        MLW = MLW + low_values[i]
		        nlows = nlows + 1
		        LTimes.append(low_dts[i])
		        LVals.append(low_values[i]-100)
		    if (low_types[i] == 0):
		        MLW = MLW + low_values[i]
		        nlows = nlows + 1
		        LTimes.append(low_dts[i])
		        LVals.append(low_values[i])
		        
		MLLW = MLLW / nllows
		MLW  = MLW  / nlows
		
	datums=[]
	datums.append(MHHW)
	datums.append(MHW)
	datums.append(MSL)
	datums.append((MHW+MLW)/2)
	datums.append((MHHW+MLLW)/2)
	datums.append(MLW)
	datums.append(MLLW)
	datums=np.array(datums)
	LTimes=np.array(LTimes)
	HTimes=np.array(HTimes)    
	LVals=np.array(LVals)    
	HVals=np.array(HVals)
	
	return datums, LTimes, HTimes, LVals, HVals, nlows, nhighs

####################################################################   

def analyze_model_tides(t_wl,wl):

	t_wl_unix = t_wl.astype('datetime64[s]').astype('i8') # convert time into an integer seconds of unix_epoch for interpolant
	
	msl = wl[t_wl_unix % 3600].mean() # CO-OPS defines MSL as average of hourly water levels
	    
	akima = Akima1DInterpolator(t_wl_unix,wl) 
	t_extrema_akima = akima.derivative(1).roots()
	akima_t_extrema = akima(t_extrema_akima)
	lows = np.array(akima.derivative(2)(t_extrema_akima)).clip(0).astype(np.bool8)
	t_mlw,mlw = t_extrema_akima[lows],akima_t_extrema[lows]
	t_mhw,mhw = t_extrema_akima[~lows],akima_t_extrema[~lows]
	
	argHigherLows,argLowerHighs = [],[] # "25 hour algorithm" to determine higher lows & lower highs
	higherlows = np.zeros(mlw.shape).astype(np.bool8)
	lowerhighs = np.zeros(mhw.shape).astype(np.bool8)
	tidalTus = 12.5 * 3600 # 25/2 hour period; nominal is (12 + 25.2/60) * 3600
	for argSubExtrema,t_mw,mw,extremaFn in ((argHigherLows,t_mlw,mlw,np.argmin),(argLowerHighs,t_mhw,mhw,np.argmax)):
	    extremaTypeDifs = np.diff(np.hstack(([0],(np.diff(t_mw) < tidalTus).astype('int8'),[0])))
	    iStarts, = np.where(extremaTypeDifs > 0)
	    iEnds, = np.where(extremaTypeDifs < 0)
	    iEnds += 1 # per slicing
	    for iStart,iEnd in zip(iStarts,iEnds):
	        mwSemiEvalArg,relSemiChecked = np.arange(iStart,iEnd),np.zeros(iEnd-iStart).astype(np.bool8)
	        while np.any(~relSemiChecked):
	            mwSemiSubEvalArg = mwSemiEvalArg[~relSemiChecked]
	            relSemiSubExtremaArg = extremaFn(mw[mwSemiSubEvalArg])
	            t_mwSemiSubExtrema = t_mw[mwSemiSubEvalArg][relSemiSubExtremaArg]
	            relSemiSubChecked = np.fabs(t_mw[mwSemiSubEvalArg] - t_mwSemiSubExtrema) < tidalTus
	            relSemiChecked[~relSemiChecked] |= relSemiSubChecked
	            argSubExtrema += list(mwSemiSubEvalArg[relSemiSubChecked])
	            argSubExtrema.remove(mwSemiSubEvalArg[relSemiSubExtremaArg])
	
	higherlows[argHigherLows] = True
	lowerhighs[argLowerHighs] = True
	
	# datetime-based tidal extrema time series 
	t_higherhighs,z_higherhighs = t_mhw[~lowerhighs],mhw[~lowerhighs]
	t_lowerlows,z_lowerlows = t_mlw[~higherlows],mlw[~higherlows]
	
	MHHW=z_higherhighs.mean()
	MHW=mhw.mean()
	MLW=mlw.mean()
	MLLW=z_lowerlows.mean()
	mhw[~lowerhighs]=mhw[~lowerhighs]+100
	mlw[~higherlows]=mlw[~higherlows]-100
	
	datums=[]
	datums.append(MHHW)
	datums.append(MHW)
	datums.append(msl)
	datums.append((MHW+MLW)/2)
	datums.append((MHHW+MLLW)/2)
	datums.append(MLW)
	datums.append(MLLW)
	datums=np.array(datums)
	LTimes=np.array(t_mlw)
	HTimes=np.array(t_mhw)    
	LVals=np.array(mlw)    
	HVals=np.array(mhw)
	nlows=len(mlw)
	nhighs=len(mhw)
	
	return datums, LTimes, HTimes, LVals, HVals, nlows, nhighs
####################################################################   

def analyze_types(h_dts,h_vals):
	
	htypes = []
	nh=len(h_dts)
	for i in range(nh):
	    htypes.append(0)
	#   Highest tide in first 13 hours, else first tide
	nc1 = 0
	if (h_dts[1]-h_dts[0]<=13.0*3600) and (h_vals[1]>h_vals[0]):
	    nc1=1
	htypes[nc1]=1    
	
	while nc1<nh:    
	    t1 = h_dts[nc1]
	    # get t2
	    t2=t1+25.0*3600  
	    nc2=NextClosest(h_dts,t2,nc1)
	    if nc2>nh-1:
	        break
	    t2=h_dts[nc2]
	# get t3
	    t3=t2+25.0*3600
	#find extrema closest to t3
	    nc3=NextClosest(h_dts,t3,nc2)
	    if nc3>nh-1:
	        break
	    t3=h_dts[nc3]
	    ta=(t1+t2)/2
	    tb=(t2+t3)/2
	    nc2 = Highest(h_dts, h_vals, ta, tb)
	    t2=h_dts[nc2]
	    tc=(t1+t2)/2
	    nc1 = Highest(h_dts, h_vals, t1, tc)
	    htypes[nc1]=1
	    nc1=nc2
	    
	if (htypes[nh-1]==0) and (htypes[nh-2]==0):
	    if h_vals[nh-2]>h_vals[nh-1]:
	        htypes[nh-2]=1
	    else:
	        htypes[nh-1]=1

	return htypes

