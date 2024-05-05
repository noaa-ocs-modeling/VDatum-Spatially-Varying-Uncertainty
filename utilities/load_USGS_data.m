%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script uploads simultanious water level records at selected USGS stations 
% given a list of the stations' ID (USGSglist_ids.txt as an example).
% The records' time window is from 2023-01-01T00:00:00.000 to 2023-11-10T23:59:59.999
% in UTC zone 5 (US Eastern Time Zone during standard time), as given by
% parameters sstart and send.
% For each station, its data (water level time histories, 
% with time counted in hours from a moment given by aa0) are saved in a netcdf file 
% in a directory USGSdata (make sure it exists), as given by a variable fdir_out. 
%
%           Elena.Tolkova@noaa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
% --- read station IDs from a text file ---
glist='USGSglist_ids.txt';
fid=fopen(glist);
nsta=0;
while 1
    sdata=textscan(fid,'%s  %*[^\n]',1);
    if isempty(sdata{1})
        break
    else
        nsta=nsta+1;
        gname{nsta}=sdata{1}{1};
    end
end
fclose(fid);

%-- set t=0 at aa0 ----
aa0=datetime('2023-01-01T00:00:00.000-05:00','InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSXXX','TimeZone','UTC');
%-- pieces to form a query ---
s1='https://nwis.waterservices.usgs.gov/nwis/iv/?sites=';
sstart='2023-01-01T00:00:00.000-05:00';
send='2023-11-10T23:59:59.999-05:00';
s2=['&parameterCd=00065&startDT=' sstart '&endDT=' send '&siteStatus=all&format=json'];
j=0; % no-data counter 

for nn=1:nsta
    ssite=gname{nn}
    clear s ss h t q 
    s=webread([s1 ssite s2]);
    if isempty(s.value.timeSeries)
        j=j+1;
        nodata{j}=ssite;
        continue
    end
    sss=[s.value.timeSeries.values];
    ss=sss(1).value;
    nt=length(ss)
    h(1:nt)=NaN;
    t(1:nt)=NaN;
    q(1:nt)=int8(ones(nt,1));
    for i=1:nt
        h(i)=str2double(ss(i).value);
        aa=datetime(ss(i).dateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSXXX','TimeZone','UTC');
        t(i)=datenum(aa-aa0);
        if ss(i).qualifiers{1}~='A'
            q(i)=0;
        end
    end
    t=24*t;     % days to hours
    h=0.3048*h; % feet to meters

    fdir_out='USGSdata/';
    fname=['gage' ssite]
    lon=s.value.timeSeries.sourceInfo.geoLocation.geogLocation.longitude;
    lat=s.value.timeSeries.sourceInfo.geoLocation.geogLocation.latitude;

    ncid=netcdf.create([fdir_out fname '.nc'],'CLOBBER');
%------------ define dims and vars -----------------
    dimid_pnt = netcdf.defDim(ncid,'point',1);
    dimid_tim = netcdf.defDim(ncid,'t',nt);

    lonid=netcdf.defVar(ncid,'longitude','NC_DOUBLE',dimid_pnt); 
    latid=netcdf.defVar(ncid,'latitude','NC_DOUBLE',dimid_pnt); 
    timid=netcdf.defVar(ncid,'time','NC_DOUBLE',dimid_tim); 
    eleid=netcdf.defVar(ncid,'hight','NC_FLOAT',dimid_tim);
    quaid=netcdf.defVar(ncid,'Data-value-qualification','NC_BYTE',dimid_tim);

    fillValue=-9999999;
    netcdf.defVarFill(ncid,timid,true,fillValue);
    netcdf.defVarFill(ncid,eleid,true,fillValue);
    netcdf.putAtt(ncid,timid,'units','hour');
    netcdf.putAtt(ncid,eleid,'units','meter');
    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid,varid,'station ID',ssite);
    netcdf.putAtt(ncid,varid,'record_start_time',[sstart ' UTC']);
    netcdf.endDef(ncid)
%----- write data ------------------------
    netcdf.putVar(ncid,lonid,lon)
    netcdf.putVar(ncid,latid,lat)
    netcdf.putVar(ncid,timid,t)
    netcdf.putVar(ncid,eleid,h)
    netcdf.putVar(ncid,quaid,q)
%-----------------------------------------------
    netcdf.close(ncid);
end

