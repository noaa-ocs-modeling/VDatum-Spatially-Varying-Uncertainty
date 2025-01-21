% reads DART data listing and creates emuation of fort.61.nc for a
% single time history, to be fed to TAD.  Input data format (sample):
% #YY  MM DD hh mm ss T   HEIGHT
% #yr  mo dy hr mn  s -      m
% 2023 06 30 23 45 00 1 4769.206
% 2023 06 30 23 30 00 1 4769.170
% 2023 06 30 23 15 00 1 4769.129
% ...
clear variables

DARTname='21415Attu_Apr_Jun2023';
fid=fopen([DARTname '.txt'],'r');
fgetl(fid);   % skip line
fgetl(fid);   % skip line 
a=fscanf(fid,'%d %d %d %d %d %d %d %g',[8,inf]); 
a=a';
fclose(fid);
%--------- remove 1 min and 15 sec data if any --------
msk=(a(:,7)>1.1);
a(msk,:)=[];
%--------- depth vector, accending time order ---------
z=flip(a(:,8));
%--------- calculate time vector ----------
tt=datetime(a(:,1:6));
EPOCH=tt(end);
ttnum=convertTo(tt,'epochtime','Epoch',EPOCH,'TicksPerSecond',1);
tm=double(flip(ttnum));
%---------- find gaps and fill in, if any ----------
clear msk
msk=(z>9000);
nz=length(z);
dmsk=msk(2:nz)-msk(1:(nz-1));
gap_start=find(dmsk==1);
gap_end=find(dmsk==-1);
gap_len=(gap_end-gap_start)*0.25;
ngp=length(gap_len);
if ngp>0
    disp(['found ' int2str(ngp) ' gaps, ' int2str(max(gap_len)), ' hr max'])
    z(msk)=[];
    tz=tm;
    tz(msk)=[];
    zz=interp1(tz,z,tm);
else
    disp('no gaps in the record')
    zz=z;
end
zmean=sum(zz)/length(zz);
zz=zz-zmean;
%-------- admire the result --------------
figure
plot(tm/3600/24,zz)
xlabel(['days, starting ' int2str(convertTo(EPOCH,'yyyymmdd'))])
ylabel('m')
%----- write emulation of fort.61.nc --------------
iwrite=1;
if iwrite
    fdir='./';
    recname=[DARTname '.nc'] 
    gid=1;
    ng=1;
    lon=0;
    lat=0;
    nt=length(tm);
    ncid = netcdf.create([fdir recname],'CLOBBER');
    ggDimId  = netcdf.defDim(ncid,'station',ng);
    tmDimId = netcdf.defDim(ncid,'time',nt);
    xgid=netcdf.defVar(ncid,'x','double',ggDimId)        
    ygid=netcdf.defVar(ncid,'y','double',ggDimId)        
    tmid=netcdf.defVar(ncid,'time','double',tmDimId)
    gageID=netcdf.defVar(ncid,'gageID','int',ggDimId)  
    zzid=netcdf.defVar(ncid,'zeta','double',[ggDimId, tmDimId]) 
    netcdf.endDef(ncid)
    netcdf.putVar(ncid,xgid,lon) 
    netcdf.putVar(ncid,ygid,lat) 
    netcdf.putVar(ncid,zzid,zz)
    netcdf.putVar(ncid,tmid,tm) 
    netcdf.putVar(ncid,gageID,gid) 

    netcdf.close(ncid);
end


