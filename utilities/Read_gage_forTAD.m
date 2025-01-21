% reads NOS gage data listing and creates emuation of fort.61.nc for a
% single time history, to be fed to TAD.  Input data format (sample):
% 9461380,W1
% Time,Value
% 06/01/2024 00:00,1.369
% 06/01/2024 00:06,1.347
% 06/01/2024 00:12,1.342
% ...
clear variables

data_name='w1_20240601_20240831_Adak9461380_plotdata';
fid=fopen([data_name '.csv'],'r');
fgetl(fid);   % skip line
fgetl(fid);   % skip line 
a=fscanf(fid,'%d/%d/%d %d:%d,%g',[6,inf]); 
a=a';
fclose(fid);
%--------- depth vector ---------
z=a(:,6);
%--------- calculate time vector ----------
ii=[3,1,2,4,5,6];
aa=a(:,ii);
aa(:,6)=0;
tt=datetime(aa);
EPOCH=tt(1);
ttnum=convertTo(tt,'epochtime','Epoch',EPOCH,'TicksPerSecond',1);
tm=double(ttnum);
%---------- find gaps and fill in, if any ----------
clear msk
msk=(abs(z)>9000);
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
iwrite=0;
if iwrite
    fdir='./';
    recname=[data_name '.nc'] 
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


