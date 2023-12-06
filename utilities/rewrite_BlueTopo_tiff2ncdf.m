% Rewrite selected tiff tiles in netcdf format, one band and one UTM zone at a time
% elena.tolkova@noaa.gov, 11/29/2023
clear variables
incdf=1;  % set to 1 to save
ires=8;
iutm=19;

fdir=['res' int2str(ires) 'mUTM' int2str(iutm)];  
ff=dir([fdir '/']);
fname={ff.name};
fbytes={ff.bytes};
nfls=length(fname)

for nfl=3:nfls
    
    ffname=fname{nfl}
    idot=strfind(ffname,'.tiff');
    if isempty(idot) | fbytes{nfl}<1111
        continue
    end
    
    clear R xx yy zz
   
    [A,R]=geotiffread([fdir '/' ffname]);
    xlm=R.XWorldLimits;
    ylm=R.YWorldLimits;
    nn=R.RasterSize;
    ny=nn(1);
    nx=nn(2);
    dx=(xlm(2)-xlm(1))/nx;
    dy=(ylm(2)-ylm(1))/ny;
    xx=xlm(1)+dx/2:dx:xlm(2);
    yy=ylm(1)+dy/2:dy:ylm(2);

    switch length(size(A))
        case 2
            AA=A;
        case 3
            AA=A(:,:,1);
        otherwise
            disp(size(A))
            continue
    end
    clear A
    AA=flip(AA);
    AA(abs(AA)>9000)=NaN; 
    zz=AA'; 
    clear AA
    
    if nfl==3
        figure
        dd=10;
        pcolor( xx(1:dd:nx), yy(1:dd:ny), (zz(1:dd:nx,1:dd:ny))' )
        shading flat
        colorbar
        keyboard
    end

    if incdf

        % subsample
        xxx=xx(1:2:end);
        yyy=yy(1:2:end);
        zzz=zz(1:2:end,1:2:end);
        nx=length(xxx);
        ny=length(yyy);

        fdir_out=[fdir '_netcdf/'];
        ncid=netcdf.create([fdir_out ffname(1:(idot-1)) '.nc'],'CLOBBER');
    %------------ define dims and vars -----------------
        dimid_lon = netcdf.defDim(ncid,'x',nx);
        dimid_lat = netcdf.defDim(ncid,'y',ny);

        lonid=netcdf.defVar(ncid,'x','NC_DOUBLE',dimid_lon); 
        latid=netcdf.defVar(ncid,'y','NC_DOUBLE',dimid_lat); 
        depid=netcdf.defVar(ncid,'z','NC_FLOAT',[dimid_lon dimid_lat]);
        fillValue=-9999999;
        netcdf.defVarFill(ncid,depid,true,fillValue);
        netcdf.endDef(ncid)
    %----- write data ------------------------
        netcdf.putVar(ncid,lonid,xxx)
        netcdf.putVar(ncid,latid,yyy)
        netcdf.putVar(ncid,depid,zzz)
    %-----------------------------------------------
        netcdf.close(ncid);
        clear xxx yyy zzz
    end
    disp(nfl)
    
end
