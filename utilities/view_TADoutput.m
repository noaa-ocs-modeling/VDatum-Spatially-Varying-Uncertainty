clear variables
%------------- read original record ----------
fname='../../w1_20240601_20240831_Adak9461380_plotdata.nc'; 
t=ncread(fname,'time');
z=ncread(fname,'zeta');
msk=(z<-999);
z(msk)=NaN;
% ----------- read TAD output ---------
xxhhll='xx_hh_ll_Adak9461380'; 
pathin='../TADoutput/';
ifile=[pathin xxhhll '.nc'];
hh=ncread(ifile,'Htime_Hval');
ll=ncread(ifile,'Ltime_Lval');
xx=ncread(ifile,'datums');
ig=ncread(ifile,'stationN');

anan=900;
thg=hh(1:240);
vhg=hh(241:480);
nanmsk=(abs(vhg)>anan);
thg(nanmsk)=[];
vhg(nanmsk)=[];
ihh=(vhg>90);
vhg(ihh)=vhg(ihh)-100;

tlw=ll(1:240);
vlw=ll(241:480);
nanmsk=(abs(vlw)>anan);
tlw(nanmsk)=[];
vlw(nanmsk)=[];
ill=(vlw<-90);
vlw(ill)=vlw(ill)+100;

%----- plot tidal extrema after TAD ------
figure
plot(t(~msk)/3600/24,z(~msk),'linewidth',0.5)
hold on
lwd=1.5;
mrk=9;

plot(thg(ihh)/3600/24,vhg(ihh),['o' 'r'],'linewidth',lwd,'markersize',mrk)
plot(thg(~ihh)/3600/24,vhg(~ihh),['x' 'r'],'linewidth',lwd,'markersize',mrk)

plot(tlw(~ill)/3600/24,vlw(~ill),['+' 'g'],'linewidth',lwd,'markersize',mrk)
plot(tlw(ill)/3600/24,vlw(ill),['o' 'g'],'linewidth',lwd,'markersize',mrk)
legend('tide','HH','LH','HL','LL','location','best')
title(['extrema after TAD, ', int2str(length(ihh)), ...
    ' highs,' int2str(length(ill)),' lows'])
set(gca,'fontsize',10)
