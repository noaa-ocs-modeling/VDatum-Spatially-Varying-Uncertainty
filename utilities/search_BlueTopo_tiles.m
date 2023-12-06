% script to search for surveys, one band and one UTM zone at a time
% puts together a list of names of the selected surveys
% elena.tolkova@noaa.gov, 11/13/2023
clear variables

P=genpath('C:\Users\elena.tolkova\Documents\OceanMesh2D');
addpath(P)

VD=m_shaperead('NBStiles'); % upload the latest
% https://noaa-ocs-nationalbathymetry-pds.s3.amazonaws.com/index.html#Test-and-Evaluation/Modeling/_Modeling_Tile_Scheme/
nb=length(VD.ncst)
msk_res=int32(zeros(nb,1));
msk_utm=int32(zeros(nb,1));

for n=1:nb
    if ~isempty(VD.dbfdata{n,6})
        s=VD.dbfdata{n,6}; %resolution
        ns=length(s);
        msk_res(n)=int32(str2double(s(1:(ns-1))));
        ss=VD.dbfdata{n,7}; % UTM zone
        msk_utm(n)=int32(str2double(ss));
        clear s ss
    end
end

A=m_shaperead('NBStiles_select_contour1');   % my work area, can use grid boundary 
xbox=A.ncst{1}(:,1);
ybox=A.ncst{1}(:,2);
clear A
figure
plot(xbox,ybox,'g')
hold on

%-- comment above to rerun with another ires ------
%---- select tiles in a given UTM zone at a given resolution ---
ires=4;  % 2,4,8,16 m
iutm=19;  % UTM zone
msk=(msk_res==ires)&(msk_utm==iutm);
ii=find(msk==1);
%---- exclude tiles outside the mesh ------
for n=1:length(ii)
    ax=VD.ncst{ii(n)}(:,1);
    ay=VD.ncst{ii(n)}(:,2);
    in=inpolygon(ax,ay,xbox,ybox);
    if sum(in)==0
        ii(n)=0;
    end
end
ii(ii==0)=[];

k=length(ii);
disp(['plotting ' int2str(k) ' data boxes in the work area'])
for n=1:k  
    ax=VD.ncst{ii(n)}(:,1);
    ay=VD.ncst{ii(n)}(:,2);
    plot(ax,ay,'r')
end
%---- write command file for uploading selected tiles ------
% fname=['tiles_res' int2str(ires) 'm_UTM' int2str(iutm)]
% fid=fopen([fname '.cmd'],'wt');
% for n=1:length(ii)  
%     s=['curl -O ' VD.dbfdata{ii(n),3}];  % load VD.dbfdata{ii(n),3}
%     fprintf(fid,'%s \n', s);
%     clear s
% end
% fclose(fid);

