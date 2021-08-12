%--------------------------------------------------------
% This program is to find nearest node to the stations (for
% computing water distance for the next step)

% Input files:  
% 1. grid file
% 2. info file for tide stations 
%       obs = 
%  struct with fields:
%
%              id: [193×1 int32]
%               y: [193×1 double]
%               x: [193×1 double]
%            mhhw: [193×1 double]
%             mhw: [193×1 double]
%             mlw: [193×1 double]
%            mllw: [193×1 double]
%          navd88: [193×1 double]
%             rms: [193×1 double]
%      station_name: {193×1 cell}

% Output file:   station_node_info.mat 
%       which contains the following variables
%
%       sinfo        original station info, 193 obs
%       pmoe_datum   stations covered by mesh and will be used for VDatum application,
%                    182 obs for the NEGOM test case
%       all          grid info in matlab format
%       bcpsi0       Station node number list

%           Version 2.0
%           Quesetions to liujuan.tang@noaa.gov
%           01/24/2020
%------------------------------------------------------------
% 
clear
% observation's dirctory and file name for tide stations
path_in='input/';
station_file='obs_negov_20170419';


path_grid='../../adcirc_output/';
gridfile='fort.14';
plotok=1;       % plot grid 1=yes; 2=no
check_inok=0;   % perform additional check whether the tide gage is within an element
excludeok=1;    % exclude stations not covered by the mesh.
%                 % 1=yes; 0= no
exdist=1000;    % excluding gagues  with distance greater than 1000m;  
                % need to check the plot to make sure the right amount of
                % gagues included; change exdist if needed.
%---------------------------------------------------------
eval(['load ' [path_in station_file]]);% load obs for tide stations 
sdata=double([obs.x obs.y obs.id]);
    
pts=sdata(:,1:3);
nsid=length(pts); % total number of tide stations 
%---------------------------------------------------------
%load ADCIRC grid into matlab format
all=f_read_bathy_adcirc_grd(path_grid,gridfile,0);

p=all.nodes;  
p(:,[1 end])=[];
q=all.elements;  
q(:,1:2)=[];
ne=length(q);

      
%---------------
%2. Find the nearest node for each tide station
    
inok=obs.x*0;       %station covered by grid: 0=no, >0 yes
snode=inok*0;       %nearest node for each tide station
dis2node=inok;
res_ele=inok;
outok=inok;
for i=1:nsid
    col='';
    iid=obs.id(i);
    fprintf(1,'%d %d\n',i,iid);
    ixy=[obs.x(i) obs.y(i)];%pts(i,1:2);
    
    %--------------------------
    % 2.1 find the closest node point
    dist = pos2dist_1(ixy(2),ixy(1),p(:,2),p(:,1),2)*1000; %  (m)
    [dist_min,eloc]=sort(dist);
    snode(i)=eloc(1);
    res_ele(i)=pos2dist_1(p(eloc(1),2),p(eloc(1),1),p(eloc(2),2),p(eloc(2),1),2)*1000;
    dis2node(i)=dist_min(1);eloc=eloc(1);
    
    fprintf(1,'Nearest node: %d; distance = %.0f m\n',eloc,dis2node(i));
    fprintf(1,'The nearest element length= %.0f m\n',res_ele(i));
  
    %--------------------------
    
    %2.2 check if point is with in elements
    if check_inok==1;
    for j=1:ne
        jloc=q(j,[1:3 1]);
        xv=p(jloc,1);
        yv=p(jloc,2);
        in = inpolygon(ixy(1),ixy(2),xv,yv);
        if in
            fprintf(1,'station %d in element %d\n',i,j);
            inok(i)=j;%col='g';  dx=dx0;
            break
        end
    end
    fprintf(1,'---------------\n');
    end
end
obs.inok=inok;
obs.snode=snode;
obs.dis2node=dis2node;
obs.res_ele=res_ele;
%------------------
nn=1:nsid;
if excludeok==1
loc_ex=find(dis2node>1000);
ns=length(loc_ex);
fprintf('%d stations with too far distance; removed\n',ns)
nn(loc_ex)=[];
end
%---------------------------------------------
% reformat station information with node and depth information
sinfo=[];

for i=1:nsid
    ii=i;%nn(i);
    sinfo(i).id=obs.id(ii);
    sinfo(i).ox=obs.x(ii) ;
    sinfo(i).oy=obs.y(ii) ;
%    sinfo(i).oxy=[obs.x(ii) obs.y(ii)];
    sinfo(i).node=snode(ii);
    sinfo(i).depth=all.nodes(snode(ii),end); %water depth at the node
    sinfo(i).element_length=dis2node(ii);
    sinfo(i).inok=inok(ii);
    sinfo(i).omhhw=obs.mhhw(ii);
    sinfo(i).omhw=obs.mhw(ii);
    sinfo(i).omllw=obs.mllw(ii);
    sinfo(i).omlw=obs.mlw(ii);
    sinfo(i).orms=obs.rms(ii);
    sinfo(i).onavd88=obs.navd88(ii);
    sinfo(i).station_name=obs.station_name(ii);
end
%----------------------
%---------------------
% check stations with rms 
orms=[sinfo.orms]';
loc=find(isnan(orms) | orms<-9);
fprintf('    %d stations has  rms. \n',length(orms)-length(loc));
pmoe_datum=sinfo;
% remove tide stations that will not be used
if excludeok==1
     pmoe_datum(loc_ex)=[];
end
bcpsi0=[pmoe_datum.node]';
fprintf('%d Original stations \n',nsid)
fprintf('%d Final used stations \n',length(pmoe_datum))

%-----------------------------------
% save to matlab data file
save station_node_info sinfo pmoe_datum all  bcpsi0
    
%-----------------------------------------
% Plot grid to check  origianl gauges locations and node locations 
if plotok==1
    clf
    set(gcf,'Paperpositionmode','auto','visible','on');
    hold on
    col=[0 0 0]+0.5;
    fv.vertices=p;
    %-------------
    fv.faces=q;
    patch(fv,'FaceColor',col,'EdgeColor',col,'FaceAlpha',0)

    hold on
    plot([pmoe_datum.ox],[pmoe_datum.oy],'g+')   
    hold on
    plot(all.nodes(bcpsi0,2),all.nodes(bcpsi0,3),'ro');
    
    hold on
    plot([sinfo(loc_ex).ox],[sinfo(loc_ex).oy],'b+') % excluded gauge in blue

    legend('Mesh','COOPS stations ','node location', 'excluded')
end
    
