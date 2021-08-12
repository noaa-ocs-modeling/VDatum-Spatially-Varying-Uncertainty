%-----------------------------------------------------------------------
% This program is to compute water distance from tide station nodes to 
% each mesh nodes in an ADCIRC grid. 

% It takes about 80-240 sec for each station for this NEGOM example.
% If the program is stopped, just rerun it again. It will continue the 
% last one station without overwrite the prevous ones.

% Input:  station_node_info.mat
%         all           ADCIRC grid
%         pmoe_datum    Tide stations info

% Output: one file for each tide station 
%      output_dist_time/OWdist8727648_2...

%                               02/14/2017
% 1. Add path output option
% 2. Add minimum travel time option for test purpose.
%                               03/06/2017
%
%        Questions to liujuan.tang@noaa.gov
%
%-----------------------------------------------------------------------
% ------- 0. Input: load  grid and tide station nodal file--------------
clear
imethod=2;                  %1 = plance distance; 2= sphereical distance
reloadok=0;                 %1 = read ADCIRC format grid file ; 
                            %0 = read grid in matlab format
                            
timeok=0;                   % minimize 1 = compute travel time; 
                            %          0 = compute surface distance
%gridfile='input/k6s4j5.grd';     % ADCIRC grid file name if reloadok==1
load station_node_info.mat    % tide station nodal file and grid file

pathout='output_dist_time/';
%-------------------------------------------------------------------
if ~exist(pathout,'dir')
    eval(['mkdir ' pathout])
end
if reloadok==1      
    fid=fopen(gridfile,'rt');
    tit=fgetl(fid);
    [nums]=fscanf(fid,'%i %i', 2);  ne=nums(1); nn=nums(2);
    A=fscanf(fid,'%f',[4,nn])';     x=A(:,2);   y=A(:,3); z=A(:,4);
    A=fscanf(fid,'%f',[5,ne])';     meshele=A(:,[3 4 5]);
    fclose(fid);                    save fort14 meshele x y z ne nn
    fprintf(1,'\n------------------------------------\nInput file: %s\n',gridfile);
else
    ne=all.line2(1); nn=all.line2(2);
    x=all.nodes(:,2);   y=all.nodes(:,3); z=all.nodes(:,4);
    meshele=all.elements(:,[3 4 5]);
end
fprintf(1,'Total elements = %d\n',ne);fprintf(1,'Total nodes =%d\n',nn);
%------- 1.0 Compute unique segements and corresponding distance-----------
A=[meshele(:,1:2);meshele(:,2:3);meshele(:,[3 1])];
[Au0] = unique(sort(A,2), 'rows','stable');
nAu=length(Au0);        fprintf(1,'Total segements = %d \n',nAu);
 
dA0 = pos2dist_1(y(Au0(:,1)),x(Au0(:,1)),y(Au0(:,2)),x(Au0(:,2)),imethod);%km
if timeok==1 
    zA=[z(Au0(:,1))+z(Au0(:,2))]/2;
    loc=find(zA<0.1);
    zA(loc)=0.1;
    v=sqrt(9.81*zA);
    dA0=dA0./v*1000/3600; % [minute]
    fprintf(1,'Compute minimum travel time [hour] \n');
else
    fprintf(1,'Compute minimum distance [km] \n');
end
dis0=[1:nn]'*nan;       
clear A meshele
%--------2.0 start loop for stations---------------------------------------
%nsid=length(sinfo);  % total number of stations
sinfo=pmoe_datum;
nsid=length(pmoe_datum);  % total number of stations

for iid=1:nsid
    if timeok==0
        outfile=['OWdist' int2str(sinfo(iid).id) '_' int2str(imethod) '.mat'];
    else
        outfile=['OWtime' int2str(sinfo(iid).id) '_' int2str(imethod) '.mat'];
    end
        if ~exist([pathout outfile],'file')
        %---2.1 Assign 0 km to station node0----------------------------
        dis=dis0; last_node=dis0;  % last_node connected to current node
        node0=sinfo(iid).node;  dis(node0)=0; last_node(node0)=0;
        t1=node0; kl=0;      
        fprintf(1,'Station %d of %d: %d; node =%d\n',iid,nsid,sinfo(iid).id, node0);
        %---2.2 loop from the node  1---------------------------------- 
        tic
        while 1
            kl=kl+1; node1=t1;t1=[];
            %---find the locations of node1 paires in segements
            [i1,i2]=find(ismember(Au0,node1));
            [i3]=sub2ind([nAu 2],i1,i2);
            iAu=Au0(i1,:);
            iAu0=Au0(i3);
            %---compute accumulated distance, get node2 connected to node1, 
            idis=dis(iAu0)+dA0(i1);
            node2=sum(iAu,2)-iAu0;
            n1=length(node2);
            %---remove duplicates, keep the shortest distance for the same node
            [idis,i1]=sort(idis); 
            node2=node2(i1);    node2_1=iAu0(i1);
            [node2u,i1,i2] = unique(node2,'stable');    
            idis=idis(i1);      node2_1=node2_1(i1);
            %---assign distance to dis
            [dis(node2u),i1]=min([dis(node2u) idis],[],2);
            i2=find(i1==2);  
            t1=node2u(i2);  node2_1=node2_1(i2);
            last_node(t1)=node2_1;
            %if (rem(kl,100)==0 | kl>2200), 
            %fprintf(1,'Layer %d: %d --> %d node for next layer\n',kl,n1,length(t1)); 
            %end
            if isempty(t1) 
                toc
                eval(['save ' pathout outfile ' dis last_node node0 ']); 
                fprintf(1,'Save to %s\n----------------------\n',outfile);  
                break;  
             end
        end
        else
            fprintf(1,'Found %s; Skipped.\n',outfile);
        end
end
%--------End of program-----------------

