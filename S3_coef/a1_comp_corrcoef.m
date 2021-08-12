% here we going to calculate the tidal datum error correlation coefficient 
% (essentially the tidal datum correlation coefficient) between COOPS data
% points and all mesh grid. The STD of the model error will be from the
% error at the station between model and data. some assumption may needed
% to spread the error over whole domain. for example, the near the open
% boundary the error will be consistent with the boundary data set error.
%                       Developed by Lei Shi
%                       Reformated by Liujuan Tang 01/29/202
%Corrected one typo in line 120 by Elena Tolkova 12/01/2020
function a1_comp_corrcoef();
runid='R58_k6s4_msl_5o2_a53_merged';
%---------------------------------------
% load the mesh grid 
fn_grd   = ['../../adcirc_output/fort.14'];% Load model grid
runlist={'mhhw';'mhw';'mlw';'mllw';'mtl';'dtl'};
ilist=2%2:6;%   Test one datum first before run all, e.g.  1:6

pathin='..\S2_lv8jmp\';
flvxx=[pathin runid '_testxx.mat'];
flvhh=[pathin runid '_testhh.mat'];
flvll=[pathin runid '_testll.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. bcpsi: node number corrosponding to COOPS' stations, size mx1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load C:\Users\Liujuan.Tang\Documents\VDatum\NEGulfofMexico\svu\uncertainty\bcpsi.dat;
load ../S1_wdist/station_node_info.mat
bcpsi=bcpsi0;%[sinfo.node]';
m=length(bcpsi);
%---------------------------
[meshele, meshnod, bathy] = xmgrid2quoddy01(fn_grd);
lnodes= length(meshnod); % lnodes=318860 for Chesapeake Bay;
leles = length(meshele);
n=lnodes;

%-------------------------------
[xeq1,xeq2,xeq3,SThhw,EThhw,STlhw,ETlhw,SThlw,EThlw,STllw,ETllw]=process_lv_output(flvxx,flvhh,flvll,bathy,meshnod,bcpsi);

rday=1;
%%%%%%%%%%%%%%%%%%%%%
for whichTS=ilist  
    coef=zeros(n,m);
    idatumname=runlist{whichTS};
    idatum=zeros(n,1);
    switch idatumname
        case 'mhhw'
            Tst01=SThhw;Ten01=EThhw;% the time unit is days
        case 'mhw'
            Tst01=min(SThhw,STlhw);Ten01=max(EThhw,ETlhw);% the time unit is days
        case 'mlw'
            Tst01=min(STllw,SThlw);Ten01=max(ETllw,EThlw);% the time unit is days
        case 'mllw'
            Tst01=STllw;Ten01=ETllw;% the time unit is days
        case 'mtl'
            Tst01=min([SThhw,STlhw,STllw,SThlw]);Ten01=max([EThhw,ETlhw,ETllw,EThlw]);% the time unit is days
        case 'dtl'
            Tst01=min(SThhw, STllw);Ten01=max(EThhw,ETllw);
    end
for k=1:m%112%1:m % for m tide stations % 115 entrance to mobile bay
    tic
    fprintf(1,'%d, %d left\n',k,m-k);
    j=bcpsi(k);
    if xeq1(j,1) == 0 & xeq1(j,2) == 0 
        coef(:,k)=0; coef(j,k)=1;
    else
        loc=[];loc2=[];
        switch idatumname
            case 'mhhw'
                loc=find(xeq3(j,:)==2);
            case 'mhw'
                loc=find(xeq3(j,:)>=1);
            case 'mlw'
                loc=find(xeq3(j,:)<=-1);
            case 'mllw'
                loc=find(xeq3(j,:)==-2);
            case 'mtl'
                loc=find(xeq3(j,:)>=1);
                loc2=find(xeq3(j,:)<=-1);
            case 'dtl'
                loc=find(xeq3(j,:)==2);
                loc2=find(xeq3(j,:)==-2);
        end
        fent1=xeq1(j,loc);fent2=xeq2(j,loc);
        breaks = linspace(min(fent1),max(fent1),round((max(fent1)-min(fent1))/2)); % about 1 break in every 2 days
        pp1 = splinefit(fent1,fent2,breaks);
        %  xx = linspace(ceil(Tst01)+2,floor(Ten01)-2,24*10*(floor(Ten01)-ceil(Tst01)-4));

        xx = linspace(ceil(Tst01)+rday,floor(Ten01)-rday,24*10*(floor(Ten01)-ceil(Tst01)-2*rday));
        yy = ppval(pp1,xx);
        if ~isempty(loc2)
            yy1=yy;
            fent12=xeq1(j,loc2);fent22=xeq2(j,loc2);
            breaks = linspace(min(fent12),max(fent12),round((max(fent12)-min(fent12))/2)); % about 1 break in every 2 days
            pp1 = splinefit(fent12,fent22,breaks);
            yy2 = ppval(pp1,xx);
            yy=(yy+yy2)/2;
        end
        
            for i=1:n%30%32474%1:n %Tif mod(i,1000)==0,i,end;
                if xeq1(i,1)==0.0 &  xeq1(i,2) == 0 
                    coef(i,k)=0.0;
                else
                    loc=[];loc2=[];
                    switch idatumname
                        case 'mhhw'
                            loc=find(xeq3(i,:)==2);
                        case 'mhw'
                            loc=find(xeq3(i,:)>=1);
                        case 'mlw'
                            loc=find(xeq3(i,:)<=-1);
                        case 'mllw'
                            loc=find(xeq3(i,:)==-2);
                        case 'mtl'
                            loc=find(xeq3(i,:)>=1);
                            loc2=find(xeq3(i,:)<=-1);
                        case 'dtl'
                            loc=find(xeq3(i,:)==2);
                            loc2=find(xeq3(i,:)==-2);
                    end
                    cent1=xeq1(i,loc);cent2=xeq2(i,loc);
                    %breaks = linspace(min(cent1),max(cent1),round((max(fent1)-min(fent1))/2)); % about 1 break in every 2 days
                    breaks = linspace(min(cent1),max(cent1),round((max(cent1)-min(cent1))/2)); % Corrected one type by Elena 2020-12-01
                    pp2 = splinefit(cent1,cent2,breaks);
                    dent2 = ppval(pp2,xx);
                    if ~isempty(loc2)
                        dent20=dent2;
                        cent12=xeq1(i,loc2);cent22=xeq2(i,loc2);
                        %breaks = linspace(min(cent12),max(cent12),round((max(fent1)-min(fent1))/2)); % about 1 break in every 2 days
                        breaks = linspace(min(cent12),max(cent12),round((max(cent1)-min(cent1))/2)); % Corrected one typo 2020-12-01
                        pp2 = splinefit(cent12,cent22,breaks);
                        dent22 = ppval(pp2,xx);
                        dent2=(dent2+dent22)/2;
                        %b1_plot
                    end

                    bent05=corrcoef(yy,dent2);
                    coef(i,k)=bent05(1,2); 
                    
                    %fprintf(1,'%6.4f\n',bent05(1,2))
                    %idatum(i)=mean(dent2);
                end
            end
    end
    fprintf(1,'%d coef max=%.4f; min=%.4f\n',k,max(coef(:,k)),min(coef(:,k)));
    toc  
end
bcpsi_coef=bcpsi;
eval(['save coef_sf_' idatumname '_' runid '.mat coef bcpsi_coef']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------

function [meshele, meshnod,bathy] = xmgrid2quoddy01(filename_xmg)

% [meshele, meshnod,bathy] = xmgrid2quoddy(filename_xmg);
% Convert the xmgredit5 mesh file into quoddy arrays.
%
fid = fopen(filename_xmg,'rt');
line = fgetl(fid);
[nums] = fscanf(fid,'%i %i',[2]);
numele=nums(1); 
numnode=nums(2);
A = fscanf(fid,'%f',[4,numnode])';
meshnod= A(:,1:3);
bathy=A(:,[1 4]);
A = fscanf(fid,'%f',[5,numele])';
meshele=A(:,[1 3 4 5]);
fclose(fid);
clear A
%---------------------------------------------------------------------
function [xeq1,xeq2,xeq3,SThhw,EThhw,STlhw,ETlhw,SThlw,EThlw,STllw,ETllw]=process_lv_output(flvxx,flvhh,flvll,bathy,meshnod,bcpsi);
%  load(flvxx,'-ascii');
%  load(flvhh,'-ascii');
%  load(flvll,'-ascii');
 textlist={'hh';'ll';'xx'};
for j=1:3
    itext=textlist{j};
    eval(['flv=flv' itext ';'])  
    eval(['load ' flv ])
    fprintf(1,'Load %s\n',flv);
end
 good_nodes=bathy(:,1);
%--------------------------------------------------------------------
% Remove land or shallow nodes since tidal wave may reach bottom
% this section is to remove land nodes or very shallow nodes, grid
% specific
%Added by Liujuan Tang on Dec. 1, 2017)
%--------------------------------------------------------------------
% rmok=1;
% if rmok==1
%     dmin=1;
%     loc=find( bathy(:,end)<0.2 | (bathy(:,end)<dmin & meshnod(:,2)>-85));
%     fprintf(1,'----\n%d nodes removed due to shallow depth\n-----\n',length(loc))
%     lia=ismember(loc,bcpsi);
%     loc(lia)=[];
%     good_nodes(loc)=[];
%     
% end
%-------------------------
fprintf(1,'process_lv_output...');

ctime=zeros(length(testhh(:,1)),480);cDT=ctime;cDTs1=ctime;cDTs2=ctime;

SThhw=-inf;EThhw=inf;
STlhw=-inf;ETlhw=inf;
SThlw=-inf;EThlw=inf;
STllw=-inf;ETllw=inf;

for ii=1:length(good_nodes)
    i=good_nodes(ii); %1:length(testhh(:,1))
    clear cth03;
    if testxx(i,9) == 0;
        ccc=testhh(i,:);
        ct01=ccc(1:240);ch01=ccc(241:480);
        ct01=ct01(abs(ct01)<995);cs01=1*ones(size(ct01));
        ch01=ch01(abs(ch01)<995);cs01(abs(ch01)>50)=2;ch01(abs(ch01)>50)=ch01(abs(ch01)>50)-100;

        cent01=ct01(cs01==2);
        if min(cent01)<10
        SThhw=max(SThhw,min(cent01));EThhw=min(EThhw,max(cent01));
        cent01=ct01(cs01==1);STlhw=max(STlhw,min(cent01));ETlhw=min(ETlhw,max(cent01));        
        
        ccc=testll(i,:);
        ct02=ccc(1:240);ch02=ccc(241:480);
        ct02=ct02(abs(ct02)<995);cs02=-1*ones(size(ct02));
        ch02=ch02(abs(ch02)<995);cs02(abs(ch02)>50)=-2;ch02(abs(ch02)>50)=ch02(abs(ch02)>50)+100;

        cent01=ct02(cs02==-1);SThlw=max(SThlw,min(cent01));EThlw=min(EThlw,max(cent01));
        cent01=ct02(cs02==-2);STllw=max(STllw,min(cent01));ETllw=min(ETllw,max(cent01));

        cth03(1,:)=[ct01 ct02];
        cth03(2,:)=[ch01 ch02];
        cth03(3,:)=[cs01 cs02];
        cth03(4,:)=min(1,max(-1,cth03(3,:)));
        cth03=cth03';
        [b,c]=sort(cth03(:,1));cth03=cth03(c,:);

%         plot(cth03(cth03(:,3)==2,1),cth03(cth03(:,3)==2,2),'r*');
%         plot(cth03(cth03(:,3)==1,1),cth03(cth03(:,3)==1,2),'b*');
%         plot(cth03(cth03(:,3)==-1,1),cth03(cth03(:,3)==-1,2),'b^');
%         plot(cth03(cth03(:,3)==-2,1),cth03(cth03(:,3)==-2,2),'r^');
        l=length(cth03(:,1));
        ctime(i,1:l)=cth03(:,1)';% this is seq1,seq2,seq3,and seq4 without adjustment of M2 phase
        cDT(i,1:l)=cth03(:,2)';
        cDTs1(i,1:l)=cth03(:,3)';
        cDTs2(i,1:l)=cth03(:,4)';
        else
            testxx(i,9)=-1;
            %fprintf('%d textxx = %d\n',i,testxx(i,9)); 
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
xeq1=ctime;xeq2=cDT;xeq3=cDTs1;xeq4=cDTs2;
%Tsave xeq1234.mat xeq1 xeq2 xeq3 xeq4;
fprintf(1,'...Done.\n--------------------\n');
