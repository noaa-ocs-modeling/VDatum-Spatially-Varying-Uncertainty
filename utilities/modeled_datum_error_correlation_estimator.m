% Example of using known modeled datum errors at gages 
% for estimating how these errors correlation coefficient  
% depends on distance (of any kind) between the gages, and
% for testing hypotheses about how this coefficient can be 
% approximated with some prescribed function(s)  
%----------------------------------------------------------
% See function correlation_estimator below the example
%----------------------------------------------------------
% variables:
% m, int32 - number of gages
% mdist, m x m double - matrix of gage2gage distances returned by wdist 
%              (measured on the mesh along the cell edges)
% gdist, m x m double - matrix of gage2gage geophysical distances
%              (measured along the Great circle)
% cid, 1 x m int32 - string of gages' IDs
% mhh,mhw,mlw,mll, 1 x m single - observed datums
% lv, 8 x m double - modeled datums at gages as in TAD output, 
%           column: 1-#, 2-mhhw,3-mhw, 4-msl,5,6-dtl/mtl,7-mlw,8-mllw  
% whatelse, m x m double - matrix of whatever factor is being tested  
%     for affecting the error correlation (for instance, coef(i,j) - 
%     correlation of tidal envelopes), set to ones(m) is none used
% eestn, 4 x m double - rows of known model errors for mhh,mhw,mlw,mll
% xmax - max distance to look for any meaninfull correlation
% nb - number of bins to divide [0 xmax] into
% xbins, 1 x nb double - bin edges on the right side
%------ output of correlation_estimator -------------
% nnn, 1 x nb int32 - number of products in a bin
% xxx, 1 x nb double - average gage2gage distance in a bin
% yyy, 4 x nb double - error correlation coefficient estimate
%                                      in a bin, per datum
% ---------------------------------------------------
% L - correlation distance as estimated from fitting yyy vs xxx
%%%%%%%%%%%%%%%   elena.tolkova@noaa.gov   %%%%%%%%%%%%

%*** gage order must be the same in all data structures ***

clear variables

path_grid='C:\Users\elena.tolkova\Documents\EastCoast\jet_run_EC1EC2DEv9\';
eval(['load ' path_grid 'pmoe_datum'])

cid=[pmoe_datum.id];     % station id
i=find(cid==8410140);  % distance map missing
pmoe_datum(i)=[];
cid(i)=[];
nodes=int32([pmoe_datum.node]);  % grid node # at station location
m=length(nodes)
% -------------------- gage 2 gage distances, wdist -----------------
mdist(1:m,1:m)=NaN;
pathin='python_dist/';
for i=1:m
    iid=cid(i);
    infile=[pathin 'station' int2str(iid) '.nc'];
    dis=ncread(infile,'dis'); % water distance from i-th station
    mdist(i,:)=dis(nodes);
end
msk=(mdist>2000);
mdist(msk)=NaN;
% -------------------- observed datums at gages -----------------
mhh= [pmoe_datum.omhhw];
mhw= [pmoe_datum.omhw];
mlw= [pmoe_datum.omlw];
mll= [pmoe_datum.omllw];
% -------------------- calculated datums at gages -----------------
ver='v21';
path_lv=['C:\Users\elena.tolkova\Documents\EastCoast\jet_run_EC1EC2DE' ver '\'];
temp=dir([path_lv 'mpdatums/']);
nfiles=length(temp)

lv(1:8,1:m)=NaN;
ngfound=0;
for n=3:nfiles
    fname=[path_lv 'mpdatums\' temp(n).name];
    nds=ncread(fname,'stationN');
    xx=ncread(fname,'datums');
    [Cnds,IA,IB]=intersect(nodes,nds);
    lv(1,IA)=cid(IA);
    lv(2:8,IA)=xx(:,IB);
    ngfound=ngfound+length(IA);
    clear xx nds fname Cnds IA IB
 end
if ngfound ~= m
    keyboard
end
%---------------- convert model datums from Model Zero to MSL ----------
lvmsl=lv(4,:);
msk=(lvmsl<-999)|isnan(lvmsl);
lv(2:8,msk)=NaN;
for i=2:8
    lv(i,:)=lv(i,:)-lvmsl;
end
% --------------------- datum errors ----------------------------
eestn(1:4,1:m)=NaN;    % model datum error at stations
eestn(1,:)=mhh-lv(2,:);
eestn(2,:)=mhw-lv(3,:);
eestn(3,:)=-(mlw-lv(7,:));
eestn(4,:)=-(mll-lv(8,:));
% ---------------- bins for averaging corr products -------------
xmax=1200; %largest distance at which correlation calculation stops
nb=31; 
xbins=logspace(0,2,nb-1); xbins=(xmax/100)*(xbins-1)+1;
xbins=[xbins xmax];
dx=(xmax-xbins(27))/4;
for i=28:nb
    xbins(i)=xbins(i-1)+dx;
end

whatelse=ones(m);

[nnn,xxx,yyy]=correlation_estimator(m,eestn,mdist,whatelse,xbins,1);

% mark bins containing fewer than nmin entries for removal in ex,ey
clear msk
minbin=1;
msk=(nnn<minbin);
xxx(msk)=NaN;
yyy(1:4,msk)=NaN;

ex=[];
ey=[];
d2fit=200; % distance within which to fit 
[tmp K]=min(abs(xxx-d2fit));

figure
hold on
for jj=1:4
    plot(xxx,yyy(jj,:),'.','markersize',20)
    ex=[ex, xxx(1:K)];
    ey=[ey, yyy(jj,1:K)];
end
clear msk
msk=isnan(ex);
ex(msk)=[];
ey(msk)=[];

xlim([0 xmax])
grid on
ax=gca;
ax.FontSize=12;
xlabel('km','fontsize',14);

% fit correlation coefficient estimate (ey vs ex) with L/(x+L) 
Ke=length(ex);
LL=0:0.1:100;
nL=length(LL);
clear r rr
r(1:nL)=0;
rr(1:Ke)=0;
for i=2:nL
    L=LL(i);
    for j=1:Ke
        rr(j)=abs(ey(j)-L/(ex(j)+L));
    end
    r(i)=sum(rr);
end
[tmp iL]=min(r(2:end));
L=ceil(LL(iL))
xr=0:0.1:500;
for j=1:length(xr)
    rr(j)=L/(xr(j)+L);
end
plot(xr,rr,'linewidth',2)
xlim([0 xxx(K+5)])
grid on
legend('mhh','mhw','mlw','mll','L/(x+L)')
title('error correlation estimate vs mdist')
%----------------------------------------------------------
% this function uses known modeled datum errors at gages 
% for estimating how these errors correlation coefficient  
% depends on distance (of any kind) between the gages, and
% for testing hypotheses about how this coefficient can be 
% approximated with some prescribed function(s)    

function [nnn,xxx,yyy]=correlation_estimator(m,eestn,dist,whatelse,xbins,iplotprod)

    nb=length(xbins);
    if iplotprod  % plot error products for m(m-1)/2 gage pairs, mhhw
        
        z=eestn(1,:);
        a=mean(z)
        estn =z-a;
        m_error = std(estn);
        estn=estn/m_error;
        ecorr=estn'*estn;

        figure
        hold on
        xlim([0 xbins(end)])
        for i=1:m
            plot(dist(i,i:m),ecorr(i,i:m),'.','markersize',9)
        end    
        xticks(xbins)
        xlb=cell(1,nb);
        for i=6:6:nb
            xlb{i}=int2str(round(xbins(i)));
        end
        ax=gca;
        ax.XTickLabels=xlb;
        xtickangle(45)
        ax.TickLength=[0.05 0.025]
        yticks([-20:10:30])
        ax.FontSize=12;
        xlabel('km','fontsize',14)
        title('error products vs mdist')
    end

%                      %  placeholders for:
    xxx(1:nb)=NaN; % average g2g distance within a bin
    nnn(1:nb)=int32(0); % number of values within a bin
    yyy(1:4,1:nb)=NaN; % mean ecorr within a bin, for each of 4 datums

    for jj=1:4
        z=eestn(jj,:);
        a=mean(z);
        estn = z-a;
        m_error = std(estn);
        estn=estn/m_error;
    % ecorr(i,j)=product of normalized jj-th datum errors at i-th and j-th gages
        ecorr=estn'*estn; 
    % divide ecorr element-wise by hypothesized correlation factor whatelse
        ecorr=ecorr./whatelse;        
    % remove duplicates and reshape into a single row
        for i=2:m
            dist(i,1:(i-1))=NaN;
            ecorr(i,1:(i-1))=NaN;
        end
        xx=reshape(dist,1,m*m);
        yy=reshape(ecorr,1,m*m);
        msk=isnan(xx);
        xx(msk)=[];
        yy(msk)=[];
    %                      %  placeholders for:
        xxx(1:nb)=NaN; % average g2g distance within a bin
        nnn(1:nb)=int32(0); % number of values within a bin
        yyy(jj,1:nb)=NaN; % sum of ecorr within a bin, jj-th datum
        for i=1:nb
            clear msk
            msk=xx<xbins(i);
            k=sum(msk);
            if k>0
                nnn(i)=k;
                xxx(i)=sum(xx(msk))/k;
                yyy(jj,i)=sum(yy(msk))/k;
                xx(msk)=[];
                yy(msk)=[];
            end
        end
    end

end

%--------------------- SCRATCH ----------------------------------------
% % ------------   gage 2 gage distances, Great circle -----------------
% gdist(1:m,1:m)=0;
% % gx=[pmoe_datum.ox]*pi/180;
% % gy=[pmoe_datum.oy]*pi/180;
% load allfort14
% p=[allgrd.nodes];
% gx=p(nodes,2)*pi/180;
% gy=p(nodes,3)*pi/180;
% clear allgrid p
% Re=111.111*180/pi;
% for i=1:m
%     ln1=gx(i);
%     lt1=gy(i);
%     for j=(i+1):m
%         ln2=gx(j);
%         lt2=gy(j);
%         alp=acos(cos(lt1)*cos(lt2)*cos(ln1-ln2)+sin(lt1)*sin(lt2));
%         d=Re*alp;
%         gdist(i,j)=d;
%         gdist(j,i)=d;
%     end
% end
% % ----------- wdist/gdist distribution ----------
% % remove duplicates and reshape into a single row
% for i=2:m
%     mdist(i,1:(i-1))=NaN;
%     gdist(i,1:(i-1))=NaN;
% end
% xx=reshape(mdist,1,m*m);
% yy=reshape(gdist,1,m*m);
% msk=isnan(xx);
% xx(msk)=[];
% yy(msk)=[];
% figure
% hold on
% plot(xx,yy,'.','markersize',3)
% grid on
% % take care of zero for calculating ratios
% clear msk
% msk=(xx<0.01);
% xx(msk)=1;
% yy(msk)=1;
% rr=yy./xx;
% keyboard
