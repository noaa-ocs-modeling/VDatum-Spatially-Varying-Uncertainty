%=======================================================================
% This program is to prepare input for spatially varying uncertainty 
%   computation
% Input files:
%   pmoe_datum.mat, allfort14, mp_datums_v2_newtad_xx.mat
%   Pmatrix_[varname]_masksNNN.mat, NNN is given by masks 
%
% Variables
%   n:      number of mesh nodes
%   m:      number of tide stations used for SVU computation
% dm:       model datum, n by 1
% covar:    n by m
% uncer:    rms for obs at tide stations, m x 1 
%     rearranged from code by Liujuan.Tang@noaa.gov
%                           06/2020-06/2023
%---------------------------------------------------------------------
% 0. Input section
clear variables
masks=[0 4 0];
i=1;
    imask_cvr=masks(i,1);
    imask_dst=masks(i,2);
    imask_dm=masks(i,3);
no_model_dm=0;  %if ~0, model dm emulates observations

adcirc_dir='C:\Users\elena.tolkova\Documents\TXmodel\jet_runTXLAMS_v1\';

isave=1;

eval(['load ' adcirc_dir 'pmoe_datum'])
bcpsi=[pmoe_datum.node]';  % grid node # at station location
m=length(bcpsi)  % number of gages used for processing
%---------------------------------------------------------------------
path_pre='pre_process/'; % directory with covar
varlist={'mhhw';'mhw';'mlw';'mllw'};

idatum=1; % select datum type here 
varname=varlist{idatum}

if imask_cvr==0
    infile=['Pmatrix_masks' int2str(imask_cvr)  int2str(imask_dst)  int2str(imask_dm) ];
else    
    infile=['Pmatrix_' varname '_masks' int2str(imask_cvr)  int2str(imask_dst)  int2str(imask_dm) ];
end
    temp=['load ' path_pre infile '.mat'];        
fprintf(1,'%s\n',temp);
eval(temp);

if rmsok~=1
    pmoe_datum(loc_norms)=[];
end
%------------------------------------------------------
    %3. Process model dm and observation ostn
    temp=['load ' adcirc_dir  'mp_datums_v2_newtad_xx.mat'];
    fprintf(1,'%s\n',temp);
    eval(temp);
    dmall=testxx;  clear testxx
    [n ntmp]=size(dmall);  % n - number of nodes

    dm=[];
    switch varname
        case 'mhhw'
            dm=dmall(:,2);
        case 'mhw'
            dm=dmall(:,3);
        case 'mlw' 
            dm=dmall(:,7);
        case 'mllw'
            dm=dmall(:,8);
    end
    dmsl=dmall(:,4); % model mean sea level
    dm=dm-dmsl;     % datums relative to MSL
%--- replace 0 on land with values interpolated from nearby waters
    dm_good=double(dm);
    dm_good(iexcl)=[];
    dm(iexcl)=NaN;
    
    eval(['load ' adcirc_dir 'allfort14.mat'])  % load grid
    p=allgrd.nodes;  
    p(:,1)=[];
    p(:,end)=[];
    q=allgrd.elements;  
    q(:,1:2)=[];
    clear allgrid
    
    x_good=p(:,1);
    xnan=x_good(iexcl);
    x_good(iexcl)=[];
    y_good=p(:,2);
    ynan=y_good(iexcl);
    y_good(iexcl)=[];
    
    dmnan=griddata(x_good,y_good,dm_good,xnan,ynan,'nearest');
    dm(iexcl)=dmnan;
fprintf(1,'max_loc=%.4f m; min_loc=%.4f m\n',max(dmnan),min(dmnan))
fprintf(1,'max=%.4f m; min=%.4f m\n',max(dm),min(dm))
clear x_good xnan y_good ynan dmnan dm_good
% --------------------------------------------------------    
    dmsave=dm;   
% --------------------------------------------------------    
    if no_model_dm
        dm(1:n,1)=0;
        ostn=dmsave(bcpsi);  
        uncer(1:m,1)=0.005;
    else
       
        switch varname
            case 'mhhw'
                ostn=[pmoe_datum.omhhw]';
            case 'mhw'
                ostn=[pmoe_datum.omhw]';
            case 'mlw' 
                ostn=[pmoe_datum.omlw]';
            case 'mllw'
                ostn=[pmoe_datum.omllw]';
        end
        
        uncer=[pmoe_datum.orms]';
        if rmsok==1
            uncer(isnan(uncer))=0.055;
        end
        uncer=max(uncer,0.005);
    end  % if no_model_dm
    %------------------------------------------------------
    % 5. prepare other input parameters
    mstn=dm(bcpsi);  % m x 1   model datum at stations
    estn = ostn-mstn;   % m x 1  model datum error at stations
    m_error = std(estn); % 1 by 1 stand. dev. of estn
    %------------------------------------------------------
    covar = (m_error*m_error)*covar;
 %---------------  moved from tidal_datums_svu  -------------------    
W01=eye(m);
mc=0;
while 1
    mc=mc+1;
    [dd,G]=stat_interp_uncertainty_jackknife(bcpsi,estn,uncer,W01,covar);
    d=dm+dd;           % corrected datum;
    d_jn=d(bcpsi);      % corrected datum at stations
    dd_jn=dd(bcpsi);   % corrections to the datum at stations
    change=zeros(m,1);
    for i=1:m
        if abs(d_jn(i)-ostn(i)) > min(0.01, uncer(i)) 
               W01(i,i)=W01(i,i)/sqrt(2.0);
               change(i)=1;
        end
    end
    
    fprintf(1,'mc=%d; sum of change is %.f \n',mc,sum(change)); 
    if sum(change) == 0 ,   break;  end
end

% here we calculate the uncertainty
diaP=m_error*m_error*ones(n,1);
diaPA=zeros(n,1);
% uncertainty estimate
R=eye(m);
for i=1:m
    R(i,i)=uncer(i)*uncer(i)*R(i,i);
end
HPH=covar(bcpsi,:); 
diaPA=sqrt(diaP-2.0*sum(G.*covar,2)+sum((G*(R+HPH)).*G,2));
diaPA_jn=diaPA(bcpsi);

if isave
    dm=dmsave;
    outfile=['d_dd_diaPA_' varname  '_masks' int2str(imask_cvr)  int2str(imask_dst)  int2str(imask_dm) ]   
    eval(['save ' path_pre outfile '   d dd diaPA dm change d_jn dd_jn diaPA_jn'])
end
%---- ANALIZE RESULT  --------------

figure
hold on
plot(ostn,'xr')      % observed datum at stations
plot(d(bcpsi),'-.g')  % model datum after corrections
plot(dm(bcpsi),'-.b') % original model datum
plot(diaPA(bcpsi),'-.r')
grid on
legend('obs dm','corrected dm','model dm', 'uncert','location','best')

% ---------------- nontidal BPs -------------------------
shpfile='C:\Users\elena.tolkova\Documents\TXmodel\BPs\nontidal_more'
C=m_shaperead(shpfile);
nc=length(C.ncst)
msk(1:n,1)=0;
for i2=1:nc
    clear cx cy
    cx=C.ncst{i2}(:,1);
    cy=C.ncst{i2}(:,2);
    msk=msk | inpolygon(p(:,1),p(:,2),cx,cy);
end
%check corrected datum - remove nontidal area from check
 switch varname(2)
            case 'h'
                loc_sgn=find(d(~msk)<0);
                fprintf(1,'\n Found %d points with negative hw\n',length(loc_sgn))
            case 'l'
                loc_sgn=find(d(~msk)>0);
                fprintf(1,'\n Found %d points with positive lw\n',length(loc_sgn))
 end

iplot=1;
if iplot

    load ../../../DEdelches01/comp_svu_v29/S4_svu/cm_uncertainty
    figure
    h=trimesh(q,p(:,1),p(:,2),d);
    view([0 90])
    switch varname
        case {'mhhw','mhw'}
            caxis([0 0.5])
            colormap(cm)
        case {'mllw','mlw'}
             caxis([-0.5 0])
             colormap(flip(cm))
    end
    colorbar
    title(varname)
    hold on
   %----- non-tidal areas ------- 
   for i2=1:nc
        clear cx cy
        cx=C.ncst{i2}(:,1);
        cy=C.ncst{i2}(:,2);
        plot3(cx,cy,100*ones(size(cx)),'r','linewidth',1)
   end
   
    figure
    h=trimesh(q,p(:,1),p(:,2),d-dm);
    view([0 90])
    switch varname
        case {'mhhw','mhw'}
            %caxis([0 0.5])
            colormap(cm)
        case {'mllw','mlw'}
             colormap(cm)
    end
    colorbar
    title(['corrected-original ',infile])
   %100*norm(d-dm)/sqrt(length(p))  % cm

end



 

