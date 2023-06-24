%=======================================================================
% This code is a variation of same-name code in S4_svn by Liujuan.Tang
% to assemble matrix 'covar' of three optional matrices - correlation of 
% tidal envelopes, distance mask, and datum ratio mask - according to masks
%
% Output:  pre_process/
%           Pmatrix_[varname]_masksNNN.mat 
% Variables
%   masks:  vector of flags for build 'covar' 
%   n:      number of mesh nodes
%   m:      number of tide stations used for SVU computation
%                           Elena.Tolkova@noaa.gov
%---------------------------------------------------------------------
% 0. Input section
clear variables
masks=[0 4 0];
i=1;
    imask_cvr=masks(i,1);
    imask_dst=masks(i,2); % e-fold length Lr=imask_dst*55.5(km)
    imask_dm=masks(i,3);

adcirc_dir='C:\Users\elena.tolkova\Documents\TXmodel\jet_runTXLAMS_v1\';
% --- do for each datum, or just once if only imask_dst is not zero ---
varname='mhw'; %{'mhhw';'mhw';'mlw';'mllw'}
%---------------------------------------------------------------------
rmsok=1; % 1 Assign rms to tide stations without rms
         % 0 Do not use tide stations without rms

path_pre='pre_process/'; % directory to save the pre_processed data
if ~exist(path_pre,'dir')
    eval(['mkdir ' path_pre])
end
%---------------------------------------------------------------------
    fprintf('\n----------------------------\n',varname);
    fprintf('Preparing SVU input file for %s\n',varname);
%----------- number of mesh nodes ------------------------------------
n=657497
    
%1. Station information; check stations with rms;    
eval(['load ' adcirc_dir 'pmoe_datum'])
    
if rmsok~=1
    loc_norms=[];
    orms=[pmoe_datum.orms]';
    loc_norms=find(isnan(orms) | abs(orms)>10);
    orms_valid=orms;
    orms_valid(loc_norms)=[];
    [maxorms,iloc]=max(orms_valid);
    fprintf('%d stations has  rms.\n max=%.4f (m); \n',...
        length(orms_valid),maxorms);
    pmoe_datum(loc_norms)=[];
end
    
bcpsi=[pmoe_datum.node]';  % grid node # at station location
m=length(bcpsi)
%---------------------------------------------------------------------
covar=ones(n,m);
%---------------------- distance mask --------------------------------
    if imask_dst
    % 2. Load water distance from stations
    %   Prepare mdist
        run4dist='../'; 
        pathin=[run4dist 'S1_wdist/output_dist_time_land_removed/'];
        cid=[pmoe_datum.id];
        mdist(1:n,1:m)=NaN;
        for i=1:m
            iid=cid(i);
            infile=['OWdist' int2str(iid) '_2.mat'];
            eval(['load ' pathin infile])
            mdist(:,i)=dis;
        end
        if imask_dst<100
            Lr = 55.5*imask_dst;
        else
            Lr=55.5/(imask_dst-100);
        end
        dsmsk = exp(-mdist/Lr);
        excl=isnan(mdist);   % mdist is NaN at dry nodes
        dsmsk(excl)=0;
        covar = covar.*dsmsk;
    end
    % -------------------- datum-based mask -----------------
    if imask_dm
        
        dmmsk=ones(n,m);
        temp=['load ' adcirc_dir  'mp_datums_v2_newtad_xx.mat'];
        fprintf(1,'%s\n',temp);
        eval(temp);
        dmall=testxx;
        clear testxx
        dm=dmall(:,3)-dmall(:,7); % mhw-mlw range
        dm(dm<0.05)=0.05;
        switch imask_dm
            case 1
                pow=1;
            case 2
                pow=0.5;
        end

        for j=1:m
            dmo=dm(bcpsi(j));  % datum at j-th gage
            for i=1:n
                a=(dm(i)/dmo)^pow;
                dmmsk(i,j)=min([a, 1/a]);
            end
        end
        clear dmall dm
        covar=covar.*dmmsk;
    end % imask_dm
    %-------------------covariance mask -----------------------------
    if imask_cvr
    %6. prepare covar  n by m
        coef(1:n,1:m)=NaN;
        corrdir=[upper(varname) '/']        
        for i=1:m
            eval(['load ' '../S3_coeff/' corrdir 'coef' varname  '_st' int2str(i) '.mat']);
            coef(:,i)=stcoef;
        end
        covar = covar.* coef;
    end
%---------------------------------------------------------------------
    if sum(sum(isnan(covar)))
        keyboard
    end
%---------------------------------------------------------------------
    %7. Save as SVU input file in .mat format
    if imask_cvr==0
        outfile=['Pmatrix_masks' int2str(imask_cvr)  int2str(imask_dst)  int2str(imask_dm) ];
    else    
        outfile=['Pmatrix_' varname '_masks' int2str(imask_cvr)  int2str(imask_dst)  int2str(imask_dm) ];
    end
     fprintf('Save to %s%s\n', path_pre,outfile)
    iexcl=find(excl(1:n)==1);
    eval(['save  ' path_pre outfile ' covar rmsok loc_norms iexcl']);

%------------------ plotting a mask for a gage #idot ----------
iplot=1;
if iplot
eval(['load ' adcirc_dir 'allfort14.mat'])  % load grid
    p=allgrd.nodes;  
    p(:,1)=[];
    p(:,end)=[];
    q=allgrd.elements;  
    q(:,1:2)=[];
    clear allgrid
idot=100;
xdot=pmoe_datum(idot).ox;
ydot=pmoe_datum(idot).oy;
d=covar(:,idot); %mdist(:,idot);  %
    figure
    h=trimesh(q,p(:,1),p(:,2),d);
    view([0 90])
    hold on
    plot3(xdot,ydot,10,'xr')
    colorbar
end    