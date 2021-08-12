%=======================================================================
% This program is to prepare input for spatially varying uncertainty 
%   computation
% Input files:
%   S1/station_node_info.mat
%   S2/*testxx.mat
%   S3_coef/saved/coef_sf_mhhw_runname.mat

% Output:  pre_process/
%           SVU_input_mhhw_m_*.mat 

% Variables
%   n:      number of mesh nodes
%   m:      number of tide stations used for SVU computation
% 
%
% dm:       model datum, n by 1
% dmsl:     model mean sea level from lv8j_mp, n x 1

% coef:     correlation coefficent, n by m
% uncer:    rms for obs at tide stations, m x 1 
%
% mstn:     model datum at station nodes , m x 1
%                           Liujuan.Tang@noaa.gov
%                           01/29/2020
%---------------------------------------------------------------------
% 0. Input section
clear
runid='R58_k6s4_msl_5o2_a53_merged';

rmsok=1; % 1 Assign rms to tide stations without rms
         % 0 Do not use tide stations without rms

path_pre='pre_process/'; % directory to save the pre_processed data

varlist={'mhhw';'mhw';'mlw';'mllw';'mtl' ;'dtl'};
rundatum=1:6%;1:6;  % 1:6 run the six datums in varlist
             % Set to 1 to test one datum first 
%---------------------------------------------------------------------

if ~exist(path_pre,'dir')
    eval(['mkdir ' path_pre])
end
%---------------------------------------------------------------------
% Start loop for datums
for idatum=rundatum 
    varname=varlist{idatum};
    fprintf('\n----------------------------\n',varname);
    fprintf('Preparing SVU input file for %s\n',varname);
    %-------------------------------------------------------
    %1. Station information; check stations with rms;
    %   select stations for svu computation 
    load ../S1_wdist/station_node_info.mat
    orms=[pmoe_datum.orms]';
    loc=find(isnan(orms) | orms<-9);
    orms_orig=orms;
    orms_orig(loc)=[];
    [maxorms,iloc]=max(orms);
    fprintf('%d stations has  rms.\n max=%.4f cm; min=%.4f\n',...
    length(orms)-length(loc),maxorms,min(orms_orig));
    
    if rmsok~=1
        pmoe_datum(loc)=[];
    end
    %-------------------------------------------------------------------
    %1.1 manually remove stations if needed
    rmstation=[
         8728171 % Outlier River station far away from NEGOM VDatum region
         %8727989 % Aucilla River
    ];
    cid0=[pmoe_datum.id];
    if ~isempty(rmstation)
        for i=1:length(rmstation)
            iid=rmstation(i);
            loc2=find(cid0==iid);
            cid0(loc2)=[];
            pmoe_datum(loc2)=[];
        end
    end
    nid=length(pmoe_datum);
    fprintf('%d stations excluded; %d station left\n',length(rmstation),nid);
    %--------------------------------------------------------
    %   1.2 Check and remove duplicate station nodes if possible
    cid=[pmoe_datum.id];
    cnode=[pmoe_datum.node];
    [c2,ia2,ic2] = unique(cnode,'stable');
    b=diff(ic2);
    loc=find(b~=1);
    b1=cnode(loc);
    rmnode=[];
    for i=1:length(b1)
        inode=b1(i);
        loc2=find(cnode==inode);
        rmnode=[rmnode loc2(2:end)];
    end
    pmoe_datum0=pmoe_datum;
    pmoe_datum(rmnode)=[];
    nid=length(pmoe_datum);
    eval(['save ' path_pre 'pmoe_datum_' int2str(nid) ' pmoe_datum']);
    fprintf('%d unique node stations. \n',nid);
    fprintf('----------------------------\n',varname);
    %------------------------------------------------------
    bcpsi=[pmoe_datum.node]';
    %--------------------------------------------------
    % 2. Load water distation from stations
    %   Prepare mdist
    pathin='../S1_wdist/output_dist_time/';
    cid=[pmoe_datum.id];
    mdist=[];
    for i=1:length(cid)
        iid=cid(i);
        infile=['OWdist' int2str(iid) '_2.mat'];
        eval(['load ' pathin infile])
        mdist(:,i)=dis;
    end
    mdist=mdist/111;
    %------------------------------------------------------
    %3. Process model dm and observation ostn
    temp=['load ../S2_lv8jmp/' runid '_testxx.mat'];
    fprintf(1,'%s\n',temp);
    eval(temp);
    dmall=testxx;

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
        case 'mtl'
            dm=dmall(:,5);
            ostn=([pmoe_datum.omhw]+[pmoe_datum.omlw])'/2;        
        case 'dtl'
            dm=dmall(:,6);
            ostn=([pmoe_datum.omhhw]+[pmoe_datum.omllw])'/2;
    end
    
    dmsl=dmall(:,4); % model mean sea level
    if idatum<5
        eval(['ostn=[pmoe_datum.o' varname ']'';']); %144x1
    end
    %b1_prep_model_data_v3
    b1_proc_model_datum
    %----------------------------------------------
    %4. Prepare observation rms for stations
    uncer=[pmoe_datum.orms]';
    if rmsok==1
        b2_proc_obs_rms
    end
    fprintf(1,'----------------------------\n');
    %----------------------------------------------------
    % 5. prepare other input parameters
    mstn=dm(bcpsi);  % m x 1
    estn = ostn-mstn;   % m x 1
    uncer=max(uncer,0.0005);
    o_error = sqrt(mean(uncer.*uncer));
    m_error = std(estn); % 1 by 1
    n=length(dm);
    m=length(bcpsi);
    nc = n;
    nstn = m;
    %------------------------------------------------------
    %6. prepare covar  n by m
    Lr = 2.0;
    datt = exp(-mdist/Lr); 
    adist=100.0*Lr;
    eval(['load ../S3_coef/saved/coef_sf_' varname '_' runid '.mat'])
    coef0=coef;
    coef=[];
    % re-combine coef based on the tide gages used
    for i=1:length(bcpsi)
        iid=bcpsi(i);
        loc=find(bcpsi_coef==iid);
        coef(:,i)=coef0(:,loc);
    end
    covar = m_error*m_error*datt.* coef;
    %----------------------------------------------------
    %7. Save as SVU input file in .mat format
    outfile=['SVU_input_' varname '_' int2str(nid) '_' runid];
    fprintf('Save to %s%s\n', path_pre,outfile)
    %mdist, coef,datt not used by tidal_datum_svu
    %eval(['save  ' path_pre outfile ' bcpsi mstn ostn estn uncer o_error m_error dm n m nc nstn coef mdist Lr datt covar varname']);
    eval(['save  ' path_pre outfile ' bcpsi mstn ostn estn uncer o_error m_error dm n m nc nstn Lr covar varname loc_model_nan']);
end