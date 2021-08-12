%------------------------------------------------------------------------
% This program is to plot the product datums and SVU
% Input files:
%   ../S1_wdist/station_node_info.mat
%   out_mat/d_dd_diaPA_xxxx*.mat  %xxxx: datum
%   pre_process/SVU_input_xxxx_181_R58_k6s4_msl_5o2_a53_merged.mat

% Ouput files:
%   out_fig*/*.png
%   out_grd/*.grd

% Suggest to test one datum (e.g. mhhw) first before going into 
% the loop for all six datums.
% Please adjust figure size, position, color scale for your case.

%                               liujuan.tang@noaa.gov
%                               Last modified: 01/29/2020
%------------------------------------------------------------------------
clear


%-----------------------------------------------------------------------
runid='181';  % number of tide stations used for SVU
runname='R58_k6s4_msl_5o2_a53_merged'; % ADCIRC run ID

grdok=0;   % Output original datum and corrected datum in ADCIRC grd format
           %1=yes; 0=no
printok=0; % save figures ? 1=yes; 2=no

path_pre='pre_process/';
path_png=['out_fig_' runid '_' runname ];
   if ~exist(path_png,'dir')
        eval(['mkdir ' path_png])
   end
if grdok==1
    path_grd=['out_grd/' ];
   if ~exist(path_grd,'dir')
        eval(['mkdir ' path_grd])
   end
end
%--------------------------------------------------
% Write to log file
delete log.out % remove exit file
diary log.out % screen output to file
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
fprintf(1,'------------------------------------------\n');
fprintf(1,'%s\n%s\n',t,runname);
%------------------------------------------------------------------------
% Setup figure size; You may need to adjust based on your grid
%text color
figwh= [ 750 500]*1.5;      %Figure width and height
gcf_po=[10   20  figwh];    %Figure position
ax_po=[0.05 0.08 0.9 0.85]; %Axis position

ta=20;tt=30;
ck=[ 0 0 0]+0;
bkc=[0 0 0]+1;
load cm_uncertainty % colormap

datumlist={'mhhw';'mhw';'mlw';'mllw';'mtl' ;'dtl'};
varlist={'d';'dd';'diaPA'}; % see vartest 
vartext={'Product' ;'Correction ';'Uncertainty'};

for ii=4 % datums 1-6
    close all
    iivarname=datumlist{ii};
    load ../S1_wdist/station_node_info.mat
    p=all.nodes;  
    p(:,1)=[];
    p(:,end)=[];
    %-------------
    q=all.elements;  
    q(:,1:2)=[];
    fprintf('----------------%s-----------------\n',iivarname)
    
    eval(['load ' path_pre 'pmoe_datum_' runid '.mat']);
    infile=['out_mat/d_dd_diaPA_' iivarname  '_' runid '_' runname];
    eval(['load(''' infile ''');']); 
    for j=1:3
        fprintf('\n%d) Plotting  %s, %s\n',j,varlist{j},vartext{j});

        figure(j)
        clf
        set(gcf,'Position',gcf_po,'visible','on','color',bkc,'paperpositionmode','auto','InvertHardcopy','off')
        axes('position',ax_po)
        
        eval(['a=' varlist{j} ';'])
        amin=min(min(a(:,:)));
        fprintf(1,'max=%.4f m; min=%.4f m\n',max(max(a(:,:))),min(min(a(:,:))))
        sf=sprintf('max=%.4f m; min=%.4f m\n',max(max(a(:,:))),min(min(a(:,:))));
        colormap(cm)

        %---------------------------
        % Check the corrected datum at nan model points

        nanfile=['pre_process/mdatum_nanloc_' iivarname '.mat' ];
        %nanfile=['pre_process/mdatum_nanloc_mhhw.mat' ];
        load(nanfile, 'loc');
        loc_mnan=loc;
        fprintf(1,'\nCheck %d no model datum points:\n',length(loc_mnan))
        fprintf(1,'max=%.4f m; min=%.4f m\n',max(a(loc_mnan)),min(a(loc_mnan)));
%         fprintf('Replace %d nan datum points with the nearest corrected datums\n',length(loc_mnan))
%         a0=a;
%         a=f_replace_model_nan(all.nodes(:,2:3),a,loc);
%         fprintf(1,'max=%.4f m; min=%.4f m\n',max(max(a(:,:))),min(min(a(:,:))))
%         sf=sprintf('max=%.4f m; min=%.4f m\n',max(max(a(:,:))),min(min(a(:,:))));

        %--------------------------
%         %check corrected datum
%         switch iivarname(2)
%             case 'h'
%                 loc=find(a<0);
%             case 'l'
%                 loc=find(a>0);
%         end
%         hold on
%         plot3(p(loc,1),p(loc,2),p(loc,1)*0+10,'+','color','k');

        %end
        %-----------------------------
        % Check  model points
        loc=find(imag(a)~=0);
        if ~isempty(loc)
            fprintf(1,'Complex number found for %d points\nCheck small value points.\n',length(loc))
            a(loc)=nan;
            hold on
            fprintf(1,'max=%.4f m; min=%.4f m\n',max(max(a(:,:))),min(min(a(:,:))))
            sf=sprintf('max=%.4f m; min=%.4f m\n',max(max(a(:,:))),min(min(a(:,:))));
        end
        %axis equal
        trisurf(q,p(:,1),p(:,2),a);   shading interp,  view([0 90])
        plot3(p(loc,1),p(loc,2),p(loc,1)*0+10,'+','color','k');
        %-----------------
        set(gca,'color',bkc)
        grid off,box on
        hb=colorbar;                    
        set(hb,'fontsize',18,'color',ck)
        axis tight
        if ii<5
            % set up color range.
            switch j
                case 1
                    %temp=sort((-1)^rem(ii+1,2)*[0.05 0.4]);
                    %  temp=[-1 1]*.4;
                    temp=[-1 1]*.4;
                    caxis(temp);
                case 2
                    caxis([-1 1]*.4);
                case 3
                    caxis([0.006 0.05]);
                    if ii==2
                        % caxis([0.006 0.16]);
                    end
            end
        else
            switch j
                case 1
                    %temp=sort((-1)^rem(ii+1,2)*[0.05 0.4]);
                    %  temp=[-1 1]*.4;
                    temp=[-1 1]*.4;
                    caxis(temp);
                case 2
                    caxis([-1 1]*.2);
                case 3
                    caxis([0.006 0.05]);
                if ii==2
                    % caxis([0.006 0.16]);
                end
            end
        end
        hold on
        xy=[pmoe_datum.ox; pmoe_datum.oy]';
        plot3(xy(:,1),xy(:,2),xy(:,1)*0+10,'o','color',[0 0 0]);
        xlabel(sf,'fontsize',20,'color',ck)
        title([vartext{j} ' ' upper(iivarname) ' ' varlist{j} ],'fontsize',tt,'color',ck)
        if printok==1
            eval(['print -dpng -r300 out_fig_' runid '_' runname '/f_SVU_' iivarname '_' varlist{j} '_' runname '.png'])
        end
        % Write to grd format for populating
        if j==1
            if grdok==1
                outfile=[path_grd 'corr_' iivarname  runid '_' runname '.grd'];
                all1=all;
                all1.nodes(:,4)=a;
                f_write_bathy_fort14(all1,outfile,0);
            end
        end
    end
%-----------------
%ip
fprintf('\n4) Plot and write  original model datum\n');
figure(4)
clf
set(gcf,'Position',gcf_po,'visible','off','color',bkc,'paperpositionmode','auto','InvertHardcopy','off')
axes('position',ax_po)
temp=['load(''' path_pre 'SVU_input_' iivarname '_' runid '_' runname ''', ''dm'')'];
fprintf('%s\n',temp)
eval(temp);
%trisurf(q,p(:,1),p(:,2),(dm*(-1)^rem(ii+1,2)));   shading interp,  view([0 90])
trisurf(q,p(:,1),p(:,2),dm);   shading interp,  view([0 90])

    hb=colorbar;                    set(hb,'fontsize',ta,'color',ck)
    axis tight
    %caxis([0.009 0.06]);
    colormap(cm)

    sf=sprintf('max=%.4f m; min=%.4f m\n',max(max(dm(:,:))),min(min(dm(:,:))));    
    fprintf(1,'max=%.4f m; min=%.4f m\n',max(max(dm(:,:))),min(min(dm(:,:))))
    %axis equal
    xlabel(sf,'fontsize',20,'color',ck)
    set(gca,'color',bkc)
    grid off,box on
    %temp=sort((-1)^rem(ii+1,2)*[0.0 0.4]);
    temp=[-1 1]*.4;
    caxis(temp);
    hold on
    plot3(xy(:,1),xy(:,2),xy(:,2)*0+10,'o','color',[0 0 0]);

    title([ upper(iivarname)],'fontsize',tt,'color',ck)
    if printok==1
        eval(['print -dpng -r300 ' path_png '/f_' iivarname '.png'])
    end
    if grdok==1
        outfile=[ path_grd  iivarname  runid '_' runname '.grd'];

        all1=all;
        all1.nodes(:,4)=dm;
        f_write_bathy_fort14(all1,outfile,0)
    end
end
diary off

