%------------------------------------------------------------------------
% this code calculates and plots as a function of L (e-fold length of
% exponential distance mask) one of 4 parameters:
%   r.m.s. deviation of the product from observed datums at the gages;
%   r.m.s. deviation of the product from the model domain-wide;
%   r.m.s. of the magnitude of the gradient of the product’s deviation 
%                                                         from the model; 
%   r.m.s. of the product uncertainty.
%
% Input: outputs of a2_tidal_datum.m for a set of 'covar' matrices used for
% data assimilation, the set is defined by an array of masks
%                   Elena.Tolkova@noaa.gov
%                          06/2023
%------------------------------------------------------------------------
clear variables
%-----------------------------------------------------------------------
path_pre='pre_process/'; % directory with Gs and products
datumlist={'mhhw';'mhw';'mlw';'mllw';'mtl' ;'dtl'};
varlist={'d';'dd';'diaPA'}; % see vartest 
vartext={'Product' ;'Correction';'Uncertainty'};
masks=[0 200 0; 0 110 0; 0 1 0; 0 2 0; 0 4 0; 0 6 0]
[im i]=size(masks);

%---------- difference between neighboring nodes -------
adcirc_dir='C:\Users\elena.tolkova\Documents\TXmodel\jet_runTXLAMS_v1\';
eval(['load ' adcirc_dir 'allfort14.mat'])
    ne=allgrd.line2(1); nn=allgrd.line2(2);
    x=allgrd.nodes(:,2);   y=allgrd.nodes(:,3); z=allgrd.nodes(:,4);
    meshele=allgrd.elements(:,[3 4 5]);
fprintf(1,'Total elements = %d\n',ne);fprintf(1,'Total nodes =%d\n',nn);
%----------------- remove elements containing dry and shallow nodes ------------------
dry=0.1;
mskdry=(z<=dry);
drynodes=int32(1:nn);
drynodes(~mskdry)=[];
[i1,i2]=find(ismember(meshele,drynodes));
i1=unique(i1);
meshele(i1,:)=[];
%------- 1.0 Compute unique segements and corresponding distance-----------
A=[meshele(:,1:2);meshele(:,2:3);meshele(:,[3 1])];
[Au0] = unique(sort(A,2), 'rows','stable');
nAu=length(Au0);        fprintf(1,'Total segements = %d \n',nAu);
dA0 = pos2dist_1(y(Au0(:,1)),x(Au0(:,1)),y(Au0(:,2)),x(Au0(:,2)),2);%km

eval(['load ' adcirc_dir 'pmoe_datum'])
bcpsi=[pmoe_datum.node]';  % grid node # at station location

figure
msz=11;
lwd=2.5;
hold off
lgd={};

for ii=1:4 % datums 1-4
    iivarname=datumlist{ii};
    fprintf('----------------%s-----------------\n',iivarname)
    lgd{ii}=iivarname;    
 switch iivarname
    case 'mhhw'
         ostn=[pmoe_datum.omhhw]';
    case 'mhw'
        ostn=[pmoe_datum.omhw]';
    case 'mlw' 
        ostn=[pmoe_datum.omlw]';
    case 'mllw'
        ostn=[pmoe_datum.omllw]';
 end

hh(1:im)=NaN;
ddstats=zeros(im,5);
Lr=zeros(im,1);
edg=0:0.01:2.0;
 for imsk=1:im
    imask_cvr=masks(imsk,1);
    imask_dst=masks(imsk,2);
    imask_dm=masks(imsk,3);
infile=[path_pre  'd_dd_diaPA_' iivarname '_masks' int2str(imask_cvr)  int2str(imask_dst)  int2str(imask_dm) '.mat'];
    eval(['load(''' infile ''');']); 
    % dd is product-model diff
    ggerr=d(bcpsi)-ostn; % product-observations diff
    ggerr(isnan(ostn))=[];
    m=length(ggerr);
    n=length(dm);

    if imask_dst<100
        Lr(imsk) = 55.5*imask_dst;
    else
        Lr(imsk)=55.5/(imask_dst-100);
    end
% --statistics--
disp([iivarname, ' imask=',int2str(imask_dst)])
[100*sqrt(dd'*dd/n),100*sum(abs(dd))/n,100*sqrt(diaPA'*diaPA/n)] %,100*sqrt(ggerr'*ggerr/m),100*max(abs(ggerr))]    

incrdd=zeros(nAu,1);
 for i=1:nAu
    i1=Au0(i,1);
    i2=Au0(i,2);
    r=dA0(i);
    incrdd(i)=1000*(dd(i1)-dd(i2))/r;% units: mm/km
 end
ddstats(imsk,1:5)=[sqrt(2*incrdd'*incrdd/nAu),100*sqrt(ggerr'*ggerr/m),100*sqrt(dd'*dd/n),100*sum(abs(dd))/n,100*sqrt(diaPA'*diaPA/n)];   

 end  % for imsk=1:im
% plot(Lr,ddstats(:,1),'o:','markersize',msz,'linewidth',lwd)
% ylabel('gradient r.m.s., mm/km','fontsize',16)
% picname='gradient_Ls';
% hold on

% plot(Lr,ddstats(:,2),'+:','markersize',msz,'linewidth',lwd)
% ylabel('product-obs, r.m.s., cm','fontsize',14)
% picname='product-obs_Ls';
% hold on

% plot(Lr,ddstats(:,3),'s:','markersize',msz,'linewidth',lwd)
% ylabel('model change, r.m.s., cm','fontsize',14)
% picname='model_change_Ls';
% hold on

% % plot(Lr,ddstats(:,4),'s:','markersize',msz,'linewidth',lwd)
% % ylabel('model change, mean abs, cm','fontsize',14)
% % picname='model_change_abs_Ls';
% % hold on

plot(Lr,ddstats(:,5),'^:','markersize',msz,'linewidth',lwd)
ylabel('uncertainty, r.m.s., cm','fontsize',14)
picname='uncer_Ls';
hold on

end
lg=legend(lgd,'location','best')
xlabel('L, km','fontsize',16)
set(gca,'fontsize',14)
grid on

