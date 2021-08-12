function tidal_datum_svu(ii,runid,runname);
% This code was developed by Lei Shi
% Liujuan Tang revised the structure. 
% Last revised on Aug. 16 2017
%      ikout1d=['d.txt'],         % final datum filed
%      ikoutdd=['dd.txt'],        % correction field
%      ikoutPA=['diaPA.txt'],     % uncertainty field
varlist={'mhhw';'mhw';'mlw';'mllw';'mtl' ;'dtl'};

path_pre='pre_process\';
path_out='out_mat\';
if ~exist(path_out,'dir')
    eval(['mkdir ' path_out])
end

eval(['load ' path_pre 'SVU_input_' varlist{ii} '_' runid '_' runname '.mat']);
W01=zeros(length(uncer),length(uncer));
for i=1:length(uncer);W01(i,i)=1;end;% 1/4/2016.
mc=0;
while 1
    mc=mc+1;
    [d, dd, diaPA, d_jn, dd_jn, diaPA_jn]=stat_interp_uncertainty_jackknife(bcpsi,estn,uncer,dm,W01,covar,mc);
    change=zeros(length(uncer),1);
    for i=1:length(uncer);
        if abs(d_jn(i)-ostn(i)) > min(0.01, uncer(i)); % we need to estabilish the one-to-one cooresponding of dd_jack and uncer!!!
        %if abs(d_jn(i)-ostn(i)) >  uncer(i); % we need to estabilish the one-to-one cooresponding of dd_jack and uncer!!!
                %fprintf(1,'%d: W01(i,i) old= %.3e; new =%.3e\n',i,W01(i,i),W01(i,i)/sqrt(2.0));
               W01(i,i)=W01(i,i)/sqrt(2.0);% W01(i)=W01(i)/2.0;
               %if W01(i,i)>1e-5
               change(i)=1;
               %end
        end;
    end;
    
    fprintf(1,'mc=%d; sum of change is %.f \n',mc,sum(change)); 
%     if rem(mc,1)==2
%         ioutfile=['output/d_dd_diaPA_' int2str(mc)];
%         eval(['save ' ioutfile ' d dd diaPA change d_jn dd_jn diaPA_jn'])
%     end
    if sum(change) == 0 ,   break;  end;
end
outfile=['d_dd_diaPA_' varname '_' runid '_' runname] 
eval(['save ' path_out outfile '   d dd diaPA dm change d_jn dd_jn diaPA_jn'])
    
function [d, dd, diaPA, d_jn, dd_jn, diaPA_jn]=stat_interp_uncertainty_jackknife(bcpsi,estn,uncer,dm,W01,covar,mc);

m_error=std(estn);% % % covar=coef.*datt;
PH=covar;% here we assume a constant model error distribution
HPH=PH(bcpsi,:);
[n, m]=size(PH);
nc=n;% warning nc is the same as n, nsta==m
go_ordinary=1;
if ~go_ordinary % assumption, node length is nc, 
    diaP=m_error*ones(nc,1);
    R=eye(length(uncer));
    for i=1:length(uncer);R(i,i)=W01(i,i)*uncer(i)*uncer(i)*R(i,i);end;
    RHPH=R+HPH;
    invRHPH=inv(R+HPH);
    PHinvRHPH=PH*invRHPH;
    G=PHinvRHPH;
else; % assumption, node length is nc, 
    diaP=m_error*ones(nc,1);    % assumption, station length is m.
    R=eye(length(uncer));
    for i=1:length(uncer);R(i,i)=W01(i,i)*uncer(i)*uncer(i)*R(i,i);end;
    RHPH=R+HPH;
    RHPH1=ones(m+1,m+1);
    RHPH1(1:m,1:m)=RHPH;
    RHPH1(m+1,m+1)=0.0;
    invRHPH1=inv(RHPH1);  

    invRHPH=inv(R+HPH);        
    PH1=ones(n,m+1);
    PH1(1:n,1:m)=PH;
    PH1invRHPH1=PH1*invRHPH1;
    PH1invRHPH1=PH1invRHPH1(1:n,1:m);
    G=PH1invRHPH1;
end;
diffom=estn;
if ~go_ordinary;
    dd=PH*(RHPH\diffom);% dd=PH*invRHPH*(do-Hdm);
else;%         size(PH1invRHPH1),size(diffom(:,2)),
    dd=PH1invRHPH1*diffom;
    %outfile=['PH1invRHPH1'];
    %eval(['save output/' outfile '_' int2str(mc) ' ' outfile])
end;
d=dm+dd;%d=dm+PH*invRHPH*(do-Hdm);
d_jn=d(bcpsi);
dd_jn=dd(bcpsi);
%----------------------
% here we calculate the uncertainty
diaP=m_error*m_error*ones(nc,1);
diaPA=zeros(nc,1);
% uncertainty estimate
cent01=G;% cent01=PH*invRHPH;
bent01=zeros(nc,1);
R=eye(length(uncer));
for i=1:length(uncer);R(i,i)=uncer(i)*uncer(i)*R(i,i);end;
% diaPA=sqrt(diaP-sum(G.*PH,2));
diaPA=sqrt(diaP-2.0*sum(G.*PH,2)+sum((G*(R+HPH)).*G,2));
diaPA_jn=diaPA(bcpsi);

