function [dd,G]=stat_interp_uncertainty_jackknife(bcpsi,estn,uncer,W01,PH)

m_error=std(estn);
HPH=PH(bcpsi,:);   % PH is covar=m_error^2*coef.*datt
[n, m]=size(PH);

%diaP=m_error*ones(n,1);
R=eye(m);
for i=1:m
    R(i,i)=W01(i,i)*uncer(i)*uncer(i)*R(i,i);
end
RHPH=R+HPH;

    RHPH1=ones(m+1,m+1);
    RHPH1(1:m,1:m)=RHPH;
    RHPH1(m+1,m+1)=0.0;
    invRHPH1=inv(RHPH1);  

    %invRHPH=inv(RHPH);        
    PH1=ones(n,m+1);
    PH1(1:n,1:m)=PH;
    PH1invRHPH1=PH1*invRHPH1;
    G=PH1invRHPH1(1:n,1:m);
    dd=G*estn; % correction to the datum
end
