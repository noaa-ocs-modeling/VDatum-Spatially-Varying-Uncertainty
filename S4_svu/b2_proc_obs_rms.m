%Find tide station without rms
%Assigned the mean of the nearest three tide stations' rms
uncer0=uncer;

loc=find(isnan(uncer) | uncer<-9);
%xy0=reshape([pmoe_datum.xy],2,[])';
xy0=[pmoe_datum.ox; pmoe_datum.oy]';
xy_good=xy0;
xy_good(loc,:)=[];
uncer_good=uncer;
uncer_good(loc,:)=[];

for i=1:length(loc)
    iloc=loc(i);
    ixy=[pmoe_datum(iloc).ox pmoe_datum(iloc).oy ];
    d=sqrt((xy_good(:,1)-ixy(1)).^2+(xy_good(:,2)-ixy(2)).^2);
    [ds, I]=sort(d);
    irms=uncer_good(I(1:3));
    fprintf('%2d: rms set to %.3f  %.3f %.3f --max %.3f\n',i,irms,max(irms));
    uncer(iloc)=max(irms);
end
    
    

