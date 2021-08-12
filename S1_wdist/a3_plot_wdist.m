clear, clf
load station_node_info.mat
load cm_dist2 % colormap
pathin='output_dist_time/';
pathout='outpng_dist/';
if ~exist(pathout,'dir')
    eval(['mkdir ' pathout])
end

reloadok=1;
if reloadok==1
    iid=1;
sid=pmoe_datum(iid).id;;
infile=[ 'OWdist' int2str(sid) '_2'];
%infile=['OWtime' int2str(sid) '_2'];
p=all.nodes;  
p(:,1)=[];
p(:,end)=[];
%-------------
q=all.elements;  
q(:,1:2)=[];
end
list=dir([pathin 'OW*_2.mat']);
list={list.name};
n=size(list,2);

r=[96.3 112];
x_lim=[min(p(:,1)) max(p(:,1))];
y_lim=[min(p(:,2)) max(p(:,2))];

fsize=round([diff(x_lim) diff(y_lim)*r(2)/r(1)]/diff(x_lim)*1200);
set(gcf,'Position',[10   10   fsize],'visible','on','color','k','paperpositionmode','auto','InvertHardcopy','off')
ck=[ 0 0 0]+.8;
for i=1:n
    fprintf(1,'%d\n',i);
    infile=list{i};
    eval(['load ' pathin infile]);
    clf,  axes('Position',[0.03 0.06 .95 .92]),  colormap(cm)
    trisurf(q,p(:,1),p(:,2),dis);   shading interp,  view([0 90])
    hb=colorbar;                    set(hb,'color',ck,'fontsize',14)
    caxis([0 700])
    box on
    xlim(x_lim);ylim(y_lim);
    hold on
    h=plot3(pmoe_datum(i).ox,pmoe_datum(i).oy,500,'r+','markersize',12);
    h2=legend(h,int2str(sinfo(i).id));
    set(h2,'fontsize',36);
    set(gca,'color','k','Xcolor',ck,'ycolor',ck)
    box on
    outname=['f' infile(1:end-3) '.png'];
    eval(['print -dpng -r150 '  pathout '/' outname])
end


