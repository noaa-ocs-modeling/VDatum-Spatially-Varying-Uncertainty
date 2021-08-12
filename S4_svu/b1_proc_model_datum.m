%reset dmsl points with dmsl=9.999 and land points
%-------------------------
dm=dm-dmsl;
loc=find(dmall(:,end)~=0 |  isnan(dm) );
eval(['save ' path_pre 'mdatum_nanloc_' varname ' loc'])

fprintf('%d nodes have no datums\n',length(loc));
loc_model_nan=loc;
dm_good=dm;
dm_good(loc)=[];

xyz=all.nodes;
xyz_good=xyz;
xyz_good(loc,:)=[];
xyznan=xyz(loc,:);

dmnan=griddata(xyz_good(:,2),xyz_good(:,3),dm_good,xyznan(:,2),xyznan(:,3),'nearest');
dm(loc)=dmnan;
%dm(loc)=0;
sf=sprintf('max=%.4f m; min=%.4f m\n',max(max(dm(:,:))),min(min(dm(:,:))));    
fprintf(1,'max=%.4f m; min=%.4f m\n',max(max(dm(:,:))),min(min(dm(:,:))))

