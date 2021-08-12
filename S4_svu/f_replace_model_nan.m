function newz=f_replace_model_nan(xy,z,loc)
xyz=[xy(:,1:2) z];
newz=z;
xyz_good=xyz;
xyz_good(loc,:)=[];
xyznan=xyz(loc,:);

z_good=z;
z_good(loc)=[];
znan=griddata(xyz_good(:,1),xyz_good(:,2),z_good,xyznan(:,1),xyznan(:,2),'nearest');
newz(loc)=znan;
