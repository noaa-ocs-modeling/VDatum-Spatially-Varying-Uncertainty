function all=f_read_bathy_adcirc_grd(path_in,infile,geneok)
% this file is to read adcirc fort.14 file
%clear

%path_in='grids\';
%infile='coastlineMesh3h.grd';
%infile='e4s4_30m60m.grd'
fprintf(1,'\n------------------------------------\n');
fprintf(1,'Input file: %s\n',infile);

fid=fopen([path_in infile],'rt');
all.line1=fgetl(fid);
[ns]=fscanf(fid,'%i %i',[1 2]); % nnode nelement
all.line2=ns;
all.nodes=fscanf(fid,'%f',[4,ns(2)])';
all.elements=fscanf(fid,'%f',[5,ns(1)])';

fprintf(1,'Total element number = %d\n',ns(1));
fprintf(1,'Total node number = %d \n',ns(2));


openb.ntype=fscanf(fid,'%i');
openb.ntype_note=fgetl(fid);

fprintf(1,'------------------------------------\n');
fprintf(1,'openb.ntype: %d %s\n',openb.ntype,openb.ntype_note);
    openb.nnode=fscanf(fid,'%i');
    openb.nnode_note=fgetl(fid);
fprintf(1,'openb.nnode: %d %s\n',openb.nnode,openb.nnode_note);
if openb.ntype>0
for i=1:openb.ntype
    open1(i).n=fscanf(fid,'%i');
    open1(i).n_note=fgetl(fid);
    open1(i).node=fscanf(fid,'%f',[1,open1(i).n])';
    fprintf(1,'open1(%d).n: %d %s\n',i,open1(i).n,open1(i).n_note);
end
else
    i=1;
    open1(i).n=[];
    open1(i).node=[];
    open1(i).n_note=[];
end   
fprintf(1,'\n------------------------------------\n');

%-----------------------------
landb.ntype=fscanf(fid,'%i');
landb.ntype_note=fgetl(fid);
landb.nnode=fscanf(fid,'%i');
landb.nnode_note=fgetl(fid);
fprintf(1,'%d %s\n',landb.ntype,landb.ntype_note);
fprintf(1,'%d %s\n',landb.nnode,landb.nnode_note);
fprintf(1,'\n------------------------------------\n');

for i=1:landb.ntype
    land(i).nname=fscanf(fid,'%i %i',[1 2]);
    land(i).nname_note=fgetl(fid);
    land(i).node=fscanf(fid,'%i',land(i).nname(1));
    fprintf(1,'%d %d %s\n',land(i).nname,land(i).nname_note);
end
if geneok==1
geneb.ntype=fscanf(fid,'%i');
geneb.ntype_note=fgetl(fid);
geneb.nnode=fscanf(fid,'%i');
geneb.nnode_note=fgetl(fid);

for i=1:geneb.ntype
    gene(i).nname=fscanf(fid,'%i %i',[1 2]);
    gene(i).nname_note=fgetl(fid);
    gene(i).node=fscanf(fid,'%i',land(i).nname(1));
end
all.geneb=geneb;
all.gene=gene;

end
fclose(fid)
all.openb=openb;
all.open1=open1;
all.landb=landb;
all.land=land;
