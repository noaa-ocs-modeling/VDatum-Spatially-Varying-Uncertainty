function f_write_bathy_fort14(all,outfile,uok)
%input all structure
%      line1: [1x81 char]
%       line2: [367019 192889]
%       nodes: [192889x4 double]
%    elements: [367019x5 double]
%       openb: [1x1 struct]
%       open1: [1x1 struct]
%       landb: [1x1 struct]
%        land: [1x38 struct]
%load allfort14

%outfile='newfort.14';
fid2=fopen(outfile,'wt');
if isempty(all.line1)
    all.line1=outfile;
end
fprintf(fid2,'%70s\n',deblank(all.line1));
fprintf(1,'%s\n',all.line1);



fprintf(fid2,'%9d  %9d \n',all.line2);
fprintf(1,'%d %d\n',all.line2);

fprintf(fid2,'%9d  %18.10f  %18.10f %12.6f  \n',all.nodes');
fprintf(1,'max depth=%.2f; min depth= %.2f m\n',max(all.nodes(:,end)),min(all.nodes(:,end)));

fprintf(fid2,'%d %d %d %d %d \n',all.elements');

% %-------
% fprintf(fid2,'%d %s\n',all.openb.ntype,all.openb.ntype_note);
% fprintf(fid2,'%d %s\n',all.openb.nnode,all.openb.nnode_note);
% if all.openb.nnode>0
% 
% fprintf(1,'%d %s\n',all.openb.ntype,all.openb.ntype_note);
% fprintf(1,'%d %s\n',all.openb.nnode,all.openb.nnode_note);
% 
% 
% fprintf(fid2,'%d %s\n',all.open1.n,all.open1.n_note);
% fprintf(fid2,'%d\n',all.open1.node);
% end
% fprintf(1,'%d %s\n',all.open1.n,all.open1.n_note);
% 
% %wrote_bnodes=all.open1.node;
% % %------------------
% % %remove duplicated nodes
% % landnode_count=0;
% % for i=1:all.landb.ntype
% %     inode=all.land(i).node;
% %     lia = ismember(inode,wrote_bnodes,'rows');
% %     na=find(lia);
% %     if ~isempty(na)
% %         fprintf(1,'\n------------\n%d %d %s\n',all.land(i).nname,all.land(i).nname_note);
% % 
% %         fprintf('%d duplicated boundary nodes found\n ',length(na));
% %         if na(1)==1 | na(end)==length(inode) 
% %             for j=1:length(na)
% %                 fprintf('Removed %d node: %d\n',na(j),inode(na(j)));
% %             end
% %             inode(na)=[];
% %             all.land(i).node=inode;
% %             all.land(i).nname(1)=length(inode);
% %         end
% %     end
% %     landnode_count=landnode_count+all.land(i).nname(1);
% %     wrote_bnodes=[wrote_bnodes;inode];
% % end
% 
% %-----------
% fprintf(fid2,'%d %s\n',all.landb.ntype,all.landb.ntype_note);
% fprintf(fid2,'%d %s\n',all.landb.nnode,all.landb.nnode_note);
% 
% fprintf(1,'%d %s\n',all.landb.ntype,all.landb.ntype_note);
% fprintf(1,'%d %s\n',all.landb.nnode,all.landb.nnode_note);
% 
% 
% for i=1:all.landb.ntype
% fprintf(fid2,'%d %d %s\n',all.land(i).nname,all.land(i).nname_note);
% %fprintf(fid2,'%d %d\n',all.land(i).nname(1),all.land(i).nname(2));
% fprintf(1,'%d %d %s\n',all.land(i).nname,all.land(i).nname_note);
% fprintf(fid2,'%d\n',all.land(i).node);
% 
% end

fprintf(1,'Wrote to %s ... Done\n',outfile);


fclose(fid2);