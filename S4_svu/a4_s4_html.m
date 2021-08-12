 clear
path_png='out_fig_181_R58_k6s4_msl_5o2_a53_merged';
outfile=[path_png '/all_plots.html'];
fid=fopen(outfile,'wt');
fprintf(fid,'<html><table>\n');

varlist={'mhhw';'mhw';'mlw';'mllw';'mtl';'dtl'};
nb=4;
na=6;

path_in=path_png;

ps=round([850 500]/850*400);
%n=size(runlist,1);
k=0;
pw=ps(1);
ph=ps(2);
for i=1:na
    fprintf(fid,'<tr>\n');
    temp=dir([path_in '/f*_' varlist{i} '*.png']);
    runlist={temp.name}';

    for j=[4 1 2 3 ]
        k=j;
       
        fprintf(fid,'<td>',pw);
        kname=runlist{k};
        %kname=['fig' int2str(k) '.png'];
        temp=['<p><a href="' kname '" target="_blank"><img src="' kname '"  width="' int2str(ps(1)) '" height="' int2str(ps(2)) '" border="0"></a></p>'];
        temp2=int2str(k);
%        fprintf(fid,'%s<p align="center">%s:  %d</p></td>\n',temp,temp2,moe_datum(k).id);
        %fprintf(fid,'%s%s</td>\n',temp,temp2);
        fprintf(fid,'%s</td>\n',temp);
    end
        
    fprintf(fid,'</tr>\n');  
    
end
  fprintf(fid,'</table></them>\n') ;
  fclose(fid)