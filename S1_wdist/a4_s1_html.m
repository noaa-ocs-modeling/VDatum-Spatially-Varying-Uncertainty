%clear
path_png=[pathout '/'];

outfile=[path_png 'all_plots.html'];
fid=fopen(outfile,'wt');
fprintf(fid,'<html> <body bgcolor="black"><table>\n');


path_in=path_png;
temp=dir([path_in '/f*.png']);
runlist={temp.name}';
a=imread([path_in runlist{1}]);

[fa,fb,fc]=size(a);

ps=round([fb fa]/fb*300);
n=size(runlist,1);
nb=6;
na=ceil(n/nb);
k=0;
pw=ps(1);
ph=ps(2);
for i=1:na
    fprintf(fid,'<tr>\n');

    for j=1:nb
        k=k+1;
        if k>n
            break
        end
        fprintf(fid,'<td>',pw);
        kname=runlist{k};
        %kname=['fig' int2str(k) '.png'];
        temp=['<p><a href="' kname '" target="_blank"><img src="' kname '"  width="' int2str(ps(1)) '" height="' int2str(ps(2)) '" border="0"></a></p>'];
        temp2=int2str(k);
        %fprintf(fid,'%s<p align="center">%s:  %s</p></td>\n',temp,temp2,kname(6:end-4));
        %fprintf(fid,'%s<p align="center">%s:  %s</p></td>\n',temp,temp2,kname(6:end-4));
        fprintf(fid,'%s</td>\n',temp);
        %fprintf(fid,'%s%s</td>\n',temp,temp2);
    end
        
    fprintf(fid,'</tr>\n');  
    
end
  fprintf(fid,'</table></body></them>\n') ;