% This program is to combine multiple files produced 
% by the parallel lv8j into one.
% 
% Input files: runid/*.*
% Each process produces 3 files containing 
% highs, lows and model datums respectively.
% The NEGOM case used 256 processors.

% Output: 3 combined matlab data files
% (1) runid_testhh.mat: highs and time for all grid nodes
% (2) runid_testll.mat: lows and time for all grid nodes
% (3) runid_testxx.mat: model datums for all grid nodes

% (1) and (2) will be used in S3 to compute correlation coefficient matrix
% (3) will be used in S4 for SVU computation

%                           Liujuan.Tang@noaa.gov
%                           Last modified 01/27/2020
%---------------Input --------------
% Directory of the inout files
clear
runid='R58_k6s4_msl_5o2_a53_merged'; %ADCIRC run ID
pathin=[ runid '/'];
%-----------------------------------
textlist={'hh';'ll';'xx'}; %hh: highs; ll: lows; xx: datums
for j=1:3
    itext=textlist{j}
    temp=dir([pathin '*' itext '*']);
    filelist={temp.name};  % input file
    n=length(filelist);
    fprintf(1,'%d files found\n',n)
    dmall=[];
    for i=1:n
        fprintf(1,'%s of %d\n',filelist{i},n);
        ifile=[pathin filelist{i}];
        temp=load(ifile);
        dmall=[dmall;temp];
    end
    eval(['test' itext '=dmall;']) 
    eval(['save ' runid '_test' itext '.mat test'  itext]) 
    eval(['test' itext '=[];']) % clear 
    fprintf(1,'Done \n--------------\n')
end