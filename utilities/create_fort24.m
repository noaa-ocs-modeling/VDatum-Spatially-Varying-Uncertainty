clearvars; clc;
addpath(genpath('utilities/'))
addpath(genpath('datasets/'))

TS        = '1992-07-01 12:00' ; % start time of simulation
TE        = '1992-08-04 12:00' ; % end time of simulation
% %CONST     = 'major8' ;
CONST     = {'K1','K2','M2','N2','O1','P1','Q1','S2'};
DT = 6;
TPXO9 = 'C:/Users/.../h_tpxo9.v1.nc';
m = msh('C:/Users/fort'); %fort.14 file
mesh1 = Make_f15( m, TS, TE, DT, 'tidal_database',TPXO9, 'const', CONST) ;
test = Make_f24(mesh1,'FES2014/SAL/','FES2014');  %SAL folder contains netcdf files
write(test,'fort','.24');

%The constituents order should be the same as your fort.15, if it's not (if fort.15 created by SMS) you
%can use the following function to reorder the fort.15
reorder_and_fixwidth_f24('fort.24','fort.24_reordered', CONST);

