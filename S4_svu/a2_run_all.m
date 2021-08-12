clear
%varlist={'mhhw';'mhw';'mlw';'mllw';'mtl' ;'dtl'};
runname='R58_k6s4_msl_5o2_a53_merged';
runid='181';

for i=1:6 % compute svu for six datums in varlist
    tidal_datums_svu(i,runid,runname);
end