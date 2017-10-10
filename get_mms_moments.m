cd /Users/Cecilia/Data/MMS/2015Aug15/
tmpDataObj = dataobj('data/mms3_des_brst_l1b_moms_20150815125500_v0.1.1.cdf');
ne = get_variable(tmpDataObj,'mms3_des_numberDensity');
ne = mms.variable2ts(ne);

peXX = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_PresXX')); 
peYY = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_PresYY')); 
peZZ = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_PresZZ')); 
pe = irf.ts_scalar(peXX.time.epochUnix,(peXX.data + peYY.data + peZZ.data)/3);
pe.units = peXX.units;
pe.userData = peXX.userData;


TeXX = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_TempXX'));
TeYY = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_TempYY'));
TeZZ = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_TempZZ'));
Te = irf.ts_scalar(TeXX.time.epochUnix,(TeXX.data + TeYY.data + TeZZ.data)/3);
Te.units = TeXX.units;
Te.userData = TeXX.userData;


veX = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_bulkX'));
veY = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_bulkY'));
veZ = mms.variable2ts(get_variable(tmpDataObj,'mms3_des_bulkZ'));
ve = irf.ts_vec_xyz(veX.time.epochUnix,[veX.data veY.data veZ.data]);
ve.units = veX.units;
ve.userData = veX.userData;

% Downsample electron moments
fs = 1/(ne.time(2)-ne.time(1));
fny = fs/2;
pe_lowres = irf_filt(pe,0,fny/2,fs,5);
Te_lowres = irf_filt(Te,0,fny/2,fs,5);
ne_lowres = irf_filt(ne,0,fny/2,fs,5);


for kk = 1:3
    h(kk) = subplot(3,1,kk);
end
%h = irf_plot(3);
isub = 1;
hca = h(isub); isub = isub + 1;
irf_plot(hca,ne)
hca = h(isub); isub = isub + 1;
irf_plot(hca,Te)
hca = h(isub); isub = isub + 1;
irf_plot(hca,pe)
