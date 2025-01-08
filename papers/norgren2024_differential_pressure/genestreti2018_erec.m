tint = irf.tint('2017-07-11T22:33:00.00Z/2017-07-11T22:34:20.00Z');
tint_e_rec = irf.tint('2017-07-11T22:34:03.00Z/2017-07-11T22:34:04.00Z');

mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db');

ic = 1:4;

c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gseE?fast = mms.get_data(''E_gse_edp_fast_l2'',tint,?)ö',ic);

%%
angle = 1;
ic = 1:4;
id_lmn = 14; % from Genestreti 2018
switch id_lmn
  case 1 % GSW
    L = [0.99986, -0.0521, 0.0052];
    M = [0.0523, 0.9980, 0.0362];
    N = [-0.0019,-0.9361,0.9993];  
    angle = 21;
  case 8 % MVA-Ve
    L = [0.9482,-0.2551,-0.1893];
    M = [0.1749,0.9168,-0.3591];
    N = [0.2651,0.3074,0.9139];  
  case 9 % MVA-E
    L = [];
    M = [];
    N = [];  
  case 14 % #14 from Genestreti 2018, and the one used by Torbert
    L = [0.9482,-0.2551,-0.1893];
    M = [0.1818,0.9245,-0.3350];
    N = [0.2604,0.2832,0.9230];  
    angle = 1;
end
lmn = [L;M;N];

c_eval("mvaE? = gseE?.tlim(tint_e_rec)*lmn';",ic)
c_eval("mvaE?fast = gseE?fast.tlim(tint_e_rec)*lmn';",ic)

%c_eval('E = mvaE?;',ic)
%c_eval('Efast = mvaE?fast;',ic)

%meanEbrst = mean(E.data);
%stdEbrst = std(E.data);
%meanEfast = mean(Efast.data);
%stdEfast = std(Efast.data);

%c_eval('mvaE?_M_nonstar = irf.ts_scalar(mvaE?.time, (mvaE?.y.data - cosd(angle)*mvaE?.y.data)/sind(angle));',ic)
%c_eval('mvaE?fast_M_nonstar = irf.ts_scalar(mvaE?fast.time, (mvaE?fast.y.data - cosd(angle)*mvaE?fast.y.data)/sind(angle));',ic)


%disp(sprintf('BRST: E = [%7.4f,%7.4f,%7.4f] +- [%7.4f,%7.4f,%7.4f], FAST: E = [%7.4f,%7.4f,%7.4f] +- [%7.4f,%7.4f,%7.4f]',meanEbrst(:),stdEbrst(:),meanEfast(:),stdEfast(:)))
disp(sprintf('id_LMN = %g',id_lmn))
c_eval("disp(sprintf('MMS ?: BRST: EM = %7.4f +- %7.4f, FAST: EM = %7.4f +- %7.4f',mean(mvaE?.data(:,2)),std(mvaE?.data(:,2)),mean(mvaE?fast.data(:,2)),std(mvaE?fast.data(:,2))))",ic)

%%
hca = subplot(1,1,1);
set(hca,'ColorOrder',mms_colors('1234'))
plot(hca,mvaE1fast.data(:,3),mvaE1fast.data(:,2),'kx',...
         mvaE2fast.data(:,3),mvaE2fast.data(:,2),'rx',...
         mvaE3fast.data(:,3),mvaE3fast.data(:,2),'gx',...
         mvaE4fast.data(:,3),mvaE4fast.data(:,2),'bx')
hca.XLim = [-20 35];
hca.YLim = [-4 8];

