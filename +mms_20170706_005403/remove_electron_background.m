mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
pathLocalUser = ['/Users/' localuser '/'];

tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
tint_fred = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:54:20.00Z');
tint_lobe = irf.tint('2017-07-06T00:54:07.00Z/2017-07-06T00:54:08.00Z');
tint_sheet = irf.tint('2017-07-06T00:54:18.70Z/2017-07-06T00:54:19.50Z');
tint_sep = irf.tint('2017-07-06T00:54:13.50Z/2017-07-06T00:54:16.50Z');
tint_phi = irf.tint('2017-07-06T00:54:14.00Z/2017-07-06T00:54:15.50Z');

ic = 1;

%% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)

% Original distribution
c_eval('eDist = ePDist?.tlim(tint_sheet);',ic)
eDist_orig = eDist;

%% Remove background
nSecondary = [0.5 1 1.5 2];
[eDist_nobg] = mms.remove_edist_background(eDist_orig);
c_eval('[eDist_nobg?] = mms.remove_edist_background(eDist_orig,''nSecondary'',nSecondary(?),''ZeroNaN'',NaN);',1:numel(nSecondary))


%% Make reduced distribution
lowerelim = 000;
vint = [-Inf Inf];
c_eval('scpot = scPot?.resample(eDist);',ic)
c_eval('dir_red = dmpaB?.resample(eDist).norm;',ic)
vgmax = 70000;
vg = -vgmax:1000:vgmax;
vg(abs(vg)>70000) = [];
nMC = 500;

tic; ef1D_orig = eDist_orig.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B
tic; ef1D_nobg = eDist_nobg.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B
c_eval('tic; ef1D_nobg? = eDist_nobg?.reduce(''1D'',dir_red,''scpot'',scpot,''lowerelim'',lowerelim,''nMC'',nMC,''vg'',vg); toc',1:numel(nSecondary)) % reduced distribution along B

fsheet_orig = ef1D_orig.tlim(tint_sheet);
fsheet_nobg = ef1D_nobg.tlim(tint_sheet);
c_eval('fsheet_nobg? = ef1D_nobg?.tlim(tint_sheet);',1:numel(nSecondary))

%% Plot
figure;
hca = subplot(1,1,1);
plot(hca,fsheet_orig.depend{1}(1,:)*1e-3,mean(fsheet_orig.data,1),...
         fsheet_nobg.depend{1}(1,:)*1e-3,mean(fsheet_nobg.data,1),...
         fsheet_nobg1.depend{1}(1,:)*1e-3,mean(fsheet_nobg1.data,1),...
         fsheet_nobg2.depend{1}(1,:)*1e-3,mean(fsheet_nobg2.data,1),...
         fsheet_nobg3.depend{1}(1,:)*1e-3,mean(fsheet_nobg3.data,1),...
         fsheet_nobg4.depend{1}(1,:)*1e-3,mean(fsheet_nobg4.data,1))
legend(hca,{'original','background removed',...
            sprintf('background removed, nSecondary = %g',nSecondary(1)),...
            sprintf('background removed, nSecondary = %g',nSecondary(2)),...
            sprintf('background removed, nSecondary = %g',nSecondary(3)),...
            sprintf('background removed, nSecondary = %g',nSecondary(4))})
%hca.XLim = [-1 1]*1e1;
