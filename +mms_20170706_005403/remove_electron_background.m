mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
pathLocalUser = ['/Users/' localuser '/'];

tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
tint_fred = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:54:20.00Z');
tint_lobe = irf.tint('2017-07-06T00:54:07.00Z/2017-07-06T00:54:08.00Z');
tint_sheet = irf.tint('2017-07-06T00:54:18.70Z/2017-07-06T00:54:19.50Z');
tint_sheet_zoom = irf.tint('2017-07-06T00:54:19.00Z/2017-07-06T00:54:19.20Z');
tint_sep = irf.tint('2017-07-06T00:54:13.50Z/2017-07-06T00:54:16.50Z');
tint_phi = irf.tint('2017-07-06T00:54:14.00Z/2017-07-06T00:54:15.50Z');

ic = 1;

%% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)

% Original distribution
c_eval('eDist = ePDist?.tlim(tint_sheet_zoom);',ic)
eDist_orig = eDist;

%% Remove background
nSecondary = [0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7];
nPhoto = 1;
[eDist_nobg] = mms.remove_edist_background(eDist_orig);
c_eval('[eDist_nobg?] = mms.remove_edist_background(eDist_orig,''nSecondary'',nSecondary(?),''Nphotoe_art'',nPhoto,''ZeroNaN'',0);',1:numel(nSecondary))


%% Make reduced distribution
lowerelim = 000;
vint = [-Inf Inf];
c_eval('scpot = scPot?.resample(eDist);',ic)
c_eval('dir_red = dmpaB?.resample(eDist).norm;',ic)
dir_red = [1 0 0];
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

%% Plot 1D 
figure(89)
hca = subplot(1,1,1);
if 0
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
else
  plot(hca,fsheet_orig.depend{1}(1,:)*1e-3,mean(fsheet_orig.data,1),...
           fsheet_nobg.depend{1}(1,:)*1e-3,mean(fsheet_nobg.data,1))
  legends = {'original','default model'};
  for inSec = 1:numel(nSecondary) % varying nSecondary
    hold(hca,'on')
    c_eval('plot(hca,fsheet_nobg?.depend{1}(1,:)*1e-3,mean(fsheet_nobg?.data,1));',inSec)         
    hold(hca,'off')
    legends{end+1} = sprintf('nSecondary = %5.2f',nSecondary(inSec));
  end
  legend(hca,legends,'location','eastoutside')
end
hca.XLabel.String = 'v_x (10^3 km/s)';
hca.YLabel.String = 'f_e (s/m^4)';
if 1
  hca.YScale = 'log';
  hca.YLim = [1e-4 2e-2];
end
hca.XGrid = 'on';
hca.YGrid = 'on';
%hca.XLim = [-1 1]*40e0;

%% Use plasma sheet as reference, because the peak at v = 0 stands out really clearly there
c_eval('scpot = scPot?.resample(eDist);',ic)
nMC = 200;
vgmax = 40000;
vg = -vgmax:1000:vgmax;
vg(abs(vg)>70000) = [];
x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

ef2Dxy_orig = eDist_orig.reduce('2D',x,y,'scpot',scpot,'nMC',nMC,'vg',vg,'base','cart');
ef2Dxy_nobg = eDist_nobg.reduce('2D',x,y,'scpot',scpot,'nMC',nMC,'vg',vg,'base','cart');
c_eval('ef2Dxy_nobg? = eDist_nobg?.reduce(''2D'',x,y,''scpot'',scpot,''nMC'',nMC,''vg'',vg,''base'',''cart'');',1:numel(nSecondary))

%%
figure(88)
it = 1:7;

%h = setup_subplots(1+numel(nSecondary)/2,2);
h = setup_subplots(4,4);

isub = 1;
if 1 % original
  hca = h(isub); isub = isub + 1;
  ef2Dxy_orig(it).plot_plane(hca);
  hca.Title.String = 'Original';
end
if 1 % default from model parameters
  hca = h(isub); isub = isub + 1;
  ef2Dxy_nobg(it).plot_plane(hca);
  hca.Title.String = 'Default model, without specifying inputs';
end
iSec = 0;
for inSec = 1:numel(nSecondary) % varying nSecondary
  iSec = iSec + 1;
  hca = h(isub); isub = isub + 1;
  c_eval('ef2Dxy_nobg?(it).plot_plane(hca);',inSec)
  hca.Title.String = sprintf('nSecondary = %g',nSecondary(iSec));
end
hlink = linkprop(h,{'XLim','YLim','CLim'});
hlink.Targets(1).XLim = 30*[-1 1];
hlink.Targets(1).YLim = 30*[-1 1];
hlink.Targets(1).CLim = [-11 -8];