ic = 1;
tint = irf.tint('2015-11-12T07:18:54.00Z/2015-11-12T07:19:45.00Z');
tint_lim = irf.tint('2015-11-12T07:19:00.00Z/2015-11-12T07:19:10.00Z');

% Load data
if 0
  c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
  c_eval('ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
else
  load /Users/cecilia/Data/ePdist_early
end
ePDist1 = ePDist1.tlim(tint_lim);
dmpaB1 = dmpaB1.tlim(tint_lim);

%% Test operations

%%
ePDist1_e64 = ePDist1.e64;
colors = mms_colors('matlab');
nrows = 1;
ncols = 1;
isub = 1;

if 0
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  plot(hca,1:size(ePDist1_e64.depend{1},2),ePDist1_e64.depend{1})
end
if 0
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  plot(hca,ePDist1_e64.depend{1}(1,:),ePDist1_e64.depend{1}(1,:))
end
hca = subplot(nrows,ncols,isub); isub = isub + 1;
lines = loglog(hca,ePDist1.depend{1}(1,:),ePDist1.depend{1}(1,:),...
                   ePDist1.depend{1}(1,:),ePDist1.depend{1}(1,:) - ePDist1.ancillary.delta_energy_minus(1,:),...
                   ePDist1.depend{1}(1,:),ePDist1.depend{1}(1,:) + ePDist1.ancillary.delta_energy_plus(1,:)...
               );
hold(hca,'on')
             %hca = subplot(nrows,ncols,isub); isub = isub + 1;
c_eval('lines(?).Color = colors(!,:);',1:3,1)
c_eval('lines(?).Marker = ''*'';',1:3,1)
lines(2).LineStyle = '--';
lines(3).LineStyle = '--';

lines = loglog(hca,ePDist1.depend{1}(2,:),ePDist1.depend{1}(2,:),...
                   ePDist1.depend{1}(2,:),ePDist1.depend{1}(2,:) - ePDist1.ancillary.delta_energy_minus(2,:),...
                   ePDist1.depend{1}(2,:),ePDist1.depend{1}(2,:) + ePDist1.ancillary.delta_energy_plus(2,:)...
               );
c_eval('lines(?).Color = colors(!,:);',1:3,2)
c_eval('lines(?).Marker = ''*'';',1:3,2)
lines(2).LineStyle = '--';
lines(3).LineStyle = '--';

if 1
  lines = loglog(hca,ePDist1_e64.depend{1}(1,:),ePDist1_e64.depend{1}(1,:),...
                     ePDist1_e64.depend{1}(1,:),ePDist1_e64.depend{1}(1,:) - ePDist1_e64.ancillary.delta_energy_minus(1,:),...
                     ePDist1_e64.depend{1}(1,:),ePDist1_e64.depend{1}(1,:) + ePDist1_e64.ancillary.delta_energy_plus(1,:)...
                 );
  c_eval('lines(?).Color = colors(!,:);',1:3,3)
  c_eval('lines(?).Marker = ''o'';',1:3,3)
  lines(2).LineStyle = '--';
  lines(3).LineStyle = '--';
end
hold(hca,'off')
         % loglog(hca,ePDist1_e64.depend{1}(1,:),ePDist1_e64.depend{1}(1,:),...
%            ePDist1_e64.depend{1}(1,:),ePDist1_e64.depend{1}(1,:) + ePDist1_e64.ancillary.delta_energy_minus(1,:),...
%            ePDist1_e64.depend{1}(1,:),ePDist1_e64.depend{1}(1,:) + ePDist1_e64.ancillary.delta_energy_plus(1,:)...
%            )

%% Yuri's operations
elim = [460 700];
h = irf_plot(2);

hca = irf_panel('e pitchangles low');
irf_spectrogram(hca,ePDist1.elim(elim).pitchangles(dmpaB1,18).specrec('pa'));

hca = irf_panel('e pitchangles low e64');
irf_spectrogram(hca,ePDist1.e64.elim(elim).pitchangles(dmpaB1,18).specrec('pa'));
