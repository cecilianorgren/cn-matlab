tint = irf.tint('2017-07-09T17:31:58.047139000Z/2017-07-09T17:34:23.048132000Z');

c_eval('scPot = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('iPDist = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic)
iPDist_counts = iPDist; iPDist_counts.data = (iPDist.data./iPDistErr.data).^2;

%% Prep data
PD_use = iPDist; % iPDist

nMovMean = 7; % Number of distributions for the moving average. Use the same for finding the energy limt and the gmm

mask_nMovMean = [7 7];
mask_Counts = 5;
run_metadata.PD.mask_nMovMean = mask_nMovMean;
run_metadata.PD.mask_nCounts = mask_Counts;

mat = iPDist_counts.find_low_counts('counts',mask_Counts,'nMovMean',mask_nMovMean,'output','mat');
PD = PD_use.mask('energy','mat',mat);  
PD = PD.movmean(nMovMean);
PD = PD(1:nMovMean:PD.length);
PDc = iPDist_counts.movmean(nMovMean);
PDc = PDc(1:nMovMean:iPDist_counts.length);

nMP = 100000; % Number of macroparticles  
MP = PD.macroparticles('ntot',nMP,'skipzero',1,'scpot',scPot.resample(PD),'v_init','v2','MParts');


%%
K = 1:4;
options = statset('MaxIter',150,'TolFun',1e-5);
GM = GMDist.fit(MP,K,'start','previous','Options',options,'SharedCovariance',false,'RegularizationValue',0.1);
      
%GM.bic
%gm{it,iK} = fitgmdist(R,K,'Start',S,'SharedCovariance',false,'Options',options,'RegularizationValue',0.1);
      

%% Diagnostics, L2 error
K = 4;
[tsL2errPartitions, tsL2NGroups] = GM.tsL2errPartitions(K,'sort','descend'); 
tsL2errBest = GM.tsL2errBest(K);

it = 40; % Example merge
M = 2; % plot maximum number of
bestPInM = tsL2errBest.userData.listPartitions{it,M};
%exGMall = GM.merge(K,{1:K});
exGMall = GM;
exGMmerg = GM.merge(K,bestPInM);

fall = exGMall.tsReduced(K,3); fall = fall(it);
for iK = 1:K
  fcomp{iK} = exGMall.tsReduced(K,3,'component',iK); fcomp{iK} = fcomp{iK}(it);
end
fmerg = exGMmerg.tsReduced(M,3); fmerg = fmerg(it);



%tsL2errPartitions.data(tsL2NGroups.data~=3) = NaN;
h = irf_plot(4);
[h,h2] = initialize_combined_plot('leftright',4,1,1,0.5,'vertical');

hca = irf_panel('omni deflux');
irf_spectrogram(hca,PD.omni.deflux.specrec,'log')

hca = irf_panel('L2err for all P');
irf_spectrogram(hca,tsL2errPartitions.specrec,'lin')

hca = irf_panel('N groups in P');
irf_spectrogram(hca,tsL2NGroups.specrec,'lin')

hca = irf_panel('Best err for nP = 1:K');
irf_spectrogram(hca,tsL2errBest.specrec,'lin')

hca = irf_panel('all L2err');
 
c_eval('h(?).LineWidth = 1;',1:numel(h))

hmark = irf_pl_mark(h,tsL2errBest.time(it),'r');

isub = 1;
hca = h2(isub); isub = isub + 1;
%hca.ColorOrder = [0 0 0; 0.3 0.3 0.3; mms_colors('matlab')];
plot(hca,fall.depend{1},fall.data,fmerg.depend{1},fmerg.data)
hold(hca,'on')
for iK = 1:K
  plot(hca,fcomp{iK}.depend{1},fcomp{iK}.data)
end
hold(hca,'off')
hca.ColorOrder = [0 0 0; 0.7 0.7 0.7; mms_colors('matlab')];
hca.XLim = [-2500 2500];
%hca.YScale = 'log'; %hca.YLim(1) = -1;
strComp = cell(1,K);
for iComp = 1:K, strComp{iComp} = num2str(iComp); end
strAll = strrep(num2str(1:K),' ','');
strMerg = strrep(strjoin(cellfun(@(x) num2str(x),bestPInM,'UniformOutput',false),', '),' ','');
irf_legend(hca,{strAll,strMerg,strComp{:}}',[0.98 0.98])
irf_legend(hca,{sprintf('L2err*(M=%g) = %g',M,tsL2errBest.data(it,M))}',[0.02 0.98])



