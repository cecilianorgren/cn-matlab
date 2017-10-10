% check_beam_solver
% in beam_solver_run_spis: run [wr,wi,wrall,wiall,wrin,wiin,param] = beam_solver(Te1,Te2/Te1,Ti/Te1,R,S,k,n,doPlot,TolFun);

wi1 = squeeze(wi(:,:,2,1,1,1)')*sqrt(1836); % in units of ion plasma frequency
wi2 = squeeze(wiall(:,:,2,1,1,1)')*sqrt(1836); % in units of ion plasma frequency
wi3 = squeeze(wiall(:,:,2,1,1,2)')*sqrt(1836); % in units of ion plasma frequency
wr1 = squeeze(wr(:,:,2,1,1,1)')*sqrt(1836); % in units of ion plasma frequency
wr2 = squeeze(wrall(:,:,2,1,1,1)')*sqrt(1836); % in units of ion plasma frequency
wr3 = squeeze(wrall(:,:,2,1,1,2)')*sqrt(1836); % in units of ion plasma frequency



wi2(wi2<-13) = NaN;
for kk = 1:6; h(kk)=subplot(2,3,kk); end
isub=1;
hca = h(isub); isub=isub+1;
surf(hca,k,S,wi1)
hca = h(isub); isub=isub+1;
surf(hca,k,S,wi2)
hca = h(isub); isub=isub+1;
surf(hca,k,S,wi3)
hca = h(isub); isub=isub+1;
surf(hca,k,S,wr1)
hca = h(isub); isub=isub+1;
surf(hca,k,S,wr2)
hca = h(isub); isub=isub+1;
surf(hca,k,S,wr3)

clim = 3*[-1 1];
for kk = 1:6; set(h(kk),'clim',clim); view(h(kk),[0 0 1]); end