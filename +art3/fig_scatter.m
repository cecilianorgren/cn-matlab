cd /Users/Cecilia/Research/EH2/Runs-solver_Buneman_kvR/
doLoad = 1;
doClearup = 0;
doPlot = 1;
doPhi = 0;
if doLoad
    loadPath = '/Users/Cecilia/Research/EH/BeamSolver/';
    matName = '2015-03-21T012941_beam_solver_67kx27Sx20Rx4Te2x3Ti.mat';
    matName = '2015-03-25T141546_beam_solver_67kx27Sx20Rx4Te2x3Ti.mat';
    matName = '2015-03-26T122805_beam_solver_67kx27Sx20Rx4Te2x3Ti.mat';
    matName = '2015-03-27T111315_beam_solver_67kx27Sx20Rx4Te2x3Ti.mat';
    matName = '2015-03-27T192411_beam_solver_67kx27Sx20Rx2Te2x2Ti.mat';
    load([loadPath matName])
end

iTe2 = 2;
vph_plot = vphmax;
vph_plot(vphmax>2) = NaN;
vph_plot(vphmax<1e-3) = NaN;
f_tot = 1;
f_vthebeam = 0;
vbeam = repmat(S',1,numel(R),nTe2,nTi)*param.vte1; % km/s
vphi = f_tot*(vbeam-f_vthebeam*param.vte2(iTe2)-vph_plot*param.vte1);
vph = vph_plot*param.vte1;

% daniels fit
pexp=1;
pcoeff=15;
dppxx = linspace(1e-4,1e2,10);
dppyy = pcoeff*dppxx.^pexp;

ppyy = pcoeff*(dppxx).^(pexp);

iS = 1:2:numel(S);
iR = 1:2:numel(R);
iTe2 = 1:numel(Te2);
iTi = 1:numel(Ti);
iTi=2;
iTe2=1:2;
%%
ppx = reshape(vph(iS,iR,iTe2,iTi),numel(vph(iS,iR,iTe2,iTi)),1)/param.vte1;
ppy = reshape(vph(iS,iR,iTe2,iTi),numel(vph(iS,iR,iTe2,iTi)),1)./reshape(vphi(iS,iR,iTe2,iTi),numel(vphi(iS,iR,iTe2,iTi)),1);
ppx2 = reshape(vph(iS,iR,1,iTi),numel(vph(iS,iR,1,iTi)),1)/param.vte1;
ppy2 = reshape(vph(iS,iR,1,iTi),numel(vph(iS,iR,1,iTi)),1)./reshape(vphi(iS,iR,2,iTi),numel(vphi(iS,iR,1,iTi)),1);
ppx3 = reshape(vph(iS,iR,2,iTi),numel(vph(iS,iR,2,iTi)),1)/param.vte1;
ppy3 = reshape(vph(iS,iR,2,iTi),numel(vph(iS,iR,2,iTi)),1)./reshape(vphi(iS,iR,2,iTi),numel(vphi(iS,iR,2,iTi)),1);

loglog(dppxx,dppyy,'k',ppx,ppy,'ob',ppx2,ppy2,'.g',ppx3,ppy3,'xr')
set(gca,'xlim',[1e-3 1e1],'ylim',[1e-3 1e1])
