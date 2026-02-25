if doLoad; load /home/cecilia/Research/BeamSolver/2015-05-12T151301_dg_solver.mat; end
doposwi = 0;
[rerunR,rerunS] = find(kmax*Ld>0.99);
nk = 400;
kvec = linspace(1e-5,2,nk)/Ld;
dk = kvec(2)-kvec(1);
nrerun = numel(rerunR);    

for kk = 1:nrerun;
    iR = rerunR(kk);
    iS = rerunS(kk);
    disp(['R=' num2str(R(iR)) ' S=' num2str(S(iS))])
    [wrtemp,witemp,kmaxtemp,wrmaxtemp,wimaxtemp,vphmaxtemp,residualtemp] = dg_solver(vde1,vde2(iS),veth1,veth2,vith,wpe,wpe1(iR),wpe2(iR),wpi,kvec,doposwi);
    wr{iR,iS} = wrtemp;
    wi{iR,iS} = witemp;
    residual{iR,iS} = residualtemp;
    wimax(iR,iS) = wimaxtemp;
    wrmax(iR,iS) = wrmaxtemp;
    kmax(iR,iS) = kmaxtemp;
    vphmax(iR,iS) = vphmaxtemp;      
end
tosave.wr = wr;
tosave.wi = wr;
tosave.residual = residual;
tosave.wimax = wimax;
tosave.wrmax = wrmax;
tosave.kmax = kmax;
tosave.vphmax = vphmax;   
eval(['tosave' num2str(doposwi) '=tosave'])

savePath = '/home/cecilia/Research/BeamSolver/';
saveName = [datestr(now,'yyyy-mm-ddTHHMMSS') '_dg_solver'];
save([savePath saveName])
disp(['saved ' savePath saveName])