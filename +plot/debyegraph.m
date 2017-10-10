

nn = 300; n = logspace(-3,7,nn);
nT = 100; T = logspace(2,9,nT); % K

ld =  @(T,n) 6.9e1*sqrt(T*1e-6/n);
debye = zeros(nn,nT);

for ii = 1:nn
    for jj = 1:nT
       debye(ii,jj) = ld(T(jj),n(ii));
    end
end
%%
units = irf_units;
TeV = T/11065;
surf(TeV,n,log(debye))
view([0 0 1])
set(gca,'xlim',TeV([1 end]),'ylim',n([1 end]),'xscale','log','yscale','log','xtick',10.^unique(round(log(TeV))))
xlabel('T [eV]')
ylabel('n [cc]')
title('Debye length')
ch = colorbar;
ylabel(ch,'log(\lambda_{De})')
shading flat