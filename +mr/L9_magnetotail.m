units = irf_units;
singleevent = 0;
if singleevent
  BL = 15e-9; % T
  nL = 0.1e6; % m-3
else
  vecBL = linspace(5,20,110)*1e-9;
  vecNL = linspace(0.01,1,100)*1e6;
  vecNL = logspace(-2,0,100)*1e6;
  [BL,nL] = meshgrid(vecBL,vecNL);
end

LLRE = 30;
LL = LLRE*units.RE; % m
LX = 0.1*LL; % m

VA = BL./sqrt(units.mu0.*nL*units.mp);

Rv = 0.1;
ER = Rv.*BL.*VA; % V/m

Phi = LX.*ER;

AR = 0.5*(0.1*LL).^2;
Flux = BL.*AR;
T = Flux./Phi;

if singleevent
  disp(sprintf(' Tail width = %.0f RE \n X line length = %.0f RE \n VA = %.0f km/s \n ER = %.3g mV/m \n X line Potential drop = %.1f kV \n Reconnection duration = %.2f s',LLRE,LX/units.RE,VA*1e-3,ER*1e3,Phi*1e-3,T))
else
  nrows = 2;
  ncols = 2;

  hca = subplot(nrows,ncols,1);
  %pcolor(hca,vecBL,vecNL,VA); shading(hca,'flat'); colorbar
  [c,h] = contour(hca,vecBL*1e9,vecNL*1e-6,VA*1e-3,[300 500:500:4000]);
  clabel(c,h)
  hca.Title.String = 'Alfv?n velocity';
  irf_legend(hca,'km/s',[0.1,0.9],'color','k','fontsize',16)


  hca = subplot(nrows,ncols,2);
  %pcolor(hca,vecBL,vecNL,ER); shading(hca,'flat'); colorbar
  [c,h] = contour(hca,vecBL*1e9,vecNL*1e-6,ER*1e3,[0.1 0.3 0.5 1:10]);
  clabel(c,h)
  hca.Title.String = 'Reconnection electric field';
  irf_legend(hca,'mV/m',[0.1,0.9],'color','k','fontsize',16)

  hca = subplot(nrows,ncols,3);
  %pcolor(hca,vecBL,vecNL,ER); shading(hca,'flat'); colorbar
  [c,h] = contour(hca,vecBL*1e9,vecNL*1e-6,Phi*1e-3,[5 10 20:20:140]);
  clabel(c,h)
  hca.Title.String = 'X line potential drop';
  irf_legend(hca,'kV',[0.1,0.9],'color','k','fontsize',16)

  hca = subplot(nrows,ncols,4);
  %pcolor(hca,vecBL,vecNL,ER); shading(hca,'flat'); colorbar
  [c,h] = contour(hca,vecBL*1e9,vecNL*1e-6,T,[10 20 30 40 50 70 100 150 200 300 400:200:600]); 
  clabel(c,h)
  hca.Title.String = 'Duration of reconnection event';
  irf_legend(hca,'s',[0.1,0.9],'color','k','fontsize',16)

  for isub = 1:4
    hca = subplot(nrows,ncols,isub); 
    hca.XLabel.String = 'B (nT)';
    hca.YLabel.String = 'n (cm^{-3})';
    hca.YScale = 'log';
    hca.FontSize = 14;
  end
end