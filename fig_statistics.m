% figure

% load Daniels data
load('/Users/Cecilia/Research/EH2/Daniel/wavespeeds.mat')
yES = wavespeeds.vESwaves./wavespeeds.vtrESwaves;
xES = wavespeeds.vESwaves./wavespeeds.vthESwaves;
yEH = wavespeeds.veholes./wavespeeds.vtreholes;
xEH = wavespeeds.veholes./wavespeeds.vtheholes;
%%

% Fit to stats
xFitStats = linspace(1e-3,1e2,10000);

% Daniels fit
pexp=1;
pcoeff=15;
yFitStatsDG = pcoeff*xFitStats.^pexp;
strFitStatsDG = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];

% My fit by eye
pexp=0.9;
pcoeff=10;
yFitStatsCN = pcoeff*xFitStats.^pexp;
strFitStatsCN = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];

% Electrostatic waves
curve = fit(xES',yES','power1');%@(y) a*x*exp(b));
pexp = curve.b*1;
pcoeff = curve.a*1;
yFitStatsES = curve(xFitStats);
strFitStatsES = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];

% Electron holes
curve = fit(xEH',yEH','power1');%@(y) a*x*exp(b));
pexp = curve.b*1;
pcoeff = curve.a*1;
yFitStatsEH = curve(xFitStats);
strFitStatsEH = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];

% Electron holes + electrostatic waves
curve = fit([xES'; xEH'],[yES'; yEH'],'power1');%@(y) a*x*exp(b));
pexp = curve.b*1;
pcoeff = curve.a*1;
yFitStatsALL = curve(xFitStats);
strFitStatsALL = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];


if 0 % illustrate fit
    loglog(xEH,yEH,'ko',xES,yES,'rx',...
        xFitStats,yFitStatsEH,'k--',...
        xFitStats,yFitStatsES,'r--',...
        xFitStats,yFitStatsDG,'g--',...
        xFitStats,yFitStatsCN,'c--',...
        xFitStats,yFitStatsALL,'b--')
    set(gca,'ylim',[1e-2 1e1],'xlim',[1e-3 0.6e0])
end
%


% Fit to model
pexp = 1;
pcoeff = 1.5;
xFitModel = linspace(1e-4,1e2,10);
yFitModel = pcoeff*xFitModel.^pexp;
strFitModel = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];



if 1
    loglog(xES,yES,'rx',xEH,yEH,'ko',xFitStats,yFitStatsCN,'k--','linewidth',2)
    set(gca,'ylim',[1e-2 1e1],'xlim',[1e-3 1e0])
    xlabel('v_{ph}/v_{te}')
    ylabel('v_{ph}/v_{T}')
else
    plot(1./xES,1./yES,'rx',1./xEH,1./yEH,'ko',1./xFitStats,1./yFitStatsCN,'k--','linewidth',2)
    set(gca,'ylim',[6e-2 3e1],'xlim',[1e0 8e2])
    xlabel('v_{te}/v_{ph}')
    ylabel('v_{T}/v_{ph}')
end

set(gcf,'defaultAxesFontSize',16);
set(gcf,'DefaultTextFontSize',16);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
