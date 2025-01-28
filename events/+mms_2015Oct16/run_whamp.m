cd /Users/Cecilia/Research/2015Oct16/Whamp

n = [6 14];
m = [0 0];
t = [05 30]*1e-3;
vd = [00 0];
d = [0.001 1];
a1 = [3 0.05];
a2 = [6 1];


h = subplot(2,2,1);
toPlot = [1];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar

h = subplot(2,2,2);
toPlot = [2];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar

h = subplot(2,2,3);
toPlot = [2 1];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar

%% try with negative densities
n = [-23 40]*1e6;
m = [0 0];
t = [24 35]*1e-3;
vd = [0 0];
d = [0.2 0.5];
a1 = [2.1 1.7];
a2 = [1.1 0];

n = [-23 40]*1e6;
m = [0 0];
t = [23 35]*1e-3;
vd = [0 0];
d = [0.2 0.5];
a1 = [2.2 1.7];
a2 = [1.1 0];
toPlot = [1 2];

vlim = 1e7*[-1 1];
h = subplot(2,2,1);
toPlot = [1];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar
h.XLim = vlim; h.YLim = h.XLim;

h = subplot(2,2,2);
toPlot = [2];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar
h.XLim = vlim; h.YLim = h.XLim;

h = subplot(2,2,3);
toPlot = [1 2  ];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar
h.XLim = vlim; h.YLim = h.XLim;

h = subplot(2,2,4);
toPlot = [2 1];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
h.YScale = 'log';
h.XScale = 'log';
h.XLim = [1e1 5e3];
h.YLim = [1e-2  1e6];
h.XTick = [1e-1 1e0 1e1 1e2 1e3];

% add real distribution
tint = irf.tint('2015-10-16T10:33:44.90Z/2015-10-16T10:33:45.00Z'); % magnetosphere-magnetosheath-magnetosphere
tint = irf.tint('2015-10-16T10:33:44.90Z/2015-10-16T10:33:44.93Z')+0.03; % magnetosphere-magnetosheath-magnetosphere
hold(h,'on')
mms.plot_pitchangles(h,ePDist1,dmpaB1,'tint',tint,'scPot',scPot1,'pitchangle',[0 90 180])

colormap('jet')


%% run whamp
distribution = 3;
switch distribution
  case 1
    n = [-23 40];
    m = [0 0];
    t = [24 35];
    vd = [0 0];
    d = [0.2 0.5];
    a1 = [2.2 1.7];
    a2 = [1.1 0];
    toPlot = [1 2];
  case 2
    n = [-24 40];
    m = [0 0];
    t = [26 35];
    vd = [0 0];
    d = [0.2 0.5];
    a1 = [2.25 1.9];
    a2 = [1.1 0];
    toPlot = [1 2];
  case 3
    n = [-23 40];%*1e6;
    m = [0 0];
    t = [23 35];%*1e-3;
    vd = [0 0];
    d = [0.2 0.5];
    a1 = [2.2 1.7];
    a2 = [1.1 0];
    toPlot = [1 2];
end


clear Species;
for ii = toPlot
  species.n = n(ii);
  species.m = m(ii);
  species.t = t(ii);
  species.vd = vd(ii);
  species.a = a1(ii);
  species.a1 = d(ii);
  species.b = a2(ii);
  Species{ii} = species;
end
  
PlasmaModel.B = 17;
PlasmaModel.Species = Species;

InputParameters.fstart = 0.000001; % Hz,  start frequency
kperp = logspace(-5,0.4,30);
InputParameters.kperp = [log10(kperp(1)) log10(kperp(2))-log10(kperp(1)) log10(kperp(end))];   % - kp (perpendicular wave vector) as scalar or, [kp_start kp_step kp_end]
kpar = logspace(-2.10,0.4,30);
InputParameters.kpar = [log10(kpar(1)) log10(kpar(2))-log10(kpar(1)) log10(kpar(end))];    % - kz (parallel wave vector) or [kz_start kz_step kz_end]
InputParameters.varyKzFirst = 0; % - 1 for each kp value, step through all kz values, - s0 for each kz value, step through all kp values
InputParameters.useLog= 1; % - 1 input log10(p) and log10(z),(default), 0 input given as p and z
InputParameters.maxIterations = 150; % - maximum number of iterations when searching for solution (default 50)


Output = whamp.run(PlasmaModel,InputParameters);

InputParameters.fstart = 0.1; % Hz,  start frequency
Output2 = whamp.run(PlasmaModel,InputParameters);


InputParameters.fstart = 0.0001; % Hz,  start frequency
kperp = logspace(-3,0.4,20);
InputParameters.kperp = [log10(kperp(1)) log10(kperp(2))-log10(kperp(1)) log10(kperp(end))];   % - kp (perpendicular wave vector) as scalar or, [kp_start kp_step kp_end]
InputParameters.kperp = [0 0.0020 0.15];   % - kp (perpendicular wave vector) as scalar or, [kp_start kp_step kp_end]
kpar = [0.002,0.0010,0.15];
InputParameters.kpar = kpar;%[log10(kpar(1)) log10(kpar(2))-log10(kpar(1)) log10(kpar(end))];    % - kz (parallel wave vector) or [kz_start kz_step kz_end]
InputParameters.useLog = 0;
Output3 = whamp.run(PlasmaModel,InputParameters);


fce = 0.0279928*1000*PlasmaModel.B; % fce in Hz

%Output = whamp.run([],[]);

h = subplot(2,2,1);

whamp.plot_f(h,n(toPlot)*1e6,m(toPlot),t(toPlot)*1e-3,vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
h.YScale = 'log';
h.XScale = 'log';

h.XLim = [1e1 5e3];
h.YLim = [1e-2  1e6];
h.XTick = [1e-1 1e0 1e1 1e2 1e3];
% add real distribution
%tint = irf.tint('2015-10-16T10:33:44.40Z',0.033*70); % magnetosphere-magnetosheath-magnetosphere
%time = tint(2);
time = tint(1);
hold(h,'on')
mms.plot_pitchangles(h,ePDist1,dmpaB1,'tint',tint,'scPot',scPot1,'pitchangle',[0 90 180])
hold(h,'off')

h = subplot(2,2,2);
%surf(h,Output.kperp,Output.kpar,(real(Output.f)),imag(Output.f))
%hold(h,'on')
%surf(h,Output2.kperp,Output2.kpar,(real(Output2.f)),imag(Output2.f))
surf(h,Output3.kperp,Output3.kpar,(real(Output3.f))/fce,imag(Output3.f))
%hold(h,'off')
shading flat

hca = h;
hca.XLabel.String = 'k_{\perp}';
hca.YLabel.String = 'k_{||}';
hca.ZLabel.String = 'f/f_{ce}';
hca.YScale = 'log';
hca.XScale = 'log';
%hca.ZLim = [76 82];
%hca.ZScale = 'log';
hcb = colorbar('peer',hca);
hca.CLim = max(abs(hca.CLim))*[-1 1]*1;
irf_colormap('poynting')
%shading(hca,'flat');

h = subplot(2,2,3);
plot(h,Output3.kpar,(real(Output3.f(:,1)))/fce,Output3.kpar,imag(Output3.f(:,1)))