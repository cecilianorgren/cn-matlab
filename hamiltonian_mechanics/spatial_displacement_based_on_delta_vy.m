%
% 
dvy = [-10000e3 -20000e3 -30000e3]; % m/s
q = -units.e;
m = units.me;
dAy = -m/q*dvy;

legends_dvy = arrayfun(@(x)sprintf('%.0f km/s',x*1e-3),dvy,'UniformOutput',false);

Ay =  @(z,B0,lz) -lz*B0*log(cosh(z/lz));
Bx = @(z,B0,lz) B0*tanh(z/lz);

B0 = 10e-9;
lz = 100e3;
zvec = linspace(-2*lz,2*lz,100);

nrows = 2;
ncols = 1;
clear h
ip = 0;
for irow = 1:nrows
  for icols = 1:ncols
    ip = ip + 1;
    h(irow,icols) = subplot(nrows,ncols,ip);
  end
end

isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,zvec*1e-3,Bx(zvec,B0,lz)*1e9)
hca.XLabel.String = 'z (km)';
hca.YLabel.String = 'B_x (nT)';

hca = h(isub); isub = isub + 1;
hA = plot(hca,zvec*1e-3,Ay(zvec,B0,lz)); % to make nTkm: *1e9*1e-3
hold(hca,'on')
hdA = plot(hca,zvec*1e-3,repmat(dAy,numel(zvec),1));
hold(hca,'off')
hca.XLabel.String = 'z (km)';
hca.YLabel.String = 'A_y (Tm)';
hleg = legend(hdA,legends_dvy,'location','best');
hleg.Title.String = '\Delta v_y';

% Formatting
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 14;',1:numel(h))
