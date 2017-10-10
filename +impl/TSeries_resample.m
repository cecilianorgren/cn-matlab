% Make TSeries
scalarT = EpochTT('2002-03-04T09:30:00Z'):50:EpochTT('2002-03-04T10:30:00Z');
t = scalarT - scalarT.start;
x = exp(0.001*t).*sin(2*pi*t/180*0.1);        % define function x(t)=exp(0.001(t-to))*sin(t-to)
scalarTS = irf.ts_scalar(scalarT,x);                   % define scalar TSeries object
scalarTS2 = irf.ts_scalar(scalarT,[x x*0.9 x*0.8 x*0.6 x*1.2]);                   % define scalar TSeries object

xyzT = EpochTT('2002-03-04T09:30:00Z'):15:EpochTT('2002-03-04T10:30:00Z');
t = xyzT - xyzT.start;
x = exp(0.001*t).*sin(2*pi*t/180*0.1);        % define function x(t)=exp(0.001(t-to))*sin(t-to)
y = exp(0.001*t).*cos(2*pi*t/180*0.1);	% z(t)=exp(0.001(t-to))*cos(t)
xyzTS = irf.ts_vec_xyz(xyzT,[x y y*0.5]);				  % B2 has two components, x & y

plus = scalarTS + 2*scalarTS;
minus = scalarTS - 2*scalarTS;
irf_plot({plus,minus,scalarTS+2},'comp')
%% Resampling a tensorOrder = 0, to the time series of tensorOrder 1

scalarTSresamp = scalarTS.resample(xyzTS); 
h = irf_plot({scalarTS,scalarTSresamp,xyzTS},'.');
h(1).YLabel.String = 'scalar';
h(2).YLabel.String = 'scalar.resample(xyz.time)';
h(3).YLabel.String = 'xyz';
% this gives error since (line 815 in TSeries):
% if obj.tensorOrder~=1, error('Not yet implemented'); end
%%
scalarTS2resamp = scalarTS2.resample(xyzTS); 
h = irf_plot({scalarTS2,scalarTS2resamp,xyzTS},'.');
h(1).YLabel.String = 'scalar';
h(2).YLabel.String = 'scalar.resample(xyz.time)';
h(3).YLabel.String = 'xyz';
%%
xyzTSresamp = xyzTS.resample(scalarTS);
h = irf_plot({scalarTS,scalarTSresamp,xyzTS,xyzTSresamp},'.');
h(1).YLabel.String = 'scalar';
h(2).YLabel.String = 'scalar.resample(xyz.time)';
h(3).YLabel.String = 'xyz';
h(4).YLabel.String = 'xyz.resample(scalar.time)';

%%
j, dmpaB1, ne1
EJxB1 = j.cross(dmpaB1)*1e-9*1e-9/e/(ni1*1e6)*1e3;
j, dmpaB2, ne2
EJxB2 = j.cross(dmpaB2)*1e-9*1e-9/e/(ni2*1e6)*1e3;