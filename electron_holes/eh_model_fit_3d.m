%% Load data, each eh interval
ic = 1:4;
units = irf_units;
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint = tint + [+5 -5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % time interval of interest
t0 = tint_zoom(1) + (tint_zoom(2)-tint_zoom(1))*0.5; % center of time interval

localuser = datastore('local','user');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
db_info = datastore('mms_db');   

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',ic);
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

% Load timeseries of phi, probably not needed
vph = -9000e3;
c_eval('[phi?,phi_progressive?,phi_ancillary?] = get_phi(gseE?par,vph,tint_zoom,tint_zoom);',1:4)
tint_phi = phi1.time([1 end]);

% Load quantities for individual EHs, obtained from parallel timing
data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
neh = numel(data_tmp.EH_properties.vel);
tmp_fields = fields(data_tmp.EH_properties);
for ieh = 1:neh
  for ifield = 1:numel(tmp_fields)
    ehprop(ieh).(tmp_fields{ifield})   = data_tmp.EH_properties.(tmp_fields{ifield})(ieh,:);  
  end
end
% obs_eh_properties = data_tmp.EH_properties;
% obs_lpp = obs_eh_properties.Lpp; % peak to peak length
% obs_potential = obs_eh_properties.potential;
% obs_potential_max = obs_eh_properties.max_potential;
% obs_velocity = obs_eh_properties.vel;
% obs_neh = numel(obs_velocity);
% c_eval('obs_t0_epoch_mms? = obs_eh_properties.time_mms?;')
% c_eval('obs_phi? = irf.ts_scalar(obs_t0_epoch_mms?,obs_potential(:,?));')
% c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')
% % charge separation assuming potential structure is single gaussian
% % Can not be correct?
% dn = units.eps0/units.e*obs_potential_max./(obs_lpp*1e3)*1e-6; % cc

c_eval('gseR? = gseR?.resample(t0);',ic);
R0 = (gseR1.data + gseR2.data + gseR3.data + gseR4.data)/4; % center of tetrahedron
c_eval('gseR?rel = gseR? - R0;',ic); % distance in km from center of tetrahedron

%% Fit to triple Gaussian model to estimate phi0, lperp, lpar, test of optimization function
% If the structure is not axissymmetric, angle, will turn the perpendicular
% plane. Implement later. Seems to work now.

% Use symbolic expressions to get derivative of phi = E
syms x y z phi0 lx ly lz angle x0 y0 z0
R = [x y z];
% Triple Gaussian potential structure
%phi = @(x,y,z,phi0,lx,ly,lz,angle) phi0*exp(-0.5*(x/lx).^2-0.5*(y/ly).^2-0.5*(y/ly).^2-0.5*(z/lz).^2);
%phi = @(x,y,z,phi0,lx,ly,lz,angle) phi0*exp(-0.5*((x*cosd(angle)-y*sind(angle))/lx).^2-0.5*((x*sind(angle)+y*cosd(angle))/ly).^2-0.5*(z/lz).^2);
%phi = @(x,y,z,phi0,lx,ly,lz,angle) phi0*exp(-0.5*((x*cosd(angle)-y*sind(angle))/lx).^2-0.5*((x*sind(angle)+y*cosd(angle))/ly).^2-0.5*(z/lz).^2);
phi = phi0*exp(-0.5*(((x-x0)*cosd(angle)-(y-y0)*sind(angle))/lx).^2-0.5*(((x-x0)*sind(angle)+(y-y0)*cosd(angle))/ly).^2-0.5*(z/lz-z0).^2);

% Electric field from potential
E = -gradient(phi,R);

% Order of the inputs defined by argument 'Vars'
mf_phi = matlabFunction(phi,'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % V
mf_Ex = matlabFunction(E(1),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)
mf_Ey = matlabFunction(E(2),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)
mf_Ez = matlabFunction(E(3),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)

%
% We have lz from oservations, peak-to-peak length
% We will do the fit at z = 0 and z0 = 0 by time-shifting the timeseries
% with the dt corresponding to the parallel speed.
% These values listed below are the starting guessing values for the 
% optimization function fminsearch.
x0 = 0; % km
y0 = 0; % km
z0 = 0; % km
lx0 = 10; % km
ly0 = 15; % km
lz0 = 10; % km, from observations later
phi0 = 100; % V
philev = 0:10:phi0;
angle = 45;

xvec = 2*linspace(-lx0,lx0,30); % km 
yvec = 2*linspace(-ly0,ly0,31); % km
zvec = 2*linspace(-lz0,lz0,32); % km 
iz0 = numel(zvec)/2;

[X,Y,Z] = meshgrid(xvec,yvec,zvec);
PHI = mf_phi(X,Y,Z,phi0,lx0,ly0,lz0,angle,x0,y0,z0);
EX = mf_Ex(X,Y,Z,phi0,lx0,ly0,lz0,angle,x0,y0,z0);
EY = mf_Ey(X,Y,Z,phi0,lx0,ly0,lz0,angle,x0,y0,z0);
EZ = mf_Ez(X,Y,Z,phi0,lx0,ly0,lz0,angle,x0,y0,z0);

% Try to make a fit to synthetic data
% r1 = [10,10,0];  % km
% r2 = [5,-7,0];   % km
% r3 = [-5,-7,0];  % km
% r4 = [-10,10,0]; % km
r1 = gseR1rel.data; r1(end) = 0; % km
r2 = gseR1re2.data; r2(end) = 0; % km
r3 = gseR1re3.data; r3(end) = 0; % km
r4 = gseR1re4.data; r4(end) = 0; % km
c_eval('ex? = mf_Ex(r?(1),r?(2),r?(3),phi0,lx0,ly0,lz0,angle,x0,y0,z0);',1:4)
c_eval('ey? = mf_Ey(r?(1),r?(2),r?(3),phi0,lx0,ly0,lz0,angle,x0,y0,z0);',1:4)
c_eval('ez? = mf_Ez(r?(1),r?(2),r?(3),phi0,lx0,ly0,lz0,angle,x0,y0,z0);',1:4)

params0 = double([phi0; lx0+1; ly0-5; lz0; angle+10; x0; y0-7; z0]);
% I need to create a cost function, which is to be evaluated. For example
% ediffx = mf_Ex(vars) - Ex_obs
% ediffy = mf_Ey(vars) - Ey_obs
% ediffz = mf_Ez(vars) - Ez_obs
% cost_function = sqrt(ediffx.^2 + ediffx.^2 + ediffz.^2)

% xdata = [R1(1),R2(1),R3(1),R4(1)];
% ydata = [R1(2),R2(2),R3(2),R4(2)];
% zdata = [R1(3),R2(3),R3(3),R4(3)];
xdata = [r1(1),r2(1),r3(1),r4(1)];
ydata = [r1(2),r2(2),r3(2),r4(2)];
zdata = [r1(3),r2(3),r3(3),r4(3)];
Ex_data = [ex1,ex2,ex3,ex4];
Ey_data = [ey1,ey2,ey3,ey4];
Ez_data = [ez1,ez2,ez3,ez4];

cost_function = @(params) eh_costfunction(params,xdata,ydata,zdata,Ex_data,Ey_data,Ez_data,mf_Ex,mf_Ey,mf_Ez);
params = fminsearch(cost_function,params0);

mf_phi_best = @(x,y,z) mf_phi(x,y,z,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));
mf_Ex_best = @(x,y,z) mf_Ex(x,y,z,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));
mf_Ey_best = @(x,y,z) mf_Ey(x,y,z,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));
mf_Ez_best = @(x,y,z) mf_Ez(x,y,z,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));

PHI_best = mf_phi_best(X,Y,Z);
c_eval('Ex?_best = mf_Ex_best(r?(1),r?(2),r?(3));',1:4)
c_eval('Ey?_best = mf_Ey_best(r?(1),r?(2),r?(3));',1:4)
c_eval('Ez?_best = mf_Ez_best(r?(1),r?(2),r?(3));',1:4)

Ex_data_best = [Ex1_best,Ex2_best,Ex3_best,Ex4_best];
Ey_data_best = [Ey1_best,Ey2_best,Ey3_best,Ey4_best];
Ez_data_best = [Ez1_best,Ez2_best,Ez3_best,Ez4_best];



% tdata = timeline-mean(timeline);
% ydata = Bobs;
% yfit = @(t,B0,B1,dt) B0 + B1*tanh(t/dt);
% fun = @(x) sseval(x,tdata,ydata); % costfuction, root mean squared
% x0 = double([theta; 1]);
% %options = optimset
% bestx = fminsearch(fun,x0);
% Bopt = yfit(tdata,bestx(1),bestx(2),bestx(3));

colors = pic_colors('matlab');

nrows = 2;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

% original data
hca = h(isub); isub = isub + 1;
contour(hca,X(:,:,iz0),Y(:,:,iz0),PHI(:,:,iz0))
axis(hca,'equal')
hca.XGrid = 'on';
hca.YGrid = 'on';
hold(hca,'on')
quiver(hca,X,Y,EX,EY,'color',[0 0 0])
c_eval('quiver(hca,r?(1),r?(2),ex?,ey?,''linewidth'',1,''color'',colors(?,:))',1:4)
c_eval('plot(hca,r?(1),r?(2),''o'',''color'',colors(?,:))',1:4)
hold(hca,'off')

% fit
hca = h(isub); isub = isub + 1;
contour(hca,X(:,:,iz0),Y(:,:,iz0),PHI_best(:,:,iz0))
axis(hca,'equal')
hca.XGrid = 'on';
hca.YGrid = 'on';
hold(hca,'on')
c_eval('quiver(hca,r?(1),r?(2),Ex?_best,Ey?_best,''linewidth'',1,''color'',colors(?,:))',1:4)
c_eval('plot(hca,r?(1),r?(2),''o'',''color'',colors(?,:))',1:4)
hold(hca,'off')

disp('Done.')

%% Recalculate time shift... Is there something wrong?
%% Apply to MMS data
% We need to time shift the series and pick Eperp from the time where the
% Epar reverses direction.
% First do without loop, just choose one EH.
ieh = 10;
lpp = ehprop(ieh).Lpp; % km
vel = ehprop(ieh).vel; % km/s
%vel = 0.67*1e4;
tpp = lpp/abs(vel); % s
T = 4*tpp; % s

c_eval('t0_eh_mms? = ehprop(ieh).time_mms?;',1:4)
c_eval('tint_eh_mms? = t0_eh_mms? + 0.5*T*[-1 1];',1:4)
c_eval('E?par = gseE?par.tlim(tint_eh_mms?);',1:4)
c_eval('E?perp = gseE?perp.tlim(tint_eh_mms?);',1:4)

dt = [t0_eh_mms1,t0_eh_mms2,t0_eh_mms3,t0_eh_mms4]-t0_eh_mms1;
%dt(2) = -0.0027;
tint_eh_all = t0_eh_mms1 + T*[-1 1];

% We need to transform the data into a field aligned coordinate system
avB = (gseB1.resample(t0_eh_mms1).data + gseB2.resample(t0_eh_mms2).data + gseB3.resample(t0_eh_mms3).data + gseB4.resample(t0_eh_mms4).data)/4;
par = avB/norm(avB);
perp1 = cross(par,cross([0 1 0],par)); perp1 = perp1/norm(perp1);
perp2 = cross(par,perp1);
lmn = [perp1;perp2;par];
c_eval('E? = gseE?.tlim(tint_eh_mms?)*lmn'';',1:4)
c_eval('R? = gseR?rel.data*lmn'';',1:4)

% Get time delay from velocity
c_eval('dt_(?) = R?(3)/vel;',1:4)
[dt dt_'-dt_(1)]*1e3 % ms

h = irf_plot(6);
hca = irf_panel('Epar');
comp = 'z';
irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp')

hca = irf_panel('Epar dt');
irf_plot(hca,{E1.z,E2.z,E3.z,E4.z},'comp','dt',dt)

hca = irf_panel('Eperp1');
irf_plot(hca,{E1.x,E2.x,E3.x,E4.x},'comp')
hca = irf_panel('Eperp1 dt');
irf_plot(hca,{E1.x,E2.x,E3.x,E4.x},'comp','dt',dt)

hca = irf_panel('Eperp2');
irf_plot(hca,{E1.y,E2.y,E3.y,E4.y},'comp')
hca = irf_panel('Eperp2 dt');
irf_plot(hca,{E1.y,E2.y,E3.y,E4.y},'comp','dt',dt)

irf_zoom(h,'x',tint_eh_all)
c_eval('irf_pl_mark(h(!),t0_eh_mms?,h(!).ColorOrder(?,:))',1:4,[1 3 5])

% hca = irf_panel('Eperp x');
% irf_plot(hca,{E1perp.x,E2perp.x,E3perp.x,E4perp.x},'comp')
% hca = irf_panel('Eperp y');
% irf_plot(hca,{E1perp.y,E2perp.y,E3perp.y,E4perp.y},'comp')
% hca = irf_panel('Eperp z');
% irf_plot(hca,{E1perp.z,E2perp.z,E3perp.z,E4perp.z},'comp')

%% Get phase speeds again
% Speeds are in locally parallel coordinates dfined by perp1, perp2, par.
% tref is approximately the center of the ESW as seen by mms1, so that
% t1 = t_ref + dt(1), etc
% read in manual timeshift into structure 'manual', i saved them
% in a different function so that they can be easily called from different
% places
manual = edi_event_manual_dt;
if 0 
  %%
  manual(1).t_ref = '2017-07-06T13:54:05.533947021Z'; % 
  manual(1).gseBav = [-22.7878    1.6498   -1.2820];
  manual(1).perp1 = [0.0720    0.9974    0.0040];
  manual(1).perp2 = [0.0562         0   -0.9984];
  manual(1).par   = [-0.9958    0.0721   -0.0560];  
  manual(1).dt  = [0.0000,-0.00115,-0.00087,-0.00132]; 
  manual(1).v  = 8.78e+03*[1.00,-0.04, 0.05];
  
  manual(2).t_ref = '2017-07-06T13:54:05.538696044Z'; % t1 = t_ref + dt(1), etc
  manual(2).gseBav = [-22.7882    1.6481   -1.2911];
  manual(2).perp1 = [0.0719    0.9974    0.0041];
  manual(2).perp2 = [0.0566   -0.0000   -0.9984];
  manual(2).par   = [-0.9958    0.0720   -0.0564];   
  manual(2).dt  = [0.0000   -0.00075   -0.0011   -0.0015]; 
  manual(2).v = 8.02e+03 * [0.94 -0.35 -0.03]; % 
  
  manual(3).t_ref = '2017-07-06T13:54:05.547171386Z'; % t1 = t_ref + dt(1), etc
  manual(3).gseBav = [-22.7855    1.6436   -1.3151];
  manual(3).perp1 = [0.0717    0.9974    0.0041];
  manual(3).perp2 = [0.0576    0.0000   -0.9983];
  manual(3).par   = [-0.9958    0.0718   -0.0575];   
  manual(3).dt  = [0.0000   -0.00162   -0.00108   -0.00125]; 
  manual(3).v = 7.62e+03 * [0.98  0.20 -0.05]; % this one is ok, clearly something wrong with EH_properties7.62e+03 * [0.98  0.20 -0.05]
  
  manual(4).t_ref = '2017-07-06T13:54:05.549964843Z'; % t1 = t_ref + dt(1), etc
  manual(4).gseBav = [-22.7882    1.6481   -1.2911];
  manual(4).perp1 = [0.0718    0.9974    0.0042];
  manual(4).perp2 = [0.0579         0   -0.9983];
  manual(4).par   = [-0.9957    0.0720   -0.0578]; 
  manual(4).gseBav = [ -22.7852    1.6469   -1.3216]; 
  manual(4).dt  = [0.0000   -0.0013   -0.00093   -0.0013];
  manual(4).v = 8.44e+03 * [1.00  0.05  0.02]; % 
  
  manual(5).t_ref = '2017-07-06T13:54:05.556753417Z'; % t1 = t_ref + dt(1), etc
  manual(5).gseBav = [-22.7891    1.6627   -1.3376];
  manual(5).perp1 = [0.0725    0.9974    0.0043];
  manual(5).perp2 = [0.0586   -0.0000   -0.9983];
  manual(5).par   = [-0.9956    0.0726   -0.0584]; 
  manual(5).gseBav = [ -22.7891    1.6627   -1.3376]; 
  manual(5).dt  = [0.0000   -0.00105   -0.00115   -0.0012]; 
  manual(5).v = 8.74e+03 * [0.98 -0.09 -0.15]; % 
  
  manual(6).t_ref = '2017-07-06T13:54:05.562840820Z'; % t1 = t_ref + dt(1), etc
  manual(6).gseBav = [-22.7938    1.6766   -1.3445];
  manual(6).perp1 = [0.0731    0.9973    0.0043];
  manual(6).perp2 = [0.0589         0   -0.9983];
  manual(6).par   = [-0.9956    0.0732   -0.0587]; 
  manual(6).dt  = [0.0000,-0.00135,-0.00098,-0.00145];  
  manual(6).v  = 7.81e+03*[1.00  0.00  0.05];
  
  manual(7).t_ref = '2017-07-06T13:54:05.567746826Z'; % t1 = t_ref + dt(1), etc
  manual(7).gseBav = [-22.8033    1.6873   -1.3534];
  manual(7).perp1 = [0.0735    0.9973    0.0044];
  manual(7).perp2 = [0.0592   -0.0000   -0.9982];
  manual(7).par   = [-0.9955    0.0737   -0.0591]; 
  manual(7).dt  = [0.0000,-0.00063,-0.00080,-0.0014];  
  manual(7).v  = 9.12e+03 * [0.92 -0.37  0.09];
  
  manual(8).t_ref = '2017-07-06T13:54:05.572459472Z'; % t1 = t_ref + dt(1), etc
  manual(8).gseBav = [-22.8114    1.6932   -1.3658];
  manual(8).perp1 = [0.0738    0.9973    0.0044];
  manual(8).perp2 = [0.0598         0   -0.9982];
  manual(8).par   = [-0.9955    0.0739   -0.0596]; 
  manual(8).dt  = [0.0000,-0.00075,-0.00120,-0.0013];  
  manual(8).v  = 8.51e+03*[0.94,-0.30,-0.16];
  
  manual(9).t_ref = '2017-07-06T13:54:05.579765136Z'; % t1 = t_ref + dt(1), etc
  manual(9).gseBav = [-22.8235    1.7000   -1.3811];
  manual(9).perp1 = [0.0740    0.9972    0.0045];
  manual(9).perp2 = [0.0604         0   -0.9982];
  manual(9).par   = [-0.9954    0.0741   -0.0602]; 
  manual(9).dt  = [0.0000,-0.00097,-0.00086,-0.0012];  
  manual(9).v  = 9.66e+03*[1.00,-0.10, 0.01];
  
  manual(10).t_ref = '2017-07-06T13:54:05.585509765Z'; % t1 = t_ref + dt(1), etc
  manual(10).gseBav = [-22.8265    1.7105   -1.3870];
  manual(10).perp1 = [0.0745    0.9972    0.0045];
  manual(10).perp2 = [0.0606   -0.0000   -0.9982];
  manual(10).par   = [-0.9954    0.0746   -0.0605]; 
  manual(10).dt = [0.0000,-0.00085,-0.00115,-0.00145]; 
  manual(10).v = 8.06e+03*[0.96,-0.28,-0.07];
  
  manual(11).t_ref ='2017-07-06T13:54:05.604416748Z'; % t1 = t_ref + dt(1), etc
  manual(11).gseBav = [ -22.8326    1.7458   -1.4273];
  manual(11).perp1 = [0.0759    0.9971    0.0047];
  manual(11).perp2 = [0.0624   -0.0000   -0.9981];
  manual(11).par   = [-0.9952    0.0761   -0.0622]; 
  manual(11).dt = [0.0000,-0.00095,-0.00075,-0.0011]; 
  manual(11).v = 1.05e+04*[1.00,-0.05, 0.04];
  
  manual(12).t_ref = '2017-07-06T13:54:05.615465576Z';
  manual(12).gseBav = [-22.8243    1.7780   -1.4472];
  manual(12).perp1 = [0.0774    0.9970    0.0049];
  manual(12).perp2 = [0.0633    0.0000   -0.9980];
  manual(12).par   = [-0.9950    0.0775   -0.0631]; 
  manual(12).dt = [nan nan nan nan]; 
  manual(12).v = nan * [nan nan nan]; % unclear timing here, only ok between 3 and 4
  
  manual(13).t_ref = '2017-07-06T13:54:05.672972656Z';
  manual(13).gseBav = [-22.8477    1.8595   -1.5054];
  manual(13).perp1 = [0.0808    0.9967    0.0053];
  manual(13).perp2 = [0.0657         0   -0.9978];
  manual(13).par   = [-0.9946    0.0809   -0.0655]; 
  manual(13).dt = [0.0000,-0.00083,-0.00095,-0.00123]; 
  manual(13).v = 9.48e+03*[0.98,-0.21,-0.05];
end

ehprop_new = [];
for ieh = 13;%2;[1 6:13];%1:numel(ehprop)
  %% Get time
  lpp = ehprop(ieh).Lpp; % km
  vel = ehprop(ieh).vel; % km/s
  tpp = lpp/abs(vel); % s
  T = 2*tpp; % s
  T = 0.006;
  
  %c_eval('t0_eh_mms? = ehprop(ieh).time_mms?;',1:4)
  c_eval('t0_eh_mms? = EpochTT(manual(ieh).t_ref);',1:4)
  c_eval('tint_eh_mms? = t0_eh_mms? + 0.5*T*[-1 1];',1:4)  
  c_eval('E?par = gseE?par.tlim(tint_eh_mms?);',1:4)
  c_eval('E?perp = gseE?perp.tlim(tint_eh_mms?);',1:4)

  % Get data from old data (which are somewhat wrong)
  
  dt_tmp = [t0_eh_mms1,t0_eh_mms2,t0_eh_mms3,t0_eh_mms4]-t0_eh_mms1;

  tint_eh_all = t0_eh_mms1 + [min(dt_tmp) max(dt_tmp)]+0.5*T*[-1 1];

  % We need to transform the data into a field aligned coordinate system
  avB = (gseB1.resample(t0_eh_mms1).data + gseB2.resample(t0_eh_mms2).data + gseB3.resample(t0_eh_mms3).data + gseB4.resample(t0_eh_mms4).data)/4;
  par = avB/norm(avB);
  perp1 = cross(par,cross([0 1 0],par)); perp1 = perp1/norm(perp1);
  perp2 = cross(par,perp1);
  lmn = [perp1;perp2;par];
  c_eval('E? = gseE?.tlim(tint_eh_mms?)*lmn'';',1:4)
  c_eval('R?rel = gseR?rel.data*lmn'';',1:4)
  pppB = avB*lmn';
  
  %tint_vicinity = EpochTT([esw_data{1}{ii}; esw_data{2}{ii}]);
  %tint_esw = EpochTT([esw_data{3}{ii}; esw_data{4}{ii}]);
  %t_center = EpochTT(esw_data{5}{ii});
  
  %% Get properties thorugh Gaussian fits
  % First, make a fit to gaussian, and the relative zero time is used to
  % get relative times of observations, and speed
  % Note: if timeseries is used to make fit, potential and lengths needs to
  % be multiplied with velocity to get the correct units.
  % Parallel component is z.
  c_eval('esw_gaussian_fit? = fit_esw_gaussian(E?.z);',1:4)
  
  t0_gaussian_fit = [esw_gaussian_fit1.x0,...
                     esw_gaussian_fit2.x0,...
                     esw_gaussian_fit3.x0,...
                     esw_gaussian_fit4.x0];
  dt_gaussian_fit = t0_gaussian_fit - t0_gaussian_fit(1);
             
  % Calculate velocity based on spacecraft position and timeshift
  % Based on mms.mms4_v, but easier to apply directly.
  DD = [R1rel; R2rel; R3rel; R4rel];
  TT = dt_gaussian_fit;
  mm = DD\TT;  
  v_gaussian_fit = mm/norm(mm)/norm(mm);

  % Collect results into structure
  fit_gaussian(ieh).gseBav = avB;
  fit_gaussian(ieh).note1 = 'ppp indicates par/perp1/perp2 coordinate system, defined by newxyz';
  fit_gaussian(ieh).pppBav = pppB;
  fit_gaussian(ieh).newxyz = lmn;
  c_eval('fit_gaussian(ieh).mms? = esw_gaussian_fit?;',1:4)  
  fit_gaussian(ieh).pppE = {E1,E2,E3,E4};
  fit_gaussian(ieh).t0 = t0_gaussian_fit;
  fit_gaussian(ieh).dt = dt_gaussian_fit;
  fit_gaussian(ieh).v = v_gaussian_fit;
  fit_gaussian(ieh).vpar = v_gaussian_fit(3);
  fit_gaussian(ieh).pitch_angle = acosd(dot(v_gaussian_fit/norm(v_gaussian_fit),pppB/norm(pppB)));
  
  %% Get properties through correlating fields with xcorr
  % Time shift is comparalbe to instrument sampling frequency, therefore,
  % we need to artifically increase the sampling frequency.
  dt_sampling_original = gseE1.time(2)-gseE1.time(1);
  upsampling_factor = 2;
  dt_sampling = dt_sampling_original*upsampling_factor;
  tint_vicinity = tint_eh_all;
  timeline = tint_vicinity(1):dt_sampling:tint_vicinity(2);
  
  % Upsample fields to same timeline and rotate into par/perp system
  c_eval('pppE? = gseE?.tlim(tint_vicinity).resample(timeline)*lmn'';',1:4)      
  
  % Cross correlate mms2-4 to mms1
  % Should cross correlate all with all, to compare
  dt_xcorr_all = zeros(4,4);
  C_xcorr = zeros(4,4);
  for ic1 = 1:4
    for ic2 = 1:4
      c_eval('[tmpC,lags] = xcorr(pppE?.z.data,pppE!.z.data,''coeff'');',ic1,ic2)  
      i_shift = find(abs(tmpC) == max(abs(tmpC)));
      C_xcorr(ic1,ic2) = tmpC(i_shift);
      di = -lags(i_shift);
      dt_xcorr_all(ic1,ic2) = di*dt_sampling;
    end
  end
   
  % Specify which reference spacecraft to use for dt = 0
  ic_ref = 1;
  dt_xcorr = tocolumn(dt_xcorr_all(ic_ref,:));
  DD = [R1rel; R2rel; R3rel; R4rel];  
  TT = dt_xcorr;
  mm = DD\TT;  
  v_xcorr = mm/norm(mm)/norm(mm);
  
  % Collect data into structure
  struct_xcorr(ieh).gseBav = avB;
  struct_xcorr(ieh).note1 = 'ppp indicates par/perp1/perp2 coordinate system, defined by newxyz';
  struct_xcorr(ieh).pppBav = pppB;
  struct_xcorr(ieh).newxyz = lmn;  
  struct_xcorr(ieh).pppE = {E1,E2,E3,E4};  
  struct_xcorr(ieh).dt_all = dt_xcorr_all;
  struct_xcorr(ieh).sc_ref = ic_ref;
  struct_xcorr(ieh).dt = dt_xcorr;
  struct_xcorr(ieh).v = v_xcorr;
  struct_xcorr(ieh).vpar = v_xcorr(3); % since third component is parallel to B
  struct_xcorr(ieh).pitch_angle = acosd(dot(v_xcorr/norm(v_xcorr),pppB/norm(pppB)));
  
  %% Get time delay manually
  dt_manual = manual(ieh).dt;
  
  %% Plot results
  if 1
    %%
  figure(78)
  h = irf_plot(5);
  hca = irf_panel('Epar');
  comp = 'z';
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp')
  hca.YLabel.String = {'E','original'};
  
  hca = irf_panel('Epar dt gaussian');
  comp = 'z';
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt_gaussian_fit)
  irf_legend(hca,sprintf('dt = [%.2f,%.2f,%.2f,%.2f] ms',dt_gaussian_fit(1)*1e3,dt_gaussian_fit(2)*1e3,dt_gaussian_fit(3)*1e3,dt_gaussian_fit(4)*1e3),[0.02 0.98])
  hca.YLabel.String = {'E','dt from','gaussian fit'};
  
  hca = irf_panel('Epar fit gaussian');
  c_eval('fitE? = irf.ts_scalar(fit_gaussian(ieh).mms?.x_ref+fit_gaussian(ieh).mms?.x,fit_gaussian(ieh).mms?.fun_E(fit_gaussian(ieh).mms?.x));',1:4)
  irf_plot(hca,{fitE1,fitE2,fitE3,fitE4},'comp','dt',dt_gaussian_fit)
  irf_legend(hca,sprintf('dt = [%.2f,%.2f,%.2f,%.2f] ms',dt_gaussian_fit(1)*1e3,dt_gaussian_fit(2)*1e3,dt_gaussian_fit(3)*1e3,dt_gaussian_fit(4)*1e3),[0.02 0.98])
  hca.YLabel.String = {'E', 'gaussian fit'};
  
  hca = irf_panel('Epar dt xcorr');
  comp = 'z';
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt_xcorr)
  irf_legend(hca,sprintf('dt = [%.2f,%.2f,%.2f,%.2f] ms',dt_xcorr(1)*1e3,dt_xcorr(2)*1e3,dt_xcorr(3)*1e3,dt_xcorr(4)*1e3),[0.02 0.98])
  hca.YLabel.String = {'E','dt from xcorr'};
  
  hca = irf_panel('Epar dt manual;');
  comp = 'z';
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt_manual)
  irf_legend(hca,sprintf('dt = [%.2f,%.2f,%.2f,%.2f] ms',dt_manual(1)*1e3,dt_manual(2)*1e3,dt_manual(3)*1e3,dt_manual(4)*1e3),[0.02 0.98])
  hca.YLabel.String = {'E','dt manual'};
  
  h(1).Title.String = sprintf('ieh = %g',ieh);
  irf_zoom(h,'x',tint_eh_all)
  %cn.print(sprintf('dt_ieh=%g_T_smaller',ieh),'path',printpath)
  pause%(0.5)
  end  
end
ehprop_new.gauss = fit_gaussian;
ehprop_new.xcorr = struct_xcorr;
ehprop_new.manual = manual;

%% Go through all EHs to get their time delay manually.
% This is done and saved in structure array manual() above.
for ieh = 1%13;[1 6:13];%1:numel(ehprop)
  %% Get approximate times
  % Get data from old data (which are somewhat wrong)  
  lpp = ehprop(ieh).Lpp; % km
  vel = ehprop(ieh).vel; % km/s
  tpp = lpp/abs(vel); % s
  T = 2*tpp; % s
  T = 0.006;
  
  c_eval('t0_eh_mms? = ehprop(ieh).time_mms?;',1:4)
  % for some ehs, the eh was not seen and is not defined, then set time to
  % mms1, who saw all of them (at least the ones who were in the list)
  c_eval('if t0_eh_mms? < tint(1), t0_eh_mms? = t0_eh_mms1; end',1:4)
  c_eval('tint_eh_mms? = t0_eh_mms? + 0.5*T*[-1 1];',1:4)  
  c_eval('E?par = gseE?par.tlim(tint_eh_mms?);',1:4)
  c_eval('E?perp = gseE?perp.tlim(tint_eh_mms?);',1:4)

  
  dt_tmp = [t0_eh_mms1,t0_eh_mms2,t0_eh_mms3,t0_eh_mms4]-t0_eh_mms1;
  tint_eh_all = t0_eh_mms1 + [min(dt_tmp) max(dt_tmp)]+0.5*T*[-1 1];

  % We need to transform the data into a field aligned coordinate system
  avB = (gseB1.resample(t0_eh_mms1).data + gseB2.resample(t0_eh_mms2).data + gseB3.resample(t0_eh_mms3).data + gseB4.resample(t0_eh_mms4).data)/4;
  par = avB/norm(avB);
  perp1 = cross(par,cross([0 1 0],par)); perp1 = perp1/norm(perp1);
  perp2 = cross(par,perp1);
  lmn = [perp1;perp2;par];
  
  % Include a slightly longer time interval, such that the structure is
  % more easily observed and identified
  c_eval('E? = gseE?.tlim(tint_eh_all)*lmn'';',1:4)
  c_eval('R?rel = gseR?rel.data*lmn'';',1:4)
  pppB = avB*lmn';
  irf_4_v_gui(E1.resample(E1).z,E2.resample(E1).z,E3.resample(E1).z,E4.resample(E1).z,gseR1.resample(E1),gseR2.resample(E1),gseR3.resample(E1),gseR4.resample(E1))
end

%% Get remaining ESW properties: length scale, potential
irf_log OFF
manual = edi_event_manual_dt;
data = [];
colors = mms_colors('1234');
% Run through to estimate peak-to-peak length scales, and possible
% potentials
for ieh = 1%:numel(manual)
  T = 0.006;
  tref = EpochTT(manual(ieh).t_ref);
  dt = manual(ieh).dt;
  tminus = EpochTT(manual(ieh).tminus);
  tplus = EpochTT(manual(ieh).tplus);
  tp1 = tref + manual(ieh).tp1;
  tp2 = tref + manual(ieh).tp2;
  tpp = tp2-tp1;
  lmn = [manual(ieh).perp1; manual(ieh).perp2; manual(ieh).par]; 
  gsev = manual(ieh).v; % velocity is given in gse coordinate system
  pppv = gsev*lmn';
  pppvnorm = pppv/norm(pppv); % unit vector of velocity
  vpar = pppv(3);  
  tint_eh_all = tref + [min(dt) max(dt)]+0.5*T*[-1 1];
  c_eval('E? = gseE?.tlim(tint_eh_all)*lmn'';',1:4) % use same time interval for all
  c_eval('R?rel = gseR?rel.data*lmn'';',1:4)
  pppB = avB*lmn';
  
  c_eval('intE? = irf_integrate(E?.z,tref+dt(?));',1:4) % (mV/m)*s
  c_eval('phi? = vpar*intE?;',1:4) % (km/s)*(mV/m)*s = V  
  
  % Try to estimate what the potential should be based on the value at
  % tminus and tplus  
  %c_eval('tminus? = tminus + dt(?);',1:4)
  %c_eval('tplus?  = tplus  + dt(?);',1:4)    
  c_eval('phiminus? = phi?.resample(tminus + dt(?));',1:4)
  c_eval('phiplus?  = phi?.resample(tplus + dt(?));',1:4)
  c_eval('phiminusplus?  = irf.ts_scalar([tminus tplus] + dt(?),[phiminus?.data,phiplus?.data]);',1:4)  
  % c_eval('phiminusplus?  = phi?.resample([tminus tplus] + dt(?));',1:4)
  % This above is not working, and I have no idea why...
  c_eval('phimean(?) = 0.4*(phiplus?.data + phiminus?.data);',1:4)  
  
  data(ieh).tpp = tpp;  
  data(ieh).lpp = tpp*vpar;
  data(ieh).phiminus = [phiminus1.data phiminus2.data phiminus3.data phiminus4.data];
  data(ieh).phiplus = [phiplus1.data phiplus2.data phiplus3.data phiplus4.data];  
  data(ieh).phimean = phimean;
  data(ieh).vpar = vpar;
  data(ieh).pitchangle = acosd(pppvnorm(3));
  
  
  % Plot
  if 0
    figure(79)
    h = irf_plot(3);

    hca = irf_panel('Epar');  
    set(hca,'ColorOrder',mms_colors('1234'))
    comp = 'z';
    irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp')
    hca.YLabel.String = {'E_{||} (mV/m)'};

    hca = irf_panel('Epar dt');  
    set(hca,'ColorOrder',mms_colors('1234'))
    comp = 'z';
    irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt)
    hca.YLabel.String = {'E_{||} (mV/m)'};
    irf_legend(hca,{sprintf('<tpp>=%.0f',tpp(1));...
                    sprintf('<tpp>=%.0f',tpp(2));...
                    sprintf('<tpp>=%.0f',tpp(3));...
                    sprintf('<tpp>=%.0f',tpp(4))},[0.02 0.02])

    hca = irf_panel('phi');    
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,{phi1,phi2,phi3,phi4},'comp','dt',dt)
    %set(hca,'ColorOrder',mms_colors('11223344'))
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,{phiminusplus1,phiminusplus2,phiminusplus3,phiminusplus4},'comp','dt',dt);
    irf_legend(hca,{sprintf('<phi1>=%.0f',phimean(1));...
                    sprintf('<phi2>=%.0f',phimean(2));...
                    sprintf('<phi3>=%.0f',phimean(3));...
                    sprintf('<phi4>=%.0f',phimean(4))},[0.02 0.02])
     %c_eval('hlines(?).Marker = ''o'';',1:2)
  %   set(hca,'ColorOrder',mms_colors('22'))
  %   hlines = irf_plot(hca,{phiminus2,phiplus2},'comp');
  % %   c_eval('hlines(?).Marker = ''o'';',1:2)   
  %   hlines = irf_plot(hca,{phiminus3,phiplus3},'comp');
  %  %  c_eval('hline.s(?).Marker = ''o'';',1:2)   
  %   hlines = irf_plot(hca,{phiminus4,phiplus4},'comp');
  %   % c_eval('hlines(?).Marker = ''o'';',1:2)
  %   hca.YLabel.String = {'\phi (V)'};

    h(1).Title.String = sprintf('ieh = %g',ieh);
    irf_zoom(h,'x',tint_eh_all)
    irf_pl_mark(h,tref,[0 0 0])
    irf_pl_mark(h,[tminus tplus],[1 1 0])
    c_eval('irf_pl_mark(h,tp1(?),colors(?,:))',1:4)
    c_eval('irf_pl_mark(h,tp2(?),colors(?,:))',1:4)

    fprintf('ieh = %g\n',ieh)

    %[tphi,levphi] = get_time(2,'epochtt');
    %fprintf('manual(ieh).tminus = ''%s'', manual(ieh).tplus = ''%s'' \n',tphi(1).utc,tphi(2).utc)

    %[tpeaks,levpeaks] = get_time(8,'epochtt');
    %tp1 = [tpeaks(1),tpeaks(3),tpeaks(5),tpeaks(7)]-tref;
    %tp2 = [tpeaks(2),tpeaks(4),tpeaks(6),tpeaks(8)]-tref;  
    %fprintf('manual(ieh).tp1 = [%8.5f,%8.5f,%8.5f,%8.5f];\n',tp1(1),tp1(2),tp1(3),tp1(4))
    %fprintf('manual(ieh).tp2 = [%8.5f,%8.5f,%8.5f,%8.5f];\n',tp2(1),tp2(2),tp2(3),tp2(4))
    %pause
  end
end
irf_log ON
