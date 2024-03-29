%% Set up database
localuser = datastore('local','user');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS'); 
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   

%% Time interval of event
tint_burst = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint_burst = tint_burst + [+5 -5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file

tint_figure = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.620Z');
t0 = tint_figure(1) + (tint_figure(2)-tint_figure(1))*0.5; % center of time interval

%% Load data
doLoad = 1;
ic = 1:4;
units = irf_units;

% New wave properties
manual = edi_event_manual_dt;

if doLoad
  c_eval('gseR? = mms.get_data(''R_gse'',tint_burst,?);',ic);
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint_burst);',ic);
  c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint_burst);',ic);
  c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)  
end

c_eval('gseR? = gseR?.resample(t0);',ic);
R0 = (gseR1.data + gseR2.data + gseR3.data + gseR4.data)/4; % center of tetrahedron
c_eval('gseR?rel = gseR? - R0;',ic); % distance in km from center of tetrahedron

%% Get remaining ESW properties that are not saved: potential
irf_log OFF

data = [];
colors = mms_colors('1234');
% Run through to estimate peak-to-peak length scales, and possible
% potentials
for ieh = 1:numel(manual)
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
  gseBav = manual(ieh).gseBav;
  pppB = gseBav*lmn';
  
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
  c_eval('phimean(?) = 0.5*(phiplus?.data + phiminus?.data);',1:4)  
  
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

%% Fit to triple Gaussian model to estimate phi0, lperp, lpar, test of optimization function
% If the structure is not axissymmetric, angle, will turn the perpendicular
% plane.

% Use symbolic expressions to get derivative of phi = E
syms x y z phi0 lx ly lz angle x0 y0 z0
R = [x y z];
% Triple Gaussian potential structure
phi = phi0*exp(-0.5*(((x-x0)*cosd(angle)-(y-y0)*sind(angle))/lx).^2-0.5*(((x-x0)*sind(angle)+(y-y0)*cosd(angle))/ly).^2-0.5*(z/lz-z0).^2);

% Electric field from potential
E = -gradient(phi,R);

% Order of the inputs defined by argument 'Vars'
mf_phi = matlabFunction(phi,'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % V
mf_Ex = matlabFunction(E(1),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)
mf_Ey = matlabFunction(E(2),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)
mf_Ez = matlabFunction(E(3),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)

% We have lz from oservations, peak-to-peak length
% We will do the fit at z = 0 and z0 = 0 by time-shifting the timeseries
% with the dt corresponding to the parallel speed.
% These values listed below are the starting guessing values for the 
% optimization function fminsearch.
x0 = 0; % km
y0 = 0; % km
z0 = 0; % km
lx0 = 5; % km
ly0 = 5; % km
lz0 = 5; % km, from observations later
phi0 = 200; % V
angle = 0;

philev = 0:10:phi0; % For plotting purposes

xvec = 2*linspace(-lx0,lx0,29); % km 
yvec = 2*linspace(-ly0,ly0,31); % km
zvec = 2*linspace(-lz0,lz0,33); % km 
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
r2 = gseR2rel.data; r2(end) = 0; % km
r3 = gseR3rel.data; r3(end) = 0; % km
r4 = gseR4rel.data; r4(end) = 0; % km
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

%% Apply to MMS data, at z = 0, i.e. do not fit Epar, only Eperp at one parallel location
ieh = 10; % Define which EH to use

% Use symbolic expressions to get derivative of phi = E
syms x y z phi0 lx ly lz angle x0 y0 z0
R = [x y z];
% Triple Gaussian potential structure
phi = phi0*exp(-0.5*(((x-x0)*cosd(angle)-(y-y0)*sind(angle))/lx).^2-0.5*(((x-x0)*sind(angle)+(y-y0)*cosd(angle))/ly).^2-0.5*(z/lz-z0).^2);

% Electric field from potential
E = -gradient(phi,R);

% Order of the inputs defined by argument 'Vars'
mf_phi = matlabFunction(phi,'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % V
mf_Ex = matlabFunction(E(1),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)
mf_Ey = matlabFunction(E(2),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)
mf_Ez = matlabFunction(E(3),'Vars',[x y z phi0 lx ly lz angle x0 y0 z0]); % mV/m (V/km)

% We have lz from oservations, peak-to-peak length
% We will do the fit at z = 0 and z0 = 0 by time-shifting the timeseries
% with the dt corresponding to the parallel speed.
% These values listed below are the starting guessing values for the 
% optimization function fminsearch.
x0 = 0; % km
y0 = 0; % km
z0 = 0; % km
lx0 = 5; % km
ly0 = 5; % km
lz0 = 5; % km, from observations later
phi0 = 200; % V
angle = 0;

philev = 0:50:phi0; % For plotting purposes

%xvec = 2*linspace(-lx0,lx0,30); % km 
%yvec = 2*linspace(-ly0,ly0,31); % km
%zvec = 2*linspace(-lz0,lz0,32); % km 
%iz0 = numel(zvec)/2;

%[X,Y,Z] = meshgrid(xvec,yvec,zvec);
%PHI = mf_phi(X,Y,Z,phi0,lx0,ly0,lz0,angle,x0,y0,z0);
%EX = mf_Ex(X,Y,Z,phi0,lx0,ly0,lz0,angle,x0,y0,z0);
%EY = mf_Ey(X,Y,Z,phi0,lx0,ly0,lz0,angle,x0,y0,z0);
%EZ = mf_Ez(X,Y,Z,phi0,lx0,ly0,lz0,angle,x0,y0,z0);


% r1 = gseR1rel.data; r1(end) = 0; % km, put z = 0
% r2 = gseR2rel.data; r2(end) = 0; % km    ...
% r3 = gseR3rel.data; r3(end) = 0; % km    ...
% r4 = gseR4rel.data; r4(end) = 0; % km    ...

% Get electric field at center of ESW for each spacecraft, only perform
% fitting of perpendicular plane for now
tref = EpochTT(manual(ieh).t_ref);
%tref = tref +- 0.0001; % check if small offset of the center of the structures matters
dt = manual(ieh).dt;
tplus = EpochTT(manual(ieh).tplus);
tminus = EpochTT(manual(ieh).tminus);
T =  tplus - tminus;
%tint_eh = [tminus,tplus];
tint_eh = tref + abs([min(dt) max(dt)]) + 0.5*T*[-1 1];
tint_eh = tref + abs(min(dt))*[-1 1] + 0.5*T*[-1 1];
%tint_eh = [tminus tplus] + [min(dt) max(dt)] + 0.001*[-1 1];

gsev = manual(ieh).v;
gseB = manual(ieh).gseBav;
pitchangle = manual(ieh).pitchangle;
v = gsev*lmn';

lpp = manual(ieh).lpp;

% Time at the center of the ESWs for each spacecraft
c_eval('t? = tref + dt(?);',1:4)

% Coordinate system: perp1, perp2, par 
lmn = [manual(ieh).perp1; manual(ieh).perp2; manual(ieh).par]; 

ffilt = 30; % for highpass filtering
c_eval('gseE?filt = gseE?.filt(ffilt,0,[],5);',1:4) % use same time interval for all
%c_eval('gseE?filt = gseE?;',1:4) % use same time interval for all

% Try instead subtracting the lowpass value
%c_eval('gseE?filt = gseE? - gseE?.filt(0,ffilt,[],5);',1:4) % use same time interval for all

% Try instead subtracting the lowpass value
%c_eval('gseE?filt = gseE? - gseE?.resample(gseB?).resample(gseE?);',1:4) % use same time interval for all

% Electric field, for plotting
%c_eval('E? = gseE?filt.tlim(tint_eh)*lmn'';',1:4) % use same time interval for all

%c_eval('ex? = e?(1);',1:4) % perp1
%c_eval('ey? = e?(2);',1:4) % perp2
%c_eval('ez? = e?(3);',1:4) % par

% Electric field, for fitting 
c_eval('e? = gseE?filt.resample(t?).data*lmn'';',1:4) % use same time interval for all
c_eval('ex? = e?(1);',1:4) % perp1
c_eval('ey? = e?(2);',1:4) % perp2
c_eval('ez? = e?(3);',1:4) % par

% Spacecraft position
c_eval('r? = gseR?rel.data*lmn'';',1:4) % km
c_eval('r?(3) = 0;',1:4) % put z = 0, because we resampled to the ESW center time

params0 = double([phi0; lx0; ly0; lz0; angle; x0; y0; z0]);
% I need to create a cost function, which is to be evaluated. For example
% ediffx = mf_Ex(vars) - Ex_obs
% ediffy = mf_Ey(vars) - Ey_obs
% ediffz = mf_Ez(vars) - Ez_obs
% cost_function = sqrt(ediffx.^2 + ediffx.^2 + ediffz.^2)

xdata = [r1(1),r2(1),r3(1),r4(1)];
ydata = [r1(2),r2(2),r3(2),r4(2)];
zdata = [r1(3),r2(3),r3(3),r4(3)];
Ex_data = [e1(1),e2(1),e3(1),e4(1)];
Ey_data = [e1(2),e2(2),e3(2),e4(2)];
Ez_data = [e1(3),e2(3),e3(3),e4(3)];

cost_function = @(params) eh_costfunction(params,xdata,ydata,zdata,Ex_data,Ey_data,Ez_data,mf_Ex,mf_Ey,mf_Ez);
params = fminsearch(cost_function,params0);

philev = linspace(0,100*round(params(1)/100),11); % For plotting purposes

xvec = 2*linspace(-params(2),params(2),30); % km 
yvec = 2*linspace(-params(3),params(3),31); % km
zvec = 2*linspace(-params(4),params(4),32); % km, z should be 0
iz0 = numel(zvec)/2;
[X,Y,Z] = meshgrid(xvec,yvec,zvec);

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
colors = mms_colors('1234');

% nrows = 2;
% ncols = 1;
% h = setup_subplots(nrows,ncols);

figure(100)
nTimePanels = 3;
nrows = 2;
ncols = 1;

[h1,h2] = initialize_combined_plot(nTimePanels,nrows,ncols,0.5,'vertical');
isub = 1;

if 1 % timeseries of data
  hca = irf_panel('E perp1');
  set(hca,'ColorOrder',mms_colors('1234'))
  comp = 'x';
  icomp = 1;
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt)
  hca.YLabel.String = 'E_{\perp,1} (mV/m)';
  irf_legend(hca,sprintf('e_{perp,1} = [%4.2f,%4.2f,%4.2f]',lmn(icomp,1),lmn(icomp,2),lmn(icomp,3)),[0.02 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % timeseries of data
  hca = irf_panel('E perp2');
  set(hca,'ColorOrder',mms_colors('1234'))
  comp = 'y';
  icomp = 2;
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt)
  hca.YLabel.String = 'E_{\perp,2} (mV/m)';
  irf_legend(hca,sprintf('e_{perp,2} = [%4.2f,%4.2f,%4.2f]',lmn(icomp,1),lmn(icomp,2),lmn(icomp,3)),[0.02 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % timeseries of data
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  comp = 'z';
  icomp = 3;
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt)
  hca.YLabel.String = 'E_{||} (mV/m)';
  irf_legend(hca,sprintf('e_{||} = [%4.2f,%4.2f,%4.2f]',lmn(icomp,1),lmn(icomp,2),lmn(icomp,3)),[0.02 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

h1(1).Title.String = sprintf('id = %g, f highpass = %.0f',ieh,ffilt);
irf_pl_mark(h1,tref,[0 0 0])
irf_zoom(h1,'x',tint_eh)

if 0 % original data, only when using synthetic data
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
end
if 1 % fit
  hca = h2(isub); isub = isub + 1;
  contour(hca,X(:,:,iz0),Y(:,:,iz0),PHI_best(:,:,iz0))
  if not(abs(params(2)) > 30 || abs(params(3)) > 30), axis(hca,'equal'); end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  c_eval('quiver(hca,r?(1),r?(2),Ex?_best,Ey?_best,''linewidth'',1,''color'',colors(?,:))',1:4)
  c_eval('plot(hca,r?(1),r?(2),''o'',''color'',colors(?,:))',1:4)
  hold(hca,'off')
  hca.XLabel.String = 'perp 1 (km)';
  hca.YLabel.String = 'perp 2 (km)';
end
if 1 % information
  hca = h2(isub); isub = isub + 1;
  hca.Visible = 'off';
  irf_legend(hca,{...
    'esw parameters:',...
    sprintf('dt = [%4.1f,%4.1f,%4.1f,%4.1f] ms',dt(1)*1e3,dt(2)*1e3,dt(3)*1e3,dt(4)*1e3),...
    sprintf('v(perp1,perp2,par) = %.0f x [%5.2f,%5.2f,%5.2f] km/s',norm(v),v(1)/norm(v),v(2)/norm(v),v(3)/norm(v)),...
    sprintf('pitchangle = %.0f deg ',pitchangle),...
    sprintf('2l_{||} = L_{pp} = dt*v_{||} = [%4.1f,%4.1f,%4.1f,%4.1f] ms',abs(lpp(1)),abs(lpp(2)),abs(lpp(3)),abs(lpp(4))),...
    }',...    
    [0.0 0.999],'color',[0 0 0])
%   irf_legend(hca,{...
%     sprintf('phi0 = %.0f',params(1)),...
%     sprintf('lperp1 = %.1f km',params(2)),...
%     sprintf('lperp2 = %.1f km',params(3)),...
%     sprintf('angle = %.0f deg',params(4)),...
%     sprintf('x0 (perp1) = %.0f km',params(5)),...
%     sprintf('y0 (perp2) = %.0f km',params(6)),...
%     sprintf('z0 (par) = %.0f km',params(7))}',...    
%     [0.02 0.02],'color',[0 0 0])
  
  hleg = irf_legend(hca,{...
    sprintf('fitting parameters: \nphi0 = %.0f V',params(1)),...
    sprintf('lperp1 = %.1f km',params(2)),...
    sprintf('lperp2 = %.1f km',params(3)),...
    sprintf('angle = %.0f deg',params(4)),...
    sprintf('x0 (perp1) = %.0f km',params(5)),...
    sprintf('y0 (perp2) = %.0f km',params(6)),...
    sprintf('z0 (par) = %.0f km',params(7))}',...    
    [0.0 0.00],'color',[0 0 0]);
  hleg(1).VerticalAlignment = 'bottom';
end

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
fontsize = 12;
c_eval('h(?).FontSize = fontsize;',1:numel(h))
disp('Done.')

%% Apply to MMS data, at z = zvec, take into account epar also
for ieh = 12:numel(manual) % Define which EH to use

% Use symbolic expressions to get derivative of phi = E
syms x y z phi0 lx ly lz azim x0 y0 z0 polar r
syms phi(x,y,z,lx,ly,lz,x0,y0,z0,azim,polar,phi0)
syms E(x,y,z,lx,ly,lz,x0,y0,z0,azim,polar,phi0)
syms Ex(x,y,z,lx,ly,lz,x0,y0,z0,azim,polar,phi0)
syms Ey(x,y,z,lx,ly,lz,x0,y0,z0,azim,polar,phi0)
syms Ez(x,y,z,lx,ly,lz,x0,y0,z0,azim,polar,phi0)
inp = [x,y,z,lx,ly,lz,x0,y0,z0,azim,polar,phi0];
% [x,y,x,lx,ly,lz,azim,polar,phi0]

R = [x y z];
%Rp = [x y z];

% Triple Gaussian potential structure. The ESWs can be aligned and travel
% with and angle with respect to the magnetic field. This means we can have
% bipolar filed variations in the components perpendicular to B.

% Rotation matrix for azimuthal rotation 'angle' and polar 'polar'
% rotR = [sin(polar)*cos(azim) sin(polar)*sin(azim)  cos(polar);...
%         cos(polar)*cos(azim) cos(polar)*sin(azim) -sin(polar);...
%           -sin(azim)             cos(azim)               0];
% frotR = symfun(rotR,[azim,polar]);    
% rotR = [1 0 0; 0 1 0; 0 0 1];
rotR = [cos(polar)*cos(azim) cos(polar)*sin(azim) -sin(polar);...
          -sin(azim)             cos(azim)               0;...
        sin(polar)*cos(azim) sin(polar)*sin(azim)  cos(polar)]; 
frotR = symfun(rotR,[azim,polar]);
R0 = [x0 y0 z0];
L = [lx ly lz];
% Rp = rotR*((R-1*R0))';
%Rp = rotR*R';
Rp = rotR*(R-R0)';
%R = transpose(rotR)*rotR;
%Rp = transpose(rotR)*R';
% So in which coordinate system do I need to define the ESW properties?

% xp = (x-x0)*cosd(angle)*sind(polar)+(y-y0)*sind(angle)*sind(polar);
% yp = -(x-x0)*sind(angle)*sind(polar)+(y-y0)*cosd(angle)*sind(polar);
% zp = (z-z0)*cosd(polar);
% 
% % Think this rotation around z worked
% xp1 =  (x-x0)*cosd(angle)+(y-y0)*sind(angle);
% yp1 = -(x-x0)*sind(angle)+(y-y0)*cosd(angle);
% zp1 =  (z-z0);
% 
% % Rotate round x
% xp2 =  xp1;
% yp2 =  (yp1)*cosd(polar)+(zp1)*sind(polar);
% zp2 = -(yp1)*sind(polar)+(zp1)*cosd(polar);
% 
% xp = xp2;
% yp = yp2;
% zp = zp2;

xp = Rp(1);
yp = Rp(2);
zp = Rp(3);

mf_xp = matlabFunction(xp,'Vars',[x y z azim x0 y0 z0 polar]); % 
mf_yp = matlabFunction(yp,'Vars',[x y z azim x0 y0 z0 polar]); %
mf_zp = matlabFunction(zp,'Vars',[x y z azim x0 y0 z0 polar]); %

phi_var = [x,y,z,lx,ly,lz,azim,polar,phi0];
phi(x,y,z,lx,ly,lz,x0,y0,z0,azim,polar,phi0) = ...
  phi0*exp(...
  -0.5*(xp/lx).^2 ...
  -0.5*(yp/ly).^2 ...
  -0.5*(zp/lz).^2 ...
  );
fphi = symfun(phi,phi_var);


% Electric field from potential
Ex = -gradient(phi,x);
Ey = -gradient(phi,y);
Ez = -gradient(phi,z);

% Order of the inputs defined by argument 'Vars'
% inp = [x y z phi0 lx ly lz azim x0 y0 z0 polar];
mf_phi = matlabFunction(phi,'Vars',inp); % V
mf_Ex = matlabFunction(Ex,'Vars',inp); % mV/m (V/km)
mf_Ey = matlabFunction(Ey,'Vars',inp); % mV/m (V/km)
mf_Ez = matlabFunction(Ez,'Vars',inp); % mV/m (V/km)

%zz = linspace(-20,20,100);
%plotyy(zz,phi(0,0,zz,15,10,3,0,0,0,0,0,100),zz,Ez(0,0,zz,15,10,3,0,0,0,0,0,100))
%grid on
%
% Plot coordinate system.
% @(x,y,z,angle,x0,y0,z0,polar)
if 0
  % Define coordinate rotation
  polar_ = 10*(pi/180);
  azim_ = 0*(pi/180);
  lx_ = 15;
  ly_ = 10;
  lz_ = 5;
  x0_ = 0;
  y0_ = 0;
  z0_ = 0;
  x_ = 0;
  y_ = 0;
  phi0_ = 100;
  zvec_ = 5*linspace(-lz_,lz_,49);  
  xvec = 5*linspace(-lx_,lx_,50); % km 
  yvec = 5*linspace(-ly_,ly_,51); % km
  % define parallel length scale from electric field time series
  zvec = zvec_;
  %zvec = 2*linspace(-params(4),params(4),33); % km, z should be 0
  iz0 = ceil(numel(zvec)/2);
  [X,Y,Z] = meshgrid(xvec,yvec,zvec);
  
  hca = subplot(2,1,1);
  X0 = 0.*[1 1 1];
  Y0 = [0 0 0];
  Z0 = [0 0 0];
  Xunit = [1 0 0];
  Yunit = [0 1 0];
  Zunit = [0 0 1];  
  XP = mf_xp(Xunit,Yunit,Zunit,azim_,X0,Y0,Z0,polar_);
  YP = mf_yp(Xunit,Yunit,Zunit,azim_,X0,Y0,Z0,polar_);
  ZP = mf_zp(Xunit,Yunit,Zunit,azim_,X0,Y0,Z0,polar_);
  % (x,y,z,phi0,lx,ly,lz,azim,x0,y0,z0,polar)
  %PHI = mf_phi(X,Y,Z,lx_,ly_,lz_,phi0_,azim_,X0(1),Y0(1),Z0(1),theta_);
  %tic
  PHI = mf_phi(X,Y,Z,lx_,ly_,lz_,X0(1),Y0(1),Z0(1),azim_,polar_,phi0_);
  %toc
  quiver3(hca,X0,Y0,Z0,Xunit*qscale,Yunit*qscale,Zunit*qscale,'linewidth',2)
  hold(hca,'on')
  qscale = 20;
  quiver3(hca,X0,Y0,Z0,XP*qscale,YP*qscale,ZP*qscale,'linewidth',2)  
  hold(hca,'off')
  hold(hca,'on')
  hiso = patch(hca,isosurface(X,Y,Z,PHI,phi0_*0.5));
  hiso.FaceAlpha = 0.3;
  hiso.FaceColor = [0.5 0 0];
  %shading(hca,'flat')
  hold(hca,'off')
  hca.XLim = 1.1*[-1 1];
  hca.YLim = 1.1*[-1 1];
  hca.ZLim = 1.1*[-1 1];
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hca.ZLabel.String = 'z';
  legend(hca,{'R','RP'})
  axis(hca,'equal')
  
  hca = subplot(2,1,2);
  % mf_Ex = @(x,y,z,phi0,lx,ly,lz,angle,x0,y0,z0,polar) ...
  plot(hca,...
       zvec_,mf_Ex(x_,y_,zvec_,phi0_,lx_,ly_,lz_,azim_,x0_,y0_,z0_,polar_),...
       zvec_,mf_Ey(x_,y_,zvec_,phi0_,lx_,ly_,lz_,azim_,x0_,y0_,z0_,polar_),...
       zvec_,mf_Ez(x_,y_,zvec_,phi0_,lx_,ly_,lz_,azim_,x0_,y0_,z0_,polar_));
  legend(hca,{'x','y','z'})
  hca.YLim = [-10 10];
end
%

% We have lz from oservations, peak-to-peak length
% We will do the fit at z = 0 and z0 = 0 by time-shifting the timeseries
% with the dt corresponding to the parallel speed.
% These values listed below are the starting guessing values for the 
% optimization function fminsearch.
x0 = 0; % km
y0 = -5; % km
z0 = 0; % km 
lx0 = 15; % km
ly0 = 10; % km
lz0 = 3; % km, from observations later
phi0 = 300; % V
azim = 0;
polar = -13*pi/180;manual(ieh).pitchangle;

philev = 0:50:phi0; % For plotting purposes

% Get electric field at center of ESW for each spacecraft, only perform
% fitting of perpendicular plane for now
tref = EpochTT(manual(ieh).t_ref);
%tref = tref +- 0.0001; % check if small offset of the center of the structures matters
dt = manual(ieh).dt;
tplus = EpochTT(manual(ieh).tplus);
tminus = EpochTT(manual(ieh).tminus);
T =  tplus - tminus;

tint_eh = tref + abs([min(dt) max(dt)]) + 0.5*T*[-1 1];
tint_eh = tref + abs(min(dt))*[-1 1] + 0.3*T*[-1 1];

% The time fitting interval is different for each spacecraft, but should
% extend the same time if we align the reference times
c_eval('tint_fit_mms? = [tminus tplus] + dt(?) + 0.2*T*[-1 1];',1:4)

gsev = manual(ieh).v;
gseB = manual(ieh).gseBav;
pitchangle = manual(ieh).pitchangle;
v = gsev*lmn';
bvec =  (gseB/norm(gseB))*lmn'; % magnetic field unit vector in new csys
gse = [1 0 0; 0 1 0; 0 0 1];
gsevec = gse*lmn'; % gse uit vectros in new coordinate system
lpp = manual(ieh).lpp;

% Time at the center of the ESWs for each spacecraft
c_eval('t? = tref + dt(?);',1:4)

% Coordinate system: perp1, perp2, par 
lmn = [manual(ieh).perp1; manual(ieh).perp2; manual(ieh).par]; 

% make new lmn coordinates which are in the wave system
par = -manual(ieh).v/norm(manual(ieh).v);
perp1 = cross(par,cross([0 1 0],par)); perp1 = perp1/norm(perp1);
perp2 = cross(par,perp1);
lmn = [perp1; perp2; par];

ffilt = 100; % for highpass filtering
c_eval('gseE?filt = gseE?.filt(ffilt,0,[],5);',1:4) % use same time interval for all
%ffilt = 0;
%c_eval('gseE?filt = gseE?;',1:4) % use same time interval for all

% Try instead subtracting the lowpass value
%c_eval('gseE?filt = gseE? - gseE?.filt(0,ffilt,[],5);',1:4) % use same time interval for all

% Try instead subtracting the lowpass value
%c_eval('gseE?filt = gseE? - gseE?.resample(gseB?).resample(gseE?);',1:4) % use same time interval for all

% Electric field, for plotting
%c_eval('E? = gseE?filt.tlim(tint_eh)*lmn'';',1:4) % use same time interval for all
c_eval('E? = gseE?filt.resample(gseE1).tlim(tint_fit_mms?)*lmn'';',1:4) % use same time interval for all
% Due to the different time intervals, it might happen that one point more
% or less falls within the interval.
minlength = min([E1.length,E2.length,E3.length,E4.length]);
c_eval('if not(E?.length == minlength), E? = E?(1:minlength); end',1:4)
  
timeline = E1.time - tref;
zline = -timeline*v(3); % km/s

% Electric field, for fitting 
c_eval('e? = E?.data;',1:4) % use same time interval for all
c_eval('ex? = e?(:,1);',1:4) % perp1
c_eval('ey? = e?(:,2);',1:4) % perp2
c_eval('ez? = e?(:,3);',1:4) % par
nt = size(ex1,1);

% Spacecraft position
c_eval('r? = gseR?rel.data*lmn'';',1:4) % km
c_eval('r? = repmat(r?,nt,1);',1:4) % km
c_eval('r?(:,3) = zline;',1:4) % put z = zline, because we resampled to the ESW center time

params0 = double([lx0; ly0; lz0; x0; y0; z0; azim; polar; phi0]);
% I need to create a cost function, which is to be evaluated. For example
% ediffx = mf_Ex(vars) - Ex_obs
% ediffy = mf_Ey(vars) - Ey_obs
% ediffz = mf_Ez(vars) - Ez_obs
% cost_function = sqrt(ediffx.^2 + ediffx.^2 + ediffz.^2)

% these needs to be changes from 1 value per spacecraft, to nt values
% size(x_data) = [nt,4]
xdata = [r1(:,1),r2(:,1),r3(:,1),r4(:,1)];
ydata = [r1(:,2),r2(:,2),r3(:,2),r4(:,2)];
zdata = [r1(:,3),r2(:,3),r3(:,3),r4(:,3)];
Ex_data = [e1(:,1),e2(:,1),e3(:,1),e4(:,1)];
Ey_data = [e1(:,2),e2(:,2),e3(:,2),e4(:,2)];
Ez_data = [e1(:,3),e2(:,3),e3(:,3),e4(:,3)];

% Cost function needs to not only include 1 point for each spacecraft, but
% many. I.e. we fit an extended time interval, including the parallel 
% electric field. We can still keep the reference level to 0 for al
% spacecraft, i.e. displace the timelines.
cost_function = @(params) eh_costfunction_finite_tilt(params,xdata,ydata,zdata,Ex_data,Ey_data,Ez_data,mf_Ex,mf_Ey,mf_Ez);
params = fminsearch(cost_function,params0);

philev = linspace(0,100*round(params(1)/100),11); % For plotting purposes

xmax = params(1)*cos(params(7)) + params(2)*sin(params(8));
ymax = -params(1)*sin(params(7)) + params(2)*cos(params(8));
xvec = 3*linspace(-xmax,xmax,59); % km 
yvec = 3*linspace(-xmax,xmax,61); % km
% define parallel length scale from electric field time series
zvec = -timeline*vpar;
%zvec = 2*linspace(-params(4),params(4),33); % km, z should be 0
iz0 = ceil(numel(zvec)/2);
[X,Y,Z] = meshgrid(xvec,yvec,zvec);

mf_phi_best = @(x,y,z) mf_phi(x,y,z,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8),params(9));
mf_Ex_best = @(x,y,z) mf_Ex(x,y,z,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8),params(9));
mf_Ey_best = @(x,y,z) mf_Ey(x,y,z,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8),params(9));
mf_Ez_best = @(x,y,z) mf_Ez(x,y,z,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8),params(9));

PHI_best = mf_phi_best(X,Y,Z);
c_eval('EX?_best = mf_Ex_best(r?(1,1),r?(1,2),zvec);',1:4)
c_eval('EY?_best = mf_Ey_best(r?(1,1),r?(1,2),zvec);',1:4)
c_eval('EZ?_best = mf_Ez_best(r?(1,1),r?(1,2),zvec);',1:4)
% c_eval('Ex?_best = mf_Ex_best(r?(1),r?(2),r?(3));',1:4)
% c_eval('Ey?_best = mf_Ey_best(r?(1),r?(2),r?(3));',1:4)
% c_eval('Ez?_best = mf_Ez_best(r?(1),r?(2),r?(3));',1:4)

c_eval('E?fit = E?.clone(E?.time,[EX?_best,EY?_best,EZ?_best]);',1:4)
Ex_data_best = [EX1_best,EX2_best,EX3_best,EX4_best];
Ey_data_best = [EY1_best,EY2_best,EY3_best,EY4_best];
Ez_data_best = [EZ1_best,EZ2_best,EZ3_best,EZ4_best];



% tdata = timeline-mean(timeline);
% ydata = Bobs;
% yfit = @(t,B0,B1,dt) B0 + B1*tanh(t/dt);
% fun = @(x) sseval(x,tdata,ydata); % costfuction, root mean squared
% x0 = double([theta; 1]);
% %options = optimset
% bestx = fminsearch(fun,x0);
% Bopt = yfit(tdata,bestx(1),bestx(2),bestx(3));
%
colors = pic_colors('matlab');
colors = mms_colors('1234');
matlab_colors = pic_colors('matlab');
% nrows = 2;
% ncols = 1;
% h = setup_subplots(nrows,ncols);
%
figure(100)
nTimePanels = 3;
nrows = 2;
ncols = 1;

[h1,h2] = initialize_combined_plot(nTimePanels,nrows,ncols,0.6,'vertical');
isub = 1;

if 1 % timeseries of data
  hca = irf_panel('E perp1');
  set(hca,'ColorOrder',mms_colors('1234'))
  comp = 'x';
  icomp = 1;
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt)
  hold(hca,'on')
  hfit = irf_plot(hca,{E1fit.(comp),E2fit.(comp),E3fit.(comp),E4fit.(comp)},'comp','dt',dt,'--');
  hold(hca,'on')
  hca.YLabel.String = 'E_{\perp,1} (mV/m)';
  irf_legend(hca,sprintf('e_{perp,1} = [%4.2f,%4.2f,%4.2f]',lmn(icomp,1),lmn(icomp,2),lmn(icomp,3)),[0.02 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % timeseries of data
  hca = irf_panel('E perp2');
  set(hca,'ColorOrder',mms_colors('1234'))
  comp = 'y';
  icomp = 2;
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt)
  hold(hca,'on')
  hfit = irf_plot(hca,{E1fit.(comp),E2fit.(comp),E3fit.(comp),E4fit.(comp)},'comp','dt',dt,'--');
  hold(hca,'on')
  hca.YLabel.String = 'E_{\perp,2} (mV/m)';
  irf_legend(hca,sprintf('e_{perp,2} = [%4.2f,%4.2f,%4.2f]',lmn(icomp,1),lmn(icomp,2),lmn(icomp,3)),[0.02 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % timeseries of data
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  comp = 'z';
  icomp = 3;
  irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt)
  hold(hca,'on')
  hfit = irf_plot(hca,{E1fit.(comp),E2fit.(comp),E3fit.(comp),E4fit.(comp)},'comp','dt',dt,'--');
  hold(hca,'on')
  hca.YLabel.String = 'E_{||} (mV/m)';
  irf_legend(hca,sprintf('e_{||} = [%4.2f,%4.2f,%4.2f]',lmn(icomp,1),lmn(icomp,2),lmn(icomp,3)),[0.02 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
  irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'}',[0.98 0.02])
end


h1(1).Title.String = sprintf('id = %g, f highpass = %.0f',ieh,ffilt);
%h_tref = irf_pl_mark(h1,tref,matlab_colors(4,:));
%c_eval('h_tref(?).LineStyle = ''-.'';',1:numel(h_tref))
hpl_z0 = irf_pl_mark(h1,E1(iz0).time,matlab_colors(7,:));
c_eval('hpl_z0(?).LineStyle = ''-.'';',1:numel(hpl_z0))
irf_zoom(h1,'x',tint_eh)
legend([hpl_z0(1)],{'z = 0'})

if 0 % original data, only when using synthetic data
  hca = h(isub); isub = isub + 1;
  [Ccont,Hcont] = contour(hca,X(:,:,iz0),Y(:,:,iz0),PHI(:,:,iz0));
  clabel(Ccont,Hcont);
  axis(hca,'equal')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  quiver(hca,X,Y,EX,EY,'color',[0 0 0])
  c_eval('quiver(hca,r?(1),r?(2),ex?,ey?,''linewidth'',1,''color'',colors(?,:))',1:4)
  c_eval('plot(hca,r?(1),r?(2),''o'',''color'',colors(?,:))',1:4)
  hold(hca,'off')
end
if 1 % fit in perpendicular plane at z = 0
  hca = h2(isub); isub = isub + 1;
  contour(hca,X(:,:,iz0),Y(:,:,iz0),PHI_best(:,:,iz0))
  if not(abs(params(2)) > 30 || abs(params(3)) > 30), axis(hca,'equal'); end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  % Observed field
  c_eval('hq_data = quiver(hca,r?(iz0,1),r?(iz0,2),ex?(iz0,:),ey?(iz0,:),''linewidth'',2,''color'',colors(?,:));',1:4)
  % Fit field
  c_eval('hf_fit = quiver(hca,r?(iz0,1),r?(iz0,2),EX?_best(iz0,:),EY?_best(iz0,:),''linewidth'',1,''color'',colors(?,:));',1:4)
  c_eval('plot(hca,r?(iz0,1),r?(iz0,2),''o'',''color'',colors(?,:))',1:4)
  hold(hca,'off')  
  legend([hq_data,hf_fit],{'MMS','fit'},'location','best','box','on')
  hca.XLabel.String = 'x, perp 1 (km)';
  hca.YLabel.String = 'y, perp 2 (km)';
  hca.XLabel.String = 'e_{perp,1} (km)';
  hca.YLabel.String = 'e_{perp,2} (km)';
  hca.Title.String = 'z = 0';
end
if 0 % isosurface in 3D space, view 1
  hca = h2(isub); isub = isub + 1;
  surf_level = 0.5*params(end); % half of max
  hiso = patch(hca,isosurface(X,Y,Z,PHI_best,surf_level));
  hiso.FaceAlpha = 0.2;
  hiso.FaceColor = mms_colors('4');
  %if not(abs(params(2)) > 30 || abs(params(3)) > 30), axis(hca,'equal'); end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  qscale = params(3)*2;
  hq_b = quiver3(hca,0,0,0,qscale*bvec(1),qscale*bvec(2),qscale*bvec(3),0,'linewidth',1,'color',[0 0 0]);
  %c_eval('quiver(hca,r?(iz0,1),r?(iz0,2),EX?_best(iz0,:),EY?_best(iz0,:),''linewidth'',1,''color'',colors(?,:))',1:4)
  %c_eval('plot(hca,r?(iz0,1),r?(iz0,2),''o'',''color'',colors(?,:))',1:4)
  legend([hq_b,hiso],{'B','phi = phi_0/2'},'location','best','box','off')
  hold(hca,'off')
  hca.XLabel.String = 'x, perp 1 (km)';
  hca.YLabel.String = 'y, perp 2 (km)';
  hca.ZLabel.String = 'z, par (km)';
  axis(hca,'equal')
  view(hca,[0 1 0])
end
if 1 % information
  hca = h2(isub); isub = isub + 1;
  hca.Visible = 'off';
  irf_legend(hca,{...
    'esw parameters:',...
    sprintf('dt = [%4.1f,%4.1f,%4.1f,%4.1f] ms',dt(1)*1e3,dt(2)*1e3,dt(3)*1e3,dt(4)*1e3),...
    sprintf('v(perp1,perp2,par) = %.0f x [%5.2f,%5.2f,%5.2f] km/s',norm(v),v(1)/norm(v),v(2)/norm(v),v(3)/norm(v)),...
    sprintf('pitchangle = %.0f deg ',pitchangle),...
    sprintf('2l_{||} = L_{pp} = dt*v_{||} = [%4.1f,%4.1f,%4.1f,%4.1f] km',abs(lpp(1)),abs(lpp(2)),abs(lpp(3)),abs(lpp(4))),...
    }',...    
    [-0.0 0.999],'color',[0 0 0])%,'horizontalalignment','left')
%   irf_legend(hca,{...
%     sprintf('phi0 = %.0f',params(1)),...
%     sprintf('lperp1 = %.1f km',params(2)),...
%     sprintf('lperp2 = %.1f km',params(3)),...
%     sprintf('angle = %.0f deg',params(4)),...
%     sprintf('x0 (perp1) = %.0f km',params(5)),...
%     sprintf('y0 (perp2) = %.0f km',params(6)),...
%     sprintf('z0 (par) = %.0f km',params(7))}',...    
%     [0.02 0.02],'color',[0 0 0])
  
  hleg = irf_legend(hca,{...
    'fitting parameters:',...
    sprintf('lperp1 = %.1f km, lperp2 = %.1f km, lpar = %.1f km',params(1),params(2),params(3)),...        
    sprintf('x0 (perp1) = %.0f km, y0 (perp2) = %.0f km, z0 (perp2) = %.0f km',params(4),params(5),params(6)),...     
    sprintf('azim = %.0f deg',params(7)*180/pi),...
    sprintf('polar = %.0f deg',params(8)*180/pi),...
    sprintf('phi0 = %.0f V ',params(9))}',...  
    [-0.0 0.00],'color',[0 0 0])%,'horizontalalignment','left');
  %hleg(1).VerticalAlignment = 'bottom';
end
if 0 % isosurface in 3D space, view 2
  hca = h2(isub); isub = isub + 1;
  hiso = patch(hca,isosurface(X,Y,Z,PHI_best));
  hiso.FaceAlpha = 0.2;
  hiso.FaceColor = mms_colors('4');
  %if not(abs(params(2)) > 30 || abs(params(3)) > 30), axis(hca,'equal'); end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  qscale = params(3)*2;
  hq_b = quiver3(hca,0,0,0,qscale*bvec(1),qscale*bvec(2),qscale*bvec(3),0,'linewidth',1,'color',[0 0 0]);
  %c_eval('quiver(hca,r?(iz0,1),r?(iz0,2),EX?_best(iz0,:),EY?_best(iz0,:),''linewidth'',1,''color'',colors(?,:))',1:4)
  %c_eval('plot(hca,r?(iz0,1),r?(iz0,2),''o'',''color'',colors(?,:))',1:4)
  %legend([hq_b],{'B'},'location','best','box','off')
  hold(hca,'off')
  hca.XLabel.String = 'x, perp 1 (km)';
  hca.YLabel.String = 'y, perp 2 (km)';
  hca.ZLabel.String = 'z, par (km)';
  axis(hca,'equal')
  view(hca,[1 0 0])
end

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
fontsize = 12;
c_eval('h(?).FontSize = fontsize;',1:numel(h))
disp('Done.')
cn.print(sprintf('3d_fit_ieh=%g_ffilt=%g',ieh,ffilt))
end

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
for ieh = 13;[1 6:13];%1:numel(ehprop)
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
for ieh = 1:numel(manual)
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
