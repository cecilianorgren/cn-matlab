%% Load data to base particles on
% Spacecraft id
ic = 1;

% Time interval of event
tint_burst = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint_burst = tint_burst + [+5 -5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file

% Time intervals for modelling the distribution
tint_model = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); 
tint_model = irf.tint('2017-07-06T13:54:05.51Z/2017-07-06T13:54:05.63Z'); 

% Time interval for figure
tint_figure = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.630Z');
tint_figure = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.620Z');

% Set up database
localuser = datastore('local','user');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS'); 
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   

% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint_burst);',ic);
%c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
%c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint_burst+[20 0]));',ic)
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint_burst);',ic);

% Bin into cartesian grid
eint = [000 40000];
vint = [-Inf Inf];
vg = (-30:1:30)*1e3;

c_eval('eDist = ePDist?.tlim(tint_figure(1)++0.03+2*[-0.015 +0.015]);',ic)
           
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist(1)).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;
%lmn = [ePara.data; ePerp1.data; ePerp2.data];

%orient = [ePara.data; ePerp1.data; ePerp2.data];
orient = [ePerp1.data; ePerp2.data; ePara.data];
lowerelim = 40;
nMC = 100e0;

%ePDistCart = eDist.elim([40 inf]).rebin('cart',{vg,vg,vg},orient);

%% Electron test particles
% To see how the electrons move through the structure, I pass them through
% a defined potential.
% Define potential as a function
% To calculate the electric field from the potential, use symbolic
% expressions

units = irf_units;
m = units.me;
q = -units.e;

% From data fit
lx = 7.4e3;  % m
ly = 12.6e3; % m
lr = sqrt(lx.^2 + ly.^2); % m
lz = 4.3e3;  % m
phi0 = 237;  % V
B0 = 30e-9;  % T

% Average properties and assumed perpendicular lengths
lx = 12.5e3;  % m
ly = 12.5e3; % m
ly = 7.4e3;  % m
lx = 12.6e3; % m
lr = sqrt(lx.^2 + ly.^2); % m
lz = 3.5e3;  % m
phi0 = 237;  % V
B0 = 20e-9;  % T

% Get particle velocities from distribution
particles = eDist.macroparticles('ntot',1,'skipzero','scpot',scpot);
vx = [particles.vx particles.vy particles.vz]*orient(1,:)'*1e3; % m/s
vy = [particles.vx particles.vy particles.vz]*orient(2,:)'*1e3; % m/s
vz = [particles.vx particles.vy particles.vz]*orient(3,:)'*1e3; % m/s

% Move vz to frame of ESW
vph = -8500e3; % m/s
% Only keep the ones going in negative vz direction 
iInwards = (vz-vph) < 0;
vx = vx(iInwards);
vy = vy(iInwards);
vz = vz(iInwards);

% Properties of electrons
T0 = [100 200 400 700]; % eV, perpendicular energy
T0 = [200]; % eV, perpendicular energy
vt = sqrt(T0*units.eV*2/units.me); % m/s
Ek0 = [10]; % eV, parallel energy, defines U
vk = sqrt(Ek0*units.eV*2/units.me); % m/s

particles = eDist.macroparticles('ntot',1,'skipzero','scpot',scpot);

% for visibility, have different groups start in different positions
% no, just make different subplots instead
% rstart = [0.25 0.5 0.75]*lr; % in terms of lr
% thetastart = [0 90 180]; 
% Since the structure is not spherically symmetric, we try different
% starting points
x0 = [0 0.25 0.5 0.75 1.0 1.25]*lx; % m
y0 = [0 0.25 0.5 0.75 1.0 1.25]*ly; % m
x0 = [0.2]*lx; % m
y0 = [0]*ly; % m
z0 = 30e3; % m
%[VX0,VY0,VZ0,X0,Y0] = ndgrid(vt,vk,x0,y0); 
[VX0,X0,Y0] = ndgrid(vx,x0,y0); 
[VY0,X0,Y0] = ndgrid(vy,x0,y0); 
[VZ0,X0,Y0] = ndgrid(vz-vph,x0,y0); 

% Symbolic expression of potential structure and electric fields
syms x y z phi rho E
R = [x,y,z];

r = sqrt(x^2 + y^2); % r = r/rs
J0 = besselj(0,r); % Bessel function of zeroth order
absJ0 = abs(besselj(0,r)); % Bessel function of zeroth order
mf_absJ0 = matlabFunction(absJ0);
l00 = fminsearch(@(x) abs(mf_absJ0(x,0)),0);
J0 = besselj(0,l00*r/lr); % redefine bessel function with 'normalized' r

% Parallel potential profile, Gaussian
phi_par = phi0.*exp(-0.5*(z./lz).^2); 

% Perpendicular potentialprofile
phi_perp = J0; % Chen 2002
phi_perp = exp(-0.5*(x./lx).^2)*exp(-0.5*(y./ly).^2); % Gaussian

% Total potential profile
phi = phi_par*phi_perp;
sphi = symfun(phi,[x,y,z]);
mf_phi = matlabFunction(phi);

% Electric field
E = -gradient(phi,R);
mf_Ex = matlabFunction(E(1));
mf_Ey = matlabFunction(E(2));
mf_Ez = matlabFunction(E(3));

% Charge density
rho = divergence(E,R)*units.mu0;
mf_rho = matlabFunction(rho,'vars',R);

% Vector potential
A = [0, B0*x, 0];
mf_Ay = matlabFunction(A(2),'vars',R); 

% Magnetic field
B = curl(A,R);
mf_Bx = matlabFunction(B(1),'vars',R);
mf_By = matlabFunction(B(2),'vars',R);
mf_Bz = matlabFunction(B(3),'vars',R);

% Set up grid
xmax = 3*lx;
ymax = 3*xmax;
zmax = 3*max([xmax ymax]);
xvec = linspace(-xmax,xmax,11);
yvec = linspace(-ymax,ymax,11);
zvec = linspace(-zmax,zmax,21);
[X,Y,Z] = meshgrid(xvec,yvec,zvec);
%xvec2 = linspace(-xmax,xmax,101);

% Prepare for integration
inpEx = @(x,y,z) mf_Ex(x,y,z);
inpEy = @(x,y,z) mf_Ey(x,y,z);
inpEz = @(x,y,z) mf_Ez(x,y,z);
inpBx = @(x,y,z) mf_Bx(x,y,z);
inpBy = @(x,y,z) mf_By(x,y,z);
inpBz = @(x,y,z) mf_Bz(x,y,z);

EoM = @(t,xyz) eom(t,xyz,inpEx,inpEy,inpEz,inpBx,inpBy,inpBz,m,q); 
boxedge = [-31 31]*1e3; % terminate integration when electron exits this region in z.
options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@(tt,xx) myEventBoxEdge(tt,xx,boxedge));
%options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,obj.xi(([1 end]))));
%options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));
%options = odeset('RelTol',1e-6);

% Initialize figure for plotting
doPlot = 1;
if doPlot
  colors = pic_colors('matlab');
  nrows = 2;
  ncols = 2;
  npanels = nrows*ncols;
  h = setup_subplots(nrows,ncols);
  c_eval('plot(h(?),nan,nan);',1:npanels)
  c_eval('hold(h(?),''on'');',1:npanels)
  c_eval('h(?).XGrid = ''on'';',1:npanels)
  c_eval('h(?).YGrid = ''on'';',1:npanels)
end

% Run through particles
S = struct([]);
for ie = 1:numel(X0)
  % Initial conditions
  % Start gyromotion radially outwards, such that R0 defines the center of
  % the gyromotion
  %az0 = atan2(Y0(ie),X0(ie));
  %vx0 = VT(ie)*cos(AZ0);
  %vy0 = VT(ie)*sin(AZ0);
  x_init = [X0(ie),Y0(ie),z0,VX0(ie),VY0(ie),VZ0(ie)]; % m, m/s
  x_init;
  T = [0 0.2];

  % Perform integration
  [t,x_sol] = ode45(EoM,T,x_init,options);
  
  % Calculate gyrocenter of particle
  r_sol = sqrt(x_sol(:,1).^2 + x_sol(:,2).^2);
  [PKS,LOCS]= findpeaks(r_sol); 
  if numel(LOCS)>1
    mean_peak_distance = diff(LOCS);
    r_center = movmean(r_sol,fix(mean(mean_peak_distance)));
  else
    r_center = r_sol*NaN;    
  end

  % Calculate total energy U = Ek + Ep
  Ep = -units.e*mf_phi(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  Ek = 0.5*units.me*(x_sol(:,4).^2 + x_sol(:,5).^2 + x_sol(:,6).^2);
  Eky = 0.5*units.me*(x_sol(:,5).^2);
  Ekpar = 0.5*units.me*(x_sol(:,6).^2);

  % Extract vector potential
  Ay = mf_Ay(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  
  % Collect output
%  S(ie).vt0 = VT(ie);
%  S(ie).T0 = units.me*VT(ie)^2/2/units.eV; % sqrt(t*units.eV*2/units.me)
  S(ie).x = x_sol(:,1);
  S(ie).y = x_sol(:,2);
  S(ie).z = x_sol(:,3);
  S(ie).z = x_sol(:,3);
  S(ie).vx = x_sol(:,4);
  S(ie).vy = x_sol(:,5);
  S(ie).vz = x_sol(:,6);
  S(ie).vperp = sqrt(x_sol(:,4).^2 + x_sol(:,5).^2);
  S(ie).vpar = x_sol(:,6);
  S(ie).pitchangle = atan2d(S(ie).vperp,S(ie).vpar);    
  S(ie).r = r_sol;
  S(ie).r_gc = r_center;
  S(ie).Ekperp = units.me*S(ie).vperp.^2/2/units.eV;
  S(ie).Ekpar = units.me*S(ie).vpar.^2/2/units.eV;
  S(ie).Ep = Ep;
  S(ie).Ek = Ek;
  S(ie).U = Ek + Ep;
  S(ie).Ex = mf_Ex(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  S(ie).Ey = mf_Ey(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  S(ie).Ez = mf_Ez(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  
  S(ie).tstart = t(1);
  S(ie).tstop = t(end);  
  S(ie).zstart = x_sol(1,3);
  S(ie).zstop = x_sol(end,3);  
  S(ie).Ekpar_start = S(ie).Ekpar(1);
  S(ie).Ekpar_stop = S(ie).Ekpar(end);
  S(ie).Ekperp_start = S(ie).Ekperp(1);
  S(ie).Ekperp_stop = S(ie).Ekperp(end);
  S(ie).pitchangle_start = S(ie).pitchangle(1);
  S(ie).pitchangle_stop = S(ie).pitchangle(end);
  
  nRef = numel(find([S.zstop]>0));
  
  if doPlot %&& S(ie).z(end)>20e3 %S(ie).z(end)<-20e3
    
    [it0,iv0,ix0,iy0] = ind2sub(size(X0),ie);
    % Work done on particles
    Exdx = cumtrapz(S(ie).x,S(ie).Ex);
    Eydy = cumtrapz(S(ie).y,S(ie).Ey);
    Ezdz = cumtrapz(S(ie).z,S(ie).Ez);
    
    phi1D = mf_phi(S(ie).x(1),S(ie).y(1),S(ie).z);
    
    isub = 1;
    
    if 1 % (z,vz)
      hca = h(isub); isub = isub + 1;   
      plot(hca,x_sol(:,3)*1e-3,x_sol(:,6)*1e-3*1e-3)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'v_z-v_{ph} (10^3 km/s)';
      hca.Title.String = sprintf('reflected = %g/%g',nRef,ie);
      hca.XLim = [-30 30];
      hca.YLim = [-25 10];
    end
    
    if 0 % trajectory in x,y plane
      hca = h(isub); isub = isub + 1;   
      plot(hca,x_sol(:,1)*1e-3,x_sol(:,2)*1e-3)
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'y (km)';
      hca.XLim = [-35 35];
      hca.YLim = [-35 35];
    end
    if 1 % trajectory in x,y,z plane
      hca = h(isub); isub = isub + 1;   
      plot3(hca,x_sol(:,1)*1e-3,x_sol(:,2)*1e-3,x_sol(:,3)*1e-3)
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'y (km)';
      hca.ZLabel.String = 'z (km)';
    end
    if 0 % zstart,zstop
      hca = h(isub); isub = isub + 1;   
      plot(hca,x_sol(1,3)*1e-3,x_sol(end,3)*1e-3,'.')
      hca.XLabel.String = 'z_{start} (km)';
      hca.YLabel.String = 'z_{stop} (km)';
    end
      
    
    if 0 % dEkpar(T0,Ekpar0)
      hca = h(isub); isub = isub + 1;   
      scatter(hca,S(ie).T0,S(ie).Ekpar(1),S(ie).Ekperp(end),S(ie).Ekperp(end)-S(ie).Ekperp(1))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = '\Delta E_{k,\perp}';
      hca.XLabel.String = 'T_0 (eV)';
      hca.YLabel.String = 'E_{kz,0} (eV)';  
    end
    if 0 % dvpar(T0,Ekpar0)
      hca = h(isub); isub = isub + 1;   
      scatter(hca,S(ie).T0,S(ie).Ekpar(1),abs(abs(S(ie).vpar(end))-abs(S(ie).vpar(1)))*1e-3,abs(abs(S(ie).vpar(end))-abs(S(ie).vpar(1)))*1e-3)
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = '\Delta |v_{par}| (km/s)';
      hca.XLabel.String = 'T_0 (eV)';
      hca.YLabel.String = 'E_{kz,0} (eV)';  
    end
    if 0 % dvpar(T0,Ekpar0)
      hca = h(isub); isub = isub + 1;   
      scatter(hca,S(ie).T0,S(ie).Ekpar(1),abs(abs(S(ie).vpar(end))-abs(S(ie).vpar(1)))*1e-3,S(ie).r(1)*1e-3)
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = 'r_{start} (km)';
      hca.XLabel.String = 'T_0 (eV)';
      hca.YLabel.String = 'E_{kz,0} (eV)';  
    end
    
    if 0 % vz(start,stop)
      hca = h(isub); isub = isub + 1;   
      plot(hca,x_sol(1,6)*1e-3,x_sol(end,6)*1e-3,'.')
      hca.XLabel.String = 'v_z(t_{start}) (km/s)';
      hca.YLabel.String = 'v_z(t_{stop}) (km/s)';
    end
    if 0 % vz_stop -vz_start vs T0
      hca = h(isub); isub = isub + 1;   
      plot(hca,S(ie).T0,(x_sol(1,6)-x_sol(end,6))*1e-3,'.')
      hca.XLabel.String = 'T_0 (eV)';
      hca.YLabel.String = 'v_z(t_{stop})-v_z(t_{start}) (km/s)';
    end
    if 0 % vperp_stop - vperp_start vs T0
      hca = h(isub); isub = isub + 1;   
      plot(hca,S(ie).T0,S(ie).vperp(end)-S(ie).vperp(1),'.','color',colors(ix0,:))
      hca.XLabel.String = 'T_0 (eV)';
      hca.YLabel.String = 'v_{\perp}(t_{stop})-v_{\perp}(t_{start}) (km/s)';
    end
    
    if 1 %  pitch angle
      hca = h(isub); isub = isub + 1;         
      plot(hca,S(ie).z*1e-3,abs(S(ie).pitchangle))
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = '|\theta| (deg)';
      hca.XLim = [-30 30];
      hca.YLim = [0 180];
    end
    if 1 %  Ezdz vs Exdx+Eydy
      hca = h(isub); isub = isub + 1;            
      plot(hca,Exdx+Eydy,Ezdz)
      hca.XLabel.String = 'Exdx+Eydy (V)';
      hca.YLabel.String = 'Ezdz (V)';
      axis(hca,'equal')
    end
    if 0 %  Ezdz
      hca = h(isub); isub = isub + 1;            
      plot(hca,S(ie).z,Ezdz)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'Ezdz (V)';
    end
    if 0 %  Exdx
      hca = h(isub); isub = isub + 1;            
      plot(hca,S(ie).z*1e-3,Exdx+Eydy)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'Exdx + Eydy (V)';
    end
    if 0 %  Ezdz-phi(z) , deviation from 1-D trajectory
      hca = h(isub); isub = isub + 1;            
      plot(hca,S(ie).z*1e-3,Ezdz+phi1D)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'Ezdz-phi(x0,y0,z) (V)';
      hca.XLim = [-30 30];
      hca.YLim = [-250 250];
    end
    
    if 0 % delta pitch angle
      hca = h(isub); isub = isub + 1;   
      plot(hca,S(ie).T0,S(ie).pitchangle(end)-S(ie).pitchangle(1),'.')
      hca.XLabel.String = 'T_0 (km)';
      hca.YLabel.String = '\Delta \theta (deg)';
    end
    
    if 0 % delta Ek, zero, so energy is conserved
      hca = h(isub); isub = isub + 1;   
      plot(hca,S(ie).T0,(Ek(end)-Ek(1))/units.eV)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = '\Delta E_k (eV)';
    end
    
    if 0 % delta Ek, zero, so energy is conserved
      hca = h(isub); isub = isub + 1;   
      plot(hca,S(ie).T0,(Ek(end)-Ek(1))/units.eV)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = '\Delta E_k (eV)';
    end
    if 0 % Ek
      hca = h(isub); isub = isub + 1;   
      plot(hca,x_sol(:,3),Ek)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'ek (..)';
    end
    if 0 % Ep
      hca = h(isub); isub = isub + 1;   
      plot(hca,x_sol(:,3),Ep)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'U (..)';
    end
    if 0 % U
      hca = h(isub); isub = isub + 1;   
      plot(hca,x_sol(:,3),Ek+Ep)
      hca.XLabel.String = 'z (km)';
      hca.YLabel.String = 'Ep (..)';    
    end
    if 0 % gyrocenter
      hca = h(isub); isub = isub + 1;
      plot(hca,x_sol(:,3)*1e-3,(r_center-r_sol(1))*1e-3)
      hca.XLabel.String = 'z (km)';
      %hca.YLabel.String = 'r_{gyrocenter}/(l_x^2 + l_y^2)^{-1/2}';
      hca.YLabel.String = '\Delta r_{gyrocenter} (km)';
      hca.YLabel.Interpreter = 'tex';
    end
    
    drawnow
    pause(0.01)
  end
end
if doPlot
  c_eval('hold(h(?),''off'');',1:npanels)
end

if doPlot*0
  %%
  figure(17)
  nrows = 2;
  ncols = 1;
  npanels = nrows*ncols;
  h = setup_subplots(nrows,ncols);
  isub = 1;
  
  % [VT,VK,X0,Y0] = ndgrid(vt,vk,x0,y0); 
  datasize = size(X0);
  x0_edges = [x0(1)-1 x0 + 1];
  y0_edges = [y0(1)-1 y0 + 1];
  vt0_edges = [vt(1)-1 vt + 1];
  vk0_edges = [vk(1)-1 vk + 1];
  if 1 % histograms for different starting parameters
    hca = h(isub); isub = isub + 1;   
    d_pitchangle = [S.pitchangle_abs_start] - [S.pitchangle_abs_stop];
    d_pitchangle = reshape(d_pitchangle,size(X0));
    %N = hist(d_pitchangle);
        
    data = permute(d_pitchangle,[1 2 3 4]); % move the dependent variable to first position    
    data = reshape(data,[datasize(1), prod(datasize(2:4))]); % reshape data into 2D matrix of size (N_depvar,prod(N_restvar))
    mean_T0 = mean(data,2); % average over remaining variables
    
    data = permute(d_pitchangle,[2 1 3 4]); % move the dependent variable to first position    
    data = reshape(data,[datasize(2), prod(datasize([1 3 4]))]); % reshape data into 2D matrix of size (N_depvar,prod(N_restvar))
    mean_Ek0 = mean(data,2); % average over remaining variables
    
    data = permute(d_pitchangle,[3 2 1 4]); % move the dependent variable to first position    
    data = reshape(data,[datasize(3), prod(datasize([1 2 4]))]); % reshape data into 2D matrix of size (N_depvar,prod(N_restvar))
    mean_x0 = mean(data,2); % average over remaining variables
    
    data = permute(d_pitchangle,[4 2 3 1]); % move the dependent variable to first position    
    data = reshape(data,[datasize(4), prod(datasize(1:3))]); % reshape data into 2D matrix of size (N_depvar,prod(N_restvar))
    mean_y0 = mean(data,2); % average over remaining variables
    plot(hca,1:4,mean_T0,1:4,mean_Ek0,1:4,mean_x0,1:4,mean_y0)

  end
  
end

function [value, isterminal, direction] = myEventBoxEdge(t, x_vect,boxedge)
  % integration is terminated when value changes sign
  % for this setup, value is initially negative
  value      = (x_vect(3)-boxedge(1))*(x_vect(3)-boxedge(2));        
  isterminal = 1;   % Stop the integration
  direction  = 0;
end
function  x_res = eom(t,x_vect,inpEx,inpEy,inpEz,inpBx,inpBy,inpBz,m,q)
  x_ = x_vect(1);
  y_ = x_vect(2);
  z_ = x_vect(3);
  vx_ = x_vect(4);
  vy_ = x_vect(5);
  vz_ = x_vect(6);      

  Ex = inpEx(x_,y_,z_);
  Ey = inpEy(x_,y_,z_);
  Ez = inpEz(x_,y_,z_);   
  Bx = inpBx(x_,y_,z_);
  By = inpBy(x_,y_,z_);
  Bz = inpBz(x_,y_,z_); 
  
%   it = size(x_sol_all,1);
%   x_sol_all(it+1,1) = Ex;
%   x_sol_all(it+1,2) = Ey;
%   x_sol_all(it+1,3) = Ez;
%   x_sol_all(it+1,4) = 0;
%   x_sol_all(it+1,5) = 0;
%   x_sol_all(it+1,6) = Bz;
%   x_sol_all(it+1,7) = t;

  % Equations to be solved
  x_res = zeros(6,1);
  x_res(1) = vx_; % dx/dt = vx;
  x_res(2) = vy_; % dy/dt = vy;
  x_res(3) = vz_; % dz/dt = vz;
  x_res(4) = (q/m)*(Ex + vy_*Bz - vz_*By);
  x_res(5) = (q/m)*(Ey + vz_*Bx - vx_*Bz);
  x_res(6) = (q/m)*(Ez + vx_*By - vy_*Bx);                                              

end      
      