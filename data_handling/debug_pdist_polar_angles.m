
dpol = 11.25;
pol_center = dpol/2:dpol:180;
pol_edges = 0:dpol:180;

dv_pol1 = sind(pol_center)*dpol*pi/180;
dv_pol2 = cosd(pol_edges(1:(end-1))) - cosd(pol_edges(2:end));

%% Makde daniels model PDist with fake data

PDist = iPDist3(1);
time = PDist.time;
Bxyz = irf.ts_vec_xyz(time,[1 0.01 0]);
SCpot = irf.ts_scalar(time,1);
ne = irf.ts_scalar(time,1);
Ve = irf.ts_vec_xyz(time,[0 0 0.1]);
T = irf.ts_tensor_xyz(time,1e3*eye(3));


modelPDist = mms.make_model_dist(PDist,Bxyz,SCpot,ne,Ve,T);

%modelPDist.data(:,:,:,1) = 0;

%c_eval('Tdsl = [tsLdsl?.resample(time).data; tsMdsl?.resample(time).data; tsNdsl?.resample(time).data];',ic)
%Ldsl = Tdsl(1,:);
%Mdsl = Tdsl(2,:);
%Ndsl = Tdsl(3,:);


Mdsl = [0 1 0]
Ndsl = [0 0 1];

Mdsl = [0 0 1]
Ndsl = [0 1 0];


vdf = modelPDist.reduce('2D',Mdsl,Ndsl);
dv = vdf.depend{1}(2) - vdf.depend{1}(1);
n = sum(vdf.data(:))*dv*1e3*dv*1e3*1e-6;

hca = subplot(1,1,1);
vdf.plot_plane(hca)
irf_legend(hca,sprintf('sum(f*dv*dv) = %g cc',n),[0.02 0.02])
axis equal
colormap(pic_colors('candy6'))

%% Current method
nMC = 10;

dpol = 11.25;
pol = dpol/5:dpol:3*dpol;
npol = numel(pol);
pol_edges = 0:dpol:3*dpol;

dazim = 11.25;
azim = dazim/2:dazim:(360-dazim/2);
azim_edges = 0:dazim:360;
naz = numel(azim);


POL = zeros(naz,npol,nMC);
AZIM = zeros(naz,npol,nMC);

for iaz = 1:naz
  az1 = azim_edges(iaz);
  az2 = azim_edges(iaz+1);
  for ipol = 1:npol    
    pol1 = pol_edges(ipol);
    pol2 = pol_edges(ipol+1);    
    POL(iaz,ipol,:) = pol1 + dpol*rand(1,nMC);
    AZIM(iaz,ipol,:) = az1 + dazim*rand(1,nMC);
  end
end


VX = cosd(AZIM).*sind(POL);
VY = sind(AZIM).*sind(POL);

scatter(VX(:),VY(:),'.')

%% Better method
nMC = 10;

dpol = 11.25;
pol = dpol/5:dpol:3*dpol;
npol = numel(pol);
pol_edges = 0:dpol:3*dpol;

dazim = 11.25;
azim = dazim/2:dazim:(360-dazim/2);
azim_edges = 0:dazim:360;
naz = numel(azim);


COS_POL = zeros(naz,npol,nMC);
AZIM = zeros(naz,npol,nMC);

for iaz = 1:naz
  az1 = azim_edges(iaz);
  az2 = azim_edges(iaz+1);
  for ipol = 1:npol    
    pol1 = pol_edges(ipol);
    pol2 = pol_edges(ipol+1);
    cos_pol1 = cosd(pol1);
    cos_pol2 = cosd(pol2);
    dcos_pol = cos_pol2 - cos_pol1;
    COS_POL(iaz,ipol,:) = cos_pol1 + dcos_pol*rand(1,nMC);
    AZIM(iaz,ipol,:) = az1 + dazim*rand(1,nMC);
  end
end
POL = acosd(COS_POL);


VX = cosd(AZIM).*sind(POL);
VY = sind(AZIM).*sind(POL);

scatter(VX(:),VY(:),'.')















