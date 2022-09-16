function out = irf_4_v_continuous(r1,r2,r3,r4,b1,b2,b3,b4)
%IRF_4_V_CONTINUOUS Calculate velocity of timeseries

% https://github.com/spedas/bleeding_edge/blob/4986e65fb9b09522326877902665e69ab548b07d/spedas_gui/panels/mdd_std_for_gui.pro

% Resample position (r) to data (d) timeseries
timeline = b1.time;
t = timeline - timeline(1);
dt = diff(t);
c_eval('r? = r?.resample(timeline).data;',1:4)
c_eval('babs? = b?.resample(timeline).abs.data; b? = b?.resample(timeline).data;',1:4)

% dA/dt = -v*grad(A)

% Initialize matrices
nt = timeline.length;
all_A = zeros(nt,3,3);
all_B = zeros(nt,1);
all_C = zeros(nt,3);
all_EXP = zeros(nt,2);
lambda1 = zeros(nt,1);
lambda1 = zeros(nt,1);
lambda2 = zeros(nt,1);
lambda3 = zeros(nt,1);
dimention = zeros(nt,3);
egvt1 = zeros(nt,3);
egvt2 = zeros(nt,3);
egvt3 = zeros(nt,3);

% Loop through time

for it = 1:nt
  % dR = [r2-r1, r3-r1, r4-r1]
  dR = [r2(it,:);r3(it,:);r4(it,:)]-[r1(it,:);r1(it,:);r1(it,:)];
  % dB = [b2-b1, b3-b1, b4-b1]
  Bx = [b1(1), b2(1), b3(1), b4(1)];
  By = [b1(2), b2(2), b3(2), b4(2)];
  Bz = [b1(3), b2(3), b3(3), b4(3)];
  
  Ar = [ dR(1,:); dR(2,:); dR(3,:)];
  
  dB = [b2(it,:);b3(it,:);b4(it,:)]-[b1(it,:);b1(it,:);b1(it,:)];
 
  dBx = [Bx(:,2)-Bx(1), Bx(:,3)-Bx(1), Bx(:,4)-Bx(1)];
  dBy = [By(:,2)-By(1), By(:,3)-By(1), By(:,4)-By(1)];
  dBz = [Bz(:,2)-Bz(1), Bz(:,3)-Bz(1), Bz(:,4)-Bz(1)];
  
  PBx = dBx*mpower(Ar,-1);
  PBy = dBy*mpower(Ar,-1);
  PBz = dBz*mpower(Ar,-1);
  % dB/dr = [dBx/dr; dBy/dr; dBz/dr]
  %       = [dBx/dx, dBx/dy, dBx/dz; 
  %          dBy/dx, dBy/dy, dBy/dz; 
  %          dBz/dx, dBz/dy, dBz/dz]
  %PB = dB/dR; % Matrix division
  %A = PB';
  
  A = [PBx(1), PBy(1), PBz(1); PBx(2), PBy(2), PBz(2); PBx(3), PBy(3), PBz(3)];
  
  all_A(it,:,:) = A;
  all_B(it) = PBx(1) + PBy(1) + PBz(3);  
  all_C(it,:) = [PBz(2)-PBy(3), PBx(3)-PBz(1), PBy(1)-PBx(2)];
  all_EXP(it,1) = abs(all_B(it))/(sqrt(all_C(it,1:3))*transpose(all_C(it,1:3)));
  all_EXP(it,2) = abs(all_B(it))/max(max(abs(all_A(it,1:3,1:3))));
  
end

for it = 1:nt
  A = squeeze(all_A(it,:,:));
  A1 = A*transpose(A);
  [vv,DD,W] = eig(A1);
end



PB = dB\dR';
dAdt = interp1(t(2:end)-0.5*dt,diff(aA,1)./[dt dt dt],t);

m = dR\dt'; % "1/v vector"

V = m/norm(m)^2;
out = V';
    
end

function dR = get_vol_ten(r1,r2,r3,r4)
% Function that returns volumetric tensor with sc1 as center. t can be
% vector or scalar.

%if length(t) == 1
%  t = [t,t,t,t];
%end

% Interpolate position
%r1=[t(1) interp1(r1(:,1),r1(:,[2 3 4]),t(1),'spline','extrap')];
%r2=[t(2) interp1(r2(:,1),r2(:,[2 3 4]),t(2),'spline','extrap')];
%r3=[t(3) interp1(r3(:,1),r3(:,[2 3 4]),t(3),'spline','extrap')];
%r4=[t(4) interp1(r4(:,1),r4(:,[2 3 4]),t(4),'spline','extrap')];

%r1 = r1(2:4);
%r2 = r2(2:4);
%r3 = r3(2:4);
%r4 = r4(2:4);

% Volumetric tensor with SC1 as center.
dR = [r2;r3;r4]-[r1;r1;r1]; %mabye good?


end