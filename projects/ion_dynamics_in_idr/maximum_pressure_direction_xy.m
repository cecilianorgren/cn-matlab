function varargout = maximum_pressure_direction_xy(tsP)


nT = tsP.length;
R_all = zeros(nT,3,3);
P_all_rot = zeros(nT,3,3);

P_all = tsP.data;

for it = 1:nT  
  P = squeeze(P_all(it,:,:));
 
  % Rotate tensor, to find the most unequal component
  % Angle 
  theta = 0.5*atan((2*P(1,2))./abs(P(2,2)-P(1,1)));

  R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
  P2 = R*(P*transpose(R));

  % Without this, Pxy is always zero. With it, sometimes it's not
  if P2(1,1) > P2(2,2)
    theta = pi-theta;
    R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    P2 = R*(P*transpose(R));
  end
  R_all(it,:,:) = R;
  P_all_rot(it,:,:) = P2;
end

varargout{1} = R_all;
varargout{2} = P_all_rot;

tsR = irf.ts_tensor_xyz(tsP.time,R_all);
tsP = irf.ts_tensor_xyz(tsP.time,P_all_rot);
tsR1 = irf.ts_vec_xyz(tsP.time,squeeze(R_all(:,1,:)));
tsR2 = irf.ts_vec_xyz(tsP.time,squeeze(R_all(:,2,:)));
tsR3 = irf.ts_vec_xyz(tsP.time,squeeze(R_all(:,3,:)));


varargout{1} = tsP;
varargout{2} = tsR;
varargout{3} = tsR1;
varargout{4} = tsR2;
varargout{5} = tsR3;
