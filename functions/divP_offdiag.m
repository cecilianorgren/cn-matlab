function out = divP(R1,R2,R3,R4,P1,P2,P3,P4)

% R are TSeries
% P are 3x3 cell array of TSeries pressure tensos components

% Extraxt pressure tensor components
if isa(P1,'cell')
  c_eval('PeXX? = P?{1,1}; PeXY? = P?{1,2}; PeXZ? = P?{1,3}; PeYY? = P?{2,2}; PeYZ? = P?{2,3}; PeZZ? = P?{3,3};')
elseif isa(P1,'TSeries')
  c_eval('PeXX? = TSeries(P?.time,P?.data(:,1,1));')
  c_eval('PeXY? = TSeries(P?.time,P?.data(:,1,2));')
  c_eval('PeXZ? = TSeries(P?.time,P?.data(:,1,3));')
  c_eval('PeYY? = TSeries(P?.time,P?.data(:,2,2));')
  c_eval('PeYZ? = TSeries(P?.time,P?.data(:,2,3));')
  c_eval('PeZZ? = TSeries(P?.time,P?.data(:,3,3));')  
end

ts = PeXX1.time;
c_eval('R?=R?.resample(ts);')
% Calculate pressure divergence
EPeXX = c_4_grad('R?','PeXX?','grad');
EPeXY = c_4_grad('R?','PeXY?','grad');
EPeXZ = c_4_grad('R?','PeXZ?','grad');
EPeYY = c_4_grad('R?','PeYY?','grad');
EPeYZ = c_4_grad('R?','PeYZ?','grad');
EPeZZ = c_4_grad('R?','PeZZ?','grad');
EPeX = (EPeXX.data(:,1)+EPeXY.data(:,2)+EPeXZ.data(:,3));
EPeY = (EPeXY.data(:,1)+EPeYY.data(:,2)+EPeYZ.data(:,3));
EPeZ = (EPeXZ.data(:,1)+EPeYZ.data(:,2)+EPeZZ.data(:,3));
EPe = irf.ts_vec_xyz(EPeXX.time,[EPeX EPeY EPeZ]);
%EPefac = irf_convert_fac(EPe,Bxyzav,Rxyz1);
%EPe.representation{2} = {'x','y','z'};
out = EPe;