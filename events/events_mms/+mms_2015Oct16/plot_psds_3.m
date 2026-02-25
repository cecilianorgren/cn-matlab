%%
  % Automatic results made from gui above
  ic = 1:4;
  tint_mva = irf.tint('2015-10-16T10:33:13.227751708Z/2015-10-16T10:33:38.076784912Z');
  c_eval('[out?,l?,mva?]=irf_minvar(dmpaB?.tlim(tint_mva));',ic)
  mva = mva3;
  % rotate around x-axis 
  turn_angle = 0; % degrees
  tmpv = mva;
  mva(2,:) = tmpv(2,:)*cosd(turn_angle) + tmpv(3,:)*sind(turn_angle);
  mva(3,:) = tmpv(3,:)*cosd(turn_angle) - tmpv(2,:)*sind(turn_angle);
  % v is 3x3 matrix:
  %v1: [0.0906 0.0896 0.9918] - first row
  %v2: [0.3510 -0.9349 0.0524] - second row
  %v3: [0.9320 0.3434 -0.1161] - third row   
  
%% Check knee distributions

tint = irf.tint('2015-10-16T10:33:45.00Z',0.2)+0*(-6);
vlim = 15*1e3;
elevlim = 90;
correctBin = 0;


c_eval('dist = desDist?;',ic)

c_eval('Vi0 = mean(vi?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatVi0 = double(irf_norm(Vi0));
c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatVe0 = double(irf_norm(Ve0));

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));

% Initialize figure
for ii = 1:8; h(ii) = subplot(2,4,ii); end

isub = 1;

% Plot projection onto a plane
hca = h(isub);  isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint,'xyz',[1 0 0; 0 0 1; 0 1 0],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 1 0; 0 0 1],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 0 1; cross(hatB0,[0 1 0])],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 0 1; hatB0],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);