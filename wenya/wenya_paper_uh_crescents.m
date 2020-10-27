%% Constant v
units = irf_units;
l = 35.4e3; % harris sheet width, m
B0 = 4.2e-9; % harris sheet asymptotic amplitude, T
AH = B0*l;

syms x y z %ah zh ap xp zp
R = [x y z];

AH = [0; -AH*log(cosh(z/l)); 0]; 
A = AH;
B = curl(A,R);
J = curl(B,R)/units.mu0; % (units.mu0*1e-9*1e-3)*

f_B = matlabFunction(B(1));
f_J = matlabFunction(J(2));

% Plot Jy vs Bx
zvec = linspace(-l,l,100)*2.7;
nrows = 3; ncols = 1; npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
if 1 % Bx, Jy, n
  hca = h(isub); isub = isub + 1;  
  plotyy(hca,zvec*1e-3,f_B(zvec)*1e9,zvec*1e-3,f_J(zvec)*1e9) 
  legend(hca,{'B_x (nT)','J_y (nA/m^2)'})
  hca.XLabel.String = 'z';
  %hca.YLabel.String = 'J_y';
end
if 1 % Bx, Jy, n
  hca = h(isub); isub = isub + 1;  
  plot(hca,abs(f_B(zvec))*1e9,f_J(zvec)*1e9) 
  hca.XLabel.String = '|B_x| (nT)';
  hca.YLabel.String = 'J_y (nA/m^2)';
end
if 1 % Bx, Jy, n
  hca = h(isub); isub = isub + 1;  
  plot(hca,abs(f_B(zvec)).^2*1e18,f_J(zvec)*1e9) 
  hca.XLabel.String = 'B_x^2 (nT)^2';
  hca.YLabel.String = 'J_y (nA/m^2)';
end

for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end
