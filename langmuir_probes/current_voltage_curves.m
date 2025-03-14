%% 
units = irf_units;
me = units.me;
mi = units.mi;
e = units.e;


Ie_Vneg = @(V,A_cm2,T_eV,n_cc) 10^11*(A_cm2*n_cc*e/4)*sqrt(8*T_eV*units.e/pi/me).*exp(V/T_eV);
Ie_Vpos = @(V,A_cm2,T_eV,n_cc) 10^11*(A_cm2*n_cc*e/4)*sqrt(8*T_eV*units.e/pi/me).*(1+V/T_eV);
Ie = @(V,A_cm2,T_eV,n) Ie_Vneg(V,A_cm2,T_eV,n).*heaviside(-V) + Ie_Vpos(V,A_cm2,T_eV,n).*heaviside(V);

Ii_Vneg = @(V,A_cm2,T_eV,n) -10^11*(A_cm2*n*e/4)*sqrt(8*T_eV*units.e/pi/mi).*(1-V/T_eV);
Ii_Vpos = @(V,A_cm2,T_eV,n) -10^11*(A_cm2*n*e/4)*sqrt(8*T_eV*units.e/pi/mi).*exp(-V/T_eV);
Ii = @(V,A_cm,Te_eV,n) Ii_Vneg(V,A_cm,Te_eV,n).*heaviside(-V) + Ii_Vpos(V,A_cm,Te_eV,n).*heaviside(V);

Jph_Vpos = @(V) (80*exp(-V/2)+3*exp(-V/7.5));%.*heaviside(V); % muA/m^2
Jph_Vneg = @(V) 80;
Iph = @(V,A) A*Jph_Vneg(V).*heaviside(-V) + A*Jph_Vpos(V).*heaviside(V); % A


Iph_Vpos = @(V) 270*(exp(-V/2)+(3/80)*exp(-V/7.5));%.*heaviside(V); % nA
Iph_Vneg = @(V) Iph_Vpos(0);
Iph = @(V,Arel) Arel*(Iph_Vneg(V).*heaviside(-V) + Iph_Vpos(V).*heaviside(V)); % nA

Itot = @(V,A_cm2,T_eV,n,Arel) -Iph(V,Arel) + Ie(V,A_cm2,T_eV,n) + Ii(V,A_cm2,T_eV,n);

%% Component currents for a single set of parameters
h = setup_subplots(5,1);
isub = 1;


V = linspace(-30,30,1000);
A_cm2 = 200*5; % m^2 : 2000 mm x 50 mm
A_m2 = A_cm2*1e-4;
A_cm2 = 50.25; % cm2, for EFW on Cluster
n_cc = 100; % cm^-3
Te_eV = 10; % eV
Te_V = Te_eV/units.eV; % V

if 1 % Photoelectron current
  hca = h(isub); isub = isub + 1;
  plot(hca,V,-Iph(V))
  hca.XLabel.String = 'V (Volt)';
  hca.YLabel.String = 'I_{ph} (n{A})';  
end
if 1 % Ambient electron current
  hca = h(isub); isub = isub + 1;
  plot(hca, V, Ie(V,A_cm2,Te_eV,n_cc),'-', V, Ie_Vpos(V,A_cm2,Te_eV,n_cc),'--', V, Ie_Vneg(V,A_cm2,Te_eV,n_cc),'--')
  hca.XLabel.String = 'V (Volt)';
  hca.YLabel.String = 'I_e (n{A})';  
  irf_legend(hca,sprintf('A = %g cm^2, T_e = %g eV, n = %g cc',A_cm2,Te_eV,n_cc),[0.02 0.98])
end
if 1 % Ambient ion current
  hca = h(isub); isub = isub + 1;
  plot(hca, V, Ii(V,A_cm2,Te_eV,n_cc),'-', V, Ii_Vpos(V,A_cm2,Te_eV,n_cc),'--', V, Ii_Vneg(V,A_cm2,Te_eV,n_cc),'--')
  hca.XLabel.String = 'V (Volt)';
  hca.YLabel.String = 'I_i (n{A})';  
  irf_legend(hca,sprintf('A = %g cm^2, T_e = %g eV, n = %g cc',A_cm2,Te_eV,n_cc),[0.02 0.98])
end
if 1 % Total current
  hca = h(isub); isub = isub + 1;
  plot(hca, V, Itot(V,A_cm2,Te_eV,n_cc),'-')
  hca.XLabel.String = 'V (Volt)';
  hca.YLabel.String = 'I_i (n{A})';  
  irf_legend(hca,sprintf('A = %g cm^2, T_e = %g eV, n = %g cc',A_cm2,Te_eV,n_cc),[0.02 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Gradient of total current, dI/dV
  hca = h(isub); isub = isub + 1;
  dIdV = gradient(Itot(V,A_cm2,Te_eV,n_cc),V);
  plot(hca, V, dIdV,'-')
  hold(hca,'on')
  dIdV2 = gradient(Itot(V,A_cm2,Te_eV,n_cc*2),V);
  dIdV3 = gradient(Itot(V,A_cm2,Te_eV,n_cc*0.5),V);
  plot(hca, V, dIdV2, V, dIdV3)
  hold(hca,'off')
  hca.XLabel.String = 'V (Volt)';
  hca.YLabel.String = 'dI_i/dV (n{A})';  
  irf_legend(hca,sprintf('A = %g cm^2, T_e = %g eV, n = %g cc',A_cm2,Te_eV,n_cc),[0.02 0.98])
end

for ip = 1:numel(h)  
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end









%% Total currents for a set of parameters
h = setup_subplots(2,1);
isub = 1;

V = linspace(-5,10,1000);
A_cm2 = 50.25; % cm2, for EFW on Cluster
n_cc = logspace(-2,1,4); % cm^-3
n_cc = [10]; % cm^-3
T_eV = logspace(1,4,4); % eV
T_eV = [100];
nN = numel(n_cc);
nT = numel(T_eV);

Arel = linspace(0,1,11);
nA = numel(Arel);


if 1 % Total current
  hca = h(isub); isub = isub + 1;
  holdon = 0;
  legs = {};
  for iN = 1:nN
    for iT = 1:nT
      for iA = 1:nA
        plot(hca, V, Itot(V,A_cm2,T_eV(iT),n_cc(iN),Arel(iA)),'-')
        legs{end+1} = sprintf('T = %g eV, n = %g cc, A/A_{max} = %.2f',T_eV(iT),n_cc(iN),Arel(iA));
        if not(holdon); hold(hca,'on'); end
      end
    end
  end
  hold(hca,'off')
  legend(hca,legs,'location','eastoutside')
  hca.XLabel.String = 'V (Volt)';
  hca.YLabel.String = 'I_i (n{A})';  
  %irf_legend(hca,sprintf('A = %g cm^2, T_e = %g eV, n = %g cc',A_cm2,Te_eV,n_cc),[0.02 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Gradient of total current, dI/dV
  hca = h(isub); isub = isub + 1;
  holdon = 0;
  legs = {};
  for iN = 1:nN
    for iT = 1:nT
      for iA = 1:nA
        dIdV = gradient(Itot(V,A_cm2,T_eV(iT),n_cc(iN),Arel(iA)),V);
        plot(hca, V, dIdV,'-')
        legs{end+1} = sprintf('T = %g eV, n = %g cc, A/A_{max} = %.2f',T_eV(iT),n_cc(iN),Arel(iA));
        if not(holdon); hold(hca,'on'); end
      end
    end
  end
  hold(hca,'off')
  legend(hca,legs,'location','eastoutside')
  hca.XLabel.String = 'V (Volt)';
  hca.YLabel.String = 'dI_i/dV (n{A}/V)';  
  %irf_legend(hca,sprintf('A = %g cm^2, T_e = %g eV, n = %g cc',A_cm2,Te_eV,n_cc),[0.02 0.98])
end

hlinks = linkprop(hca,{'XLim'});
for ip = 1:numel(h)  
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).FontSize = 16;
end








