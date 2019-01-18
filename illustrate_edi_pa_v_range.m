% illustrate 
units = irf_units;

% EDI energy and corresponding velocity
E_edi = 500; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 25; % eV

E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
v_edi_plusminus = v_edi_plus-v_edi_minus;
dv_edi_minus = v_edi_minus - v_edi;
dv_edi_plus = v_edi_plus - v_edi;
dv_edi = dv_edi_plus - dv_edi_minus; % m/s

v_edi_edges = [v_edi_minus,v_edi_plus]*1e-3;
pa_edi_edges = -45:11.25:45;
pa_edi_centers = [(-45+11.25):11.25:0 0:11.25:(45-11.25)];
az_edges = 0:10:360;
[V,PA,AZ] = meshgrid(v_edi_edges,pa_edi_edges,az_edges);

VX = V.*sind(PA).*cosd(AZ);
VY = V.*sind(PA).*sind(AZ);
VZ = V.*cosd(PA);

hca = subplot(1,1,1);
C = PA;
surf(hca,VX,VY,VZ)
%%
view(hca,[0 0 1])
axis(hca,'equal')
hca.YLim(1) = 0;%[0*min(X(:)) max(X(:))];
hca.XLabel.String = 'v_{||} (km/s)';
hca.YLabel.String = 'v_{\perp} (km/s)';
%polar(PA,V);
