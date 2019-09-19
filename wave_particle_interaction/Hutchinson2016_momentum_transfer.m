% Hutchinson 2016. Momentum change

units = irf_units;
fun_f_max = @(n,v,vd,vt) n.*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
fun_phi = @(phi0,lphi,x) phi0.*exp(-x.^2/2./lphi.^2);
fun_v = @(v0,vph,phi0,lphi,x,m) vph + sign(v0-vph).*((v0-vph).^2 + 2*units.e*fun_phi(phi0,lphi,x)/m).^0.5;

% What are reasoanly expected acceleration ranges? 
Eperp = 20e-3; % V/m
Epar = sqrt(units.me/units.mp).*Eperp; % V/m
a = units.e.*Epar/units.me*1e3; % F = qE = ma -> a = -e*E/me

% EH properties;
phi0 = 3000; % V
vph = 27000e3; % m/s
lphi = 60e3; % m

% properties of beam and sheet electron distributions
vd1 = 0e3; % m/s
n1 = 0.08e6; % m^-3
T1 = 800; % eV
vt1 = sqrt(2*units.e*T1./units.me); % m/s

vd2 = 35000e3; % m/s
n2 = 0.02e6; % m^-3
T2 = 200; % eV
vt2 = sqrt(2*units.e*T2./units.me); % m/s

vd3 = 00000e3; % m/s
n3 = n1+n2; % m^-3
T3 = 5000; % eV
vt3 = sqrt(2*units.e*T3./units.mp); % m/s

% transform speeds into EH frame
vd1 = vd1 - vph;
vd2 = vd2 - vph;
vd3 = vd3 - vph;
vph = vph - vph;

fun_Pocb_integrand = @(x,phi0,vph,Udot,m,v1,n1,vt1,vd1) m.*Udot.*(-2+3*v1./fun_v(v1,vph,phi0,lphi,x,m)-(v1./fun_v(v1,vph,phi0,lphi,x,m)).^3).*fun_f_max(n1,v1,vd1,vt1);
fun_Pocbt_integrand = @(x,phi0,vph,Udot,m,v1,n1,vt1,vd1) m.*Udot.*(-1+2*v1./fun_v(v1,vph,phi0,lphi,x,m)-(v1./fun_v(v1,vph,phi0,lphi,x,m)).^3).*fun_f_max(n1,v1,vd1,vt1);


nx = 1001;
xmax = 4*lphi;
xmin = -xmax;
x = linspace(xmin,xmax,nx);

nv = 1000; 
v1_min = -60000e3; % m/s
v1_max = 60000e3; % m/s
v1 = linspace(v1_min,v1_max,nv);

% Make 2D representation of integrand
% I guess the U dot 'sign' problemis taken vare of by v1. Although
% Hutchinson noted that v should always be larger than v1. That's not how
% it is now.
% Maybe divide the two electron distributions into larger than and smaller
% than vph?
Udot1 = a; % some value
Udot2 = a; % some value, opposite direction
Udot3 = a; % some value
[V1,X] = meshgrid(v1,x);
Poch1_INTEGRAND = fun_Pocbt_integrand(X,phi0,vph,Udot1,units.me,V1,n1,vt1,vd1);
Poch2_INTEGRAND = fun_Pocbt_integrand(X,phi0,vph,Udot2,units.me,V1,n2,vt2,vd2);
Poch3_INTEGRAND = fun_Pocb_integrand(X,phi0,vph,Udot3,units.mp,V1,n3,vt3,vd3);
Pocb12_top = Poch1_INTEGRAND + Poch2_INTEGRAND;
Pocb12_top(V1<0) = NaN;

Pocb12_bottom = Poch1_INTEGRAND + Poch2_INTEGRAND;
Pocb12_bottom(V1>0) = NaN;


h = setup_subplots(2,2); isub = 1;

if 1 % phi
  hca = h(isub); isub = isub + 1;
  plot(hca,x,fun_phi(phi0,lphi,x))
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = '\phi (V)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % distributions in EH frame
  hca = h(isub); isub = isub + 1;
  f_12_vpos = fun_f_max(n1,v1,vd1,vt1) + fun_f_max(n2,v1,vd2,vt2);
  f_12_vpos(v1<0) = NaN;
  f_12_vneg = fun_f_max(n1,v1,vd1,vt1) + fun_f_max(n2,v1,vd2,vt2);
  f_12_vneg(v1>0) = NaN;
  plotyy(hca,v1,[fun_f_max(n1,v1,vd1,vt1)' fun_f_max(n2,v1,vd2,vt2)'],... 
             v1,fun_f_max(n3,v1,vd3,vt3)) % NB: v1 does not refer to distribution, but to velocities outside of EH
  hold(hca,'on')
  plot(hca,v1,[f_12_vpos' f_12_vneg'],'linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'v_1 (m/s)';
  hca.YLabel.String = 'f (s/m^4)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  legend(hca,{'f_1','f_2','f_3','f_1+f_2 (v_1>0)','f_1+f_2 (v_1<0)'},'Box','off','location','northwest')
end
if 1 % integrands summed over v1
  hca = h(isub); isub = isub + 1;
  plot(hca,x,nansum(Poch1_INTEGRAND,2),x,nansum(Poch2_INTEGRAND,2),x,nansum(Poch3_INTEGRAND,2),x,nansum(Pocb12_top,2),x,nansum(Pocb12_bottom,2))    
  hca.XLabel.String = 'integrand sum over v_1';
  hca.YLabel.String = 'x (m)';  
  hca.XLim = x([1 end]);
  legend(hca,{'f_1','f_2','f_3','f_1+f_2 (v_1>0)','f_1+f_2 (v_1<0)'},'Box','off')
end
if 0 % integrand over v1,x space
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V1,Poch1_INTEGRAND)
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'v_1 (m/s)';
end
if 0 % integrand over v1,x space
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V1,Poch2_INTEGRAND)
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'v_1 (m/s)';
end
if 0 % integrand over v1,x space
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V1,Poch3_INTEGRAND)
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'v_1 (m/s)';
end
if 0 % integrand over v1,x space, v1>0
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V1,Pocb12_top)
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'v_1 (m/s)';
end
if 0 % integrand over v1,x space, v1<0
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V1,Pocb12_bottom)
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'v_1 (m/s)';
end
if 1 % integrand over v1,x space, v1<0 and -v1>0
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V1,nansum(cat(3,cat(3,Pocb12_top,-Pocb12_bottom),-Poch3_INTEGRAND),3))
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colormap(hca,pic_colors('blue_red'))
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'v_1 (m/s)';
end
if 0 % difference in integrand over v1,x space
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V1,Poch2_INTEGRAND-Poch1_INTEGRAND)
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'v_1 (m/s)';
end
if 0 % integrands summed over v1 + added up with different signs of Udot, depending on if v1>0 or v1<0 
  hca = h(isub); isub = isub + 1;
  plot(hca,x,nansum(Pocb12_top,2)-nansum(Pocb12_bottom,2)-nansum(Poch3_INTEGRAND,2))    
  hca.XLabel.String = 'integrand sum over v_1';
  hca.YLabel.String = 'x (m)';  
  hca.XLim = x([1 end]);
  legend(hca,{'f_1','f_2','f_3','f_1+f_2 (v_1>0)','f_1+f_2 (v_1<0)'},'Box','off')
end