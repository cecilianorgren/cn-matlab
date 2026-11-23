
A0 = @(z,B0) -B0*abs(z);
A0 = @(z,B0) -B0*log(cosh(z));

A = @(z,t,B0,E0) A0(z,B0) - E0*t;
Bx = @(z,t,B0,E0) -gradient(A(z,t,B0,E0),z);

z = linspace(-5,5,1000);
x = linspace(-5,5,2);

[X,Z] = ndgrid(x,z);


B0 = 1;
E0 = 1;
t1 = 0;
t2 = 0.5;
minA = floor(min([A(z(:),t1,B0,E0);A(z(:),t2,B0,E0)]));
maxA = ceil(max([A(z(:),t1,B0,E0);A(z(:),t2,B0,E0)]));
levA = linspace(minA,maxA,maxA-minA+1);
traceA = -3;

h = setup_subplots(2,2); 
isub = 1;

hca = h(isub); isub = isub + 1;
%plot(hca,z,B(z,t1,B0,E0),z,B(z,t2,B0,E0))
plot(hca,z,Bx(z,t1,B0,E0))
hca.YLim = hca.YLim*1.2;
hca.XLabel.String = 'z';
hca.YLabel.String = 'B_x';
%irf_legend(hca,{"B_x = " + ''},[0.02 0.98])

hca = h(isub); isub = isub + 1;
plot(hca,z,A(z,t1,B0,E0),z,A(z,t2,B0,E0))
hold(hca,'on')
[~,iz] = min(abs(A(z,t1,B0,E0)-traceA)); plot(hca,z(iz),A(z(iz),t1,B0,E0),'k.'); iz1 = iz;
[~,iz] = min(abs(A(z,t2,B0,E0)-traceA)); plot(hca,z(iz),A(z(iz),t2,B0,E0),'k.'); iz2 = iz;
%quiver(hca,[z(iz1),traceA],[z(iz2),traceA],0)
q = quiver(hca,z(iz1),traceA,z(iz2)-z(iz1),0,0);
q.LineWidth = 0.5;
q.MaxHeadSize = 1;
q.Color = 'k';
hold(hca,'off')
hca.XLabel.String = 'z';
hca.YLabel.String = 'A_y';
hca.Title.String = {"A_y = A_y(z,t=0)-E_0t"};
irf_legend(hca,{'A_y(t_1)','A_y(t_2)'}',[0.98 0.98])
%irf_legend(hca,{"A_y = A_y(z,t=0)-E_0t"},[0.02 0.02],'color','k')

hca = h(isub); isub = isub + 1;
contourf(hca,X,Z,A(Z,t1,B0,E0),levA,"ShowText","on")
hca.XLabel.String = 'x';
hca.YLabel.String = 'z';
hca.Title.String = 'A_y(t_1)';

hca = h(isub); isub = isub + 1;
contourf(hca,X,Z,A(Z,t2,B0,E0),levA,"ShowText","on")
hca.Title.String = 'A_y(t_2)';

hlinks = linkprop(h(3:4),{'CLim'});


%% SI units
units = irf_units;
L = 700e3; % m
A0 = @(z,B0) -B0*L.*abs(z/L);
A0 = @(z,B0) -B0*L.*log(cosh(z/L));

A = @(z,t,B0,E0) A0(z,B0) - E0*t;
Bx = @(z,t,B0,E0) -gradient(A(z,t,B0,E0),z);
Bz = @(z,t,B0,E0) gradient(A(z,t,B0,E0),x);

z = 3*L*linspace(-1,1,1000);
x = linspace(-1,1,2);

[X,Z] = ndgrid(x,z);


B0 = 5e-9;
E0 = 3e-3;
t1 = 0;
t2 = 0.8; % s
m = units.mp;
q = units.e;
minA = (min([A(z(:),t1,B0,E0);A(z(:),t2,B0,E0)])); % Wb
maxA = (max([A(z(:),t1,B0,E0);A(z(:),t2,B0,E0)]));
levA = linspace(minA,maxA,15);
showText = "off";
if 1
minA = floor(min([A(z(:),t1,B0,E0);A(z(:),t2,B0,E0)])*1e3)*1e-3;
maxA = ceil(max([A(z(:),t1,B0,E0);A(z(:),t2,B0,E0)])*1e3)*1e-3;
levA = linspace(minA,maxA,(maxA-minA)*1e3+1);
showText = "on";
end
traceA = (maxA+minA)/2;
traceA = -3*1e-3;

h = setup_subplots(2,2); 
isub = 1;

hca = h(isub); isub = isub + 1;
%plot(hca,z,B(z,t1,B0,E0),z,B(z,t2,B0,E0))
plot(hca,z*1e-3,Bx(z,t1,B0,E0)*1e9)
hca.YLim = hca.YLim*1.2;
hca.XLabel.String = 'z (km)';
hca.YLabel.String = 'B_x (nT)';
%irf_legend(hca,{"B_x = " + ''},[0.02 0.98])
hca.XLim = z([1 end])*1e-3;

hca = h(isub); isub = isub + 1;
plot(hca,z*1e-3,A(z,t1,B0,E0)*1e3,z*1e-3,A(z,t2,B0,E0)*1e3)
if 1
hold(hca,'on')
[~,iz] = min(abs(A(z,t1,B0,E0)-traceA)); plot(hca,z(iz)*1e-3,A(z(iz),t1,B0,E0)*1e3,'k.'); iz1 = iz;
[~,iz] = min(abs(A(z,t2,B0,E0)-traceA)); plot(hca,z(iz)*1e-3,A(z(iz),t2,B0,E0)*1e3,'k.'); iz2 = iz;

%q = quiver(hca,z(iz1)*1e-3,traceA,(z(iz2)-z(iz1))*1e-3,0,0);
%q.LineWidth = 0.5;
%q.MaxHeadSize = 1;
%q.Color = 'k';
hold(hca,'off')
end
hca.XLabel.String = 'z (km)';
hca.YLabel.String = 'A_y (mWb)';
hca.Title.String = {"A_y = A_y(z,t=0)-E_0t"};
%irf_legend(hca,{'A_y(t_1)','A_y(t_2)'}',[0.98 0.98])
%irf_legend(hca,{"A_y = A_y(z,t=0)-E_0t"},[0.02 0.02],'color','k')
hca.XLim = z([1 end])*1e-3;
ax = axes('position',hca.Position);
ax.YAxisLocation = "right";
ax.YLim = hca.YLim*q/m*1e-3*1e-3;
ax.Color = "none";
ax.YLabel.String = '-(m_p/e)\Delta{v_y} (km/s)';
ax.XTick = [];
%yyaxis(hca,'right')
%ylim(ylim(gca)*units.e/units.mp)
%ylabel('Scaled quantity')

hca = h(isub); isub = isub + 1;
contourf(hca,X,Z*1e-3,A(Z,t1,B0,E0)*1e3,levA*1e3,"ShowText",showText)
hca.XLabel.String = 'x';
hca.YLabel.String = 'z';
hca.Title.String = 'A_y(t_1)';

hca = h(isub); isub = isub + 1;
contourf(hca,X,Z*1e-3,A(Z,t2,B0,E0)*1e3,levA*1e3,"ShowText",showText)
hca.Title.String = 'A_y(t_2)';

hlinks = linkprop(h(3:4),{'CLim'});
colormap(flipdim(irf_colormap(gca,"Blues"),1))