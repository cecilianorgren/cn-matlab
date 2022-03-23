% Constants
units = irf_units;
G = units.G;
R = units.R_Sun;
M = units.M_Sun;
m = units.mp;
%T = 1e6; % Kelvin
T_examp = [1e6 4e6 7e6 8e6];
kB = units.kB;

vs_iso = @(T) sqrt(2*kB*T/m);
vs = vs_iso;
rc = @(T) G*M./(2*vs(T).^2);
v_esc = sqrt(2*G*M/R)

Tvec = logspace(5.5,7,100);

% Solar wind solution, just for illustrative purpose, so use rc = 1, vs = 1
vs_ = 1;
rc_ = 1;
fun = @(r,v,c) (v/vs_).^2 - log((v/vs_).^2) - 4*r/rc_ + 4*log(r/rc_) - c;
fun_ = @(r,v,rc,vs,c) (v/vs).^2 - log((v/vs).^2) - 4*r/rc + 4*log(r/rc) - c;
v_ = vs_*logspace(-2,log10(2.5),300);
r_ = rc_*logspace(-2,log10(5),300);
[R_,V_] = meshgrid(r_,v_);
colors = pic_colors('matlab');

cvec = -9:2:3;

% Plot
nrows = 4; ncols = 2;
h = setup_subplots(nrows,ncols,'vertical');

isub = 1;


hca = h(isub); isub = isub + 1; %hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,Tvec,vs_iso(Tvec))
hca.XLabel.String = 'T [K]';
hca.YLabel.String = 'v_s(T) [m/s]';
irf_legend(hca,'v_s = (2k_BT/m_i)^{1/2}  -  the sound speed increases with temperature',[0.02 0.98]);


hca = h(isub); isub = isub + 1; %hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,Tvec,rc(Tvec)) % ,Tvec,Tvec*0+R
hca.XLabel.String = 'T [K]';
hca.YLabel.String = 'r_c(T) [m]';
irf_legend(hca,'r_c = GM/2v_s^2  -  the critical radius decreases with temperature',[0.02 0.98]);

hca = h(isub); isub = isub + 1; %hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,Tvec,rc(Tvec)/R,Tvec,Tvec*0+1)
hold(hca,'on')
h_t = gobjects();
for TT = T_examp
  plot(hca,TT,rc(TT)/R,'o')
  h_t(end+1) = text(hca,TT,rc(TT)/R,sprintf(' T = %.0f K',TT),'verticalalignment','bottom');
end
h_t(end-1).VerticalAlignment = 'top';
h_t(end-1).HorizontalAlignment = 'right';
hold(hca,'off')
hca.XLabel.String = 'T [K]';
hca.YLabel.String = 'r_c(T)/R';
hca.YLim = [0 10];
txt = {'When the temperature is high enough\n'...
       'the critical radius will become\n'...
       'smaller than the solar radius.\n'...
       '',...
       'This means you only have access to the',...
       'solar wind solutions in the right part',...
       'of the graph. If you start out as subsonic,',...
       'the wind speed can only decrease.'};
txt = {'When the temperature is high enough the critical radius will become',...
       'smaller than the solar radius. Since the particle r is always bigger',...
       'than R_{Sun}, the particle r will always be larger than the critical radius.',...
       'r_c<R_{Sun}<r'};
irf_legend(hca,txt',[0.98 0.98],'color',[0 0 0]);
drawnow

hca = h(isub); isub = isub + 1; %hca = subplot(nrows,ncols,isub); isub = isub + 1;
for ic = 1:numel(cvec)
  c = cvec(ic);    
  vpow = 1;  
  hc = contour(hca,R_,V_,fun(R_,V_,c),[0 0],'color',colors(ic,:),'linewidth',1,'tag',num2str(c));  
  
  legs{ic} = sprintf('C = %4.0f',c);
  if ic == 1, hold(hca,'on'); end
  if ic == numel(cvec), hold(hca,'off'); end
end
hca.XLabel.String = 'r/r_c';
hca.YLabel.String = 'v/v_s';
hca.Title.String = {'Solar wind solutions: f(r,v) = (v/v_s)^2 - log((v/v_s)^2)-4r/r_c + 4log(r/r_c) - C = 0'};
legend(hca,legs,'location','east')
grid(hca,'on')

for T = T_examp
  hca = h(isub); isub = isub + 1; %hca = subplot(nrows,ncols,isub); isub = isub + 1;
  
  vs_tmp = vs(T);
  rc_tmp = rc(T);
  r_tmp = logspace(log10(1),log10(6e9),300); % m
  v_tmp = logspace(log10(1),log10(10e5),300); % m/s
  [R_tmp,V_tmp] = meshgrid(r_tmp,v_tmp);
  
  
  for ic = 1:numel(cvec)
    c = cvec(ic);    
    vpow = 1;  
    hc = contour(hca,R_tmp,V_tmp,fun_(R_tmp,V_tmp,rc_tmp,vs_tmp,c),[0 0],'color',colors(ic,:),'linewidth',1,'tag',num2str(c));  

    legs{ic} = sprintf('C = %4.0f',c);
    if ic == 1, hold(hca,'on'); end
    if ic == numel(cvec), hold(hca,'off'); end
  end
  % Add vs, rc, Rsun in plot
  hold(hca,'on')
  h_sun = plot(hca,v_tmp*0 + R,v_tmp,'k','linewidth',2); % Rsun
  h_rc = plot(hca,v_tmp*0 + rc_tmp,v_tmp,'r','linewidth',1); % rc
  h_vs = plot(hca,r_tmp,r_tmp*0 + vs_tmp,'g','linewidth',1); % vs
  h_rc_patch = patch(hca,[0,R,R,0],[0,0,v_tmp(end),v_tmp(end)],0.5*[1 1 1],'facealpha',0.5);
  hold(hca,'off')
  legend([h_sun,h_rc,h_vs],{'R_{sun}','r_c','v_s'},'location','eastoutside')
  
  hca.XLabel.String = 'r [m]';
  hca.YLabel.String = 'v [m/s]';
  %hca.Title.String = {'Solar wind solutions: f(r,v) = (v/v_s)^2 - log((v/v_s)^2)-4r/r_c + 4log(r/r_c) - C = 0'};
  hca.Title.String = sprintf('T = %.0f K, rc(T) = %.0f m, vs(T) = %.0f m/s',T,rc_tmp,vs_tmp);
  %legend(hca,legs,'location','east')
  grid(hca,'on')
end

% Add some text to specific axis
txt = {'Here the solution is plotted as a function of the radius and speed in SI units.',...
  'The solar radius (black vertical line) is always at the same r-location.',...
  'Depending on T, the rc and the vs changes, and therefore the solution changes.',...
  'However, relative to rc and vs, the solution is always the same.',...
  'The solar wind is always outside of the solar surface, so we can only',...
  'consider solutions that are to the right of the black vertical line.'};
irf_legend(h(5),txt',[0.13 0.98],'color',[0 0 0]);

txt = {'For this T, the critical radius',...
  'is larger than the radius of the Sun.',...
  'This means we can start from a',...
  'subsonic flow and reach a supersonic flow.'};
irf_legend(h(6),txt',[0.4 0.7],'color',[0 0 0],'horizontalalignment','left');

txt = {'For this T, rc < Rsun.',...
  'So to ever have a supersonic flow, we',...
  'need to start from a supersonic flow.',...
  'If we start from a subsonic flow,',...
  'the flow speed can only decrease as',...
  'we move away from the Sun (increasing r).'};
irf_legend(h(7),txt',[0.4 0.7],'color',[0 0 0],'horizontalalignment','left');




