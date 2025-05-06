%no02m = PIC('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
no02m = PIC('/Users/cecilianorgren/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%%
twpe = 19000;
pic = no02m.twpelim(twpe).xlim([75 125]).zlim([-7 7]);

varstrs = {'Ez','By','intEz','Ax','intEz./Ax','vx(3)','vx(3)./vz(3)'}';
cmap_br = pic_colors('blue_red');
cmaps = {cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br};
clims = {[-1 1],[-1 1],[-1 1],[-1 1],2*[-1 1],[-1 1],[-1 1],0.5*[-1 1]};

pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'sep')


h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 17;',1:numel(h))

%%
twpe = 19000;
pic = no02m.twpelim(twpe).xlim([92 93]).zlim([1 5]);

varstr_mod_v2div2 = '(vx(3).^2+vz(3).^2)/2';
varstr_mod_vx = '-(Ax./intEz).*(vx(3).^2+(vz(3)+0.).^2)/2';
varstr_mod_vz = '-(vx(3).^2+vz(3).^2).^(0.5).*(1-1.*(Ax./intEz).^2.*(vx(3).^2+vz(3).^2)/4).^0.5';

%varstr_mod_vz = '(1-4.*(vx(3).^2+vz(3).^2).^(-1).*(intEz./Ax).^2).^0.5';
varstrs = {{'Ez','By'},{varstr_mod_v2div2,'intEz',[varstr_mod_v2div2,'-intEz']},{'intEz','Ax'},{'intEz./Ax'},{'vx(3)','vz(3)'},{'vx(3)',varstr_mod_vx},{'vz(3)',varstr_mod_vz},{'vx(3)./vz(3)'},{'atand(vx(3)./vz(3))'}}';
cmap_br = pic_colors('blue_red');
cmaps = {cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br};
clims = {[-1 1],[-1 1],[-1 1],[-1 1],2*[-1 1],[-1 1],[-1 1],0.5*[-1 1]};

%pic.plot_line('z',varstrs,'A',1,'clim',clims,'cmap','sep')
pic.plot_line('z',varstrs)


h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 17;',1:numel(h))

%% Do everything manually
fontsize = 16;

twpe = 19000;
xlim = [92 93];
xlim = 98+0.25*[-1 1]-7;

iSpecies = 3;
zlim = [2 5];
zref = 5.0;

if 0
iSpecies = 5;
zlim = [-5 -1];
zref = -4.5;
end

lim_nfrac = 0.5;


pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
pic_map = no02m.twpelim(twpe).xlim(xlim+[-20 20]).zlim(zlim+[-2 2]);

z = pic.zi;
dz = z(2) - z(1);
iz_ref = pic.ind_from_lim(z,zref);

ntot = mean(pic.ni,1);
n = mean(pic.n(iSpecies),1);
nfrac = n./ntot;
% Find limit where there is a mixture of densities
if zref > 0
  iz_nlim = find(nfrac<lim_nfrac,1,'last');
  if isempty(iz_nlim), iz_nlim = 1; end
  bad_n_region = [zlim(1) z(iz_nlim) z(iz_nlim) zlim(1)];
elseif zref < 0
  iz_nlim = find(nfrac<lim_nfrac,1,'first');
  if isempty(iz_nlim), iz_nlim = 1; end
  bad_n_region = [zlim(2) z(iz_nlim) z(iz_nlim) zlim(2)];
end
%if numel(bad_n_region) < 4, bad_n_region = [NaN NaN NaN NaN]; end

vx = mean(pic.vx(iSpecies),1);
vy = mean(pic.vy(iSpecies),1);
vz = mean(pic.vz(iSpecies),1);
%v = sqrt(vx.^2 + vz.^2);

vx = vx - vx(iz_ref);
vz = vz - vz(iz_ref);
v = sqrt(vx.^2 + vz.^2);
vxyz = sqrt(vx.^2 + vy.^2 + vz.^2);


Epar = mean(pic.Epar,1);
Eperpx = mean(pic.Eperpx,1);
Eperpz = mean(pic.Eperpz,1);
Ex = mean(pic.Ex,1);
Ez = mean(pic.Ez,1);
By = mean(pic.By,1);

Ez_map = pic_map.Ez;

if 1
intEz = cumtrapz(z,Ez); %intEz = intEz - intEz(iz_ref);
intBy = cumtrapz(z,By); %intBy = intBy - intBy(iz_ref);

phi = -intEz; phi = phi - phi(iz_ref);
Ax = intBy;  Ax = Ax - Ax(iz_ref);
else
intEz = cumtrapz(z,Ez); intEz = intEz - intEz(iz_ref);
intBy = cumtrapz(z,By); intBy = intBy - intBy(iz_ref);

phi = -intEz; 
Ax = intBy;  

end
phiAx = phi./Ax; 
phiAx(abs(phiAx)>5) = NaN;
phiAx(abs(phiAx)<0.3) = NaN;

mod_vx_from_px = -Ax;
mod_vx = -(1./phiAx).*v.^2/2;
mod_vz_from_phi = sign(phi).*sqrt(2*abs(phi)); % not correct according to model, but testing out
mod_vz = -v.*(1-(1./phiAx).^2.*v.^2/4).^0.5;

nRows = 8;
nCols = 2;
h = setup_subplots(nRows,nCols,'horizontal');
isub = 1;

hca = h(isub); isub = isub + 1;
pic_map.plot_map(hca,{'Ez'},'A',1,'sep','cmap',pic_colors('blue_red'))
hmap = findobj(hca.Children,'Type','Image');
hca.CLim = prctile(abs(hmap.CData(:)),99)*[-1 1];
hold(hca,'on')
plot(hca,xlim(1)*[1 1],hca.YLim,'k',xlim(2)*[1 1],hca.YLim,'k')
hold(hca,'off')
hca.XGrid = "on"; hca.YGrid = "on";

hca = h(isub); isub = isub + 1;
pic_map.plot_map(hca,{['vx(' num2str(iSpecies) ')']},'A',1,'sep','cmap',pic_colors('blue_red'))
hmap = findobj(hca.Children,'Type','Image');
hca.CLim = prctile(abs(hmap.CData(:)),99)*[-1 1];
hold(hca,'on')
plot(hca,xlim(1)*[1 1],hca.YLim,'k',xlim(2)*[1 1],hca.YLim,'k')
hold(hca,'off')
hca.XGrid = "on"; hca.YGrid = "on";

hca = h(isub); isub = isub + 1;
pic_map.plot_map(hca,{['vz(' num2str(iSpecies) ')']},'A',1,'sep','cmap',pic_colors('blue_red'))
hmap = findobj(hca.Children,'Type','Image');
hca.CLim = prctile(abs(hmap.CData(:)),99)*[-1 1];
hold(hca,'on')
plot(hca,xlim(1)*[1 1],hca.YLim,'k',xlim(2)*[1 1],hca.YLim,'k')
hold(hca,'off')
hca.XGrid = "on"; hca.YGrid = "on";

hca = h(isub); isub = isub + 1;
pic_map.plot_map(hca,{['atand(vx(' num2str(iSpecies) ')./vz(' num2str(iSpecies) '))']},'A',1,'sep','cmap',pic_colors('blue_red'))
hmap = findobj(hca.Children,'Type','Image');
hca.CLim = prctile(abs(hmap.CData(:)),99)*[-1 1];
hca.CLim = [0 90];
hold(hca,'on')
plot(hca,xlim(1)*[1 1],hca.YLim,'k',xlim(2)*[1 1],hca.YLim,'k')
hold(hca,'off')
hca.XGrid = "on"; hca.YGrid = "on";

isub_lim = isub;

hca = h(isub); isub = isub + 1;
plot(hca,z,ntot,z,n)
irf_legend(hca,{'n_{i,tot}','n_i'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,Ez,z,By,z,Ex)
irf_legend(hca,{'E_z','B_y','E_x'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,Eperpz,z,Eperpx,z,Epar)
irf_legend(hca,{'E_{\perp,z}','E_{\perp,x}','E_{||}'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,phi,z,Ax)
irf_legend(hca,{'\phi','A_x'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,phiAx)
irf_legend(hca,{'\phi/A_x'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,phi,z,v.^2/2,z,phi+v.^2/2,z,vxyz.^2/2)
irf_legend(hca,{'\phi','v^2/2','\phi+v^2/2','v_{xyz}^2'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,vx,z,mod_vx,z,mod_vx_from_px)
irf_legend(hca,{'v_x','mod v_x','v_x=-A_x'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,vz,z,mod_vz,z,mod_vz_from_phi)
irf_legend(hca,{'v_z','mod v_z','sign(\phi)(2|\phi|)^{1/2}'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,vx,z,vz,z,-v,z,vy)
irf_legend(hca,{'v_x','v_z','-v','v_y'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
plot(hca,z,mod_vx,z,mod_vz)
irf_legend(hca,{'mod v_x',' mod v_z'},[0.98 0.98],'fontsize',fontsize)

hca = h(isub); isub = isub + 1;
%plot(hca,z,mod_vx,z,mod_vz,z,mod_vz./v)
%irf_legend(hca,{'mod v_x',' mod v_z',' mod v_z/v'},[0.98 0.98],'fontsize',fontsize)
plot(hca,z,mod_vx,z,mod_vz)
irf_legend(hca,{'mod v_x',' mod v_z'},[0.98 0.98],'fontsize',fontsize)

if 1
hca = h(isub); isub = isub + 1;
plot(hca,z,atand(vx./vz),z,atand(mod_vx./mod_vz))
irf_legend(hca,{'tan^{-1}v_x/v_z','mod tan^{-1} v_x/v_z'},[0.98 0.98],'fontsize',fontsize)
end


linkprop(h(1:(isub_lim-1)),{"XLim","YLim"})
compact_panels(h(isub_lim:end),0.01)
drawnow

for ip = isub_lim:numel(h)
  h(ip).XLim = zlim;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).FontSize = fontsize;
  h(ip).LineWidth = 1.5;
  % Plot reference z
  hold(h(ip),'on')
  plot(h(ip),zref*[1 1],h(ip).YLim,'color',[0.0 0.0 0.0],'LineStyle','--')
  hold(h(ip),'off')
  % Plot reference where n/ntot>xxx
  hold(h(ip),'on')
  %plot(h(ip),z(iz_nlim)*[1 1],h(ip).YLim,'color',[0.5 0.5 0.5],'LineStyle','--')  
  patch(h(ip),bad_n_region,h(ip).YLim([1 1 2 2]),'k','FaceAlpha',0.1,'EdgeColor','none') 
  hold(h(ip),'off')
end
text(h(1),zref,h(1).YLim(1),sprintf('z_{ref} = %3.1f ',zref),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',fontsize)
text(h(1),z(iz_nlim),h(1).YLim(2),sprintf('n/n_{tot} = %g ',lim_nfrac),'VerticalAlignment','top','HorizontalAlignment','right','FontSize',fontsize)


h(end).XLabel.String = 'z (d_i)';
h(1).Title.String = sprintf('twpe = %g, x/d_i = [%5.1f, %5.1f]',twpe,xlim(1), xlim(2));
hl = findobj(h,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))

%%
twpe = 19000;
pic = no02m.twpelim(twpe).xlim([94 95]).zlim([-3.5 -1]);

varstr_mod_v2div2 = '(vx(5).^2+vz(5).^2)/2';
varstr_mod_vx = '-(Ax./intEz).*(vx(5).^2+(vz(5)+0.).^2)/2';
varstr_mod_vz = '-(vx(5).^2+vz(5).^2).^(0.5).*(1-(Ax./intEz).^2.*(vx(5).^2+vz(5).^2)/4).^0.5';

%varstr_mod_vz = '(1-4.*(vx(3).^2+vz(3).^2).^(-1).*(intEz./Ax).^2).^0.5';
varstrs = {{'Ez','By'},{varstr_mod_v2div2,'intEz',[varstr_mod_v2div2,'-intEz']},{'intEz','Ax'},{'intEz./Ax'},{'vx(5)','vz(5)'},{'vx(5)',varstr_mod_vx},{'vz(5)',varstr_mod_vz},{'vx(5)./vz(5)'},{'atand(vx(5)./vz(5))'}}';
cmap_br = pic_colors('blue_red');
cmaps = {cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br};
clims = {[-1 1],[-1 1],[-1 1],[-1 1],2*[-1 1],[-1 1],[-1 1],0.5*[-1 1]};

%pic.plot_line('z',varstrs,'A',1,'clim',clims,'cmap','sep')
pic.plot_line('z',varstrs)


h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 17;',1:numel(h))
%%
twpe = 20000;
pic = no02m.twpelim(20000).xlim([75 125]).zlim([-7 7]);

varstrs = {'Ez','vx(3)','vz(3)','vx(3)./vz(3)'; 'By','vx(5)','vz(5)','vx(5)./vz(5)'}';
clims = {[-1 1],[-1 1],[-1 1],0.5*[-1 1],[-1 1],[-1 1],[-1 1],0.5*[-1 1]};
pic.plot_map(varstrs,'A',1,'clim',clims)

%%
twpe = 20000;
pic = no02m.twpelim(twpe).xlim([75 125]).zlim([-7 7]);

varstrs = {'veperpx','viperpx','vExBx'}';
cmap_br = pic_colors('blue_red');
cmaps = {cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br};
clims = {[-1 1],5*[-1 1],5*[-1 1],[-1 1],[-1 1],2*[-1 1],[-1 1],[-1 1],0.5*[-1 1]};
clims = {2*[-1 1],2*[-1 1],2*[-1 1],0.5*[-1 1],[-1 1],[-1 1],[-1 1],0.5*[-1 1]};
pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'sep')

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 17;',1:numel(h))
%%
twpe = 18000;
pic = no02m.twpelim(twpe).xlim([75 125]).zlim([-7 7]);

varstrs = {'Epar','vex','vepar'}';
cmap_br = pic_colors('blue_red');
cmaps = {cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br, cmap_br};
clims = {[-1 1],5*[-1 1],5*[-1 1],[-1 1],[-1 1],2*[-1 1],[-1 1],[-1 1],0.5*[-1 1]};

pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'sep')

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 17;',1:numel(h))
%%
twpe = 20000;
pic = no02m.twpelim(20000).xlim(100+0.5*[-1 1]).zlim([0 1]);

varstrs = {{'Ez','By'},...
  {'cumsum(Ez,2)','cumsum(By,2)'},...
  {'cumsum(Ez,2)./cumsum(By,2)'},...
  {'vx([3])','vz([3])'},...
  {'vx([3 5])','vx([3])','vx([5])'},...
  {'vz([3 5])','vz([3])','vz([5])'},...
  ...%{'vx([3 5]).^2 + vz([3 5]).^2'},...
  {'(vx([3]).^2 + vz([3]).^2)'},...
  {'vx(3)./vz(3)'}}';
clims = {[-1 1],[-1 1],[-1 1],0.5*[-1 1],[-1 1],[-1 1],[-1 1],0.5*[-1 1]};
pic.plot_line('z',varstrs)

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
h(end).YLim = [0 6];

%% Check ratios in a region limted to inside separatrix
%varstrs = {'ti','te','Babs'};
%cmaps = {pic_colors('thermal'),pic_colors('thermal'),pic_colors('thermal')};

pic = no02m.twpelim(21000).xlim([60 140]).zlim([-10 10]);
vx3 = pic.vx(3);
vx5 = pic.vx(5);
vz3 = pic.vz(3);
vz5 = pic.vz(5);
Ay = pic.A;

[Ainds,Avals] = saddle(Ay,'sort','minpeakprominence',0.001);
%Ay_xline = Avals(1);
Ay_xline = pic.Axline;
Ay_presheet = 10;
[X,Z] = ndgrid(pic.xi,pic.zi);


for iA = 1
%iA = -1;
Amin = Ay_xline + 0.5*(iA-1);
Amax = Ay_xline + 0.5*iA;

Amin = 2.9;
Amax = 3.4;

Amin = -3.0;
Amax = -2.3;

%Amin = -inf;
%Amax = inf;

% Coming from the top
vx3(Ay<Amin) = NaN;
vx3(Ay>Amax) = NaN;
vx3(Z<0) = NaN; % we want to check the incoming side
vz3(Ay<Amin) = NaN;
vz3(Ay>Amax) = NaN;
vz3(Z<0) = NaN; % we want to check the incoming side


% Coming from the bottom
vx5(Ay<Amin) = NaN;
vx5(Ay>Amax) = NaN;
vx5(Z>0) = NaN; % we want to check the incoming side
vz5(Ay<Amin) = NaN;
vz5(Ay>Amax) = NaN;
vz5(Z>0) = NaN; % we want to check the incoming side


%ti(ti<0) = 0;
%te(te<0) = 0;

Alevels = Ay_xline + -20:0.5:20;


v_edges = linspace(-1,1,51);
%v_edges = linspace(1,1,100);



%ti_edges = linspace(0,0.15,100); % for pressure
%te_edges = linspace(0,0.15,100);

%histcn_plot([ti(:), te(:)], t_edges,t_edges)
[count3 edges3 mid3 loc3] = histcn([vx3(:), vz3(:)], v_edges,v_edges(2:end));
[count5 edges5 mid5 loc5] = histcn([vx5(:), vz5(:)], v_edges,v_edges(2:end));

hca = subplot(2,4,[1 2]);
pcolor(hca,pic.xi,pic.zi,vx3')
shading(hca,'flat')
hca.XLabel.String = 'x (c/\omega_{pi})';
hca.YLabel.String = 'z (c/\omega_{pi})';
hcb = colorbar(hca);
hcb.YLabel.String = 'v_{ix}';
hold(hca,'on')
hca.CLim = [-1 1];
clim = hca.CLim;
contour(hca,pic.xi,pic.zi,Ay',Alevels,'k')
hca.CLim = clim;
hold(hca,'off')
colormap(hca,pic_colors('blue_red'))

hca = subplot(2,4,[5 6]);
pcolor(hca,pic.xi,pic.zi,vz3')
shading(hca,'flat')
hca.XLabel.String = 'x (c/\omega_{pi})';
hca.YLabel.String = 'z (c/\omega_{pi})';
hcb = colorbar(hca);
hcb.YLabel.String = 'v_z';
hold(hca,'on')
hca.CLim = [-1 1];
clim = hca.CLim;
contour(hca,pic.xi,pic.zi,Ay',Alevels,'k')
hca.CLim = clim;
hold(hca,'off')
colormap(hca,pic_colors('blue_red'))

hca = subplot(2,4,[3 4 7 8]);
hca.Position = [0.57 0.15 0.3 0.7];

pcolor(hca,mid3{1},mid3{2},log10(count3)')
shading(hca,'flat')
%hca.CLim = prctile(count(:),[0 99.7]);
hcb = colorbar(hca);
hcb.YLabel.String = 'log_{10} Counts';
hca.Box = 'on';
hca.Layer = 'top';
hca.FontSize = 16;
hca.XLabel.String = 'v_x';
hca.YLabel.String = 'v_z';

colormap(hca,pic_colors('thermal'))

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 16;',1:numel(h))
drawnow
pause(0.1)
%cn.print(sprintf('no02m_twpe_22000_ne_te_exhaust_A%g',iA))
end

%%
pic = no02m.twpelim(21000).xlim([60 140]).zlim([-10 10]);
Ez = pic.Ez;
By = pic.By;
%Ay = pic.A;
%%
no02m.twpelim(23000).plot_map({'Ez','By','abs(Ez./By)','abs(vz(3)./vx(3))'}','clim',{[-1 1],0.3*[-1 1],5*[0 1],5*[0 1]},'smooth',2)




%% On observed data, data can be loaded with paper_ion_dynamics_inside_idr_hall_fields.m
time_xline_ion = irf_time('2017-07-11T22:34:02.00Z','utc>EpochTT');
tref = time_xline_ion +- 140;
intEz = irf_integrate(mvaE3.z,tref);
intJxBz = irf_integrate(mvaJxBne3_mVm.z,tref);
intBy = irf_integrate(mvaB3.y,tref);
aa = ne3.data*1e6;
bb = mvaE3.z.resample(ne3).data*1e-3;

cc = mvaJLBM3.data*1e-18;
dd = mvaJ3.y.data*1e-9;

Bx_Hall_data = -1e9*(units.e*aa.*bb - cc)./dd;
Bx_Hall_data(abs(Bx_Hall_data) > prctile(abs(Bx_Hall_data),90)) = NaN;
Bx_Hall = irf.ts_vec_xyz(mvaJLBM3.time, [Bx_Hall_data 0*Bx_Hall_data 0*Bx_Hall_data]);
Bx_Hall = Bx_Hall.x;

ee = sqrt(mvaB3.x.resample(Bx_Hall).data.^2 + Bx_Hall.data.^2);
Bxy = irf.ts_scalar(Bx_Hall.time, ee);

intBxy = irf_integrate(Bxy,tref);
%%


h = irf_plot(7);

hca = irf_panel('Vi');
irf_plot(hca,mvaVi3);
hca.ColorOrder = mms_colors('xyz');
hca.YLabel.String = 'v_i (km/s)';

hca = irf_panel('B');
irf_plot(hca,mvaB3);
hca.ColorOrder = mms_colors('xyz');
hca.YLabel.String = 'B (nT)';

hca = irf_panel('Bx Hall');
hca.ColorOrder = mms_colors('1234');
irf_plot(hca,{Bx_Hall.x.filt(0,5,[],5),mvaB3.x},'comp');
hca.YLabel.String = 'B_x (nT)';
hca.YLabel.Interpreter = 'tex';
hca.ColorOrder = mms_colors('1234');
irf_legend(hca,{'B_x^{Hall}','B_x'},[0.02 0.98])

hca = irf_panel('E');
irf_plot(hca,mvaE3.resample(mvaVi3));
hca.ColorOrder = mms_colors('xyz');
hca.YLabel.String = 'E (mV/m)';

hca = irf_panel('JxB/ne');
irf_plot(hca,mvaJxBne3_mVm.resample(mvaVi3));
hca.ColorOrder = mms_colors('xyz');
hca.YLabel.String = 'JxB/ne (mV/m)';

hca = irf_panel('int By, Ez, JxB/ne');
hca.ColorOrder = mms_colors('1234');
irf_plot(hca,{intEz, intJxBz, intBy, intBxy},'comp');
hca.YLabel.String = '\int X dt (nTs, s mV/m)';
hca.YLabel.Interpreter = 'tex';
hca.ColorOrder = mms_colors('1234');
irf_legend(hca,{'\int Ez dt', '\int JxB/ne dt', '\int By dt', '\int Bxy dt'},[0.02 0.98])

if 0
  hca = irf_panel('int By');
  irf_plot(hca,intBy);
  hca.YLabel.String = '\int By dt (nTs)';
  hca.YLabel.Interpreter = 'tex';
  
  hca = irf_panel('int Ez');
  irf_plot(hca,intEz);
  hca.YLabel.String = '\int Ez dt (s mV/m)';
  hca.YLabel.Interpreter = 'tex';
  
  hca = irf_panel('int JxB/ne');
  irf_plot(hca,intJxBz);
  hca.YLabel.String = '\int JxB/ne dt (s mV/m)';
  hca.YLabel.Interpreter = 'tex';
end


hca = irf_panel('int Ez / int By');
hca.ColorOrder = mms_colors('1234');
toplot1 = 1e-3*(intEz*1e-3)/(intBy*1e-9);
toplot2 = 1e-3*(intJxBz*1e-3)/(intBy*1e-9);
toplot3 = 1e-3*(intEz*1e-3)/(intBxy*1e-9);
irf_plot(hca,{toplot1, toplot2, toplot3},'comp');
hca.YLabel.String = '\int Ez dt/\int By dt (km/s)';
hca.ColorOrder = mms_colors('1234');
irf_legend(hca,{'\int Ez dt/\int By dt', '\int JxB/ne dt/\int By dt','\int Ez dt/\int Bxy dt'},[0.02 0.98])
hca.YLabel.Interpreter = 'tex';
hca.YLim = prctile([toplot1.data; toplot2.data; toplot3.data],[3 97]);



hl = findobj(h,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))

irf_zoom(h,'x',tint)