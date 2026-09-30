pdist = PD_use(1:50);
times = pdist.time;
nt = times.length;
n = irf.ts_scalar(times,ones(nt,1));
vx = linspace(10,2000,nt)';
v = irf.ts_vec_xyz(times,[vx*1.1 vx*0 vx*0.1]);
t_in = 200; % eV
t = irf.ts_tensor_xyz(times,permute(repmat(t_in*eye(3,3),[1 1 nt]),[3 1 2]));
t_diag = irf.ts_vec_xyz(times,ones(nt,3)*t_in);
b = irf.ts_vec_xyz(times,rand(nt,3)+[vx*0+10 vx*0+2 vx*0]);
scpot = irf.ts_scalar(times,0*ones(nt,1));

Tfac = mms.rotate_tensor(t,'fac',b,'pp'); 


pd_mod_bim = mms.make_model_dist(pdist,b,scpot,n,v,t);
pd_mod_gen = pdist_generalized_maxwellian(pdist,n,v,t); 
pd_mod_pdi = pdist.generate_dist(n,t,v,'max','species','ions'); 
pd_mod_gen_ups = pd_mod_gen.interpn(0+[62 62 32]);

pd_mod = pd_mod_gen_ups;
pd_mod = pd_mod_gen;
it = 2; 
semilogy(pd_mod_gen.depend{1}(it,:),pd_mod_gen.data(it,:,1,1),'*', ...
     pd_mod_gen_ups.depend{1}(it,:),pd_mod_gen_ups.data(it,:,1,1),'*')

%loglog(pd_mod_gen.depend{1}(it,:),pd_mod_gen.omni.data(it,:),'*', ...
%     pd_mod_gen_ups.depend{1}(it,:),pd_mod_gen_ups.omni.data(it,:),'*')
if 0
  %%
irf_plot({b,n,v,t, ...
  pd_mod_bim.deflux.omni.specrec, ...
  pd_mod_gen.deflux.omni.specrec, ...
  pd_mod_pdi.deflux.omni.specrec, ...
  pd_mod_ups_gen.deflux.omni.specrec})
h = findobj(gcf,'type','axes'); h = h(end:-1:1);

%h(end-1).YScale = 'log';
hlinks = linkprop(h(5:end),{'CLim','YScale'});
h(end).YScale = 'log';
colorbar
irf_plot_axis_align
end
%% GMM'ing
nMP = 1e6;
allMP = pd_mod.macroparticles('ntot',nMP,'skipzero',1);
%allMP = pd_mod.macroparticles('ntot',nMP*10,'skipzero',1);

vecK = 1;
nK = numel(vecK);
gm = cell(nt,nK);
for iK = 1:nK
  K = vecK(iK); 
  for it = 1:nt
    MP = allMP(it);
    ntot(it) = sum(MP.df.*MP.dv);
    options = statset('MaxIter',150,'TolFun',1e-5);
        
    R = [MP.vx, MP.vy, MP.vz];
    %gm{it,nGroups} = fitgmdist(X,nGroups,'Start',S,'SharedCovariance',false,'Options',options);
    gm{it,iK} = fitgmdist(R,K,'SharedCovariance',false,'Options',options,'RegularizationValue',0.1);
  end
end

% Reconstructing PDist from gmm
%moms_orig = pdist.moments;
it = 1;
%iK = 1;
n_kk = n.data*1e-3; % What are these units?
w = cellfun(@(x)x.ComponentProportion',gm,'UniformOutput',false); w = cat(2,w{:});
mu = cellfun(@(x)x.mu,gm,'UniformOutput',false); mu = cat(3,mu{:});
sigma = cellfun(@(x)x.Sigma,gm,'UniformOutput',false); sigma = cat(4,sigma{:});
PD_gmm = pdist_generalized_maxwellian(pd_mod,w.*repmat(n_kk',[gm{1}.NumComponents 1]),mu,sigma);

moms_mod = pd_mod.moments;
moms_gmm_ = PD_gmm.moments;
moms_gmm.n = PD_gmm.n; 
moms_gmm.V = PD_gmm.vel; 
moms_gmm.T = PD_gmm.T; 
% First component
moms_gm.n = irf.ts_scalar(pdist.time,n.data.*w(1,:));
moms_gm.V = irf.ts_vec_xyz(pdist.time,permute(mu(1,:,:),[3 2 1]));
v2_gm = permute(sigma(:,:,1,:),[4 1 2 3]);
T_gm = v2_gm*1e6*units.mp/units.eV;
moms_gm.T = irf.ts_tensor_xyz(pdist.time,T_gm);

nref = ntot*1e-12; % cc
%
%h = irf_plot(9);
[h,h2] = initialize_combined_plot('leftright',9,3,1,0.6,'vertical');
fontsize = 12;


if 1 % dEFlux input
  hca = irf_panel('omni deflux inp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,pd_mod.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
irf_legend(hca,{'PDist from input (gen. from n, v, T)'}',[0.99 0.1],'fontsize',fontsize+2,'backgroundcolor','w');
end
if 1 % dEFlux model
  hca = irf_panel('omni deflux mod');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,PD_gmm.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
irf_legend(hca,{'PDist from GMM'}',[0.99 0.1],'fontsize',fontsize+2,'backgroundcolor','w');
end

hca = irf_panel('n');
hca.ColorOrder = mms_colors('1234');
irf_plot(hca,{n,moms_mod.n,moms_gm.n,moms_gmm.n},'comp')
hca.YLabel.String = 'n';
hca.ColorOrder = mms_colors('1234');
irf_legend(hca,{'input TS','PDist from input','from GMM output','PDist from GMM from input'}',[0.02 0.98],'fontsize',fontsize+2);

for comp = ["x","y","z"]
  hca = irf_panel(char(["v" + comp]));
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{v.(comp),moms_mod.V.(comp),moms_gm.V.(comp),moms_gmm.V.(comp)},'comp')
  hca.YLabel.String = ['v_' + comp];
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'input TS','PDist from input','from GMM output','PDist from GMM from input'}',[0.02 0.98],'fontsize',fontsize+2);
end

for comp = ["xx","yy","zz"]
  hca = irf_panel(char(["T" + comp]));
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{t.(comp),moms_mod.T.(comp),moms_gm.T.(comp),moms_gmm.T.(comp)},'comp')
  hca.YLabel.String = ['T_{' + comp + '}'];
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'input TS','PDist from input','from GMM output','PDist from GMM from input'}',[0.02 0.98],'fontsize',fontsize+2);

end

h(end).XTickLabelRotation = 0;
hlinks = linkprop(h(1:2),{'CLim'});
hlinks.Targets(1).CLim = [0 10];
colormap([flipdim(irf_colormap('Spectral'),1)])
irf_plot_axis_align
h(1).Title.String = sprintf('N_{MP} = %g, nK = %g',nMP,numel(vecK));

% Other plots
Ek_inp = units.mp*(v.abs2.data*1e6)/2/units.eV;
ET_inp = t.trace/3;

isub = 1;

hca = h2(isub); isub = isub + 1;
semilogx(hca,Ek_inp,n.data, ...
             Ek_inp,moms_mod.n.data, ...
             Ek_inp,moms_gm.n.data, ...
             Ek_inp,moms_gmm.n.data)
hca.ColorOrder = mms_colors('1234b');
hca.XLabel.String = '(m/2)v^2_{inp} (eV)';
hca.YLabel.String = 'Density (cc)';
irf_legend(hca,{'Input','PDist from input','from GMM output','PDist from GMM from input'}',[0.02 0.98],'fontsize',fontsize+2);
hca.YLim = [0.5 1.5];


hca = h2(isub); isub = isub + 1;
ebinwidth = Ek_inp*0.15;
semilogx(hca,Ek_inp,t.trace.data/3, ...
             Ek_inp,moms_mod.T.trace.data/3, ...
             Ek_inp,moms_gm.T.trace.data/3, ...
             Ek_inp,moms_gmm.T.trace.data/3,...
             Ek_inp,ebinwidth)
hca.ColorOrder = mms_colors('1234b');
hca.XLabel.String = '(m/2)v^2_{inp} (eV)';
hca.YLabel.String = 'trace(T)/3 (eV)';
irf_legend(hca,{'Input','PDist from input','from GMM output','PDist from GMM from input','0.15*(m/2)v_{inp}^2'}',[0.02 0.98],'fontsize',fontsize+2);
hca.YLim = [0 5]*t_in;

%% Check different T

pdist = PD_use(1:100);
pdist = pdist.interpn([62 62 32]);
times = pdist.time;
nt = times.length;
n = irf.ts_scalar(times,ones(nt,1));
vabs = linspace(10,2000,nt)';
vdir = [0.1 0 1.1]; vdir = vdir/norm(vdir);
v = irf.ts_vec_xyz(times,vabs*vdir);

ts_eye = irf.ts_tensor_xyz(times,permute(repmat(eye(3,3),[1 1 nt]),[3 1 2]));
b = irf.ts_vec_xyz(times,rand(nt,3)+[vabs*0+10 vabs*0 vabs*0]);
scpot = irf.ts_scalar(times,0*ones(nt,1));


t_in = logspace(1,3.5,12);
%t_in = 100; % eV

nT = numel(t_in);
pd_mod = cell(1,nT);
for iT = 1:nT
  t = ts_eye*t_in(iT);
  pd_mod{iT} = pdist_generalized_maxwellian(pdist,n,v,t); 
end

%
n_all_ts = cellfun(@(x){x.n},pd_mod);
v_all_ts = cellfun(@(x){x.vel},pd_mod);
T_all_ts = cellfun(@(x){x.T},pd_mod);
specrec_all_ts = cellfun(@(x){x.deflux.omni.specrec},pd_mod);

n_all_mat = cellfun(@(x)x.data ,n_all_ts ,'UniformOutput' ,false);

vx_all_mat = cellfun(@(x)x.x.data ,v_all_ts ,'UniformOutput' ,false);
vy_all_mat = cellfun(@(x)x.y.data ,v_all_ts ,'UniformOutput' ,false);
vz_all_mat = cellfun(@(x)x.z.data ,v_all_ts ,'UniformOutput' ,false);

txx_all_mat = cellfun(@(x)x.xx.data ,T_all_ts ,'UniformOutput' ,false);
tyy_all_mat = cellfun(@(x)x.yy.data ,T_all_ts ,'UniformOutput' ,false);
tzz_all_mat = cellfun(@(x)x.zz.data ,T_all_ts ,'UniformOutput' ,false);

  
legs = arrayfun(@(x)sprintf('%.0f eV',x),t_in,'UniformOutput',false);
cmap = irf_colormap('gem','interp',nT);
%cmap = irf_colormap('Blues');

%irf_plot(specrec_all_ts)
%h = findobj(gcf,'type','axes'); h = h(end:-1:1); c_eval('h(?).YScale = ''log'';',1:numel(h))
% Plot moments as colormap
fontsize = 14;

h = setup_subplots(2,2,'vertical'); isub = 1;

vabs = v.abs.data;

Ek_inp = units.mp*(v.abs2.data*1e6)/2/units.eV;
Ebinwidth = Ek_inp*0.15;
plotx = Ek_inp;
ploty = t_in;
[PX,PY] = ndgrid(plotx,ploty);
n0 = repmat(n.data,[1,nT]);
vx0 = repmat(v.x.data,[1,nT]);
vy0 = repmat(v.y.data,[1,nT]);
vz0 = repmat(v.z.data,[1,nT]);
T0 = repmat(t_in,[nt,1]);

x_label = 'mv_{beam}^2/2 (eV)';
y_label = 'T_{beam} (eV)';


hca = h(isub); isub = isub + 1;
plotc = cat(2,n_all_mat{:})./n0;
pcolor(hca,PX,PY,plotc)
shading(hca,'flat')
hcb = colorbar(hca);
hcb.YLabel.String = 'n/n_0';
hca.CLim = [0 2];
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;

hca = h(isub); isub = isub + 1;
plotc = cat(2,txx_all_mat{:})./T0;
pcolor(hca,PX,PY,plotc)
shading(hca,'flat')
hcb = colorbar(hca);
hcb.YLabel.String = 'T_{xx}/T_0';
hca.CLim = [0 2];
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;

hca = h(isub); isub = isub + 1;
plotc = cat(2,tyy_all_mat{:})./T0;
pcolor(hca,PX,PY,plotc)
shading(hca,'flat')
hcb = colorbar(hca);
hcb.YLabel.String = 'T_{yy}/T_0';
hca.CLim = [0 2];
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;

hca = h(isub); isub = isub + 1;
plotc = cat(2,tzz_all_mat{:})./T0;
pcolor(hca,PX,PY,plotc)
shading(hca,'flat')
hcb = colorbar(hca);
hcb.YLabel.String = 'T_{zz}/T_0';
hca.CLim = [0 2];
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;

for ip = 1:4
  hca = h(ip);
  hold(hca,'on')
  plot(hca,plotx,plotx*0.15,'k--')
  hold(hca,'on')
end

hlinks = linkprop(h,{'XLim','YLim','CLim'});
colormap([flipdim(irf_colormap('Spectral'),1)])
c_eval('h(?).XScale = ''log'';',1:numel(h))
c_eval('h(?).YScale = ''log'';',1:numel(h))
c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))

h(1).Title.String = sprintf('v_{dir} = [%.2f, %.2f, %.2f]',vdir(1),vdir(2),vdir(3));

%% Plot moments as lines

fontsize = 14;

h = setup_subplots(3,3,'vertical'); isub = 1;

vabs = v.abs.data;

Ek_inp = units.mp*(v.abs2.data*1e6)/2/units.eV;
Ebinwidth = Ek_inp*0.15;
plotx = Ek_inp;
% ploty = t_in;
% [PX,PY] = ndgrid(plotx,ploty);
% x_label = 'mv_{drift}^2/2 (eV)';
% y_label = 'T_{in} (eV)';
% n0 = repmat(n.data,[1 nT]);
% vx0 = repmat(v.x.data,[1 nT]);
% vy0 = repmat(v.y.data,[1 nT]);
% vz0 = repmat(v.z.data,[1 nT]);
% T0 = repmat(t_in,[nt 1]);


hca = h(isub); isub = isub + 1;
%plotc = cat(2,vx_all_mat{:});
%plotc = plotc./vx0;
%pcolor(hca,PX,PY,plotc)
ploty = cat(2,vx_all_mat{:});
plot(hca,plotx,ploty)
hca.ColorOrder = cmap;
hca.XLabel.String = x_label;
hca.YLabel.String = 'v_x (cc)';
%irf_legend(hca,legs',[1.02 0.1])

hca = h(isub); isub = isub + 1;
ploty = cat(2,vy_all_mat{:});
plot(hca,plotx,ploty)
hca.ColorOrder = cmap;
hca.XLabel.String = x_label;
hca.YLabel.String = 'v_y (cc)';
irf_legend(hca,legs',[1.02 0.1])

hca = h(isub); isub = isub + 1;
ploty = cat(2,vz_all_mat{:});
plot(hca,plotx,ploty)
hca.ColorOrder = cmap;
hca.XLabel.String = x_label;
hca.YLabel.String = 'v_z (cc)';
irf_legend(hca,legs',[1.02 0.1])

hca = h(isub); isub = isub + 1;
ploty = cat(2,txx_all_mat{:});
plot(hca,plotx,ploty)
hca.ColorOrder = cmap;
hold(hca,'on')
plot(hca,plotx,Ebinwidth,'k--')
irf_legend(hca,['- - 0.15*' x_label],[0.02 0.1],'color','k','fontsize',fontsize)
hold(hca,'off')
hca.XLabel.String = x_label;
hca.YLabel.String = 'T_{xx} (cc)';
irf_legend(hca,legs',[1.02 0.1])
hca.YScale = 'log';

hca = h(isub); isub = isub + 1;
ploty = cat(2,tyy_all_mat{:});
plot(hca,plotx,ploty)
hca.ColorOrder = cmap;
hold(hca,'on')
plot(hca,plotx,Ebinwidth,'k--')
irf_legend(hca,['- - 0.15*' x_label],[0.02 0.1],'color','k','fontsize',fontsize)
hold(hca,'off')
hca.XLabel.String = x_label;
hca.YLabel.String = 'T_{yy} (cc)';
irf_legend(hca,legs',[1.02 0.1])
hca.YScale = 'log';

hca = h(isub); isub = isub + 1;
ploty = cat(2,tzz_all_mat{:});
plot(hca,plotx,ploty)
hca.ColorOrder = cmap;
hold(hca,'on')
plot(hca,plotx,Ebinwidth,'k--')
irf_legend(hca,['- - 0.15*' x_label],[0.02 0.1],'color','k','fontsize',fontsize)
hold(hca,'off')
hca.XLabel.String = x_label;
hca.YLabel.String = 'T_{zz} (cc)';
irf_legend(hca,legs',[1.02 0.1])
hca.YScale = 'log';


hca = h(isub); isub = isub + 1;
ploty = cat(2,n_all_mat{:});
plot(hca,plotx,ploty)
hca.ColorOrder = cmap;
hca.XLabel.String = x_label;
hca.YLabel.String = 'density (cc)';
irf_legend(hca,legs',[1.02 0.1])

c_eval('h(?).XScale = ''log'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))


