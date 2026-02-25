localuser = 'cecilia';
%db_table = readtable(['/Users/' localuser '/Data/MMS/DB_Lalti/database_and_overview_plots/SDB_10-Mar-2022_V1.0.csv']);

units = irf_units;

x = db_table.pos_x;
y = db_table.pos_y;
z = db_table.pos_z;
xRE = x*1e3/units.RE;
yRE = y*1e3/units.RE;
zRE = z*1e3/units.RE;

%varstrs = {'beta_i_us','Pdyn_us','thBn','Vx_us','Ni_us','Ti_us','B_us_abs','MA','','','','','',''};
%varstrs = {'beta_i_us','Pdyn_us','thBn','Vx_us','Ni_us','Ti_us','B_us_abs','MA','Mf','B_jump','Ni_jump','Te_jump'};


cmap = irf_colormap('bluered');
% %cmap = ('waterfall2');
cmap = pic_colors('bluegrayred');
color_ksr = [0 0 0];
color_bg = [1 1 1];
color_text = [0 0 0];
fontsize = 12;


Dp = 2;
Bz = 0;
[x_mp,y_mp,~] = irf_magnetosphere('mp_shue1998',Dp,Bz);
[x_bs,y_bs,~] = irf_magnetosphere('bs',Dp,Bz);
x_mp = [x_mp(end:-1:1) x_mp];
y_mp = [y_mp(end:-1:1) -y_mp];
%ir_mp = find(abs(y_mp) > 20);
%x_mp(ir_mp) = [];
%y_mp(ir_mp) = [];
x_bs = [x_bs(end:-1:1) x_bs];
y_bs = [y_bs(end:-1:1) -y_bs];
%ir_bs = find(abs(y_bs) > 20);
%x_bs(ir_bs) = [];
%y_bs(ir_bs) = [];

xlim = [-4.99 18];
ylim = [-25 25];

nrows = 1;
ncols = 2;
h = setup_subplots(nrows,ncols);
isub = 1;


if 1 % shock normal angle
  hca = h(isub); isub = isub + 1;    
  var = db_table.('thBn');
  %var(var < -1e20) = NaN;
  scatter(hca,xRE,yRE,15,var,'filled')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'Shock normal angle';
  hcb.XLabel.Interpreter = 'none';
  hcb.Location = 'west'; 
  hcb.YLabel.Color = color_text;
  hcb.Color = color_text;
  hca.XLabel.String = 'x_{GSE} (R_E)';
  hca.YLabel.String = 'y_{GSE} (R_E)';
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  axis(hca,'equal')
  colormap(hca,cmap)
  hca.Color = color_bg;
  hca.CLim = prctile(var,[5 95]);
  hca.FontSize = fontsize;
  hcb.FontSize = fontsize;
  
  if 1 % add Shue model    
    hold(hca,'on')
    xlim_tmp = hca.XLim;
    ylim_tmp = hca.YLim;
    plot(hca,x_mp,y_mp,'color','k','linewidth',1)
    plot(hca,x_bs,y_bs,'color','k','linewidth',1)
    hca.XLim = xlim_tmp;
    hca.YLim = ylim_tmp;
    hold(hca,'off')
  end
end

if 1 % dynamoic pressure
  hca = h(isub); isub = isub + 1;    
  var = db_table.('Pdyn_us');
  var(var < -1e20) = NaN;
  scatter(hca,xRE,yRE,15,var,'filled')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'Dynamic pressure';
  hcb.XLabel.Interpreter = 'none';
  hcb.Location = 'west'; 
  hcb.YLabel.Color = color_text;
  hcb.Color = color_text;
  hca.XLabel.String = 'x_{GSE} (R_E)';
  hca.YLabel.String = 'y_{GSE} (R_E)';
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  axis(hca,'equal')
  %colormap(hca,pic_colors('candy_gray'))
  colormap(hca,irf_colormap('waterfall'))
  hca.Color = color_bg;
  hca.CLim = prctile(var,[5 95]);
  hca.FontSize = fontsize;
  hcb.FontSize = fontsize;
  
  if 1 % add Shue model    
    hold(hca,'on')
    xlim_tmp = hca.XLim;
    ylim_tmp = hca.YLim;
    plot(hca,x_mp,y_mp,'color','k','linewidth',1)
    plot(hca,x_bs,y_bs,'color','k','linewidth',1)
    hca.XLim = xlim_tmp;
    hca.YLim = ylim_tmp;
    hold(hca,'off')
  end
end

c_eval('h(?).XLim = xlim; h(?).YLim = ylim;',1:numel(h))

compact_panels(h,0.02,0.01)

if 1 % Shrink colorbar
    hb = findobj(gcf,'type','colorbar');
    for ihb = 1:numel(hb)
      hb(ihb).Position(2) = hb(ihb).Position(2) + hb(ihb).Position(4)*0.25;
      hb(ihb).Position(4) = hb(ihb).Position(4)*0.5;  
      hb(ihb).Position(3) = hb(ihb).Position(3)*0.75;
    end
  end
%%
for ivar = 1:numel(varstrs)
  varstr = varstrs{ivar};
  var = db_table.(varstr);
  var(var < -1e20) = NaN;

  hca = h(isub); isub = isub + 1;
  scatter(hca,xRE,yRE,10,var,'filled')
  hcb = colorbar(hca);
  hcb.YLabel.String = varstr;
  hcb.XLabel.Interpreter = 'none';
  hcb.Location = 'west'; 
  hcb.YLabel.Color = color_text;
  hcb.Color = color_text;
  hca.XLabel.String = 'x_{GSE} (R_E)';
  hca.YLabel.String = 'y_{GSE} (R_E)';
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  axis(hca,'equal')
  colormap(hca,cmap)
  hca.Color = color_bg;
  hca.CLim = prctile(var,[5 95]);
  
  if 1 % Add KSR boxes
    Dp = 2;
    Bz = 0;
    fun_alpha = @(Bz,Dp) (0.58-0.007*Bz)*(1+0.024*log(Dp));
    alpha = fun_alpha(Bz,Dp);

    r_fun = @(rzero,theta) rzero.*(2./(1+cos(theta))).^alpha;
    %r_fun_scaling = @(theta) (2./(1+cos(theta))).^alpha;    

    % Binning specifications for each region
    deltaR = 0;
    % Magnetopause
    r_MP = [9.5 10.5] + deltaR;
    ca_MP = [-1 45]*pi/180;   

    % Magnetosheath
    r_MSH = [10.5 12.5] + deltaR;
    ca_MSH = [-1 45]*pi/180;

    % Bowshock
    r_BS = [12.5 13.5] + deltaR;
    ca_BS = [-1 45]*pi/180;

    % Foreshock
    r_FS = [13.5 17.0] + deltaR;
    ca_FS = [-1 45]*pi/180;
    
    hold(hca,'on')
    angle = [-45:45, 45:-1:-45 -45];
    r_ = r_MP;    
    rr = [repmat(r_(1),1,91), repmat(r_(2),1,91) r_(1)];
    rr_ = r_fun(rr,angle*pi/180);
    xx = rr_.*cosd(angle); yy = rr_.*sind(angle);
    plot(hca,xx,yy,'color',color_ksr)
    r_ = r_MSH;    
    rr = [repmat(r_(1),1,91), repmat(r_(2),1,91) r_(1)];
    rr_ = r_fun(rr,angle*pi/180);
    xx = rr_.*cosd(angle); yy = rr_.*sind(angle);
    plot(hca,xx,yy,'color',color_ksr)
    r_ = r_BS;    
    rr = [repmat(r_(1),1,91), repmat(r_(2),1,91) r_(1)];
    rr_ = r_fun(rr,angle*pi/180);
    xx = rr_.*cosd(angle); yy = rr_.*sind(angle);
    plot(hca,xx,yy,'color',color_ksr)
    r_ = r_FS;    
    rr = [repmat(r_(1),1,91), repmat(r_(2),1,91) r_(1)];
    rr_ = r_fun(rr,angle*pi/180);
    xx = rr_.*cosd(angle); yy = rr_.*sind(angle);
    plot(hca,xx,yy,'color',color_ksr)
    
    hold(hca,'off')
  end
  if 1 % Add circle at PO apogee
    hold(hca,'on')
    r_PO = 17;
    plot(hca,r_PO*sind(0:360),r_PO*cosd(0:360),'linestyle','-.','color',[0.5 0.5 0.5],'linewidth',1)
    hold(hca,'off')
  end
  
end

compact_panels(0.02,0.02)
linkprop(h,{'XLim','YLim'})
h(1).XLim = [-7 18];
h(1).YLim = [-1 1]*25*0.99;

ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
c_eval('ht(?).Color = color_text;',1:numel(ht))



if 0 % Shrink colorbar
  hb = findobj(gcf,'type','colorbar');
  for ihb = 1:numel(hb)
    hb(ihb).Position(2) = hb(ihb).Position(2) + hb(ihb).Position(4)*0.25;
    hb(ihb).Position(4) = hb(ihb).Position(4)*0.5;  
    hb(ihb).Position(3) = hb(ihb).Position(3)*0.5;
  end
end
