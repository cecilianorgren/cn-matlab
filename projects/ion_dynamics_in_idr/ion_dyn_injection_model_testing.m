no02m = PIC('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%%
twpe = 20000;
pic = no02m.twpelim(20000).xlim([75 125]).zlim([-7 7]);

varstrs = {'Ez','vx(3)','vz(3)','vx(3)./vz(3)'; 'By','vx(5)','vz(5)','vx(5)./vz(5)'}';
clims = {[-1 1],[-1 1],[-1 1],0.5*[-1 1],[-1 1],[-1 1],[-1 1],0.5*[-1 1]};
pic.plot_map(varstrs,'A',1,'clim',clims)

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












