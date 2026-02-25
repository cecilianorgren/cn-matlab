h = setup_subplots(1,1);
npanels = numel(h);

xlim = [00 25]; %xlim = [-5 40];%[x(1) x(end)] + 150*[1 -1];
zlim = [-5 -0]; %zlim = [-10 10];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;
    
% Flux function
doAx = 1; % plot separatrix
doA = 1;
cA = 0*[0.8 0.8 0.8];
nA = 20;
nA = [0:-0.5:min(A(:))];
ipxA = ipx1:20:ipx2;
ipzA = ipz1:20:ipz2;

isub = 1;
if 1 % E.par
  hca = h(isub); isub = isub + 1;
  %variable = eval('E.par');  
  variable = eval('log10(ne12)');  
  [X,Z] = ndgrid(x(ipx),z(ipz));
  himag = surf(hca,X,X*0,Z,squeeze(variable(ipx,ipz)));
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'log_{10}(n_e/n_0)';
  hcb(isub-1) = hb;
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'y (d_i)';
  hca.ZLabel.String = 'z (d_i)';
  hca.CLim = [-1 1];
  colormap(hca,pic_colors('candy'))
  hca.CLim = [-2 0];
  
  hca.Visible = 'off';  
  hca.XAxis.Visible = 'on';
  hca.YAxis.Visible = 'off';
  hca.ZAxis.Visible = 'on';
  hca.XAxisLocation = 'top';
  hca.XAxis.FirstCrossoverValue  = hca.YLim(2);
  hca.XAxis.SecondCrossoverValue = hca.XLim(1);
  hca.Position(1) = 0.15;  
  hca.Position(3) = 0.75;  
  hca.Position(4) = 0.75;  
  axis(hca,'equal')
  set(hca,'view',[-41.9000,6.8000])
  
  if doA
    hold(hca,'on')
    if 1
      [cline_x,cline_y] = contourc_split(x(ipx),z(ipz),squeeze(A(ipx,ipz))',nA);
      nLines = numel(cline_x);
      for iLine = 1:nLines
        plot3(hca,cline_x{iLine},cline_x{iLine}*0,cline_y{iLine},'k')
      end
    else
      cline = contourc(x(ipx),z(ipz),squeeze(A(ipx,ipz))',nA); % ,'color',cA,'linewidth',0.5,'displayname','A'    
      cline_tmp = cline;
      for iA = 1:numel(nA)
        valA = cline_tmp(1,1);
        nSegA = cline_tmp(2,1);
        x_tmp = cline_tmp(1,1+(1:nSegA));
        y_tmp = x_tmp*0;
        z_tmp = cline_tmp(2,1+(1:nSegA));
        plot3(hca,x_tmp,y_tmp,z_tmp,'k')
        cline_tmp(:,1:(nSegA+1)) = [];
        size(cline_tmp,2);
      end
    end
    hold(hca,'off')  
  end
  if doAx
    hold(hca,'on')
    [saddle_locations,saddle_values] = saddle(A,'sort');
    if 1
      [cline_x,cline_y] = contourc_split(x(ipx),z(ipz),squeeze(A(ipx,ipz))',saddle_values(1)*[1 1]);
      nLines = numel(cline_x);
      linewidth_tmp = 1.0;
      dyA = 4;
      for iLine = 1:nLines        
        plot3(hca,cline_x{iLine},cline_x{iLine}*0,cline_y{iLine},'color','k','linewidth',2)
        plot3(hca,cline_x{iLine},cline_x{iLine}*0 - dyA*0.1,cline_y{iLine},'color','k','linewidth',linewidth_tmp)
        plot3(hca,cline_x{iLine},cline_x{iLine}*0 - dyA*0.2,cline_y{iLine},'color','k','linewidth',linewidth_tmp)
        plot3(hca,cline_x{iLine},cline_x{iLine}*0 - dyA*0.3,cline_y{iLine},'color','k','linewidth',linewidth_tmp)
        plot3(hca,cline_x{iLine},cline_x{iLine}*0 - dyA*0.4,cline_y{iLine},'color','k','linewidth',linewidth_tmp)        
      end      
    else
      hcont = contour3(hca,X,X*0,Z,squeeze(A(ipx,ipz)),saddle_values(1)*[1 1],'color',cA,'linewidth',2,'displayname','A_X','linestyle','-'); 
    end
    hold(hca,'off')  
  end
end

for ipanel = 1:npanels  
  h(ipanel).FontSize = 14;
  h(ipanel).YDir = 'normal';
  h(ipanel).XDir = 'reverse';'normal';
  h(ipanel).XLim = xlim;
  h(ipanel).YLim = [-2 0];
  h(ipanel).ZLim = zlim;  
end
