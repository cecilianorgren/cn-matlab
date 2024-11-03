time = time_xline + -6;
nMovMean = 7;
vel = mean(gseVi3.tlim(time + nMovMean*0.15*0.5*[-1 1]).data,1);
pressure = gsePi3.tlim(time + 0.15*0.5*[-1 1]);
pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).tlim(time + 0.15*0.5*[-1 1]).elim([300 Inf]);

T = squeeze(pressure.data);

% Rotate tensor twice, to find the most unequal component
theta1 = 0.5*atan((2*T(2,3))./(T(3,3)-T(2,2)));
theta2 = 0.5*atan((2*T(1,2))./(T(2,2)-T(1,1)));
%theta2 = 0.5*atan((2*T(2,1))./(T(1,1)-T(2,2)));


R1 = [1 0 0; 0 cos(theta1) -sin(theta1); 0 sin(theta1) cos(theta1)];
T1 = R1*(T*transpose(R1));

R2 = [cos(theta2) -sin(theta2) 0; sin(theta2) cos(theta2) 0; 0 0 1];
R2 = [cos(theta2) -sin(theta2) 0; sin(theta2) cos(theta2) 0; 0 0 1];
T2 = R2*(T*transpose(R2));



pdist_rot1 = pdist.shift(squeeze(vel), 1000, eye(3), 'mms');
pdist_rot2 = pdist.shift(squeeze(vel), 1000, R2, 'mms');
%function varargout = shift(pdist,v_nf,nMC,orient,sc,varargin)

nRows = 3;
nCols = 2;
h = setup_subplots(nRows,nCols,'vertical');
isub = 1;

if 1 % Original distribution
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',[1 0 0],[0 1 0]);
  vdf.plot_plane(hca);
  if 1
    hold(hca,'on')
    quiver(hca,2000*cos(-theta2),2000*sin(-theta2),'color','k','LineWidth',2)
    hold(hca,'off')
  end
  if 1
    hold(hca,'on')
    TT = T/trace(T);
    quiver(hca,0, 0, 5000*TT(1,1), 0, 'color','r','LineWidth',2)
    quiver(hca,0, 0, 0, 5000*TT(2,2), 'color','r','LineWidth',2)
    hold(hca,'off')
  end
  axis(hca,'equal')
  axis(hca,'square')
end

if 1 % Original distribution, shifted to bulk speed frame
  hca = h(isub); isub = isub + 1;
  vdf = pdist_rot1.reduce('2D',[1 0 0],[0 1 0]);
  vdf.plot_plane(hca);
  if 1
    hold(hca,'on')
    quiver(hca,0,0,2000*cos(-theta2),2000*sin(-theta2),'color','k','LineWidth',2)
    hold(hca,'off')
  end
  axis(hca,'equal')
  axis(hca,'square')
end

if 1 % Rotated distribution
  hca = h(isub); isub = isub + 1;
  vdf = pdist_rot2.reduce('2D',[1 0 0],[0 1 0]);
  vdf.plot_plane(hca);
  if 1
    hold(hca,'on')
    quiver(hca,0,0,2000*cos(-theta2+theta2),2000*sin(-theta2+theta2),'color','k','LineWidth',2)
    hold(hca,'off')
  end
  if 1
    hold(hca,'on')
    TT = T2/trace(T2);
    quiver(hca,0, 0, 5000*TT(1,1), 0, 'color','r','LineWidth',2)
    quiver(hca,0, 0, 0, 5000*TT(2,2), 'color','r','LineWidth',2)
    hold(hca,'off')
  end
  axis(hca,'equal')
  axis(hca,'square')
end

if 1 % Original distribution
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('1D',[0 1 0]);
  plot(hca,vdf.depend{1},vdf.data);  
  hca.XGrid = 'on'; hca.YGrid = 'on';
end
if 1 % Original distribution
  hca = h(isub); isub = isub + 1;
  vdf = pdist_rot1.reduce('1D',[0 1 0]);
  plot(hca,vdf.depend{1},vdf.data);  
  hca.XGrid = 'on'; hca.YGrid = 'on';
end
if 1 % Original distribution
  hca = h(isub); isub = isub + 1;
  vdf = pdist_rot2.reduce('1D',[0 1 0]);
  plot(hca,vdf.depend{1},vdf.data);  
  hca.XGrid = 'on'; hca.YGrid = 'on';
end

hlinks = linkprop(h,{'CLim'});







