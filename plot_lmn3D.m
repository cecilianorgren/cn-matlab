function plot_lmn3D(varargin)

[ax,args,nargs] = axescheck(varargin{:});
delete(ax.Children)
R1 = args{1};
R2 = args{2};
R3 = args{3};
R4 = args{4};
xyz = args{5}; X = xyz(1,:); Y = xyz(2,:); Z = xyz(3,:);
coord_system = args{6};
args = args(7:end);

have_options = 0;
have_vectors = 0;

if nargs > 6, have_options = 1; end
while have_options
  l = 1;
  switch(lower(args{1}))   
    case 'vectors'
      l = 2;
      vectors = args{2};
      have_vectors = 1;
  end
  args = args(l+1:end);  
  if isempty(args), break, end    
end

time = R1.time(1);

%c_eval('mvaR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(L).data gseR?.dot(M).data gseR?.dot(N).data]);')
%c_eval('lmnR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(L).data gseR?.dot(M).data gseR?.dot(N).data]); coord_system = {''L'',''M'',''N''};')
%c_eval('nmlR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(N).data gseR?.dot(-M).data gseR?.dot(L).data]); coord_system = {''N'',''-M'',''L''};')
%c_eval('mvaR? = nmlR?;')
%c_eval('mvaR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(L).data gseR?.dot(M).data gseR?.dot(N).data]);')

mms_marker={{'ks','markersize',12},{'rd','markersize',12},...
	{'go','markersize',12,'color',[0 0.6 0]},{'bv','markersize',12}};
mms_marker_small={{'ks','markersize',8},{'rd','markersize',8},...
	{'go','markersize',8,'color',[0 0.6 0]},{'bv','markersize',8}};
mms_marker_shaded={{'ks','color',[0.3 0.3 0.3]},...
	{'rd','color',[1 0.3 0.3]},{'go','color',[.3 1 .3]},{'bv','color',[.3 .3 1]}};
sc_list = 1:4;

%x = {mvaR1.resample(time).data,mvaR2.resample(time).data,mvaR3.resample(time).data,mvaR4.resample(time).data};

drref=0;
R0 = 0;
for ic = sc_list
  c_eval('rr{ic}=R?.resample(time).data;',ic)
  R0=R0+rr{ic}/4;
end
for ic = sc_list
  x{ic}=rr{ic}-R0;
  %x{ic}(1)=data.t;
  
  absx=irf_abs(x{ic});
  absx = absx(end);
  drref=max([drref absx]);				
end
for ic = sc_list
  xOld = x;
  x{ic} = [x{ic}*X' x{ic}*Y' x{ic}*Z'];
end
  

%delete h
if isempty(ax)
  h(1)=axes('position',[0.15  0.16 0.7 0.7]); % [x y dx dy]  
else
  h = ax;
end

for ic=sc_list
  % put Cluster markers
  plot3(h(1),x{ic}(1),x{ic}(2),x{ic}(3),mms_marker{ic}{:});
  hold(h(1),'on');
  % lines from sc to projection planes
  lineProperties = {'parent',h(1),'linestyle',':','linewidth',0.6};
  line([x{ic}(1) -drref],[x{ic}(2) x{ic}(2)],[x{ic}(3) x{ic}(3)],lineProperties{:});
  line([x{ic}(1) x{ic}(1)],[x{ic}(2) -drref],[x{ic}(3) x{ic}(3)],lineProperties{:});
  line([x{ic}(1) x{ic}(1)],[x{ic}(2) x{ic}(2)],[x{ic}(3) -drref],lineProperties{:});
  % put Cluster markers on projections
  plot3(h(1),-drref,x{ic}(2),x{ic}(3),mms_marker_shaded{ic}{:});
  plot3(h(1),x{ic}(1),-drref,x{ic}(3),mms_marker_shaded{ic}{:});
  plot3(h(1),x{ic}(1),x{ic}(2),-drref,mms_marker_shaded{ic}{:});
end

axis(h(1),[-drref drref -drref drref -drref drref ]);
for ii=1:4,
  for jj=ii+1:4,
    if any(find(sc_list==ii)) && any(find(sc_list==jj)),
      line([x{ii}(1) x{jj}(1)],...
        [x{ii}(2) x{jj}(2)],...
        [x{ii}(3) x{jj}(3)],...
        'parent',h(1),'linewidth',2,'linestyle','-',...
        'color',[0.6 0.6 0.6]);
    end
  end
end
%

if 0 % draw origo axis
  line([-drref drref],[-drref -drref],[-drref -drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
  line([-drref -drref],[-drref drref],[-drref -drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
  line([-drref -drref],[-drref -drref],[-drref drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
end
text(0.1,1,0,time.utc,'parent',h(1),'units','normalized','horizontalalignment','center','fontsize',9);
xlabel(h(1),['{\Delta}' coord_system{1} ' [km] ' ]);
ylabel(h(1),['{\Delta}' coord_system{2} ' [km] ' ]);
zlabel(h(1),['{\Delta}' coord_system{3} ' [km] ' ]);
set(h(1),'xdir','reverse');
set(h(1),'ydir','reverse');
grid(h(1),'on');
axis(h(1),[-drref drref -drref drref]);
hold(h(1),'off');
axis equal

% plot vectors
% Plot vectors
xLims = h.XLim; yLims = h.YLim; zLims = h.ZLim;
vectorColors = {[0 0 0],mms_colors('b'),mms_colors('2'),mms_colors('3'),mms_colors('4')};
while have_vectors  
  hold(h,'on');
  if isempty(vectors), break, end 
  vecnorm = vectors{1,3};
  for ic = 1:4;
    vector = vectors{1,2}{ic};
    if isa(vector,'TSeries')
      vector = vector.data;
    end
    vector = [vector*X' vector*Y' vector*Z']/vecnorm;
    vecTxt = vectors{1,1};
    quiver3(h,x{ic}(1),x{ic}(2),x{ic}(3),vector(1),vector(2),vector(3),'linewidth',1,'linestyle','-','color',vectorColors{1})
  end     
  vectorColors = vectorColors(2:end);
  vectors = vectors(2:end,:);
  hold(h,'off');
end

h.XLim = xLims; h.YLim = yLims; h.ZLim = zLims;

%rotate(h,[180 0],00)

%irf_legend(h,{'MMS1','MMS2','MMS3','MMS4'})
mainAxesPosition = h.Position;
if 0
  %%
labelAxisPosition = mainAxesPosition+[0.15  0.16 0.7 0.7]./[0.5 0.8 0.5 0.2];

  
hca = axes('Position',labelAxesPosition);
hold(hca,'on');
axis(hca,[0 1 0 1]);
yy=.5;dxx=.15;xs=.3;
plot(hca,xs,yy,'ks',xs+1*dxx,yy,'rd',xs+2*dxx,yy,'go',xs+3*dxx,yy,'bv','LineWidth',1.5);
text(xs+0.03,yy,'C1','parent',hca);
text(xs+1*dxx+0.03,yy,'C2','parent',hca);
text(xs+2*dxx+0.03,yy,'C3','parent',hca);
text(xs+3*dxx+0.03,yy,'C4','parent',hca);
axis(hca,'off');    
end