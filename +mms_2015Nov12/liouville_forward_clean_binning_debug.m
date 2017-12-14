%% Binning
% Decide binning, match with real observed bins if possible.
times = obsPDist.tlim(tintObs).time;
posN = (times-times(1) +- (times(end)-times(1))/2)*CS_normal_velocity;
posN(abs(posN)>1.1*limN*1e-3) = [];
dN = posN(2)-posN(1);
edgesN = [posN(1)-dN/2; posN+dN/2];
%edgesE = [0 f_center.depend{1}];
edgesE = [0 min(ePDist1.depend{1}(1:2,:))]; % bin according to lowest energy levels since most particles are there

edgesPolarAngle = 0:11.25:180;
edgesAzimuthalAngle = 0:11.25:360;

f_binned_all = zeros(numel(edgesN)-1,numel(edgesE)-1,numel(edgesAzimuthalAngle)-1,numel(edgesPolarAngle)-1);
f_occupied_all = zeros(numel(edgesN)-1,numel(edgesE)-1,numel(edgesAzimuthalAngle)-1,numel(edgesPolarAngle)-1);

allPa = nan(numel(saveParticle),1);
allF = nan(numel(saveParticle),1);
allE0 = zeros(numel(saveParticle),1);
allEend = zeros(numel(saveParticle),1);
allv0 = zeros(numel(saveParticle),1);
allTend = zeros(numel(saveParticle),1);
allzend = zeros(numel(saveParticle),1);
alla = zeros(numel(saveParticle),1);

disp('Binning particles')
fprintf('iP = %5.0f/%5.0f\n',0,nP)
tic;
for iP = (nP/2+1):nP %:numel(saveParticle) % step through particles
  if mod(iP,100) == 0, fprintf([repmat('\b', 1, 12) '%5.0f/%5.0f\n'],iP,nP); end % display progress  
  %if mod(iP,100) == 0, disp(sprintf('iParticle = %g/%g',iP,nParticles)); end
  % Pick out the data
  ind = 1:numel(saveParticle{iP}.t);
  z = saveParticle{iP}.r(ind,3); % N  
  allzend(iP) = z(end);
  allTend(iP) = saveParticle{iP}.t(end);
  pa = saveParticle{iP}.pa(ind); % pitch angle  
  
  
  v = saveParticle{iP}.v(ind,:); % m/s
  energy = saveParticle{iP}.energy;
  allE0(iP) = energy(1);
  allEend(iP) = energy(end);
  
  % make selections on which particles to include
  if 0%energy(1)>2000
    continue;
  end
  if z(1)<0
    continue;
  end
  % [count edges mid loc] = histcn(X, edge1, edge2, ..., edgeN)
  zbin = [25 30]*1e3;
  pabin = [75 105];
  [a b c d] = histcn([z, pa],zbin,pabin);
  if any(a>0)
    continue;
  end
  if 0%any(a>0)
    nrows = 4;
    ncols = 1;
    npanels = nrows*ncols;    
    for ipanel = 1:npanels
      h(ipanel) = subplot(nrows,ncols,ipanel);
    end
    difft =  diff(saveParticle{iP}.t);
    diffe =  diff(saveParticle{iP}.energy);
    isub = 1;
    if 1
      hca = h(isub); isub = isub + 1;
      hs = scatter(hca,z(2:end)*1e-3,pa(2:end),z(2:end)*0+5,diffe./difft);
      hca.Title.String = sprintf('iP = %g',iP);    
      hca.XLim = [-40 40];
      hca.YLim = [0 180];
      hca.XLabel.String = 'km';
      hold(hca,'on')
      hp = patch(hca,[zbin zbin(2) zbin(1)]*1e-3,[pabin(1) pabin(1) pabin(2) pabin(2)],mms_colors('x'));
      hp.FaceAlpha = 0.5;
      hold(hca,'off')
      hcb = colorbar('peer',hca);
      colormap(hca,'parula');
      hcb.YLabel.String = 'diff(energy)/diff(t)';
    end
    if 1
      hca = h(isub); isub = isub + 1;
      scatter(hca,saveParticle{iP}.t(2:end),saveParticle{iP}.energy(2:end),saveParticle{iP}.t(2:end)*0+2,diffe) 
      hca.Title.String = sprintf('iP = %g',iP);    
      hca.XLabel.String = 't';
    end
    pause
  end
  %alla(iP) = a;

  vtot = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2); % m/s
  allv0(iP) = vtot(1);
  normv = [v(:,1)./vtot v(:,2)./vtot v(:,3)./vtot];
  
  lmnxyz = saveParticle{iP}.r;
  xyzB = [Bx(lmnxyz(:,1),lmnxyz(:,2),lmnxyz(:,3)) By(lmnxyz(:,1),lmnxyz(:,2),lmnxyz(:,3)) Bz(lmnxyz(:,1),lmnxyz(:,2),lmnxyz(:,3))];
  xyzBnorm = xyzB./repmat(sqrt(sum(xyzB.^2,2)),1,3);
  oldxyz = saveParticle{iP}.r*inv(lmn'); % from LMN to dsl/gse xyz
  oldx = oldxyz(:,1);
  oldy = oldxyz(:,2);
  oldz = oldxyz(:,3);
  oldvxyz = normv*inv(lmn');  % from LMN to dsl/gse
  oldvx = oldvxyz(:,1);
  oldvy = oldvxyz(:,2);
  oldvz = oldvxyz(:,3); 
    
  polarAngle = acosd(-oldvz); % angle from z, oldvz is normalized and ind gse/dsl coordinates
  azimuthalAngle = atan2d(-oldvy,-oldvx); % angle is 0 if vx>0 (vy=0)
  % atan gives angles
  negangles = find(azimuthalAngle<0);
  azimuthalAngle(negangles) = 360+azimuthalAngle(negangles);
  
  % Sort into N and E bins
  [bins_occupied,~,mid,loc] = histcn([z energy azimuthalAngle polarAngle],edgesN*1e3,edgesE,edgesAzimuthalAngle,edgesPolarAngle);
  
  if 1%z(1)>0
  f_occupied = bins_occupied(1:(numel(edgesN)-1),1:32,1:32,1:16);
  f_occupied(f_occupied~=0) = 1;
  f_occupied_all = f_occupied_all + f_occupied;
  f_binned = bins_occupied; f_binned(f_binned~=0) = saveParticle{iP}.f;
  f_binned_all = f_binned_all + f_binned(1:42,1:32,1:32,1:16);
  
  allF(iP) = saveParticle{iP}.f;
  end
  if 0 % diagnostic plotting
    pitchangle = acosd(sum(xyzBnorm.*normv,2));
    
    figure(40)
    nrows = 6;
    ncols = 1;
    npanels = nrows*ncols;
    for ip = 1:npanels
      h(ip) = subplot(nrows,ncols,ip);
    end
    isub = 1;

    colors = mms_colors('matlab');
    if 1
      hca = h(isub); isub = isub + 1;  
      scatter(hca,lmnxyz(:,3)*1e-3,pitchangle)
      hca.XLabel.String = 'N';
      hca.YLabel.String = 'PA';      
      hca.YLim = [0 180];
    end    
    if 1
      hca = h(isub); isub = isub + 1;  
      ax = plotyy(hca,z,polarAngle,z,azimuthalAngle);      
      hca.XLabel.String = 'N';
      hca.YLabel.String = 'Angle';      
      legend('polar','azimuthal')
      ax(1).YLim = [0 180];
      ax(2).YLim = [0 360];
    end   
    if 0
      hca = h(isub); isub = isub + 1;  
      plot3(hca,lmnxyz(:,1)*1e-3,lmnxyz(:,2)*1e-3,lmnxyz(:,3)*1e-3)
      hold(hca,'on')
      plot3(hca,lmnxyz(1,1)*1e-3,lmnxyz(1,2)*1e-3,lmnxyz(1,3)*1e-3,'go')
      plot3(hca,lmnxyz(end,1)*1e-3,lmnxyz(end,2)*1e-3,lmnxyz(end,3)*1e-3,'rx')
      hold(hca,'off')      
      hca.XLabel.String = 'L';
      hca.YLabel.String = 'M';
      hca.ZLabel.String = 'N';  
      view(hca,[0 1 0]); 
      axis(hca,'equal')
    end  
    if 0
      hca = h(isub); isub = isub + 1;  
      plot3(hca,oldx*1e-3,oldy*1e-3,oldz*1e-3)
      hold(hca,'on')
      plot3(hca,oldx(1)*1e-3,oldy(1)*1e-3,oldz(1)*1e-3,'go')
      plot3(hca,oldx(end)*1e-3,oldy(end)*1e-3,oldz(end)*1e-3,'rx')
      hold(hca,'off')      
      hca.XLabel.String = 'X';
      hca.YLabel.String = 'Y';
      hca.ZLabel.String = 'Z';   
      view(hca,double(lmn(2,:)));
      axis(hca,'equal')
    end
    if 1 % pitch angle from tsPDist, psd
      midN_ = mid{1};
      midTime_ = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5;
      timesMap_ = midTime_ + midN_*1e-3/CS_normal_velocity;
      timesMap_ = timesMap_([1 end]);%+0.015*[-1 1];
      tsFmap_ = obsPDist.tlim(timesMap_([1 end])+0.015*[-1 1]); % need to add one sampling period to get the right number
      tsFmap_.data = f_binned;
      
      hca = h(isub); isub = isub + 1;
      elim = [0 1000];
      irf_spectrogram(hca,tsFmap_.pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa'),'log'); %hca.CLim = h1(isub-2).CLim;
    end
    if 1 % pitch angle from tsPDist, deflux
      midN_ = mid{1};
      midTime_ = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5;
      timesMap_ = midTime_ + midN_*1e-3/CS_normal_velocity;
      timesMap_ = timesMap_([1 end]);%+0.015*[-1 1];
      tsFmap_ = obsPDist.tlim(timesMap_([1 end])+0.015*[-1 1]); % need to add one sampling period to get the right number
      tsFmap_.data = f_binned;
      
      hca = h(isub); isub = isub + 1;
      elim = [0 1000];
      irf_spectrogram(hca,tsFmap_.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa'),'log'); %hca.CLim = h1(isub-2).CLim;
    end
    if 1      
      hca = h(isub); isub = isub + 1;      
      edgesPA = linspace(0,180,16);
      [bins_occupied,~,mid,loc] = histcn([z energy pitchangle],edgesN*1e3,edgesE,edgesPA);
      xx = mid{1};
      yy = mid{3};
      cc = squeeze(sum(bins_occupied,2))*saveParticle{iP}.f;
      pcolor(hca,xx,yy,log10(cc)');      
      hb = colorbar('peer',hca);
    end
  end
end
f_binned_all = f_binned_all./f_occupied_all;
toc
%f_binned_all(f_binned_all==0) = NaN;

% Making PDist from test particles
midN = mid{1};
midTime = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5;
timesMap = midTime + midN*1e-3/CS_normal_velocity;
timesMap_ = timesMap([1 end])+0.015*[-1 1];
%tsFmap = obsPDist.tlim(timesMap([1 end])+0.015*[-1 1]); % need to add one sampling period to get the right number
tsFmap.data = f_binned_all;
tsFmap = obsPDist.clone(timesMap,f_binned_all);
tsFmap.depend{1} = repmat(mid{2},tsFmap.length,1); % energy
tsFmap.depend{2} = repmat(mid{3},tsFmap.length,1); % azimuthal angle
tsFmap.depend{3} = mid{4};                         % polar angle
tsFmap.ancillary.energy0 = mid{2};
tsFmap.ancillary.energy1 = mid{2};
tsFmap.ancillary.energy = tsFmap.depend{1};
tsFmap.ancillary.esteptable = ones(tsFmap.length,1);

