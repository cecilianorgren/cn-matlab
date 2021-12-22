function out = funFitVDF(vdf,nPop)
% Fit Maxwellians to velocity distribution function.

% Collect input (for now, only input is PDist)
nTimes = vdf.length;
datasize = vdf.datasize;
nDim = numel(datasize)-1;
f = vdf.data; % matrix of size [nTimes datasize]
v = vdf.depend; % cell array, km/s
for iDim = 1:nDim
  v{iDim} = v{iDim}*1e3; % m/s
end

% Hardcoded parameers that should be available as input
% weight of cost function, i.e. how much weight should be put on the
% different moments of f: v^0*f, v^1*f, v^2*f
% For 1D, weight = [a b c], means cost function is
% CF = a*CF(v^0*f) + b*CF(v^1*f) + c*CF(v^2*f)
weight = ones(1,nPop*(1+nDim*2));
doPlot = 1;
iterPlot = 100;

% How should we treat the initial guess? For now, maybe just choose some
% reasonable values and see how it works.
units = irf_units;
n0 = 0.01e6; % m^-3
vd0 = 1000e3; % m/s
T0 = 100; % eV
vt0 = sqrt(T0*units.eV*2/units.me); % m/s


% Set initial guess
switch nDim
  case 1 % 1D
    X0 = [n0 vd0 vt0];    
    % Apply to X so that we can return data in more conventional units, cc,
    % km, eV.
    X_units = repmat({'m-6','m/s','m/s'},1,nPop);
    funits = @(X) reshape([X(1:3:end)*1e-6; X(2:3:end)*1e-3; X(3:3:end).^2*units.me/2/units.eV],1,numel(X));
    Xu_units = repmat({'cm-3','km/s','eV'},1,nPop);
  case 2 % 2D
    X0 = [n0 vd0 vd0 vt0 vt0];
    X_units = repmat({'m-6','m/s','m/s','m/s','m/s'},1,nPop);
    funits = @(X) reshape([X(1:5:end)*1e-6; X(2:5:end)*1e-3; X(3:5:end)*1e-3; X(4:5:end).^2*units.me/2/units.eV X(5:5:end).^2*units.me/2/units.eV],1,numel(X));
    Xu_units = repmat({'cm-3','km/s','km/s','eV','eV'},1,nPop);
  case 3 % 3D
    X0 = [n0 vd0 vd0 vd0 vt0 vt0 vt0];
end

% Repeat X0 for each population: [n1, vd1, vt1, ..., nN, vdN vtN]
X0 = repmat(X0,1,nPop);

% Initialize arrays to save the data
data = struct([]);
X_all = cell(nTimes,1);
Q_all = cell(nTimes,1); % quality parameters

% fminsearch uses the Nelder-Mead simplex (direct search) method.
% I'm not sure if the large difference between parameters affects the 
% method. For example, vt/n = 2e2. If the method varies each parameter 
% which the same "step", vt will not be varied much...
options = optimset('OutputFcn', @myoutput);

% Initial guess, leter one, the initial guess is the previous/neighbouring function
X = X0;
% Step through all the distribution functions
for iTime = 1:nTimes
  history = [];
  % Option to plot output. 2-D plots for each time step is consuming quite a
  % lot of resources and is only good or initial checks. 
  switch nDim
    case 1
      ftmp = f(iTime,:);
      vtmp = {v{1}(iTime,:)}; % use cellfun?
    case 2
      ftmp = squeeze(reshape(f(iTime,:),datasize(2),datasize(3)));
      vtmp = {v{1}(iTime,:),v{2}(iTime,:)};
    case 3 % not implemented      
  end
  
  cost_function = @(X) costfunction_maxwellian(X,vtmp,ftmp,nPop);
  try
  [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(cost_function,X,options);
  catch
    1;
  end
  data(iTime).nDim = nDim;
  data(iTime).nPop = nPop;
  data(iTime).X_units = X_units;
  data(iTime).X = X;  
  data(iTime).Xu_units = Xu_units;
  data(iTime).Xu = funits(X);
  data(iTime).history = history;
end
% Collect output
out = data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Auxiliary/help functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stop = myoutput(x,optimvalues,state)
  % Collect data
  stop = false;
  if isequal(state,'iter')
    history = [history; x];    
  end
end 
function out = costfunction_maxwellian(X,v,f_obs,nPop)
% Cost function for bi-Maxwellian
% Params should be ordered

% Can params be 2D array? If yes, then it could have ndim2 = # maxwellians.
% What do we need. 
% Specify grid? Yes, vx, vy, vz are 1-D arrays which corresponds to the
% different bins in f_obs.
% What we need to compare is then datafit = f(vx,vy,vz) and data_obs

% Specify fit function, this is affects the number of inputs and is 
% specified in calling function.
% For each separate population nInp = nPop*(1 + nDim*2).
% E.g. nDim = 1: n, v,                    T
%      nDim = 2: n, vpar, vperp,          Tpar, Tperp
%      nDim = 3: n, vpar, vperp1, vperp2, Tpar, Tperp1, Tperp2

units = irf_units;
nInp = numel(X);
nDim = numel(v);

fmax = maxwellian(nDim);

f_fit = zeros(size(f_obs));
% Evaluate the fit function for current parameters
switch nDim
  case 1
    f_fit_separate = zeros(nPop,numel(f_obs));
    for iPop = 1:nPop
      iRef = nInp/nPop*(iPop-1);
      f_tmp = fmax(v{1},X(iRef+1),X(iRef+2),X(iRef+3));
      f_fit_separate(iPop,:) = f_tmp;
      f_fit = f_fit + f_tmp;
    end
  case 2
    [VX,VY] = meshgrid(v{1},v{2});
    f_fit_separate = zeros([nPop,size(f_obs)]);
    for iPop = 1:nPop
      iRef = nInp/nPop*(iPop-1);
      f_tmp = fmax(VX,VY,X(iRef+1),X(iRef+2),X(iRef+3),X(iRef+4),X(iRef+5));
      f_fit_separate(iPop,:,:) = f_tmp;
      f_fit = f_fit + f_tmp;
    end
  case 3 % NOT implemented
    for iInp = 1:nInp
      f_fit = fmax(v{1},X(iInp+0),X(iInp+2),X(iInp+3));
    end
end

% Evaluate f with current parameters. For each search step, this evaluation
% is done for different n, vd, T etc.
switch nDim
  case 1
    % proportional to density
    data_fit_mom0 = tocolumn(v{1}).^0.*tocolumn(f_fit(:));
    data_obs_mom0 = tocolumn(v{1}).^0.*tocolumn(f_obs(:));
    % proportional to speed
    data_fit_mom1 = tocolumn(v{1}).^1.*tocolumn(f_fit(:));
    data_obs_mom1 = tocolumn(v{1}).^1.*tocolumn(f_obs(:));
    % proportional to temperature
    data_fit_mom2 = tocolumn(v{1}).^2.*tocolumn(f_fit(:));
    data_obs_mom2 = tocolumn(v{1}).^2.*tocolumn(f_obs(:));
    
    out_mom0 = sum((data_fit_mom0-data_obs_mom0).^2)/sum(data_obs_mom0.^2);    
    out_mom1 = sum((data_fit_mom1-data_obs_mom1).^2)/sum(data_obs_mom1.^2);    
    out_mom2 = sum((data_fit_mom2-data_obs_mom2).^2)/sum(data_obs_mom2.^2);    
    
    % In case of negative density, is it working?
    out_mom0(data_fit_mom0 < 0) = Inf;
    
    % Final costfunction can be based on one or a combination of the data
    % above.
    %weight = [1 1 1];
    out = weight(1)*out_mom0 + weight(2)*out_mom1 + weight(3)*out_mom2;
    
    if doPlot && mod(iterPlot,size(history,1)) == 0
      nrows = 3; ncols = 1; isub = 1;
      
      % Make string to display parameters
      switch nPop
        case 1
          str_param = sprintf('n = %g cc, vd = %g km/s, vt = %g km/s',X(1)*1e-6,X(2)*1e-3,X(3)*1e-3);  
        case 2
          str_param = sprintf('n = [%g, %g] cc, vd = [%g, %g] km/s, vt = [%g, %g] km/s',X([1 4])*1e-6,X([2 5])*1e-3,X([3 6])*1e-3);
        case 3
          str_param = sprintf('n = [%g, %g, %g] cc, vd = [%g, %g, %g] km/s, vt = [%g, %g, %g] km/s',X([1 4 7])*1e-6,X([2 5 8])*1e-3,X([3 6 9])*1e-3);
        otherwise
          str_param = '';
      end
      
      % Plot
      hca = subplot(nrows,ncols,isub); isub = isub + 1;    
      plot(hca,v{1},data_obs_mom0,v{1},data_fit_mom0,v{:},f_fit_separate)
      hca.Title.String = {...
        str_param,...
        'Cost function (to be minimized):',...
        sprintf(' CF = %g x CF0 + %g x CF1 + %g x CF2 = %g',weight(1),weight(2),weight(3),out),...
        sprintf('(f) CF0 = %g',out_mom0)
        };

      hca = subplot(nrows,ncols,isub); isub = isub + 1;    
      plot(hca,v{1},data_obs_mom1,v{1},data_fit_mom1)
      hca.Title.String = {sprintf('(vf) CF1 = %g',out_mom1)};

      hca = subplot(nrows,ncols,isub); isub = isub + 1;    
      plot(hca,v{1},data_obs_mom2,v{1},data_fit_mom2)
      hca.Title.String = {sprintf('(v^2f) CF2 = %g',out_mom2)};  

      drawnow
      %pause(0.1)
    end

  case 2
    V = sqrt(VX.^2 + VY.^2);
    % proportional to density
    data_fit_mom0 = V.^0.*f_fit;
    data_obs_mom0 = V.^0.*f_obs;
    % proportional to speed
    data_fit_mom1x = VX.^1.*f_fit;
    data_obs_mom1x = VX.^1.*f_obs;
    data_fit_mom1y = VY.^1.*f_fit;
    data_obs_mom1y = VY.^1.*f_obs;
    % proportional to temperature
    data_fit_mom2x = VX.^2.*f_fit;
    data_obs_mom2x = VX.^2.*f_obs;
    data_fit_mom2y = VY.^2.*f_fit;
    data_obs_mom2y = VY.^2.*f_obs;
    
    out_mom0 = sum((data_fit_mom0(:)-data_obs_mom0(:)).^2)/sum(data_obs_mom0(:).^2);
    out_mom1x = sum((data_fit_mom1x(:)-data_obs_mom1x(:)).^2)/sum(data_obs_mom1x(:).^2);
    out_mom1y = sum((data_fit_mom1y(:)-data_obs_mom1y(:)).^2)/sum(data_obs_mom1y(:).^2);
    out_mom2x = sum((data_fit_mom2x(:)-data_obs_mom2x(:)).^2)/sum(data_obs_mom2x(:).^2);
    out_mom2y = sum((data_fit_mom2y(:)-data_obs_mom2y(:)).^2)/sum(data_obs_mom2y(:).^2);
    
    % In case of negative density
    out_mom0(data_fit_mom0 < 0) = Inf;
    
    % Final costfunction can be based on one or a combination of the data
    % above.
    %weight = [1 1 1 1 1];
    out = weight(1)*out_mom0 + weight(2)*out_mom1x + weight(3)*out_mom1y + weight(4)*out_mom2x + weight(5)*out_mom2y;
    
    if doPlot  
      switch nPop
        case 1
          str_param = sprintf('n = %g cc, vd1 = %g km/s, vd2 = %g km/s, vt1 = %g km/s, vt2 = %g km/s',X(1)*1e-6,X(2)*1e-3,X(3)*1e-3,X(4)*1e-3,X(5)*1e-3);  
        case 2
          str_param = sprintf('n = [%g, %g] cc, vd = [%g, %g] km/s, vd1 = [%g, %g] km/s, vt1 = [%g, %g] km/s, vt2 = [%g, %g] km/s',X([1 6])*1e-6,X([2 7])*1e-3,X([3 8])*1e-3,X([4 9])*1e-3,X([5 10])*1e-3);
          str_param = sprintf('n = [%g, %g] cc \nvd = [%g, %g] km/s \nvd1 = [%g, %g] km/s \nvt1 = [%g, %g] km/s \nvt2 = [%g, %g] km/s',X([1 6])*1e-6,X([2 7])*1e-3,X([3 8])*1e-3,X([4 9])*1e-3,X([5 10])*1e-3);
        case 3
          str_param = sprintf('n = [%g, %g, %g] cc, vd = [%g, %g, %g] km/s, vt = [%g, %g, %g] km/s',X([1 4 7])*1e-6,X([2 5 8])*1e-3,X([3 6 9])*1e-3);        
        otherwise
          str_param = '';
      end
      
      nrows = 3; ncols = 2+1; isub = 1;
      h = setup_subplots(nrows,ncols,'vertical');

      hca = h(isub); isub = isub + 1;
      plot(hca,nan,nan)
      str_print = sprintf('Cost function (to be minimized): \nCF = %g x CF0 + \n%g x CF1x + %g x CF1y + \n%g x CF2x + %g x CF2y \n= %g',weight(1),weight(2),weight(3),weight(4),weight(5),out);
      hleg = irf_legend(hca,str_print,[0.98 0.98]);
      hleg.HorizontalAlignment = 'right';
      hca.Visible = 'off';
      
      hca = h(isub); isub = isub + 1;
      plot(hca,nan,nan)
      irf_legend(hca,str_param,[0.02 0.98])  
      hca.Visible = 'off';
      
      hca = h(isub); isub = isub + 1;
      hca.Visible = 'off';
      
      hca = h(isub); isub = isub + 1;
      pcolor(hca,VX,VY,data_obs_mom0)
      shading(hca,'flat')
      hca.Title.String = {sprintf('(f) CF0 = %g',out_mom0)};
      hcb = colorbar('peer',hca);
      
      hca = h(isub); isub = isub + 1;
      pcolor(hca,VX,VY,data_obs_mom1y)
      shading(hca,'flat')      
      hca.Title.String = {sprintf('(vf) CF1 = %g',out_mom1y)};  
      hcb = colorbar('peer',hca);
      
      hca = h(isub); isub = isub + 1;
      pcolor(hca,VX,VY,data_obs_mom2y)
      shading(hca,'flat')      
      hca.Title.String = {sprintf('(v^2f) CF2 = %g',out_mom2y)};  
      hcb = colorbar('peer',hca);
      

      hca = h(isub); isub = isub + 1;
      pcolor(hca,VX,VY,data_fit_mom0)
      shading(hca,'flat')
      hca.Title.String = {sprintf('(f) CF0 = %g',out_mom0)};
      hcb = colorbar('peer',hca);
      hca.CLim = h(3+1).CLim;
           
      hca = h(isub); isub = isub + 1;
      pcolor(hca,VX,VY,data_fit_mom1y)
      shading(hca,'flat')      
      hca.Title.String = {sprintf('(vf) CF1 = %g',out_mom1y)};  
      hcb = colorbar('peer',hca);
      hca.CLim = h(3+2).CLim;
            
      hca = h(isub); isub = isub + 1;
      pcolor(hca,VX,VY,data_fit_mom2y)
      shading(hca,'flat')      
      hca.Title.String = {sprintf('(v^2f) CF2 = %g',out_mom2y)};  
      hcb = colorbar('peer',hca);
      hca.CLim = h(3+3).CLim;

      drawnow
      pause(0.1)
    end

end

%doPlot = 1;
end
function out = maxwellian(nDim)
  % Return Maxwellian function, 1D, 2D
  switch nDim
    case 1
      % 1-D Maxwellian function
      fmax = @(v,n,vd,vt) n.*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
    case 2
      % 2-D Maxwellian function
      fmax = @(vx,vy,n,vdx,vdy,vtx,vty) (1/n).*fmax1D(vx,n,vdx,vtx).*fmax1D(vy,n,vdy,vty);
    case 3 % not implemented
  end   
  out = fmax;
end
end