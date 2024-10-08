function varargout = ph_fit_kappa(vdf,varargin)
% Fit Maxwellians to velocity distribution function.
% Currently only "fully" implemented for 1D distributions.
%   data = funFitVDF(f1D);
%   [data,timeseries] = funFitVDF(f1D);
% 
%   Input:
%     vdf - PDist object of type 'line (reduced)', or '1DCart'
%     'nPop',__ - Number of Maxwellians used for fitting, default is 2.
%     'plot',0/1 - Track progress and fitting.
%     'plotIter',__ - Option to plot progress/state every __th iteration.
%                     Default is 2000, based on maxIter = 2000*nInp.
%     'weight',__ - Specify weight for cost function, must have
%                   nPop*(2*nDim+1) number of elements.
%     'maxIter',__ - Maximum number of iterations per distribution, passed
%                    as input argument to fminsearch.
%     'maxFunEvals',__ - Maximum number of function evaluations per 
%                        distribution (iteration?), passed as input 
%                        argument to fminsearch.
%     'iSortPop',__ - Option to order output by certain parameter. For
%                     example iSortPop = 1 would order them be ascending
%                     density, iSortPop = 2 would order them by ascending
%                     speed vx, etc. This prevents populations form
%                     arbitrarily switch between each other. Default is 2,
%                     sort by first speed component.
%     'X0',__ - Initial guess for fminsearch, must have nPop*(2*nDim+1) 
%               number of elements. Units are ,m-3 for density and m/s for
%               speed and thermal speed, e.g.
%               X0 = repmat([0.02e6 1000e3 10000e3],1,nPop);
%     'guessprevious',0/1 - Use previous result as input. This speeds up
%                           the process but might result in some population
%                           "loosing its way", for example always having
%                           zero density.
%
%   Output:
%     data - structure array with vdf.length entries and field:
%       nRedo      - 
%       iterations - number of iterations
%       exitflag   - exitflag, 0 if max number of iterations or function
%                    evaluations wer reached (i think)
%       funcCount  - number of function evaluations
%       cost       - final value of costfunction
%       nDim       - dimensions of data
%       nPop       - number of populations used for fitting, the number of
%                    inputs required for weight and intial guess is
%                    calculated as nInp = nPop*(1+2*nDim)
%       X_units    - units of output X (m-3, m/s ,m/s)
%       X          - final values of fit parameters, main result
%       Xu_units   - units of output X (cm-3, km/s, eV)
%       Xu         - final values of fit parameters, main result in easier units                    
%       history    - fit parameters for all iterations
%       f          - resulting fit in same dimensions as input data
%       f_separate - resulting fit in same dimensions as input data,
%                    divided into each (nPop) separate population 
%     
%     timeseries - structure with each element being a timeseries
%            n: density of each population (cm-3)
%           vd: bulkd drift speed of each population (km/s)
%            T: temperature of each population (eV)
%            f: PDist of best fit vdf
%           cf: cost function
%     exitflag: exitflag
%         iter: number of iterations
%        feval: number of function evaluations
%
%   Example:
%     % Default options
%     [fitdata,ts] = funFitVDF(vdf);
%     
%     % Some options given
%     nPop = 2;
%     X0 = repmat([0.02e6 1000e3 10000e3],1,nPop);
%     [fitdata,ts] = funFitVDF(vdf,'nPop',nPop,'plot',1,'guessprevious',1);
%
%
%     h = irf_plot(5);
%     hca = irf_panel('vdf_obs');
%     irf_spectrogram(hca,vdf.specrec('velocity'),'lin')
%     hca = irf_panel('vdf_fit');
%     irf_spectrogram(hca,ts.f.specrec('velocity'),'lin')
%     hold(hca,'on'); irf_plot(hca,ts.vd); hold(hca,'off'); 
%     hca = irf_panel('n_fit');
%     irf_plot(hca,ts.n)
%     hca = irf_panel('vd_fit');
%     irf_plot(hca,ts.vd)
%     hca = irf_panel('T_fit');
%     irf_plot(hca,ts.T)
%     irf_plot_axis_align(h)
%     hlinks = linkprop([irf_panel('vdf_obs'),irf_panel('vdf_fit')],{'CLim','YLim'});
%
% See also FMINSEARCH, OPTIMSET

units = irf_units;

% Collect input (for now, only input is PDist)
nTimes = vdf.length;
datasize = vdf.datasize;
nDim = numel(datasize)-1;
f = vdf.data; % matrix of size [nTimes datasize]
v = vdf.depend{1}; % array, km/s
v = v; % m/s

if strcmp(vdf.species,'ions'); m = units.mp; else, m = units.me; end

% Hardcoded parameers that should be available as input
% weight of cost function, i.e. how much weight should be put on the
% different moments of f: v^0*f, v^1*f, v^2*f
% For 1D, weight = [a b c], means cost function is
% CF = a*CF(v^0*f) + b*CF(v^1*f) + c*CF(v^2*f)
doGuessPrevious = 1;
weight = [];
doPlot = 0;
iterPlot = 2000; % Plot every nth iteration
nPop = 2; % default
maxFunEvals = [];
maxIter = [];
X0 = [];
X = [1 1 1e3];
method = 'fminsearch'; % fminbnd
maxRedo = 0; % if EXITFLAG == 0, we redo fminsearch with current X
doSort = 0;
iSortPop = 2; % sort output by this fit parameter, e.g. 2 = vd (1D) or vdx (2/3D)
              % this prevents data from jumping around between populations              
nmin  = 0.0001e6; nmax  = 200e6; % min/max values for initial guess
vdmin = -60000e3; vdmax = 60000e3;
%if strcmp(vdf.species,'ions'); vtmin = 4000e3/1836;  vtmax = 100000e3/1836; else, vtmin = 4000e3;   vtmax = 100000e3; end
vtmin = 4000e3;   vtmax = 100000e3;

% % If we use fminbnd, we need to give upper (hi) and lower (lo) bound
% % fminbnd only works for single variable functions. Options to make intial
% % fit to density for starting guess?
% nlo  = 0;    nhi  = 1e6; % for magnetotail
% vdlo = -1e9; vdhi = 1e9;
% vtlo = 0;    vthi = 1e9;


% Read input provided by the user
have_options = 0;
nargs = numel(varargin);
if nargs > 0, have_options = 1; args = varargin(:); end

while have_options
  l = 1;
  switch(lower(args{1}))
    case 'weight'
      l = 2;
      weight = args{2};
    case 'npop'
      l = 2;
      nPop = args{2};
    case 'plot'
      l = 2;
      doPlot = args{2};
    case 'iterplot'
      l = 2;
      iterPlot = args{2};
    case 'weight'
    case 'maxfunevals'
      l = 2;
      maxFunEvals = args{2};
    case 'maxiter'
      l = 2;
      maxIter = args{2};
    case 'x0'
      l = 2;
      X0 = args{2};
    case 'isortpop'
      l = 2;
      doSort = 1;
      iSortPop = args{2};
    case 'guessprevious'
      l = 2;
      doGuessPrevious = args{2};
    otherwise
      irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
      args = args(l+1:end);    
  end
  args = args(l+1:end);
  if isempty(args), break, end
end

% Some things needs to be assigned after, when we have other parameters,
% such as nPop
nInp = nPop*(1+2*nDim);
if isempty(maxIter); maxIter = 2000*nInp; end % default 200*length(x0)
if isempty(maxFunEvals); maxFunEvals = 2000*nInp; end % default 200*length(x0)
if isempty(weight); weight = ones(1,nPop*(1+nDim*2)); end

% Check validity of certain parameters
%if numel(weight) < nPop*(1+nDim*2)
%  warning(sprintf('Provided weight array has %g elements, while the minimum required is %g. Using weight = ones(...) instead.',numel(weight),nPop*(1+nDim*2)))
%end
        
% Set initial guess, unless provided
if isempty(X0)
  % How should we treat the initial guess? For now, maybe just choose some
  % reasonable values and see how it works.
  n0 = 0.01e6; % m^-3
  vd0 = 1000e3; % m/s
  T0 = 200; % eV
  vt0 = sqrt(T0*units.eV*2/m); % m/s
  switch nDim
    case 1 % 1D
      X0 = [n0 vd0 vt0]; 
    case 2 % 2D
      X0 = [n0 vd0 vd0 vt0 vt0];
    case 3 % 3D
      X0 = [n0 vd0 vd0 vd0 vt0 vt0 vt0];
  end
  % Repeat X0 for each population: [n1, vd1, vt1, ..., nN, vdN vtN]
  X0  = repmat(X0, 1,nPop);
end
% Apply funits to X so that we can return data in more conventional units, 
% cc, km, eV.
switch nDim
  case 1
    X_units = repmat({'m-6','m/s','m/s'},1,nPop);
    funits = @(X) reshape([X(1:3:end)*1e-6; X(2:3:end)*1e-3; X(3:3:end).^2*m/2/units.eV],1,numel(X));
    Xu_units = repmat({'cm-3','km/s','eV'},1,nPop);
  case 2
    X_units = repmat({'m-6','m/s','m/s','m/s','m/s'},1,nPop);
    funits = @(X) reshape([X(1:5:end)*1e-6; X(2:5:end)*1e-3; X(3:5:end)*1e-3; X(4:5:end).^2*m/2/units.eV; X(5:5:end).^2*m/2/units.eV],1,numel(X));
    Xu_units = repmat({'cm-3','km/s','km/s','eV','eV'},1,nPop);
  case 3 % not implemented
end

% Initialize arrays to save the data
data = struct([]);
%X_all = cell(nTimes,1);
%Q_all = cell(nTimes,1); % quality parameters

% fminsearch uses the Nelder-Mead simplex (direct search) method.
% I'm not sure if the large difference between parameters affects the 
% method. For example, vt/n = 2e2. If the method varies each parameter 
% which the same "step", vt will not be varied much...
options = optimset('OutputFcn', @myoutput,'MaxIter',maxIter,'MaxFunEvals',maxFunEvals);

% Initial guess, leter (?) one, the initial guess is the previous/neighbouring function
X = X0;


% Step through all the distribution functions
disp(sprintf('Fitting %g Maxwellians to distribution. Costfunction weight is [%g,%g,%g].',nPop,weight(1),weight(2),weight(3)))
fprintf('iTime = %4.0f/%4.0f\n',0,nTimes) % display progress       
for iTime = 1:nTimes
  if mod(iTime,1) == 0, fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],iTime,nTimes); end % display progress
  if not(doGuessPrevious)    
    X = X0; % Comment/remove this to use previous results as initial guess
            % for next step.  
  end
  str_time = irf_time(vdf(iTime).time,'EpochTT>utc_yyyy-mm-dd HH:MM:SS:mmm');
  history = [];
  nRedo = 0;
  EXITFLAG = 0;
  nIter = 0;
  % Option to plot output. 2-D plots for each time step is consuming quite a
  % lot of resources and is only good or initial checks. 
  switch nDim
    case 1
      ftmp = f(iTime,:);
      vtmp = v(iTime,:); % use cellfun?
    case 2
      ftmp = squeeze(reshape(f(iTime,:),datasize(2),datasize(3)));
      vtmp = {v{1}(iTime,:),v{2}(iTime,:)};
    case 3 % not implemented      
  end
  % If the density is zero, or perhaps if the thermal spread is zero, the 
  % function has a hard time finding it's way back, so try to put it at
  % least so some small but finite value.
 cost_function = @(X) costfunction_maxwellian(X,vtmp,ftmp);
  %X
  while EXITFLAG == 0 && nRedo <= maxRedo
    [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(cost_function,X,options);  
    nIter = nIter + OUTPUT.iterations;
    nRedo = nRedo + 1;
  end
  f_final = fit_function(X,vtmp);
  
  
  % Sort output by some parameter, to keep population to jump between
  % values too much.
  if doSort
    iSort = iSortPop:(2*nDim+1):numel(X);
    [~,iSortedPop] = sort(X(iSort));
    Xreshape = reshape(X,nPop,2*nDim+1);
    Xreshape = Xreshape(iSortedPop,:);
    X = reshape(Xreshape,1,numel(X));   
  end
  %catch
  %  1;
  %end
  % only for single variable
  % [X,FVAL,EXITFLAG,OUTPUT] = fminbnd(cost_function,Xlo,Xhi,options); %
  
  if nDim > 1 % This is to allow to use cat later, but for individual it makes it worse
    f_final = reshape(f_final,[1,size(f_final)]);
    f_final_separate = reshape(f_final_separate,[1,size(f_final_separate)]);
  end
  
  data(iTime).nRedo = nRedo - 1;
  data(iTime).iterations = nIter;
  data(iTime).exitflag = EXITFLAG;
  data(iTime).funcCount = OUTPUT.funcCount; % not sure what this is...
  data(iTime).cost = FVAL;
  data(iTime).nDim = nDim;
  data(iTime).nPop = nPop;
  data(iTime).X_units = X_units;
  data(iTime).X = X;  
  data(iTime).Xu_units = Xu_units;
  data(iTime).Xu = funits(X);
  data(iTime).history = history;
  data(iTime).f = f_final;
end

% Collect output
varargout{1} = data;
if nargout == 2
  data_all = cat(1,data.Xu);
  f_all = cat(1,data.f);   
  f_all_separate = cat(1,data.f_separate);
  switch nDim
    case 1
      % Make TSeries
      ts_cf = irf.ts_scalar(vdf.time,cat(1,data.cost)); ts_cf.name = 'cost function';
      ts_exitflag = irf.ts_scalar(vdf.time,cat(1,data.exitflag)); ts_exitflag.name = 'exitflag';
      ts_iter = irf.ts_scalar(vdf.time,cat(1,data.iterations)); ts_iter.name = 'iterations';
      ts_feval = irf.ts_scalar(vdf.time,cat(1,data.funcCount)); ts_feval.name = 'function evaluations';
      ts_n = irf.ts_scalar(vdf.time,data_all(:,1:(2*nDim+1):end)); ts_n.units = Xu_units{1}; ts_n.name = 'n';
      ts_vd = irf.ts_scalar(vdf.time,data_all(:,2:(2*nDim+1):end)); ts_vd.units = Xu_units{2}; ts_vd.name = 'vd';
      ts_T = irf.ts_scalar(vdf.time,data_all(:,3:(2*nDim+1):end)); ts_T.units = Xu_units{3}; ts_T.name = 'T'; 
      ts_f =  PDist(vdf.time,f_all,'1Dcart',v{1}*1e-3);
      
      % Apply limits to data, data outside of limits are put to NaN
      ts_T.data(ts_T.data>30e4) = NaN; % 30000 eV
      ts_n.data(ts_n.data>200) = NaN; % 200 cc
      ts_vd.data(ts_vd.data>1e5) = NaN; % 100000 km/s
      
      % Collect TSeries into structure
      ts.n = ts_n;
      ts.vd = ts_vd;
      ts.T = ts_T;
      ts.f = ts_f;
      ts.cf = ts_cf;
      ts.exitflag = ts_exitflag;
      ts.iter = ts_iter;
      ts.feval = ts_feval;
      
      % Add to output
      varargout{2} = ts;      
    case 2 % Working on implementing            
      % Make TSeries
      ts_cf = irf.ts_scalar(vdf.time,cat(1,data.cost)); ts_cf.name = 'cost function';
      ts_exitflag = irf.ts_scalar(vdf.time,cat(1,data.exitflag)); ts_exitflag.name = 'exitflag';
      ts_iter = irf.ts_scalar(vdf.time,cat(1,data.iterations)); ts_iter.name = 'iterations';
      ts_feval = irf.ts_scalar(vdf.time,cat(1,data.funcCount)); ts_feval.name = 'function evaluations';
      ts_n = irf.ts_scalar(vdf.time,data_all(:,1:(2*nDim+1):end)); ts_n.units = Xu_units{1}; ts_n.name = 'n';
      ts_vd1 = irf.ts_scalar(vdf.time,data_all(:,2:(2*nDim+1):end)); ts_vd1.units = Xu_units{2}; ts_vd1.name = 'vd, dir1';
      ts_vd2 = irf.ts_scalar(vdf.time,data_all(:,3:(2*nDim+1):end)); ts_vd2.units = Xu_units{2}; ts_vd2.name = 'vd, dir2';
      ts_T1 = irf.ts_scalar(vdf.time,data_all(:,4:(2*nDim+1):end)); ts_T1.units = Xu_units{3}; ts_T1.name = 'T, dir1'; 
      ts_T2 = irf.ts_scalar(vdf.time,data_all(:,5:(2*nDim+1):end)); ts_T2.units = Xu_units{3}; ts_T1.name = 'T, dir2'; 
      
      %if vdf.length == 1, f_all = reshape(f_all,[1 size(f_all)]);end      
      ts_f =  PDist(vdf.time,f_all,'2Dcart',v{1}*1e-3,v{2}*1e-3);

      % Also pass the different population separately
      for iPop = 1:nPop
        ts_f_sep{iPop} =  PDist(vdf.time,f_all_separate(:,:,:,iPop),'2Dcart',v{1}*1e-3,v{2}*1e-3);
      end
      
      % Apply limits to data, data outside of limits are put to NaN
      ts_T1.data(ts_T1.data>30e4) = NaN; % 30000 eV
      ts_T2.data(ts_T2.data>30e4) = NaN; % 30000 eV
      ts_n.data(ts_n.data>200) = NaN; % 200 cc
      ts_vd1.data(ts_vd1.data>1e5) = NaN; % 100000 km/s
      ts_vd2.data(ts_vd2.data>1e5) = NaN; % 100000 km/s

      % Collect TSeries into structure
      ts.n = ts_n;
      ts.vd1 = ts_vd1;
      ts.vd2 = ts_vd2;
      ts.T1 = ts_T1;
      ts.T2 = ts_T2;
      ts.f = ts_f;
      ts.cf = ts_cf;
      ts.exitflag = ts_exitflag;
      ts.iter = ts_iter;
      ts.feval = ts_feval;
      ts.f_sep = ts_f_sep;

      % Add to output
      varargout{2} = ts;     
    case 3 % not implemented
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Auxiliary/help functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stop = myoutput(x,optimvalues,state)
  % Collect data during evaluation of fminsearch
  stop = false;
  if isequal(state,'iter')
    history = [history; x];    
  end
end 
function out = costfunction_maxwellian(X,v,f_obs)
% Cost function for N Maxwellian distributions
% Params should be ordered as follows for each population
%   1-D: n1, vd1_dim1, vt1_dim1, n2, vd2_dim1, vt2_dim1, ...
%   2-D: n1, vd1_dim1, vd1_dim2, vt1_dim1, vt1_dim2, ...
%        n2, vd2_dim1, vd2_dim2, vt2_dim1, vt2_dim2, ...


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
%nInp = numel(X);
%nDim = numel(v);

% Evaluate the function for current parameters
% Better do do this in subfunction, then one can also call it from the main
% function and return the data with the output.
f_fit = fit_function(X,v); 

% Here we do the cost function, i.e. the difference between the function we
% want to fit to (for the current parameters) and the data.

% proportional to density
data_fit_mom0 = tocolumn(v).^0.*tocolumn(f_fit(:));
data_obs_mom0 = tocolumn(v).^0.*tocolumn(f_obs(:));
% proportional to speed
data_fit_mom1 = tocolumn(v).^1.*tocolumn(f_fit(:));
data_obs_mom1 = tocolumn(v).^1.*tocolumn(f_obs(:));
% proportional to temperature
data_fit_mom2 = tocolumn(v).^2.*tocolumn(f_fit(:));
data_obs_mom2 = tocolumn(v).^2.*tocolumn(f_obs(:));


% cost functions for different components and moments
out_mom0 = sum((data_fit_mom0-data_obs_mom0).^2)/sum(data_obs_mom0.^2);    
out_mom1 = sum((data_fit_mom1-data_obs_mom1).^2)/sum(data_obs_mom1.^2);    
out_mom2 = sum((data_fit_mom2-data_obs_mom2).^2)/sum(data_obs_mom2.^2);    
        
if any(f_fit(:) < 0)
  out_mom0 = Inf;
end

% Total costfunction is a combination (specified by weight) of the 
% separate costfunctions above.
out = weight(1)*out_mom0 + weight(2)*out_mom1 + weight(3)*out_mom2;

if not(isreal(out))
  1;
end
disp(out)

if doPlot && mod(iterPlot,size(history,1)) == 0
  nrows = 3; ncols = 1; isub = 1;
  
  % Make string to display parameters
  str_param = sprintf('N = %g, k = %g, theta = %g',X(1),X(2),X(3));  
  
  % Plot
  if max(v) > 1e6
    vscale = 1e-6;
    vstr = 'v (10^3 km/s)';
  else
    vscale = 1e-3;
    vstr = 'v (km/s)';
  end
  hca = subplot(nrows,ncols,isub); isub = isub + 1;    
  plot(hca,v*vscale,data_obs_mom0,v*vscale,data_fit_mom0)
  hca.Title.String = {...
    'Cost function (to be minimized):',...
    sprintf(' CF = %g x CF0 + %g x CF1 + %g x CF2 = %g',weight(1),weight(2),weight(3),out),...
    sprintf('(f) CF0 = %g',out_mom0)
    };
  irf_legend(hca,str_param,[0.02 0.98],'color',[0 0 0])
  irf_legend(hca,str_time,[0.98 0.98],'color',[0 0 0])
%       hca.XLabel.String = vstr;
%       hca.YLabel.String = 'f';
%       hca.XLim = [min(v{:}) max(v{:})]*vscale;
  %legend(hca,{'obs','fit'},'location','northeast','box','off')

  hca = subplot(nrows,ncols,isub); isub = isub + 1;    
  plot(hca,v*vscale,data_obs_mom1,v*vscale,data_fit_mom1)
  hca.Title.String = {sprintf('(vf) CF1 = %g',out_mom1)};
%       hca.XLabel.String = vstr;
%       hca.YLabel.String = 'vf';
%       hca.XLim = [min(v{:}) max(v{:})]*vscale;

  hca = subplot(nrows,ncols,isub); isub = isub + 1;    
  plot(hca,v*vscale,data_obs_mom2,v*vscale,data_fit_mom2)
  hca.Title.String = {sprintf('(v^2f) CF2 = %g',out_mom2)};  
%       hca.XLabel.String = vstr;
%       hca.YLabel.String = 'v^2f';
%       hca.XLim = [min(v{:}) max(v{:})]*vscale;

  drawnow
  pause(0.1)
  %pause
end

end

function f_fit = fit_function(X,v)
  %%
  A = @(N,k,theta) N*pi*k*theta^(-3/2)*gamma(k+1)/gamma(k-1/2);
  fit_fun = @(v,N,k,theta) A(N,k,theta).*(1 + v.^2./k./theta.^2).^(-k-1);
  
  % flat top
  
  A = @(N,k,theta) (3/2/pi)*N*theta^(-3)*gamma(1/k)/gamma(1+3/2/k)/abs(gamma(-1/2/k));
  fit_fun = @(v,N,k,theta) A(N,k,theta).*(1 + (v./theta).^(2*k)).^((-k-1)/k);
  
  %fit_fun = zeros(size(f_obs));
  
  % Evaluate the fit function for current parameters
  f_fit = fit_fun(v,X(1),X(2),X(3));
end
end