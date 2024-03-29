function varargout = fit_esw_gaussian(x_orig,E,varargin)
% fit_params = fit_esw_gaussian(x,E,guess);
% [fit_params,fit_fun] = fit_esw_gaussian(x,E,guess);
%
% Based on the electric field, makes a fit to a gaussian potential:
%   phi = phi0*exp(-x^2/2/lx^2)
%   Ex = (phi0*z*exp(-x^2/(2*lx^2)))/lx^2
%
% Input:
%   x - length or time
%   E - electric field data
%   guess - optional - intial guess of fitting parameters [phi0,x0,lx]
%
% Output:
%   Structure with following fields>
%     fit_params - best fit of [phi0, x0, lx]
%     fit_fun - @(x) Ex(x,phi0,x0,lx) - matlab function - Efit = fit_fun(x)
%
% Example:
%   esw_fit = fit_esw_gaussian(E3);
%   plot(esw_fit.x,esw_fit.E,esw_fit.x,esw_fit.fun_E(esw_fit.x))  
% 

  if isnumeric(x_orig) && isnumeric(E)
    x_ref = x_orig(1);
    x = x_orig - x_ref;    
    E = E;
  elseif isa(x_orig,'TSeries')
    %x = x_orig.time - x_orig.time(1);
    x_ref = x_orig.time(1);
    x = x_orig.time - x_ref;
    E = x_orig.data;
  end
  
  % Estimate some parameters
  % Value and index of max(E) and min(E)
  [Emin,i_Emin] = min(E);
  [Emax,i_Emax] = max(E);
  % Location of max(E) and min(E)
  x_Emin = x(i_Emin);
  x_Emax = x(i_Emax); 
  % Locations of 1st and 2nd E peak
  x_peak1 = x(min([i_Emin,i_Emax]));
  x_peak2 = x(max([i_Emin,i_Emax]));
  
  % Estimate potential parameters based on amplitudes and locations of peak
  % electric field
  lx_est   = 0.5*abs(x_Emin-x_Emax);   % ls = lpp/2
  x0_est   = mean([x_peak1, x_peak2]); % location of center of structure
  phi0_est = (Emax-Emin)*lx_est;       % dE/dx
  
  % Function we will fit, defined here for plotting purposes
  mf_phi = @(z,phi0,lz,z0) phi0.*exp(-z.^2./2./lz.^2);
  mf_E = @(z,phi0,lz,z0) (phi0.*exp(-(z - z0).^2./(2*lz^2)).*(2*z - 2*z0))./(2*lz.^2);
  
  % Option to supply initial guesses for parameters
  if not(isempty(varargin)) % initial guess of parameters phi, z0, lz
    params0 = varargin{1};
  else % make gaussion fit   
    guess_params = [phi0_est; lx_est; x0_est];    
    fun = @(params) costfunction(params,x,E); % costfuction, root mean squared
    params0 = double(guess_params);
    
    % Save value of costfunction.
    history = [];
    history_param = [];
    options = optimset('OutputFcn', @myoutput);
    [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(fun,params0,options);
    params = X;

    if 0
      figure(77)      
      plot(x,E,...
           x,mf_E(x,params(1),params(2),params(3)),...
           x,mf_E(x,params0(1),params0(2),params0(3)),...
           x(i_Emin),E(i_Emin),'*',...
           x(i_Emax),E(i_Emax),'*',...
           params0(3)*[1 1],[max(E) min(E)])
      legend({'Data','best fit','initial guess','min(E)','max(E)','x0'},'location','southeast')
      pause(0.1)
    end
  end
  
  out.x = x;
  out.x_ref = x_ref;
  out.E = E;
  out.phi0 = params(1);
  out.x0 = x_ref + params(3);
  out.lx = params(2);
  out.note = 'reference function starts at x_ref';
  out.fun_phi = eval(sprintf('@(x) mf_phi(x,%g,%g,%g);',params(1),params(2),params(3)));
  out.fun_phi_ = @(x) mf_phi(x,params(1),params(2),params(3));
  out.fun_E = eval(sprintf('@(x) mf_E(x,%g,%g,%g);',params(1),params(2),params(3)));
  out.fun_E_ = @(x) mf_E(x,params(1),params(2),params(3));
  out.history.phi0 = history(1,:);
  out.history.lx = history(2,:);
  out.history.x0 = history(3,:);
 
  varargout{1} = out;


  function stop = myoutput(x,optimvalues,state)
    stop = false;
    if isequal(state,'iter')
      history = [history, x];
      %history_param = [history_param x]
    end
  end
  function sse = costfunction(params,tdata,Edata)
    % tdata or zdata
    z = tdata;
    % fit parameters
    phi0 = params(1);
    lz   = params(2);
    z0   = params(3);
    %sse = sum((ydata - (B0 + B1*tanh(tdata/dt))).^2);
    sse = sum((Edata - (phi0.*exp(-(z - z0).^2./(2*lz^2)).*(2*z - 2*z0))./(2*lz.^2)).^2);
    if phi0<0
      %sse = inf;
    end
  end
  function sse = constraint(params,tdata,Edata)
  end
end