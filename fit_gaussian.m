function varargout = fit_gaussian(tsE,tsL,varargin)
  % Returns the 

  % Electric field from Gaussian
  % phi = phi0*exp(-z^2/2/lz^2)
  % Ez = (phi0*z*exp(-z^2/(2*lz^2)))/lz^2
 

  if isnumeric(tsL) && numel(tsL == 1) % tsL is a speed
    L = (tsE.time - tsE.time(1))*tsL;
  elseif isa(tsL,'TSeries')
    tsL = tsL.resample(tsE);
    L = tsL.data;
  end
  L = L - mean(L);
  E = tsE.data;
  
  % Estimate some parameters
  % Location of max(E) and min(E)
  [Emin,i_Emin] = min(E);
  [Emax,i_Emax] = max(E);
  L_Emin = L(i_Emin);
  L_Emax = L(i_Emax);
  L_peak1 = L(min([i_Emin,i_Emax]));
  L_peak2 = L(max([i_Emin,i_Emax]));
  lz_est = 0.5*abs(L_Emin-L_Emax); % ls = lpp/2
  z0_est = mean([L_peak1, L_peak2]);
  phi0_est = (Emax-Emin)*lz_est;
  
  mf_E = @(z,phi0,lz,z0)(phi0.*exp(-(z - z0).^2./(2*lz^2)).*(2*z - 2*z0))./(2*lz.^2);
  
  if isempty(varargin) % initial guess of half width
    imax = find(E == max(E));
    imin = find(E == min(E));
    guessL = 0.5*abs(L(imax)-L(imin));
  else % make gaussion fit
    imax = find(E == max(E));
    imin = find(E == min(E));
    guess_lz = 0.5*(L(min([imax imin]))-L(max([imax imin])));
    guess_z0 = mean(L(max([imax imin]))-L(min([imax imin])));
    guess_phi0 = 2*(E(min([imax imin]))-E(max([imax imin])))/guess_lz;
    
    guess_params = varargin{1};
    guess_params(3) = guess_z0;
    %guess_params = [guess_phi0; guess_lz; guess_z0];
    guess_params = [phi0_est; lz_est; z0_est];
    guess_params
    fun = @(params) costfunction(params,L,E); % costfuction, root mean squared
    params0 = double(guess_params);
    %options = optimset
    % params = [phi0,lz,z0]'
    if 1
      bestx = fminsearch(fun,params0);
    else        
      ub = [100 100 100];
      lb = [0 0 -100];
      bestx = fmincon(fun,params0,[],[],[],[],lb,ub);
    end
    guessL = bestx;
    if 1
      figure(77)
      
      plot(L,E,...
           L,mf_E(L,bestx(1),bestx(2),bestx(3)),...
           L,mf_E(L,guess_params(1),guess_params(2),guess_params(3)),...
           L(imin),E(imin),'*',...
           L(imax),E(imax),'*',...
           guess_params(3)*[1 1],[max(E) min(E)])
      pause(0.1)
    end
  end

  guessL;
  varargout = {guessL,@(z) mf_E(z,bestx(1),bestx(2),bestx(3))};


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