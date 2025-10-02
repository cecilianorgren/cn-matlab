function varargout = fit_esw(E)
% FIT_ESW Fits a bipolar function to electric field data
%   Returns center time and half width
% phi = phi0*exp(-t^2/tp2p^2);
% E = -phi0*(2t/tp2p^2)*exp(-t^2/tp2p^2);

units = irf_units;

E_obs = E.data*1e-3;

t = E.time - E.time.start;
T = E.time.stop - E.time.start;
t0 = 0.5*T;

f_esw = @(t,A,tp2p,t0,offset) -A.*(2*(t-t0)./tp2p^2).*exp(-((t-t0)./tp2p).^2) + offset;
%f_esw = @(t,A,tp2p,t0) (t-t0).*A.*exp(-((t-t0)./tp2p).^2);

cost_function = @(X) costfunction(X,t,E_obs,f_esw);
%A0 = 10*units.eV; % V
A0 = -max(E_obs)/2*0.5*T;
X0 = [A0*2/T^2 0.25*T t0 0]; % phi0, tp2p, t0, offset
X0 = [A0 0.25*T t0 0]; % phi0, tp2p, t0, offset
X0(1);
%X0 = [0 0 t0]; % phi0, tp2p, t0
[X,FVAL,EXITFLAG,OUTPUT] = fminsearch(cost_function,X0);  


%plot(t,E_obs,t,f_esw(t,X(1),X(2),X(3)),t,f_esw(t,X0(1),X0(2),X0(3)))

tsESW = irf.ts_scalar(E.time,1e3*f_esw(t,X(1),X(2),X(3),X(4)));

varargout{1} = X;
varargout{2} = tsESW;
varargout{3} = EXITFLAG;


function out = costfunction(X,t,E_obs,f_esw)
% Evaluates sum(E_mod-E_obs)^2/sum(E_obs^2)

  E_mod = f_esw(t,X(1),X(2),X(3),X(4));
  out = sum(E_mod-E_obs).^2/sum(E_obs.^2);
  out = norm(E_mod-E_obs);

  %plot(t,E_obs,t,E_mod,t,E_mod-E_obs)
  %pause(0.1)