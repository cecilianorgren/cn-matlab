function out = whamp_string(n,t,vd,a1)
% DISTR.WHAMP_STRING Makes string of input parameters.
%   output = DISTR.WHAMP_STRING(n,t,vd,a1)
%       output - cell array with as many components as n entries.
%            n - density in m^-3
%            t - temperature in keV
%           vd - v_d/vt
%           a1 - T_perp/T_par
%   Example use:
%       text(2e1,1e-4,distr.whamp_string(n,t,vd,a1))

n=n*1e-6; % m^3 -> cc

n_comp=numel(n);

for ii = 1:n_comp
    out{ii} = ['n=',num2str(n(ii)),'cc  ',...
               'T=',num2str(t(ii)),'keV  ',...
               'T_{perp}/T_{par}=',num2str(a1(ii)),'  ',...
               'v_{d}/v_{t}=',num2str(vd(ii))];          
end