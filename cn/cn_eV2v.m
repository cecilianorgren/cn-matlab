function res=cn_eV2v(data,which)
% CN_EV2V Transforms between eV and v for electrons.
%   res=cn_eV2v(data,which);
%       data - energy in eV or velocity in km/s       
%       which - 'eV' if input data is an energy
%               'v'  if input data is an velocity

%Units=irf_units;
Units.eV = 1.6022e-19;
Units.me = 9.109399999999e-31;


% E=mv^2/2
%lower(which)
switch lower(which)    
    case 'ev'
        res=sqrt(data*Units.eV*2/Units.me)/1000;
    case 'v'
        res=Units.me*data.^2*10^6/2/Units.eV;
end