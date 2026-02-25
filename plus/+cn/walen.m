function out = walen(v0,b0,n0,tperp0,tpar0,tw)
% CN.WALEN Performes the Walen-test.
%   The theoretical velocity jump should be:
%       +- sqrt((1-alfa1)/mu0/rho)*(B2t(1-alfa2)-B1t(1-alfa1));
%        alfa = mu0*(ppar-pper)/B^2;

if size(b0,2)==4
    b0=irf_abs(b0);
end
if size(v0,2)==4
    v0=irf_abs(v0);
end

%b=irf_resamp(b0,v);
%n=irf_resamp(n0,v);
%tperp=irf_resamp(tperp0,v);
%tpar=irf_resamp(tpar0,v);

units=irf_units;

for ii=1:2;
    tperp_prel = irf_tlim(tperp0,tw{ii}); tperp{ii} = mean(tperp_prel(:,2));
    tpar_prel = irf_tlim(tpar0,tw{ii});   tpar{ii} = mean(tpar_prel(:,2));
    n_prel = irf_tlim(n0,tw{ii});         n{ii} = mean(n_prel(:,2));
    b_prel = irf_tlim(b0,tw{ii});         b{ii} = mean(b_prel(:,2:5));
    v_prel = irf_tlim(v0,tw{ii});         v{ii} = mean(v_prel(:,2:5));
    
    alpha{ii} = units.mu0*n{ii}*(tpar{ii}-tperp{ii})/b{ii}(4);
end

%alpha{} = (tpar(:,2)-tperp(:,2)).*(units.mu0*n(:,2)./b(:,5).^2)

%dv = sqrt(1-)
