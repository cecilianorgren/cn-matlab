
% All quantities are given in km and km/s

vthe = 3e3;
n = 3.5;
vExB = 200*0+1;
fce = 400; % Hz
n0 = sqrt(2*vExB/fce);

titleStr = sprintf('n_0=%s km, n = %s km, v_{ExB} = %s km/s, f_{ce} = %s Hz, v_{the} = %s km/s',...
             num2str(n0,'%g'),num2str(n,'%g'),num2str(vExB,'%g'),num2str(fce,'%g'),num2str(vthe,'%g'));
fStr = 'f = e^{-P_M^2/v_{the}^2}e^{-W_N/v_{the}^2}(1+tanh(W_N/v_{the}^2))/2';
vm = linspace(-10e3,10e3,50);
vn = linspace(-10e3,10e3,50);

[VN,VM] = meshgrid(vn,vm);

wn = @(vn,vm) vn.^2 - 2*(vExB-vm).*vExB.*n.^2./n0^2+vExB^2*n.^4./n0.^4; % km^2/s^2
pm = @(vn,vm) vm-vExB*n.^2./n0^2; % km/s

%WN = wn(VN,VM);
%PM = pm(VN,VM);

f = @(vn,vm) exp(-pm(VN,VM).^2/vthe^2).*exp(-wn(VN,VM)/vthe^2).*(1+tanh(wn(VN,VM)/vthe^2))/2;
th = @(vn,vm) 1+tanh(wn(VN,VM)/vthe^2);

nRows = 2;
nCols = 2;
for ii = 1:(nRows*nCols); h(ii) = subplot(nRows,nCols,ii); end
isub  = 1;

vnorm = 1e3;

hca = h(isub); isub = isub + 1;
pcolor(hca,VN/vnorm,VM/vnorm,wn(VN,VM))
shading(hca,'flat')
hca.XLabel.String = 'v_N';
hca.YLabel.String = 'v_M';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'W_N (km^2/s^2)';

hca = h(isub); isub = isub + 1;
pcolor(hca,VN/vnorm,VM/vnorm,pm(VN,VM))
shading(hca,'flat')
hca.XLabel.String = 'v_N';
hca.YLabel.String = 'v_M';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_M (km/s)';

hca = h(isub); isub = isub + 1;
pcolor(hca,VN/vnorm,VM/vnorm,th(VN,VM))
shading(hca,'flat')
hca.XLabel.String = 'v_N';
hca.YLabel.String = 'v_M';
hcb = colorbar('peer',hca);
hca.Title.String = '1+tanh(W_N/v_{the}^2)';

hca = h(isub); isub = isub + 1;
pcolor(hca,VN/vnorm,VM/vnorm,f(VN,VM))
shading(hca,'flat')
hca.XLabel.String = 'v_N';
hca.YLabel.String = 'v_M';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'f(v_N,v_M)';
hca.Title.String = fStr;

h(1).Title.String = titleStr;