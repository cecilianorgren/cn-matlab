nrows = 2;
m = 0:5;  
n = 0:5;

m_edges = [m(1)-0.5 m+0.5];
n_edges = [n(1)-0.5 n+0.5];

a = 2.28*1e-2; % m
b = 1.01*1e-2; % m

c = 299792458; % m/s
fd = 1.70*1e10; % Hz
wd = fd*2*pi; % w = 2*pi*f

wmn = @(a,b,m,n) c*pi*sqrt((m/a).^2 + (n/b).^2);

[M,N] = meshgrid(m,n);
[M_edges,N_edges] = meshgrid(m_edges,n_edges);

WMN = wmn(a,b,M,N);
WMN(1,1) = NaN;
%plot()
hca = subplot(nrows,1,1);
htext = findobj(hca,'type','text');
delete(htext)
%pcolor(hca,M,N,WMN)
surf(hca,M_edges,N_edges,zeros(size(N_edges)),WMN)
view(hca,[0 0 1])
hcb = colorbar('peer',hca);
hcb.XLabel.String = '\omega_{mn} (Hz)';
hca.XLabel.String = 'm'; 
hca.YLabel.String = 'n'; 
hca.FontSize = 16;
for imode = 1:numel(WMN)  
  text(hca,N(imode),M(imode),1,sprintf('  %.2e',WMN(imode)),'horizontalalignment','center')
end

hca = subplot(nrows,1,2);
htext = findobj(hca,'type','text');
delete(htext)
%pcolor(hca,M,N,WMN-wd)
surf(hca,M_edges,N_edges,zeros(size(N_edges)),(WMN-wd))
view(hca,[0 0 1])
hcb = colorbar('peer',hca);
hcb.XLabel.String = '\omega_{mn}-\omega_d (Hz)';
%hca.Title.String = sprintf('f_d = %g, wd = 2pi*f_d = %g',fd,wd);
hca.XLabel.String = 'm'; 
hca.YLabel.String = 'n'; 
hca.FontSize = 16;
colormap(pic_colors('blue_red'))
hca.CLim = max(hca.CLim)*[-1 1];
for imode = 1:numel(WMN)  
  text(hca,M(imode),N(imode),1,sprintf('  %.2e',(WMN(imode)-wd)),'horizontalalignment','center')
end
%hca.CLim(1) = 0;
%hcb.YLim(1) = 0;
htext = findobj(gcf,'type','text');
c_eval('htext(?).FontSize = 16;',1:numel(htext))
%%
hca = subplot(1,1,1);
%pcolor(hca,M,N,WMN-wd)
plot(nan,nan)
delete(htext)
[WMN_sorted,isort] = sort(WMN(:));
M_sorted = M(isort);
N_sorted = N(isort);
hold(hca,'on')
for imode = 1:numel(WMN)
  plot(hca,imode,WMN_sorted(imode)/(1),'o')
  text(hca,imode,WMN_sorted(imode)/(1),sprintf('  %g%g',M_sorted(imode),N_sorted(imode)))
end
hca.XLabel.String = 'different modes, sorted by cutoff frequency'; 
hca.YLabel.String = '\omega_{mn}'; 
hold(hca,'off')
hold(hca,'on')
plot(hca.XLim,wd*[1 1],'k--')
text(hca,mean(hca.XLim),wd,sprintf('\omega_d = %g',wd),'verticalalignment','bottom')
hold(hca,'off')
grid(hca,'on')
hca.FontSize = 16;

