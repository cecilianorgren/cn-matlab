cPs = @(w) w./(w.^2+1);
cP = @(wii,wee) cPs(wii) + cPs(wee);

cHs = @(w) 1./(w.^2+1);
cH = @(wii,wee) cHs(wii) - cHs(wee);



wii = logspace(-4,4,200);
wee = logspace(-4,4,201);

[Wee,Wii] = ndgrid(wii,wee);

CH = cH(Wii,Wee);
CP = cP(Wii,Wee);
%CH = cHs(Wii);

ipanel = 0;
nrows = 2;
ncols = 1;
h = gobjects([nrows,ncols]);
for irow = 1:nrows, for icol = 1:ncols, ipanel = ipanel + 1; h(irow,icol) = subplot(nrows,ncols,ipanel); end, end
isub = 1;


hca = h(isub); isub = isub + 1;
pcolor(hca,Wee,Wii,CH)
shading(hca,'flat')
hca.XScale = 'log';
hca.YScale = 'log';
hca.XLabel.String = '\nu_{en}/\Omega_{e}';
hca.YLabel.String = '\nu_{in}/\Omega_{i}';
hcb = colorbar(hca);
hcb.Label.String = '\sigma_P/(ne/B)';
colormap(hca,pic_colors('blue_red'))
hca.CLim = [-1 1];
hold(hca,'on')
%contour(hca,Wee,Wii,CH,[0 0],'k')
plot(hca,wee([1 end]),wii([1 end]),'k-')
hold(hca,'off')
hca.Title.String = 'Hall conductivity';



hca = h(isub); isub = isub + 1;
pcolor(hca,Wee,Wii,CP)
shading(hca,'flat')
hca.XScale = 'log';
hca.YScale = 'log';
hca.XLabel.String = '\nu_{en}/\Omega_{e}';
hca.YLabel.String = '\nu_{in}/\Omega_{i}';
hcb = colorbar(hca);
hcb.Label.String = '\sigma_H/(ne/B)';
colormap(hca,pic_colors('blue_red'))
hca.CLim = [-1 1];
hold(hca,'on')
%[Cc,hc] = contour(hca,Wee,Wii,CP,[0 0],'k');
plot(hca,wee([1 end]),wii([1 end]),'k-')
hold(hca,'off')
hca.Title.String = 'Pedersen conductivity';

%text_handles=clabel(Cc,hc,0,'LabelSpacing',100);
%[m n]=size(V);
%for i=1:n,
%    a=get(text_handles(i),'UserData');
%    y=char(num2str(a,'%5.2e'));
%    set(text_handles(i),'String', y, 'FontSize',14);% ,'Interpreter','LaTex');
%end


c_eval('h(?).FontSize = 16;',1:numel(h))
c_eval('h(?).Box = ''on'';',1:numel(h))
c_eval('h(?).LineWidth = 1.5;',1:numel(h))

