%syms T fH f f0
f0 = 1;
fH = 1;
T = @(f,fH,f0) f0*fH./sqrt(f)./(fH-f).^1.5;

f = logspace(-3,log10(3),100);

colors = pic_colors('matlab');
isub = 1;
nrows = 1;
if 0
  isub = 1;
  hca = subplot(nrows,1,isub);
  plot(hca,T(f,1,1),f,T(f,1.5,1),f,T(f,2,1),f)
  hca.Title.String = 'Varying f_H, constant f_0';
  hca.XLabel.String = 'Delay time, dT (arb. units)';
  hca.YLabel.String = 'f (arb. units)';
  legs1 = arrayfun(@(x)(sprintf('f_H = %3.1f, f_0 = 1.0',x)), [1 1.5 2], 'UniformOutput', false);
  legend(hca,legs1,'location','southeast','box','off')
end

if 0
  isub = 1;
  hca = subplot(nrows,1,isub);
  hl = plot(hca,T(f,1,0.5),f,T(f,1,1),f,T(f,1,2),f);
  c_eval('hl(?).Color = colors(?+3,:);',1:3)
  hca.Title.String = 'Constant f_H, varying f_0';
  hca.XLabel.String = 'Delay time, dT (arb. units)';
  hca.YLabel.String = 'f (arb. units)';
  legs2 = arrayfun(@(x)(sprintf('f_H = 1.0, f_0 = %3.1f',x)), [0.5 1 2], 'UniformOutput', false);
  legend(hca,legs2,'location','southeast','box','off')
end
if 0
isub = 1;
hca = subplot(nrows,1,isub);
%plot(hca,T(f,1,1),f,T(f,1.5,1),f,T(f,2,1),f,T(f,1,0.5),f,T(f,1,1),f,T(f,1,2),f)
fh0 = [1 1.5 2];
f0 = [1 1 1];
plot(hca,T(f,fh0(1),f0(1)),f,     T(f,fh0(2),f0(2)),f,     T(f,fh0(3),f0(3)),f)
hold(hca,'on')
fh0 = [1 1.5 2];
f0 = [1.5 1.5 1.5];
plot(hca,T(f,fh0(1),f0(1)),f,'--',T(f,fh0(2),f0(2)),f,'--',     T(f,fh0(3),f0(3)),f,'--')
hold(hca,'off')
hca.XLabel.String = 'Delay time, T (arb. units)';
hca.YLabel.String = 'f (arb. units)';

legs3a = arrayfun(@(x)(sprintf('f_H = %3.1f, f_0 = 1.0',x)), [1 1.5 2], 'UniformOutput', false);
legs3b = arrayfun(@(x)(sprintf('f_H = %3.1f, f_0 = 1.5',x)), [1 1.5 2], 'UniformOutput', false);
legend(hca,{legs3a{:},legs3b{:}},'location','southeast','box','off')
end

hca = subplot(nrows,1,isub);
%plot(hca,T(f,1,1),f,T(f,1.5,1),f,T(f,2,1),f,T(f,1,0.5),f,T(f,1,1),f,T(f,1,2),f)
fH = logspace(log10(5),log10(0.2),8);
f0 = [2 2 2 2 0.2 0.2 0.2 0.2];
%f0 = [2 2 2 2 2 2 2 2];

Tall = zeros(numel(f0),numel(f));
for ii = 1:8%numel(f0)
  Ttmp = T(f,fH(ii),f0(ii));
  Ttmp(f>fH(ii)) = NaN;
  Tall(ii,:) = Ttmp;
  
  leginp(ii).f0 = f0(ii);
  leginp(ii).fH = fH(ii);
end

plot(hca,Tall,f)
hca.XLabel.String = 'Delay time, T (arb. units)';
hca.YLabel.String = 'f (arb. units)';
legs = arrayfun(@(x)(sprintf('f_H = %3.1f, f_0 = %3.1f',x.fH,x.f0)), leginp, 'UniformOutput', false);
legend(hca,legs,'location','eastoutside','box','off')

hl = findobj(hca,'type','line');
c_eval('hl(?).LineStyle = ''--'';',1:4)


h = findobj(gcf,'type','axes'); h = h(end:-1:1);
linkprop(h,{'XLim','YLim'});
h(1).XLim = [0 22];
h(1).YLim = [0 0.7];

c_eval('h(?).LineWidth = 1;',1:nrows)
c_eval('h(?).FontSize = 13;',1:nrows)
c_eval('h(?).Position(2) = 0.15;',1:nrows)
c_eval('h(?).Position(3) = 0.6;',1:nrows)

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
