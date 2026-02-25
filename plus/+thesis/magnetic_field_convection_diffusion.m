
B = @(x,x0,d) exp(-(x-x0).^2/d^2)/d;

xlim = 4;


nRows = 1; nCols = 3;
isub = 1;

h(1) = subplot(nRows,nCols,isub); isub = isub +1; hca = h(1);
xlim = [-2 3];
x = linspace(xlim(1),xlim(2),100);
lines = plot(hca,x,B(x,0,1),'--',x,B(x,1,1));
lines(1).Color = [0 0 0];
lines(2).Color = [0 0 0];
hca.Title.String = 'Convection';

h(2) = subplot(nRows,nCols,isub); isub = isub +1; hca = h(2);
xlim = [-4 4];
x = linspace(xlim(1),xlim(2),100);
lines = plot(hca,x,B(x,0,1),'--',x,B(x,0,2));
lines(1).Color = [0 0 0];
lines(2).Color = [0 0 0];
hca.Title.String = 'Diffusion';

h(3) = subplot(nRows,nCols,isub); isub = isub +1; hca = h(3);
xlim = [-2.5 7];
x = linspace(xlim(1),xlim(2),100);
lines = plot(hca,x,B(x,0,1),'--',x,B(x,3,2));
lines(1).Color = [0 0 0];
lines(2).Color = [0 0 0];
hca.Title.String = 'Convection and diffusion';


for ii = 1:3
  %hca = subplot(nRows,nCols,ii);
  h(ii).Box = 'off';
  h(ii).Box = 'off';
  h(ii).FontSize = 14;
  axis(h(ii),'tight');
  axis(h(ii),'off');
  h(ii).Title.Position(2) = -0.2;
  h(ii).Position(2) = 0.2;
  h(ii).Position(4) = 0.7;
end

  
%%

lines = plot(x,B(x,0,1),'--',x,B(x,1,1));
lines(1).Color = [0 0 0];
lines(2).Color = [0 0 0];
%%

lines = plot(x,B(x,0,1),'--',x,B(x,0,2));
lines(1).Color = [0 0 0];
lines(2).Color = [0 0 0];

%%

lines = plot(x,B(x,0,1),'--',x,B(x,3,2));
lines(1).Color = [0 0 0];
lines(2).Color = [0 0 0];

