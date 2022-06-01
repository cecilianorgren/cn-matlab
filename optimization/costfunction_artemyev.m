function out = costfunction_artemyev(params,x,y,doPlot)
% COSTFUNCTION_ARTEMYEV(params,x,y,optionPlot)
%   Cost function 
%     CF = sum(sqrt((yfun-y)^2))
%   for 
%     yfun = y0*(1-a*x^2)   
%   where y is n or T, and x is Bx (Artemyev 2017).
% 
%   Params = [y0, a]
%   optionPlot = 1, to plot fmin search

f = @(x,y0,a) y0*(1-a*x.^2);

y0 = params(1);
a = params(2);

out = sum(sqrt((f(x,y0,a) - y).^2));

% if a < 0 || a > 1 || y0 > 10000
%   out = Inf;
%   return;
% end

if doPlot
  hca = subplot(1,1,1);
  plot(hca,x,y,'.',x,f(x,y0,a),'-')  
  hca.Title.String = sprintf('y0 = %g, a = %g',y0,a);
  drawnow
  disp(sprintf('y0 = %g, a = %g',y0,a))
  pause(0.1)
end

