function varargout = identify_esws(data,time_window,window_step,Emin)


% pot = phi0*exp(-x^2/L^2)
% E = -dpot/dx = phi0*(2*x/L^2)*exp(-x^2/L^2) = (E0/L^2)*exp(-x^2/L^2)
%f_esw = @(x,E0,L) (E0./L.^2)*exp(-x^2/L^2);
%f_esw = @(t,E0,t0) (E0./t0.^2)*exp(-t^2/t0^2);
f_esw = @(t,E0,t0) E0*exp(-t^2/t0^2);


tstep = window_step; % step of moving window

% Length of data in seconds
tEpoch0 = B.time(1);
Tdata = B.time(end)-B.time(1);

% Plot B field
if 0 
  irf_plot(B)
  hleg = legend({'B_x','B_y','B_z'});
  ylabel('B (GSM)')
  drawnow
end
% Save info as table
data = [];

% Step through time series
t1 = 0;
t2 = twindow;
it = 1;
while t2 < Tdata  
  % Resample and limit data to current time window
  tEpoch1 = tEpoch0 + t1;
  tEpoch2 = tEpoch0 + t2;
  Btmp = B.z.tlim([tEpoch1 tEpoch2]);  
  timeline = Btmp.time - Btmp.time(1); % s
  Bobs = Btmp.data; 
  
  % Linear regression
  % Find B0 and B1: Bfit = B0 + B1*tanh(t/dt)
  theta0 = [1; 10];
  alpha = 0.01;
  iterations = 1000;
  m = numel(Bobs);
  X = [ones(m,1) tanh((timeline-mean(timeline))/dt)];
  y = Bobs;
  [theta, J_history] = gradientDescent(X, y, theta0, alpha, iterations);
  Breg = X*theta; 
  
  % Optimization: Search parameters again, including dt this time.  
  % Use regression results as initial values.
  tdata = timeline-mean(timeline);
  ydata = Bobs;
  yfit = @(t,B0,B1,dt) B0 + B1*tanh(t/dt);
  fun = @(x) sseval(x,tdata,ydata); % costfuction, root mean squared
  x0 = double([theta; 1]);
  %options = optimset
  bestx = fminsearch(fun,x0);
  Bopt = yfit(tdata,bestx(1),bestx(2),bestx(3));
    
  % Save data in table
  data(it).it = it;
  data(it).t1 = tEpoch1;
  data(it).t2 = tEpoch2;
  data(it).Bts = Btmp;
  data(it).reg_guess = theta0;
  data(it).reg_best = theta;
  data(it).reg_RMS = rms(Breg-Bobs);  
  data(it).reg_RMS_J = J_history(end)*2;  
  data(it).reg_xcorr = xcorr(Breg,Bobs,0,'coeff');
  data(it).reg_B = Breg;   
  data(it).opt_guess = x0;
  data(it).opt_best = bestx;
  data(it).opt_RMS = rms(Bopt-Bobs);  
  data(it).opt_xcorr = xcorr(Bopt,Bobs,0,'coeff');
  data(it).opt_B = Bopt;
  
  if 0 % Plot results  
    colors =[ ...
         0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
  
    hca = subplot(2,1,1);
    plot(hca,timeline,Bobs,...
             timeline,Breg,...
             timeline,yfit_data) % timeline,Bobs,
    hca.Title.String = [tEpoch1.utc ' - ' tEpoch1.utc];
    irf_legend(hca,{sprintf('B = %.2f + %.2f*tanh(t/%.2f)',theta(1),theta(2),dt);sprintf('C_{RMS} = %.2f',J_history(end))},[1.03 0.98],'color',[0 0 0])
    legend(hca,{'B_{obs}',...
      sprintf('B = %.1f + %.1f*tanh(t/%.2f)',theta(1),theta(2),dt),...
      sprintf('B = %.1f + %.1f*tanh(t/%.2f)',bestx(1),bestx(2),bestx(3))
      },'location','eastoutside')

    hca = subplot(2,1,2);
    plot(hca,timeline,cBfit,timeline,cBobs)  

    irf_legend(hca,{sprintf('C = %.2f',C),sprintf('RMS = %.2f',dRMS)},[0.02 0.98],'color',[0 0 0])
    %pause(0.1)
  end
  
  % Go to next time window
  it = it + 1;
  t1 = t1 + tstep; 
  t2 = t1 + twindow;
end

out = data;

end
% Fitting/regression functions
function J = computeCost(X,y,theta)
  m = length(y); % number of training examples
  h = X*theta;
  J = (1/2/m)*sum((h-y).^2);
end
function [theta, J_history] = gradientDescent(X, y, theta, alpha, num_iters)
  m = length(y); % number of training examples
  J_history = zeros(num_iters, 1);
  for iter = 1:num_iters      
      h = X*theta;    
      hy = (h-y); 
      Xhy = X'*hy;
      theta = theta - (alpha/m)*Xhy;
      % Save the cost J in every iteration    
      J_history(iter) = computeCost(X, y, theta);
  end
end
function sse = sseval(x,tdata,ydata)
  B0 = x(1);
  B1 = x(2);
  dt = x(3);
  sse = sum((ydata - (B0 + B1*tanh(tdata/dt))).^2);
end

