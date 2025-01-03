function y = irf_derivative(timeseries,order,accuracy)
%IRF_DERIVATIVE Differenciate time series
%
% y = irf_derivative(x,derivative,accuracy)
%

% Finite central difference coefficients
% https://en.wikipedia.org/wiki/Finite_difference_coefficient
switch order % derivative order
  case 1
    coeff = [1/2     0       0       0;
             2/3     -1/12   0       0; 
             3/4     -3/20   1/60    0;
             4/5     -1/5    4/105   -1/280];
  case 2
    coeff = [-2      1       0       0       0;
             -5/2    4/3     -1/12   0       0;
             -49/18  3/2     -3/20   1/90    0;
             -205/72 8/5     -1/5    8/315   -1/560];
end

all_data = timeseries.data;
n = timeseries.length;
dt = diff(timeseries.time - timeseries.time.start);
dt = median(dt);
all_derivative = all_data*0;


pfa = accuracy/2+1:n-accuracy/2;    
pra = [1:accuracy/2 n-accuracy/2+1:n];
lcns = pfa;

for icomp = 1:size(all_data,2)
  var = all_data(:,icomp);
  derivative = zeros(length(lcns),size(var,2)); % initialize size of derivative
  if order == 1 % for first order derivative
  %  du/dx =  A(1,1)*( u(n+1,:) - u(n-1) )/dx % Second order accurate
  %  du/dx =  A(2,1)*( u(n+1,:) - u(n-1) )/dx + A(2,2)*( u(n+2,:) - u(n-2) )/dx % Fourth order accurate
  %  and so on...
    for term_cnt = 1:accuracy/2
      derivative = derivative + coeff(accuracy/2,term_cnt)*(var(lcns+term_cnt,:) - var(lcns-term_cnt,:))/dt;
    end
  elseif order == 2
  %  d2u/dx2 =  (A(1,2)*(u(n-1,:) + A(1,1)*(u(n,:) + A(1,2)*(u(n+1,:))/dx^2 % Second order accurate
    for term_cnt = -accuracy/2:accuracy/2
      derivative = derivative + coeff(accuracy/2,abs(term_cnt)+1)*(var(lcns+term_cnt,:))/dt^2;
    end
  end
  
  derivative = [nan(accuracy/2,1); derivative; nan(accuracy/2,1)];
  all_derivative(:,icomp) = derivative; 
end
y = timeseries.clone(timeseries.time,all_derivative);

if 0
  [~,ind_min]=min(time_steps);
  time_steps(ind_min)=[]; % remove the smallest time step in case some problems
  time_step=min(time_steps);


if size(x,2)==1, data_columns=1; else, data_columns=2:size(x,2);end

data_gaps=find(dt>3*time_step);
dt(data_gaps)=0;
y=x;
if data_gaps
  y(1:data_gaps(1),data_columns)=gradient(x(1:data_gaps(1),data_columns)',time_step)';
  for j=1:(length(data_gaps)-1)
    y((data_gaps(j)+1):data_gaps(j+1),data_columns)=gradient(x((data_gaps(j)+1):data_gaps(j+1),data_columns)',time_step)';
  end
  y(data_gaps(end):end,data_columns)=gradient(x(data_gaps(end):end,data_columns)',time_step)';
else
  y(:,data_columns)=gradient(x(:,data_columns)',time_step)';
end
end
 