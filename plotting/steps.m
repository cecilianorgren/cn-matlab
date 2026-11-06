function h = steps(ax,x,y,varargin)

if numel(x)-numel(y) == 1 % x is edges  
  x = [x(1:end-1)' x(2:end)'];
else
  x = [x(:) x(:)];
end

y = [y(:) y(:)];

x = reshape(x',[numel(x),1]);
y = reshape(y',[numel(y),1]);

plot(ax,x,y,varargin{:})

