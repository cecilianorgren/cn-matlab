a = rand(10,1);
a(a < 0.5) = 0;
%%
na = numel(a);
ava = 0;
for ii = 1:(na)
  ava = (ava*(ii-1) + a(ii))/(ii);
  disp([ava])
  
end



% mean(a(1:4)) = (a(1) + a(2) + a(3) + a(4)/4;aca()
% mean(a(1:2)) = (a(1) + a(2))/2
