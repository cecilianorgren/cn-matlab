function v = randv(N,vt,veh)
% Generate random velocities from distribution v/vt*exp(v^2/vt^2).
% v = ExB.randv(N,vt);
%       N - number of velocities
%       vt - thermal velocity

% x = v/vt;
%pdf = @(x) x.*exp(-(x.^2-x0.^2));
%cdf =  @(x) 1-1*exp(-x.^2);

cdfy = rand(N,1);
cdfx = @(cdfy) sqrt(-log(1-cdfy));

v = cdfx(cdfy).*vt;

%pdf = @(x) exp(-v.^2/vt^2)/sqrt(pi)/vt;
%cdf =  @(x) 1-1*exp(-x.^2);


%pdf = @(v,vt,veh) (v-veh)/vt.*exp(-(v.^2-veh.^2)/vt^2);
v = -50000:200000;
%plot(v,abs(pdf(v,vt,veh)))
cdfy = rand(N,1);
cdf = exp(-(v-veh).^2/vt.^2);

v = vt*1e3*sqrt(-log(2*cdfy/vt/1e3))+veh*1e3;

%syms x
%for ii = 1:N    
%    v(ii) = solve(0.5*(vt*1e3)^2*exp(-(x-veh*1e3)^2/(vt*1e3)^2)+0.5*vt*1e3*sqrt(pi)*1i*erf(1i*(veh*1e3-x)/(vt*1e3)) == cdfy(ii), x);    
%end


