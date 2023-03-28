AU = units.AU;
RE = units.RE;

rbow = 15*RE; % m
t =  linspace(0,60*60,1000); % s


v = 400 + 200 + 200*(tanh((t-mean(t))/(60*1))); % km/s

L1 = 0.01*AU;
L = (L1 - rbow)*1e-3; % km



T = L./v;

for ii = 1:numel(v)
  T(ii) = L/v(ii);
end

tnew = t + T;

hca = subplot(4,1,1);
plot(hca,t,v)
hca.YLabel.String  = 'v_{SW} (km/s)';
hca.XLabel.String  = 't';
hca.YLim = [000 1000];

hca = subplot(4,1,2);
plot(hca,t,T)
hca.YLabel.String  = 'T (s)';
hca.YLim(1) = 0;
hca.XLabel.String  = 't';

hca = subplot(4,1,3);
plot(hca,t,tnew)
hca.XLabel.String  = 't';
hca.YLabel.String  = 't+T';
hca.Title.String = 'Arrival time at bow shock: t'' = t + T';
%plot(tnew,v)

hca = subplot(4,1,4);
[~,isort] = sort(tnew);
%plot(hca,tnew(isort),v(isort))
plot(hca,tnew,v,'.')
hca.YLabel.String  = 'v_{SW} (km/s)';
hca.XLabel.String  = 't + T';
hca.YLim = [000 1000];