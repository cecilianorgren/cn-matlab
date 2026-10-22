ffun = @(x)sind(x); % spin tone

ph_ex = 0:5:360; % some phase binning

f = 1*ffun(ph_ex) + 0.2*ffun(ph_ex*2)  + 0.1*ffun(ph_ex*4); % some "spin tone data"

ph_unwrap = linspace(0,2000,3000); % some phase unwrapped
ph_wrap = mod(ph_unwrap,360); % phase wrapped
t = 1:numel(ph); % "time index"

% vq = interp1(x,v,xq)
ft = interp1(ph_ex,f,ph_wrap); % get spin tone as a function of phase
ft_ = ffun(ph_wrap); % evaluate "original function"


nrows = 5;
hca = subplot(nrows,1,1);
plot(ph_ex,f,'*',ph_ex,ffun(ph_ex))
hca.YLabel.String = 'f';

hca = subplot(nrows,1,2);
plot(t,ph_unwrap,'*')
hca.YLabel.String = 'phase unwrapped';

hca = subplot(nrows,1,3);
plot(t,ph_wrap,'*')
hca.YLabel.String = 'phase wrapped';


hca = subplot(nrows,1,4);
plot(t,ft,'*')
hca.YLabel.String = 'f(t)';


hca = subplot(nrows,1,5);
plot(t,ft_,'*')
hca.YLabel.String = 'f(t)';
