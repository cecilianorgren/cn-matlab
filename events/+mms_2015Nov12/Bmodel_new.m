%% Observed data
CS_normal_velocity = 70; % km/s
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z')+-0.01;
ic = 1;
c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE = obsE.resample(obsB);'...
],ic)
zObsB = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsE = (obsE.time.epochUnix-mean(obsE.time.epochUnix))*CS_normal_velocity;

%% Model
% BM, two gussians
% BL, two Harris
z = linspace(-30,30,1000)';
x = z*0;
y = z*0;

Bn = 1.5;
B0 = 5;
Bm = 5;
Bg = 5;

Bx = @(x,y,z) x*0 + y*0 - 4.5*tanh((z+11)/5)- 6*tanh((z-12)/5) -1.5; % L
By = @(x,y,z) x*0 + y*0 + Bg + 6*exp(-(z+10).^2/(5^2)) - 4.5*exp(-(z-15).^2/(8^2)); % M
Bz = @(x,y,z) x*0 + y*0 + z*0 + 1.5; % N

% EL = @(x,y,z) x*0 + y*0 + z*0 + 0*EL*(-exp(-(z.^2)/d^2));
% EM = @(x,y,z) x*0 + y*0 + z*0 + 1*(Er*(exp(-((z-8e3).^2)/g^2)) + Ey_inflow);
% EN = @(x,y,z) x*0 + y*0 + z*0 + 0*(E0*sin((z-dzE)/dE).*exp(-(z-dzE).^2/(2*b^2)) + 1*Enas*exp(-(z-offsEnas).^2./d.^2) + -0.8*Enas*exp(-(z+offsEnas).^2./2/d.^2));

BL = Bx(z,y,z);
BM = By(z,y,z);
BN = Bz(z,y,z);
Babs = sqrt(BL.^2 + BM.^2 + BN.^2);

hca = subplot(1,1,1); 
set(hca,'ColorOrder',mms_colors('xyza'))
plot(hca,z,[BL BM BN Babs])
hold(hca,'on')
plot(hca,zObsB,[obsB.data obsB.abs.data],'--')
hold(hca,'off')
