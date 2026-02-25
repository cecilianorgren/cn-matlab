tsl = 729600;

FID = fopen('mms3_scb_20150828_125300t.dat');
Bscmt = fread(FID,729600,'int64');
Bscmt = int64(Bscmt);
Bscmt =  EpochTT(Bscmt);

FID = fopen('mms3_scb_20150828_125300Bxyz.dat');
B = fread(FID,729600*3,'float');
B = [B(1:tsl) B((1+tsl):2*tsl) B((1+2*tsl):3*tsl)];

Bscm = TSeries(Bscmt,B,'to',1);


%load brst mode interval
tmpDataObj = dataobj('data/mms3_edp_brst_ql_dce2d_20150828125314_v0.2.0.cdf');
Exyz = get_variable(tmpDataObj,'mms3_edp_dce_xyz_dsl');
Exyz = mms.variable2ts(Exyz);
Exyz.data(:,3) = Exyz.data(:,3)*1.5;
Exyz.data(find(abs(Exyz.data) > 100)) = NaN;

tmpDataObj = dataobj('data/mms3_dfg_srvy_ql_20150828_v0.0.3.cdf');
Bxyz = get_variable(tmpDataObj,'mms3_dfg_srvy_dmpa');
Bxyz = mms.variable2ts(Bxyz);

tlimit = irf.tint(Exyz.time.start.utc,Exyz.time.stop.utc);
tlimit = irf.tint('2015-08-28T12:54:25.0Z/2015-08-28T12:54:28.0Z');

Exyz = Exyz.tlim(tlimit);
Bscm = Bscm.tlim(tlimit);
tlimitl = tlimit+[-1 1];
Bxyz = Bxyz.tlim(tlimitl);
Bxyz = Bxyz.resample(Exyz);
Bscm = Bscm.resample(Exyz);

Bscm = irf_filt(Bscm,10,1000,8192,5);
Exyz = irf_filt(Exyz,10,1000,8192,5);

SCpos = [0 1 0];

Bmag = Bxyz.abs.data;
Rpar = Bxyz.data./[Bmag Bmag Bmag];
Rperpy = irf_cross(Rpar,SCpos);
Rmag   = irf_abs(Rperpy,1);
Rperpy = Rperpy./[Rmag Rmag Rmag];
Rperpx = irf_cross(Rperpy, Rpar);
Rmag   = irf_abs(Rperpx,1);
Rperpx = Rperpx./[Rmag Rmag Rmag];

Epar = dot(Rpar,Exyz.data,2);
Eperp = dot(Rperpx,Exyz.data,2);
Eperp2 = dot(Rperpy,Exyz.data,2);

Bpar = dot(Rpar,Bscm.data,2);
Bperp = dot(Rperpx,Bscm.data,2);
Bperp2 = dot(Rperpy,Bscm.data,2);

Efac = TSeries(Exyz.time,[Eperp Eperp2 Epar],'to',1);
Bfac = TSeries(Bscm.time,[Bperp Bperp2 Bpar],'to',1);
Bpar = TSeries(Bscm.time,Bpar,'to',1);

ne = 1*1e6;
qe = 1.6e-19;
mu = 4*pi*1.e-7;
Bmag = Bxyz.abs.data;
phiB = (Bpar.data).*Bmag*1e-18/(ne*qe*mu);



[Emva,El,Ev] = irf_minvar(Efac);
Emax = TSeries(Emva.time,Emva.data(:,1),'to',1);
t = Emax.time.epochUnix;
Emax = Emax.data/1000;
Emax = [t Emax];

vphvec = [1e4:1e4:1e6];
phiE = irf_integrate(Emax);
corr = zeros(1,length(vphvec));


for ii=1:length(vphvec);
    phiEtemp = phiE(:,2)*vphvec(ii);
    corr(ii)=sum(abs(phiEtemp-phiB).^2);
end

[maxcorr,idx] = min(corr);
vphbest = vphvec(idx)
phiEbest = vphbest*phiE(:,2);
irf_plot([t phiB phiEbest])

xcorr(phiB,phiEbest,0,'coeff')
