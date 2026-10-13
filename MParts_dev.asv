PDc = iPDist_counts.movmean(nMovMean);
PD_ = PD_use(1:2);
PDc_ = PDc(1:2);
MParr = PD_.macroparticles('ntot',1e6,'skipzero',1,'v_init','v3','counts',PDc_,'Mparts');
%%
farr = arrayfun(@(x){x.df},MParr);
dcarr = arrayfun(@(x){x.dc},MParr);
dnarr = arrayfun(@(x){x.df.*x.dv},MParr);
vxarr = arrayfun(@(x){x.vx},MParr);
vyarr = arrayfun(@(x){x.vy},MParr);
vzarr = arrayfun(@(x){x.vz},MParr);
dvarr = arrayfun(@(x){x.dv},MParr);
idep1arr = arrayfun(@(x){x.iDep1},MParr);
idep2arr = arrayfun(@(x){x.iDep2},MParr);
idep3arr = arrayfun(@(x){x.iDep3},MParr);

%MP = MParts(PD.time,farr,'','3D-cart',vxarr,vyarr,vzarr,'dv',dvarr,'idepparent',{idep1arr,idep2arr,idep3arr});



MP = MParts(PD_.time,'3D-cart',vxarr,vyarr,vzarr,'f',farr,'dV',dvarr,'dc',dcarr,'dn',dnarr,'idepparent',{idep1arr,idep2arr,idep3arr});

%%
PD.length