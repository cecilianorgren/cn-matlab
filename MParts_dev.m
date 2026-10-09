PD = PD_use(1:50);
MParr = PD.macroparticles('ntot',1e6,'skipzero',1,'v_init','E');
%%
farr = arrayfun(@(x){x.df},MParr);
vxarr = arrayfun(@(x){x.vx},MParr);
vyarr = arrayfun(@(x){x.vy},MParr);
vzarr = arrayfun(@(x){x.vz},MParr);
dvarr = arrayfun(@(x){x.dv},MParr);
idep1arr = arrayfun(@(x){x.iDep1},MParr);
idep2arr = arrayfun(@(x){x.iDep2},MParr);
idep3arr = arrayfun(@(x){x.iDep3},MParr);

%MP = MParts(PD.time,farr,'','3D-cart',vxarr,vyarr,vzarr,'dv',dvarr,'idepparent',{idep1arr,idep2arr,idep3arr});



MP = MParts(PD.time,farr,'','3D-cart',vxarr,vyarr,vzarr,'dv',dvarr,'idepparent',{idep1arr,idep2arr,idep3arr});

%%
PD.length