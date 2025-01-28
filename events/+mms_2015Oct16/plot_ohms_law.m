%% Get pressure gradient from velocity
% Take only the diagonal pressure term
% Rotate tensor into frame aligned with a trajectory
% Pick out normal direction
% mms_2015Oct16.coordinate_system
tint = irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z');
[out,l,v] = irf_minvar(dmpaB1.tlim(tint));

L = v(1,:);
M = v(2,:); 
N = v(3,:);

units = irf_units;
velocity = 10;
c_eval('Peall? = TSeries(Pe?.time,[Pe?.data(:,1,1) Pe?.data(:,1,3) Pe?.data(:,1,3) Pe?.data(:,2,2) Pe?.data(:,2,3) Pe?.data(:,3,3)]);');
%c_eval('[PeXXp?,PeXYp?,PeXZp?,PeYYp?,PeYZp?,PeZZp?] = mms.rotate_tensor(Peall?,''rot'',N);')
c_eval('[Pe?rot] = mms.rotate_tensor(Peall?,''rot'',N,L,M); Pe?rot = irf.ts_tensor_xyz(Pe?rot.time,Pe?rot.data);')

c_eval('Teall? = TSeries(Te?.time,[Te?.data(:,1,1) Te?.data(:,1,3) Te?.data(:,1,3) Te?.data(:,2,2) Te?.data(:,2,3) Te?.data(:,3,3)]);');
c_eval('[Te?rot] = mms.rotate_tensor(Teall?,''rot'',N,L,M); Te?rot = irf.ts_tensor_xyz(Te?rot.time,Te?rot.data);')

c_eval('dx? = velocity*(PeXXp?.time(2)-PeXXp?.time(1));')
c_eval('dP? = TSeries(Pe?rot.time(1:(Pe?rot.length-1)),Pe?rot.xx.data(2:Pe?rot.length,:)-Pe?rot.xx.data(1:(Pe?rot.length-1),:));')
c_eval('gradPe? = dP?/dx?;')

c_eval('dPdx? = dP?/dx?;')
%%
h = irf_plot(8); isub = 1;

hca = h(isub); isub = isub + 1;
irf_plot(hca,'ne?','comp')
ylabel(hca,{['n_e'],' (cm^{-3})'})

hca = h(isub); isub = isub + 1;
irf_plot(hca,{Te1rot.xx,Te2rot.xx,Te3rot.xx,Te4rot.xx},'comp') 
ylabel(hca,'T_{NN} (eV)')

hca = h(isub); isub = isub + 1;
irf_plot(hca,{Pe1rot.xx,Pe2rot.xx,Pe3rot.xx,Pe4rot.xx},'comp')
ylabel(hca,'P_{NN} (nPa)')

hca = h(isub); isub = isub + 1;
irf_plot(hca,gradPe)
ylabel(hca,{'grad(P)','(nPa/km)'})

hca = h(isub); isub = isub + 1;
irf_plot(hca,'dP?','comp')
ylabel(hca,'diff(P_{NN}) (nPa)')

hca = h(isub); isub = isub + 1;
irf_plot(hca,'dPdx?','comp')
ylabel(hca,{['dP_{NN}/dx'],' (nPa/km)'})
irf_legend(hca,['dx = ' num2str(dx1,'%.2f') ' km'],[0.01 0.95])

hca = h(isub); isub = isub + 1;
irf_plot(hca,{dPdx1,dPdx2,dPdx3,dPdx4,(dPdx1+dPdx2.resample(dPdx1.time)+dPdx3.resample(dPdx1.time)+dPdx4.resample(dPdx1.time))/4,gradPe},'comp')
ylabel(hca,{['grad(P)'],' (nPa/km)'})
irf_legend(hca,{'\Delta P_1/\Delta x','\Delta P_2/\Delta x','\Delta P_3/\Delta x','\Delta P_4/\Delta x','\Delta P/\Delta x (av)','\nabla\cdot P (4sc)'},[0.01 0.95])

hca = h(isub); isub = isub + 1;
irf_plot(hca,{(dPdx1+dPdx2.resample(dPdx1.time)+dPdx3.resample(dPdx1.time)+dPdx4.resample(dPdx1.time))/4,gradPe},'comp')
ylabel(hca,{['grad(P)'],' (nPa/km)'})
irf_legend(hca,{'\Delta P/\Delta x (av)','\nabla\cdot P (4sc)'},[0.01 0.95])
