%% Calculate som derived quantities
Ld3=irf_plasma_calc(diB3,peaNe3,irf_tappl(peaNe3,'*0'),eVTe3,codTi3eV,'Ld');
fpe3=irf_plasma_calc(diB3,peaNe3,irf_tappl(peaNe3,'*0'),eVTe3,codTi3eV,'Fpe');
fce3=irf_plasma_calc(diB3,peaNe3,irf_tappl(peaNe3,'*0'),eVTe3,codTi3eV,'Fce');
fpefce3=irf_multiply(1,fpe3,1,fce3,-1);
units=irf_units;
gamma=irf_tappl(fpe3,'*(9.1e-31/1.67e-27)^(1/3)');
%% make overview figure

h=irf_plot(9);

isub=1;
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diB3)
    ylabel(hca,'B_{ISR2} [nT]')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diE3)
    ylabel(hca,'E_{ISR2} [mV/m]')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,eVTe3)
    ylabel(hca,'T3_e [eV]')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,peaNe3)
    ylabel(hca,'N3_e [cc]')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,Ld3)
    ylabel(hca,'L3_{de} [m]')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,fpe3)
    ylabel(hca,'f_{pe} [Hz]')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,gamma)
    ylabel(hca,'gamma [Hz]')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,fce3)
    ylabel(hca,'f_{ce} [Hz]')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 1
    hca=h(isub); isub=isub+1;
    irf_plot(hca,fpefce3)
    ylabel(hca,'f_{pe}/f_{ce}')
    set(hca,'ylim',[0 5])
    irf_zoom('y')
    irf_legend(hca,{'C3'},[0.02 0.9])
end
if 0
    hca=h(isub); isub=isub+1;
    irf_plot(hca,fcefpe3)
    ylabel(hca,'f_{ce}/f_{pe}')
    irf_legend(hca,{'C3'},[0.02 0.9])
end

irf_plot_axis_align
t_zoom=toepoch([2007 08 31 10 17 00;2007 08 31 10 18 30])';
irf_zoom(h,'x',t_zoom)
for k=1:numel(h); grid(h(k),'off'); end
tints=[2007 08 31 10 17 38.125;2007 08 31 10 17 39.625;2007 08 31 10 17 42.375;2007 08 31 10 17 44.000;2007 08 31 10 17 46.625;2007 08 31 10 17 48.125];
tints_ep=toepoch(tints)';
for k=1:numel(tints_ep); irf_pl_mark(h,[tints_ep(k)-0.06,tints_ep(k)+0.06]); end