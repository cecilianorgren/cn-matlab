% 2014-02-04
% 1. Plots simultaneously B_FGM and B_FGMSTAFF in order to compare them.
% 2. Creates average fac system to check closer on wave B_par.

%% Compare fgm to staff
t1 = [2007 08 31 10 12 0.0]; t2 = [2007 08 31 10 26 0.0]; tint = toepoch([t1;t2])';
h=irf_plot(5); isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(diE3,tint)); irf_legend(hca,{'x','y','z'},[0.95 0.92]); ylabel(hca,'E_{ISR2}')
hca=h(isub);isub=isub+1;irf_plot(hca,B3staff); irf_legend(hca,{'x','y','z'},[0.95 0.92]); ylabel(hca,'B_{STAFF,GSM}')
hca=h(isub);isub=isub+1;irf_plot(hca,{gsmB3(:,[1 2]),gsmBfs3(:,[1 2])},'comp'); irf_legend(hca,{'FGM','FGM+STAFF'},[0.95 0.92]); ylabel(hca,'B_{x,GSM}') 
hca=h(isub);isub=isub+1;irf_plot(hca,{gsmB3(:,[1 3]),gsmBfs3(:,[1 3])},'comp'); irf_legend(hca,{'FGM','FGM+STAFF'},[0.95 0.92]); ylabel(hca,'B_{y,GSM}')
hca=h(isub);isub=isub+1;irf_plot(hca,{gsmB3(:,[1 4]),gsmBfs3(:,[1 4])},'comp'); irf_legend(hca,{'FGM','FGM+STAFF'},[0.95 0.92]); ylabel(hca,'B_{z,GSM}')
art2.markbeam(h)
irf_zoom(h,'x',tint)

%% Create fac staff and fac combined
% vectors in gse-system:
% first - in spin plane, perp to B  
% second - BxSP
% third - par to B, average between C3 and C4

c_eval('b1? = irf_dot(first,gseBfs?);',3:4)  % in spin plane, perp to B  
c_eval('b2? = irf_dot(second,gseBfs?);',3:4) % BxSP
c_eval('b3? = irf_dot(third,gseBfs?);',3:4)  % par to B
c_eval('facBfs? = [b1? b2?(:,2) b3?(:,2)];',3:4); 

c_eval('b1? = irf_dot(first,gseB?staff);',3:4)  % in spin plane, perp to B  
c_eval('b2? = irf_dot(second,gseB?staff);',3:4) % BxSP
c_eval('b3? = irf_dot(third,gseB?staff);',3:4)  % par to B
c_eval('facB?staff = [b1? b2?(:,2) b3?(:,2)];',3:4); 

c_eval('b1? = irf_dot(first,gseB?);',3:4)  % in spin plane, perp to B  
c_eval('b2? = irf_dot(second,gseB?);',3:4) % BxSP
c_eval('b3? = irf_dot(third,gseB?);',3:4)  % par to B
c_eval('facB? = [b1? b2?(:,2) b3?(:,2)];',3:4); 

clear b13 b23 b33 b14 b24 b34

%% Plot fac's to compare
t1 = [2007 08 31 10 12 0.0]; t2 = [2007 08 31 10 26 0.0]; tint = toepoch([t1;t2])';
h=irf_plot(5); isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(diE3,tint)); irf_legend(hca,{'x','y','z'},[0.95 0.92]); ylabel(hca,'E_{ISR2}')
hca=h(isub);isub=isub+1;irf_plot(hca,facB3staff); irf_legend(hca,{'x','y','z'},[0.95 0.92]); ylabel(hca,'B_{STAFF,FAC}')
hca=h(isub);isub=isub+1;irf_plot(hca,{facB3(:,[1 2]),facBfs3(:,[1 2])},'comp'); irf_legend(hca,{'FGM','FGM+STAFF'},[0.95 0.92]); ylabel(hca,'B_{x,FAC}') 
hca=h(isub);isub=isub+1;irf_plot(hca,{facB3(:,[1 3]),facBfs3(:,[1 3])},'comp'); irf_legend(hca,{'FGM','FGM+STAFF'},[0.95 0.92]); ylabel(hca,'B_{y,FAC}')
hca=h(isub);isub=isub+1;irf_plot(hca,{facB3(:,[1 4]),facBfs3(:,[1 4])},'comp'); irf_legend(hca,{'FGM','FGM+STAFF'},[0.95 0.92]); ylabel(hca,'B_{z,FAC}')
art2.markbeam(h)
irf_zoom(h,'x',tint)
title(h(1),'x - sp, perp to B; y - BxSP; z - par to B')

%% Compare STAFF fac to STAFF gse
t1 = [2007 08 31 10 12 0.0]; t2 = [2007 08 31 10 26 0.0]; tint = toepoch([t1;t2])';
h=irf_plot(3); isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,irf_tlim(diE3,tint)); irf_legend(hca,{'x','y','z'},[0.95 0.92]); ylabel(hca,'E_{ISR2}')
hca=h(isub);isub=isub+1;irf_plot(hca,facB3staff); irf_legend(hca,{'x','y','z'},[0.95 0.92]); ylabel(hca,'B_{STAFF,FAC}')
hca=h(isub);isub=isub+1;irf_plot(hca,gseB3staff); irf_legend(hca,{'x','y','z'},[0.95 0.92]); ylabel(hca,'B_{STAFF,GSE}')
art2.markbeam(h)
irf_zoom(h,'x',tint)





