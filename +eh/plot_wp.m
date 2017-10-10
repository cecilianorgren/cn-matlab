% plots e field, locations of holes and 
% NOT TRUE: products Ni3, Ni4, diE3, diE4 should already be loaded and ready
% First calculate wpe, wpi
tic
tint_ov=eh.load_tint_ov;
[tints, quality, comments, color]=eh.load_tint;
cd /Users/Cecilia/Data/BM/20070831;
load matlabdiEB; % loads electric and magnetic field in ISR2 coordinates
load matlabSC; %loads spacecraft potential and position
load matlabN;
toc
% Calculate plasma frequencies
units=irf_units;
c_eval('wpe?=[peaNe?(:,1) sqrt(units.e^2*peaNe?(:,2:end)/units.mu0/units.me)];',3:4);
c_eval('wpi?=[peaNe?(:,1) sqrt(units.e^2*peaNe?(:,2:end)/units.mu0/units.mp)];',3:4);
%%
nPanels=4;
h=irf_plot(nPanels);
iSub=1;
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,irf_tlim(diE3,tint_ov));
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,irf_tlim(diE4,tint_ov));
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,{wpi3,wpi4},'comp');
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,{wpe3,wpe4},'comp');
end
%if 1
%    hca=h(iSub); iSub=iSub+1;
%    irf_plot(hca,P4);
%end
irf_zoom(h,'x',tint_ov)
if 1 % mark eh's
% mark all eh's, and especially the presented one.
% prio 1/2/3/4/5/rubbish/dl/presented interval
    for p = 1:size(tints,1)        
        %cin=scncobs(p);        
        irf_pl_mark(h,tints{p},color{quality(p)})               
    end    
end