% Compares the locations og eh's to density changes
tic
tint_ov=eh.load_tint_ov;
[tints, quality, comments, color]=eh.load_tint;
cd /Users/Cecilia/Data/BM/20070831;
load matlabdiEB; % loads electric and magnetic field in ISR2 coordinates
load matlabSC; %loads spacecraft potential and position
toc
%%
% Construct parallel E
tic
c_eval('[diE?,angle?]=irf_edb(diE?,diB?,90,''Epar'');',3:4);
toc
c_eval('Bplane?=irf_cross(diB?,[diB?(:,[1 2 3]) 0*diB?(:,1)]);',3:4)
toc
%c_eval('Ep?=irf_newxyz(diE?,Bplane?,0,diB?(:,1:4))',3:4); % [diB?(:,1) diB?(:,2:4)*0]
%c_eval('Ep?=cn_xyz(diE?,Bplane?,0,diB?)',3:4); % [diB?(:,1) diB?(:,2:4)*0]
%toc
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
    irf_plot(hca,P3);
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,irf_tlim(diE4,tint_ov));
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,P4);
end
irf_zoom(h,'x',tint_ov)
if 1 % mark eh's
% mark all eh's, and especially the presented one.
% prio 1/2/3/4/5/rubbish/dl/presented interval
    for p = 1:size(tints,1)        
        %cin=scncobs(p);        
        irf_pl_mark(h,tints{p},color{quality(p)})               
    end    
end
