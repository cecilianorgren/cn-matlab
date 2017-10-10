% hump script
cd /Users/Cecilia/Data/BM/20070831/
%% Get B in dsi coordinates
c_eval('diB?=c_coord_trans(''gsm'',''dsi'',gsmB?,''cl_id'',?);',3:4);
%% Get Epar if E.B=0 and E.B~0
c_eval('[diE?t,angle?t]=irf_edb(diE?,diB?,0);',3:4);
c_eval('[diE?p,angle?p]=irf_edb(diE?,diB?,90,''Epar'');',3:4);

%% Plot fields and angle, B, E and sc separately
t1=[2007 08 31 10 18 00.00];  
t2=[2007 08 31 10 20 00.00]; 
tint=[toepoch(t1) toepoch(t2)];
figure;
h=irf_plot(7);
isub=1;

if 1 % diB3
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diB3,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diB3')
end
if 1 % diB4
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diB4,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diB4')
end
if 1 % diE3t
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE3t,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diE3_{\perp}')
end
if 1 % diE4t
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE4t,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diE4_{\perp}')
end
if 1 % diE3p
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE3p,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diE3_{||}')
end
if 1 % diE4p
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE4p,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diE4_{||}')
end
if 1 % angle
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim([diE3t(:,1),angle3],tint),'k');hold(hca,'on');
    irf_plot(hca,irf_tlim([diE4t(:,1),angle4],tint),'r');hold(hca,'on'); 
    irf_legend(hca,{'C3','C4'},[0.95 0.95]);
    ylabel('Angle (BtoSpinPlane)')
end

irf_zoom(h,'x',tint)
%% Plot fields and angle, comparing between sc
t1=[2007 08 31 10 18 00.00];  
t2=[2007 08 31 10 20 00.00]; 
tint=[toepoch(t1) toepoch(t2)];

h=irf_plot(7);
isub=1;

if 1 % diB3
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diB3,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diB3')
end
if 1 % diB4
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diB4,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diB4')
end
if 1 % diE3t
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE3t,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diE3_{\perp}')
end
if 1 % diE4t
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE4t,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diE4_{\perp}')
end
if 1 % diE3p
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE3p,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diE3_{||}')
end
if 1 % diE4p
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE4p,tint));
    irf_legend(hca,{'x','y','z'},[0.95 0.95]);
    ylabel(hca,'diE4_{||}')
end
if 1 % angle
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim([diE3t(:,1),angle3],tint),'k');hold(hca,'on');
    irf_plot(hca,irf_tlim([diE4t(:,1),angle4],tint),'r');hold(hca,'on'); 
    irf_legend(hca,{'C3','C4'},[0.95 0.95]);
    ylabel('Angle (BtoSpinPlane)')
end

irf_zoom(h,'x',tint)

%% Localized MVA on hump
diMv3=ud.v3;
