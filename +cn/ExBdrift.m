% Calculates the ExB drift in the field aligned coordinate system
% Assumes diB3, diB4 are loaded.

t1 = [2007 08 31 10 17 00];
t2 = [2007 08 31 10 18 30];
tint = toepoch([t1;t2])';
sclist = 3:4;

if ~exist('diExB3','var');
    c_eval('diExB?=c_caa_var_get(''v_drift_ISR2__C4_CP_EFW_L2_V3D_INERT'',''mat'');',sclist);
end

% So, should I make an average coordinate system or an ever changing one.
% I think average is best, especially since its approximately constant
% during the electron holes time period.

mB=cn.mean([diB3(:,1:4);diB4(:,1:4)],1);
x=[1 0 0];

z=irf_norm(mB);
y=irf_norm(cross(z,x));
x=cross(y,z);

if ~exist('facExB3','var')
c_eval('facExB?=irf_newxyz(diExB?,x,y,z);',sclist);
tref=toepoch([2007 08 31 10 10 00]);
c_eval('facDExB?=irf_integrate(facExB?);',sclist);
c_eval('facDExBabs?=[facDExB?(:,1) sqrt(facDExB?(:,2).*facDExB?(:,2)+facDExB?(:,3).*facDExB?(:,3))];',sclist);
end
%%
%load matlabdiB
%%
h=irf_plot(10); 
iSub=1;

if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diB3); hold(hca,'on'); 
    yLabelString={'B_{ISR2}','[nT] C3'};
    ylabel(hca,yLabelString)
    irf_legend(hca,{'x','y','z','abs'},[0.2 0.1])
end
if 0
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,wletEx3); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
    text(hca,1,fce(1,2)*1e-3,'f_{ce}');
    set(hca,'yscale','log')
    ylabel(hca,'f [Hz/kHz?] (Ex3)')
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,irf_tlim(diE3,tint)); hold(hca,'on'); 
    yLabelString={'E_{ISR2}','[mV/m] C3'}; 
    ylabel(hca,yLabelString)
    irf_legend(hca,{'x','y'},[0.9 0.1])
end
if 0
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diE4); hold(hca,'on'); 
    ylabel(hca,'E [mV/m] C4')
end

if 1 % 3 panels
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,irf_tlim(facDExB3(:,[1 2]),tint)); hold(hca,'on'); 
    yLabelString={'int(v_{X,ExB})','[km] C3'};
    ylabel(hca,yLabelString)
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,irf_tlim(facDExB3(:,[1 3]),tint)); hold(hca,'on'); 
    yLabelString={'int(v_{Y,ExB})','[km] C3'};
    ylabel(hca,yLabelString)
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,irf_tlim(facDExB3(:,[1 4]),tint)); hold(hca,'on'); 
    yLabelString={'int(v_{Z,ExB})','[km] C3'};
    ylabel(hca,yLabelString)
end
if 0
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diE4); hold(hca,'on'); 
    ylabel(hca,'E [mV/m] C4')
end
if 1 % panel to delete afterwards
    iDelete=iSub;
    iZoom=iSub-1;
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,[1 1 1 1]); hold(hca,'on'); 
end
if 1 % 3 panels
    hca=h(iSub); iSub=iSub+1;
    hh(1)=irf_plot(hca,irf_tlim(facDExB3(:,[1 2]),tint_z)); hold(hca,'on'); 
    yLabelString={'int(v_{X,ExB})','[km] C3'};
    ylabel(hca,yLabelString)
    hca=h(iSub); iSub=iSub+1;
    hh(2)=irf_plot(hca,irf_tlim(facDExB3(:,[1 3]),tint_z)); hold(hca,'on'); 
    yLabelString={'int(v_{Y,ExB})','[km] C3'};
    ylabel(hca,yLabelString')
    hca=h(iSub); iSub=iSub+1;
    hh(3)=irf_plot(hca,irf_tlim(facDExB3(:,[1 4]),tint_z)); hold(hca,'on'); 
    yLabelString={'int(v_{Z,ExB})','[km] C3'};
    ylabel(hca,yLabelString)
end
if 1 % add absolute value of perpendicular distance
    hca=h(iSub); iSub=iSub+1;
    hh(4)=irf_plot(hca,irf_tlim(facDExBabs3,tint_z)); hold(hca,'on'); 
    yLabelString={'int(v_{ExB}) [km]','(X^2+Y^2)^{1/2} C3'};
    ylabel(hca,yLabelString)
end

if 1 % markings
    [tint_eh nscobs comments] = eh_tint; % load time intervals
    ax=h(:); % plot to be marked is
    % seen on one/seen on 2/developing seen on one/developing seen on two
    color=cn.colors(7);
    for k = 1:size(tint_eh,1)
        color{nscobs(k)};
        irf_pl_mark(ax,tint_eh{k},color{nscobs(k)})
    end
end
irf_zoom(h(1:5),'x',tint)
zoomin=1;
if zoomin
    t1=[2007 08 31 10 17 35];
    t2=[2007 08 31 10 17 50];
    tint_z=toepoch([t1;t2])';
    irf_zoom(hh,'x',tint_z)
end
irf_plot_zoomin_lines_between_panels(h(iDelete-1),h(iDelete+1))
delete(h(iDelete));

title(h(1),['X_{FAC}=',cn.vector_to_string(x,1),...
          '  Y_{FAC}=',cn.vector_to_string(y,1),...
          '  Z_{FAC}=',cn.vector_to_string(z,1)])
irf_plot_axis_align
set(gcf,'paperpositionmode','auto');

cn.print('integrated_ExB_fac_abs')














