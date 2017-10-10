% make wavelet of "forming" and formed electron holes to compare
t1=[2007 08 31 10 17 00];
t2=[2007 08 31 10 18 30];
tint=toepoch([t1;t2])';
sclist=3:4;
nf=50;
c_eval('wletEx?=irf_wavelet(irf_tlim(diE?(:,[1 2]),tint));',sclist)
c_eval('wletEy?=irf_wavelet(irf_tlim(diE?(:,[1 3]),tint));',sclist)
%% make plasma parameters
B=irf_add(0.5,irf_tlim(diB3,tint),0.5,irf_tlim(diB4,tint));

fce=irf_plasma_calc(B,0.1,0,100,100,'Fce');
fpe=irf_plasma_calc(B,0.1,0,100,100,'Fpe');
fcp=irf_plasma_calc(B,0.1,0,100,100,'Fcp');
fpp=irf_plasma_calc(B,0.1,0,100,100,'Fpp');
flh=irf_plasma_calc(B,0.1,0,100,100,'Flh');

%% Check if there is no large scale ExB drift
fhigh=5;
c_eval('diE?lf=irf_filt(diE?,0,fhigh);',sclist)
c_eval('diE?llf=irf_filt(diE?,0,2);',sclist)
%c_eval('Babs?=diB?; Babs(:,2:4)=repmat(diB?(:,5),1,3);',sclist);
%c_eval('diExB?=irf_cross(diE?llf,diB?);'
c_eval('diExB?=c_caa_var_get(''v_drift_ISR2__C4_CP_EFW_L2_V3D_INERT'',''mat'');',sclist);
%%
c_eval('diDr?=irf_integrate(diExB?);',sclist);
c_eval('diDrz?=irf_tlim(diDr?,tint_z);',sclist)
hh=irf_plot({diDr3,diDr4},'comp');
irf_zoom(hh,'x',tint_z)
%% make field aligned coordinate system
%% plot
h=irf_plot(6);
iSub=1;
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
    irf_plot(hca,wletEy3); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
    text(hca,1,fce(1,2)*1e-3,'f_{ce}');
    set(hca,'yscale','log')
    ylabel(hca,'f [Hz/kHz?] (Ey3)')
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,wletEx4); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
    text(hca,1,fce(1,2)*1e-3,'f_{ce}');
    set(hca,'yscale','log')
    ylabel(hca,'f [Hz/kHz?] (Ex4)')
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,wletEy4); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
    text(hca,1,fce(1,2)*1e-3,'f_{ce}');
    set(hca,'yscale','log')
    ylabel(hca,'f [Hz/kHz?] (Ey4)')
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diE3); hold(hca,'on'); 
    ylabel('E [mV/m] C3') 
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diE4); hold(hca,'on'); 
    ylabel('E [mV/m] C4')
end

if 1 % markings
    [tint_eh nscobs comments] = eh_tint; % load time intervals
    ax=h(5:6); % plot to be marked is
    % seen on one/seen on 2/developing seen on one/developing seen on two
    color=cn.colors(7);
    for k = 1:size(tint_eh,1)
        color{nscobs(k)};
        irf_pl_mark(ax,tint_eh{k},color{nscobs(k)})
    end
end
irf_zoom(h,'x',tint)
irf_plot_axis_align

%% plot exB drift
h=irf_plot(4);
iSub=1;
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diE3); hold(hca,'on'); 
    ylabel('E [mV/m] C3') 
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diE4); hold(hca,'on'); 
    ylabel('E [mV/m] C4')
end

if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diE3); hold(hca,'on'); 
    ylabel('E [mV/m] C3') 
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,diE4); hold(hca,'on'); 
    ylabel('E [mV/m] C4')
end

if 1 % markings
    [tint_eh nscobs comments] = eh_tint; % load time intervals
    ax=h(5:6); % plot to be marked is
    % seen on one/seen on 2/developing seen on one/developing seen on two
    color=cn.colors(7);
    for k = 1:size(tint_eh,1)
        1
        color{nscobs(k)};
        irf_pl_mark(ax,tint_eh{k},color{nscobs(k)})
    end
end
irf_zoom(h,'x',tint)
irf_plot_axis_align
