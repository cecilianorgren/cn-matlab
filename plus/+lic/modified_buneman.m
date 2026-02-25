t3 = toepoch([2007 08 31 10 17 38.14]); 
t4 = toepoch([2007 08 31 10 17 39.70]);

close
%n = (0.035)

%[h,f,vp,vz]=whamp.plot_f(n,m,t,vd,d,a1,a2);

for kk=1:2; h(kk)=subplot(1,2,kk); end

isub = 1;
if 1
    hca = h(isub); isub = isub + 1; hold(hca,'on');
    n = [0.07 0.002]*1e6;
    m = [0 0];
    T = [2500 60]*1e-3;
    vd = [0.0 -4.1];
    d = [1 1];
    a1 = [1 1];
    a2 = [0 0];
    axes(hca)
    whamp.plot_f(n,m,T,vd,d,a1,a2,[0 90 180],1);
    c_caa_plot_distribution_function(hca,'tint',t3,'cross-section',peaPSD3);
    hold(hca,'off')
end
if 1
    hca = h(isub); isub = isub + 1; hold(hca,'on');
    n = [0.07 0.003]*1e6;
    m = [0 0];
    T = [2000 60]*1e-3;
    vd = [0.0 -3.3];
    d = [1 1];
    a1 = [1 1];
    a2 = [0 0];
    axes(hca)
    whamp.plot_f(n,m,T,vd,d,a1,a2,[0 90 180],1);
    c_caa_plot_distribution_function(hca,'tint',t4,'cross-section',peaPSD4,'emin',1.5e2);   
end
for  kk=1:2; set(h(kk),'xlim',[0 0.5]*1e4,'ylim',[0 7],'xscale','lin','yscale','lin'); end

%%
close
nPanels = 1;
for kk=nPanels; h(kk)=subplot(1,nPanels,kk); end
isub = 1;

if 1
    hca = h(isub); isub = isub + 1;
    n = [0.07 0.003]*1e6;
    m = [0 0];
    T = [2000 60]*1e-3;
    vd = [0.0 -3.3];
    d = [1 1];
    a1 = [1 1];
    a2 = [0 0];
    axes(hca)
    whamp.plot_f(n,m,T,vd,d,a1,a2,[0 180]);
    hold(gca,'on')
end
if 1
    isub = isub - 1;
    hca = h(isub); isub = isub + 1;
    n = [0.07 0.003]*1e6;
    m = [0 0];
    T = [2000 60]*1e-3;
    vd = [0.0 -6];
    d = [1 1];
    a1 = [1 1];
    a2 = [0 0];
    axes(hca)
    whamp.plot_f(n,m,T,vd,d,a1,a2,[0 180]);    
end
for  kk=1; 
    set(h(kk),'xlim',[0 0.6]*1e5,'ylim',[0 7],'xscale','lin','yscale','lin'); 
    grid(h(kk),'off')
    box(h(kk),'on')
    title(h(kk),' ')
    xlabel(h(kk),'v [10^3 km/s]')
    xticks = tocolumn(get(h(kk),'xtick'));
    xticklabels = num2str(xticks*1e-3);
    set(h(kk),'xticklabel',xticklabels);
end

axesObjs = get(h(1), 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes  
%
% axesObjs(4) - slower beam, if 0 90 180
% axesObjs(1) - faster beam, if 0 90 180
%
if 0
    delete(axesObjs(3))
    delete(axesObjs(2))
    set(axesObjs(1),'linestyle','--')
else
    %delete(axesObjs(2))
    %delete(axesObjs(4))
    set(axesObjs(1),'linestyle','--')
end