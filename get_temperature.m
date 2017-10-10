if ~exist('peaPSD3','var'); load peaPSD; end

set(0,'DefaultAxesLineStyleOrder',{'*',':'});
t = toepoch([2007 08 31 10 17 51.75]);
plot_type = 'cross-section';
pitchangles = [0 75 105 135 165];
[ax ax_cb] = c_caa_plot_distribution_function('tint',t,plot_type,...
            peaPSD4,'pitchangle',pitchangles);
        
hold(ax,'on');

%
% Model
    
    w_n=  [0.08 0]*1e6;       % density m^3
    w_t=  [2 0];            % temperature [keV]
    w_vd= [0 0];         % V_drift/V_term
    w_a1= [1 1];           % T_perp/T_par
    w_m=  [0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
    w_d=  [1 1];            % loss cone parameter, 1.0 = no loss cone)
    w_a2= [0 0];           % ??? (use 0)
    w_pa= [0 75 105 135 165];    % pitch angeles
    plotoption=[1];
    titleoption=[0];
    cmp=[1 ];


% Plot
set(gca,'NextPlot','add','LineStyleOrder',{'-',':'});   
plotoption = 1;
[h,f,Etot]=whamp.plot_f(w_n(cmp),w_m(cmp),w_t(cmp),w_vd(cmp),...
                            w_d(cmp),w_a1(cmp),w_a2(cmp),w_pa,...
                            plotoption);

info = {['n1 = ' num2str(w_n(cmp(1))*1e-6) ' cc'],...
        ['T1 = ' num2str(w_t(cmp(1))*1e3) ' eV'],...    
        ['vte1 = ' num2str(cn_eV2v(w_t(cmp(1))*1e3,'eV'),'%.0f') 'km/s'],...    
        ['v1 = ' num2str(w_vd(cmp(1))) ' vte1'],...
       % ['n2 = ' num2str(w_n(cmp(2))*1e-6) ' cc'],...
       % ['T2 = ' num2str(w_t(cmp(2))*1e3) ' eV'],...
       % ['Tperp2/Tpar2 = ' num2str(w_a1(cmp(2)))],...
       % ['vte2 = ' num2str(cn_eV2v(w_t(cmp(2))*1e3,'eV'),'%.0f') 'km/s'],...    
       % ['v2 = ' num2str(w_vd(cmp(2))) ' vte2'],...        
       % ['v2 = ' num2str(w_vd(cmp(2))*cn_eV2v(w_t(cmp(2)),'eV')/cn_eV2v(w_t(cmp(1)),'eV'),'%.1f') ' vte1'],...    
       % ['n2/(n1+n2) = ' num2str(w_n(cmp(2))/(w_n(cmp(1))+w_n(cmp(2))),'%.2f')],...    
};
    text(5e1,1e-2,info)
    title(irf_time(t,'epoch2iso'))
    grid(ax,'off')
hold(ax,'off') 

%% make surf plot
figure(2)
[h,f,vp,vz]=whamp.plot_f(w_n(cmp),w_m(cmp),w_t(cmp),w_vd(cmp),...
                         w_d(cmp),w_a1(cmp),w_a2(cmp));
hcb = colorbar;

%% 
surf(vp,vz,f)
hcb = colorbar;
shading flat;

%% integrate
vp = 10000; % km/s
% vz(1,:) are identical.

for ii = 1:121
    curve = fit(vz(ii,1),f(ii,1),'exp2');
end
    
    
    
    


