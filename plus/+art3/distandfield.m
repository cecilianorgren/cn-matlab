% Distribution with fit + electric field with sweep marked yellow
cd /Users/Cecilia/Data/BM/20070831

% Select event.
event = 6;
switch event
    case 1; 
        time = [2007 08 31 10 17 38.158]; % center of beam ?
        % Electron distribution model       
        w_n=  [0.07 0.0075]*1e6; % density m^3
        w_t=  [1.6 0.06];       % temperature [keV]
        w_vd= [0 -4.1];           % V_drift/V_term
        w_a1= [1 1];            % T_perp/T_par        
    case 2; 
        time = [2007 08 31 10 17 39.675]; % center of beam
        % Electron distribution model       
        w_n=  [0.055 0.0075]*1e6; % density m^3
        w_t=  [1.6 0.060];       % temperature [keV]
        w_vd= [0 -3.5];           % V_drift/V_term
        w_a1= [1 1];            % T_perp/T_par  
    case 3; 
        time = [2007 08 31 10 17 42.375]; % center of smear
        % Electron distribution model       
        w_n=  [0.06 0.008]*1e6; % density m^3
        w_t=  [1.3 0.25];       % temperature [keV]
        w_vd= [0 -0.7];           % V_drift/V_term
        w_a1= [1 0.2];            % T_perp/T_par  
    case 4; 
        time = [2007 08 31 10 17 44.000]; % center of smear
        % Electron distribution model       
        w_n=  [0.055 0.012]*1e6; % density m^3
        w_t=  [1.6 0.18];       % temperature [keV]
        w_vd= [0 -1.7];           % V_drift/V_term
        w_a1= [1 0.2];            % T_perp/T_par  
    case 5; 
        time = [2007 08 31 10 17 46.400]; % center of smear ?
        % Electron distribution model       
        w_n=  [0.065 0.002]*1e6; % density m^3
        w_t=  [1.6 0.08];       % temperature [keV]
        w_vd= [0 -2];           % V_drift/V_term
        w_a1= [1 0.4];            % T_perp/T_par  
    case 6; 
        time = [2007 08 31 10 17 48.200]; % center of smear ? 
        % Electron distribution model       
        w_n=  [0.07 0.006]*1e6; % density m^3
        w_t=  [2.5 0.17];       % temperature [keV]
        w_vd= [0 -0.8];           % V_drift/V_term
        w_a1= [1 0.4];            % T_perp/T_par  
end
t = toepoch(time);

% Common parameters for model
w_m  = [0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
w_d  = [1 1];            % loss cone parameter, 1.0 = no loss cone)
w_a2 = [0 0];           % ??? (use 0)
w_pa = [15 45 75 105 135 165];    % pitch angeles
plotoption=[1];
titleoption=0;
cmp=[1 2];

% Distribution info
switch numel(cmp)
    case 1
        info = {['n1 = ' num2str(w_n(cmp(1))*1e-6) ' cc'],...
                ['T1 = ' num2str(w_t(cmp(1))*1e3) ' eV'],...    
                ['vte1 = ' num2str(cn_eV2v(w_t(cmp(1))*1e3,'eV'),'%.0f') 'km/s'],...    
                ['v1 = ' num2str(w_vd(cmp(1))) ' vte1'],...                       
                };
    case 2
        info = {['n1 = ' num2str(w_n(cmp(1))*1e-6) ' cc'],...
                ['T1 = ' num2str(w_t(cmp(1))*1e3) ' eV'],...    
                ['vte1 = ' num2str(cn_eV2v(w_t(cmp(1))*1e3,'eV'),'%.0f') 'km/s'],...    
                ['v1 = ' num2str(w_vd(cmp(1))) ' vte1'],...
                ['n2 = ' num2str(w_n(cmp(2))*1e-6) ' cc'],...
                ['T2 = ' num2str(w_t(cmp(2))*1e3) ' eV'],...
                ['Tperp2/Tpar2 = ' num2str(w_a1(cmp(2)))],...
                ['vte2 = ' num2str(cn_eV2v(w_t(cmp(2))*1e3,'eV'),'%.0f') 'km/s'],...    
                ['v2 = ' num2str(w_vd(cmp(2))) ' vte2'],...        
                ['v2 = ' num2str(w_vd(cmp(2))*cn_eV2v(w_t(cmp(2)),'eV')/cn_eV2v(w_t(cmp(1)),'eV'),'%.1f') ' vte1'],...    
                ['n2/(n1+n2) = ' num2str(w_n(cmp(2))/(w_n(cmp(1))+w_n(cmp(2))),'%.2f')],...    
                };
end            

% Set up figure.
figure(37)
set(gcf,'position',[10 400 450 550]);
h(1) = subplot(2,1,1);
h(2) = subplot(2,1,2);


set(h(1),'position',[0.13 0.37 0.775 0.58]);
set(h(2),'position',[0.13 0.12 0.775 0.17]);
%h(1) = axes('position',[0.13 0.37 0.775 0.58]); 
%h(2) = axes('position',[0.13 0.12 0.775 0.17]);

% Plot distribution.
if ~exist('peaPSD3','var'); load peaPSD; end

% Pick out time index
[tt,ind_t{1}]=irf_tlim(peaPSD3.t,t+[-1 1]);
ind_min=find(abs((tt-t))==min(abs((tt-t))));
ind_t{1}=ind_t{1}(ind_min); 
[tt,ind_t{2}]=irf_tlim(peaPSD4.t,t+[-1 1]);
ind_min=find(abs((tt-t))==min(abs((tt-t))));
ind_t{2}=ind_t{2}(ind_min); 
clear tt ind_min

% Pick out angles
pitchangles = peaPSD3.f_cs;
n_pa=length(pitchangles);  

% Pick energy levels
emin = min(min([peaPSD3.en_cs peaPSD3.en_cs]));
ind_en{1}=find(peaPSD3.en_cs>emin); % -1 in order to include emin
ind_en{2}=find(peaPSD4.en_cs>emin); % -1 in order to include emin

% Take away nans and squeeze
for ii = 1:n_pa;
    psd3{ii} = squeeze(mean(irf.nanmean(peaPSD3.p(ind_t{1},ii,ind_en{1}),1),2));
    psd4{ii} = squeeze(mean(irf.nanmean(peaPSD4.p(ind_t{2},ii,ind_en{2}),1),2));
    pa_legends{ii} = num2str(pitchangles(ii),'%.0f');
end
bg3 = squeeze(irf.nanmean(irf.nanmean(peaPSD3.p_bg(ind_t{1},:,ind_en{1}),1),2));
bg4 = squeeze(irf.nanmean(irf.nanmean(peaPSD4.p_bg(ind_t{2},:,ind_en{2}),1),2));
%pa_legends{n_pa+1}='Bg';

% Define colors
colors = {[0 0.8 0],[0 0 1],[0.8 0 0],[1 0.5 0],[0 1 1],[1 0 1]};
colorsarr=[]; for ii = 1:numel(colors); colorsarr = [colorsarr;colors{ii}]; end
set(0,'DefaultAxesColorOrder',colorsarr);

% Plot distributions
hold(h(1),'on')

for ii = 1:n_pa % the pitch angles
    loglog(h(1),peaPSD3.en_cs(ind_en{1}),psd3{ii},'color',colors{ii},'linestyle','*')   
end
for ii = 1:n_pa % the pitch angles    
    loglog(h(1),peaPSD4.en_cs(ind_en{2}),psd4{ii},'color',colors{ii},'linestyle','o')
end

loglog(h(1),peaPSD3.en_cs(ind_en{1}),bg3,'color',[0 0 0],'linestyle','--')
loglog(h(1),peaPSD4.en_cs(ind_en{2}),bg4,'color',[0 0 0],'linestyle','--')
set(h(1),'yscale','log','xscale','log')
legend(h(1),pa_legends,'edgecolor','w') 
ylabel(h(1),peaPSD3.p_label)
xlabel(h(1),'Energy  [eV]')

% Plot fit to distribution.
axes=h(1);
plotoption = 1;
[h_whamp,f,Etot]=art3.plot_f(h(1),...
                        w_n(cmp),...
                        w_m(cmp),...
                        w_t(cmp),...
                        w_vd(cmp),...
                        w_d(cmp),...
                        w_a1(cmp),...
                        w_a2(cmp),...
                        w_pa,...
                        plotoption,...
                        titleoption);
set(h(1),'ylim',[5e-5 5e1],'xlim',[7e1 3e4],'xtick',[1e0 1e1 1e2 1e3 1e4 1e5])
%h_info = text(h(1),1e3,1e-2,info,'fontsize',12);
h_info = text(80,0.1,1e-1,info,'fontsize',12,'verticalalignment','top','Parent',h(1));
%set(h_info,'position',[85 0.08 17],'verticalalignment','top')
h2_info = text(4200,15,{'* C3','o C4'},'fontsize',12,'Parent',h(1));
%set(h2_info,'position',[4200 15],'verticalalignment','top')
title(h(1),' ')

% Plot electric field.
%set(0,'DefaultAxesColorOrder',[0 0.8 0;0 0 0.8]);
irf_plot(h(2),{irf_tlim(Epar3,t+[-1 1]),irf_tlim(Epar4,t+[-1 1])},'comp')
irf_pl_mark(h(2),peaPSD3.t(ind_t{1})+[-0.125/2 0.125/2])
irf_legend(h(2),{'C3','C4'},[0.98 0.95])
irf_zoom(h(2),'x',t+[-1 1])



grid(h(1),'off'); grid(h(2),'off');
hold(h(1),'off'); hold(h(2),'off');
colormap('default');
