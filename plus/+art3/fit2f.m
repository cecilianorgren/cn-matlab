% Distribution with fit + electric field with sweep marked yellow
cd /Users/Cecilia/Data/Cluster/20070831

% Common parameters for model
w_m  = [0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
w_d  = [1 1];            % loss cone parameter, 1.0 = no loss cone)
w_a2 = [0 0];           % ??? (use 0)
w_pa = [15 45 75 105 135 165];    % pitch angeles
w_pa = [0 30 60 90 120 150 180];    % pitch angeles
w_pa = [15 45 75 105 135 165];    % pitch angeles
plotoption=[1];
titleoption=0;


% Set up figure.
figure(37)
set(gcf,'position',[64         289        1463         613]);
ndists=2;
for ii = 1:ndists; h(ii) = subplot(1,ndists,ii); end
       
% Load distributions.
if ~exist('peaPSD3','var'); load peaPSD; end

% Loop between distributions.
cmp=[1 2];
for ii = 1:ndists;
hca = h(ii);
% Select event and specify rest of model parameters
event = ii;
switch event
    case 4; 
        time = [2007 08 31 10 17 38.158]; % center of beam 
        % Electron distribution model       
        w_n=  [0.08 0.0031]*1e6; % density m^3
        w_t=  [1.7 0.075];       % temperature [keV]
        w_vd= [0 -3.5];           % V_drift/V_term
        w_a1= [1 1];            % T_perp/T_par        
    case 3; 
        time = [2007 08 31 10 17 39.689]; % center of beam
        % Electron distribution model       
        w_n=  [0.06 0.002]*1e6; % density m^3
        w_t=  [1.6 0.070];       % temperature [keV]
        w_vd= [0 -3.1];           % V_drift/V_term
        w_a1= [1 0.5];            % T_perp/T_par  
    case 2; 
        time = [2007 08 31 10 15 35.0]; % center of beam
        % Electron distribution model       
        w_n=  [0.06 0.015]*1e6; % density m^3
        w_t=  [1.6 0.70];       % temperature [keV]
        w_vd= [0 1.1];           % V_drift/V_term
        w_a1= [1 0.5];            % T_perp/T_par      
    case 5; 
        time = [2007 08 31 10 20 26.0]; % center of beam
        % Electron distribution model       
        w_n=  [0.08 0.02]*1e6; % density m^3
        w_t=  [3.2 0.70];       % temperature [keV]
        w_vd= [0 0.1];           % V_drift/V_term
        w_a1= [1  0.5];            % T_perp/T_par 
    case 1; 
        time = [2007 08 31 10 16 17.0]; % center of beam
        % Electron distribution model       
        w_n=  [0.08 0.02]*1e6; % density m^3
        w_t=  [3.2 0.70];       % temperature [keV]
        w_vd= [0 0.1];           % V_drift/V_term
        w_a1= [1  0.5];            % T_perp/T_par 
    
end
t = toepoch(time);

% Distribution info
switch_info = 2;
switch switch_info; %numel(cmp)
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
    case 3
        info = {['T1 = ' num2str(w_t(cmp(1))*1e3) ' eV'],...                                            
                ['T2 = ' num2str(w_t(cmp(2))*1e3) ' eV'],...                                
                ['n2/(n1+n2) = ' num2str(w_n(cmp(2))/(w_n(cmp(1))+w_n(cmp(2))),'%.2f')],...    
                };
end            

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
for jj = 1:n_pa;
    psd3{jj} = squeeze(mean(irf.nanmean(peaPSD3.p(ind_t{1},jj,ind_en{1}),1),2));
    psd4{jj} = squeeze(mean(irf.nanmean(peaPSD4.p(ind_t{2},jj,ind_en{2}),1),2));
    pa_legends{jj} = num2str(pitchangles(jj),'%.0f');
end
bg3 = squeeze(irf.nanmean(irf.nanmean(peaPSD3.p_bg(ind_t{1},:,ind_en{1}),1),2));
bg4 = squeeze(irf.nanmean(irf.nanmean(peaPSD4.p_bg(ind_t{2},:,ind_en{2}),1),2));
%pa_legends{n_pa+1}='Bg';

% Define colors
colors = {[0 0.8 0],[0 0 1],[0.8 0 0],[1 0.5 0],[0 1 1],[1 0 1]};
colorsarr=[]; for jj = 1:numel(colors); colorsarr = [colorsarr;colors{jj}]; end
set(0,'DefaultAxesColorOrder',colorsarr);

% Plot distributions
hold(hca,'on')

for jj = 1:n_pa % the pitch angles
    loglog(hca,peaPSD3.en_cs(ind_en{1}),psd3{jj},'color',colors{jj},'linestyle','*')   
end
for jj = 1:n_pa % the pitch angles    
    loglog(hca,peaPSD4.en_cs(ind_en{2}),psd4{jj},'color',colors{jj},'linestyle','o')
end

loglog(hca,peaPSD3.en_cs(ind_en{1}),bg3,'color',[0 0 0],'linestyle','--')
loglog(hca,peaPSD4.en_cs(ind_en{2}),bg4,'color',[0 0 0],'linestyle','--')
set(hca,'yscale','log','xscale','log')
%legend(hca,pa_legends,'edgecolor','w') 
ylabel(hca,peaPSD3.p_label)
xlabel(hca,'Energy  [eV]')

% Plot fit to distribution.
axes=hca;
plotoption = 1;
[h_whamp,f,Etot]=art3.plot_f(hca,...
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

if 0
    set(hca,'ylim',[0 7],'xlim',[7e1 3e4],...
        'xtick',[1e0 1e1 1e2 1e3 1e4 1e5],...
        'ytick',[0 1 2 3 4 5 6 7 8],...
        'yscale','lin')
    h_info = text(4e3,5,0,info,'fontsize',12,'verticalalignment','top','Parent',hca);
else
    set(hca,'ylim',[5e-5 5e1],'xlim',[7e1 3e4],...
        'xtick',[1e0 1e1 1e2 1e3 1e4 1e5],...
        'ytick',[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],...
        'yscale','log')
    h_info = text(80,1e-4,0,info,'fontsize',12,'verticalalignment','bottom','Parent',hca);
end

%h_info = text(80,0.1,2e-1,info,'fontsize',12,'verticalalignment','top','Parent',hca);

%h2_info = text(2000,15,{'* C3','o C4'},'fontsize',12,'Parent',hca);
h2_info = text(1e4,15,{'* C3','o C4'},'fontsize',12,'Parent',hca);
%title(hca,' ')

box(hca,'on')
% new legends
pa_legends = {'0-30','30-60','60-90','90-120','120-150','150-180'};
end
%legend(hca,pa_legends,'edgecolor','w') 
% new legends
pa_legends = {'0-30','30-60','60-90','90-120','120-150','150-180'};
legend(h(9),pa_legends,'edgecolor','w');



for pp = 1:ndists
    grid(h(pp),'off'); 
    hold(h(pp),'off'); 
    box(h(pp),'on')
end
%colormap('default');
