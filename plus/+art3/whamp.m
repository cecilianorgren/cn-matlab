%% 1. Run whamp with quite distinct beam and try to find beam mode
R=0.5;
ni=0.06;
S=0.75; % vd/vtbg;
n   = [ni*(1-R) ni*R ni 0 0 0 0 0 0 0]*1e6;        % density
t   = [1.6 0.06 2.2 0 0 0 0 0 0 0];            % temperature [keV]
d   = [1 1 1 1 1 1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
a   = [1 1 1 1 1 1 1 1 1 1];           % T_perp/T_par
b   = [0.0 0.0 0 0 0 0 0 0 0 0];        % loss cone parameter
mass= [0 0 1 0 0 0 0 0 0 0];        % mass
vd  = [0 S*sqrt(t(1)/t(2)) 0 0 0 0 0 0 0 0];         % V_drift/V_term
fce = 0.7; % kHz
pzl = 1;  
pitchangle = [0 45 90 135 180];
plotoption = 0; % f vs E
%excl = []; 
%n(excl)=0;
incl=1:3;
[h,f,v] = whamp.plot_f(n(incl),mass(incl),t(incl),vd(incl),d(incl),a(incl),b(incl),pitchangle,plotoption);
close
txtModel = ['eh2-R' num2str(R*100,'%03.0f') '-S' num2str(S*1000,'%04.0f') '.txt'];
%% Write model to text file
myformat= [repmat('%.6f ',1,10) '\n'];

% open the file with permission to append
%txtModel = '20130831-beam.txt';
%fid = fopen(['/Users/Cecilia/Research/EH2/whamp/models/' txtModel],'w');
fid = fopen(['/Users/Cecilia/Whamp/Models/' txtModel],'w');

% write values at end of file
fprintf(fid, myformat,n);
fprintf(fid, myformat,t);
fprintf(fid, myformat,d);
fprintf(fid, myformat,a);
fprintf(fid, myformat,b);
fprintf(fid, myformat,mass);
fprintf(fid, myformat,vd);
fprintf(fid, '%.3f \n',fce);
fprintf(fid, '%1.0f \n',pzl);
fprintf(fid, myformat,vd);

% close the file 
fclose(fid);
disp(['Model file saved to: /Users/Cecilia/Whamp/Models/' txtModel])

% Create textfile for results if it doesnt already exist.
folder = '';
fname1  = ['/Users/Cecilia/Research/EH2/whamp/results/' folder txtModel];
fopen(fname1,'a');
disp(['Save WHAMP results to: /Users/Cecilia/Whamp/Models/' txtModel])

%% Now run whamp and save data
% where to look for solutions
% fr = 1-10 fpe
fpe = irf_plasma_calc(25,n(1)*1e-6,0,1.6,2.2,'Fpe');
xf = 1;
fn = fpe/(fce*1e3);
fr = xf*fn;
p = 0; % parallel propagation
vth1 = irf_plasma_calc(25,0.06,0,t(1)*1e3,2.2,'Vte'); % m/s
lamDe = irf_plasma_calc(25,0.06,0,t(1)*1e3,2.2,'Ld'); % m/s
xl = 0.4; % klamDe = xl
zn = vth1/2/pi/lamDe/(fce*1e3);
z = xl*zn;

%%
%cd      = '/Users/Cecilia/irfu-matlab/+whamp/';
%txtFile = 'eh2-vd0500_2.txt';
%txtFile = txtModel;
%folder = 'Sx_R20/';
%fname1  = ['/Users/Cecilia/Research/EH2/whamp/results/' folder txtFile];
file1   = load(fname1);

fig=figure(2);
set(fig,'defaultAxesFontUnits','pixels');
set(fig,'defaultAxesFontSize',14);
set(fig,'defaultTextFontSize',14);

[p,z,fr,fim]=m2xyz(file1); 
vp=(fr/fn)./(z'/zn)/sqrt(2);
vt = cn_eV2v(t(1)*1e3,'eV');
vreal = vp*xvd*vt;
v_fimax = mean(vreal(fim==max(fim)));

% Buneman phase velocity
units=irf_units;
vph_Bun = (units.me/16/units.mp)^(1/3)*S; % normalized to vtbg


% make figure
nPlots=3;
for kk = 1:nPlots; h(kk)=subplot(1,nPlots,kk); end
isub=1;
if 1 % psd
    hca=h(isub); isub=isub+1;
    plot(hca,v_fimax*[1 1],[min(min(f)) max(max(f))],'--b',...
             sqrt(2)*v_fimax*[1 1],[min(min(f)) max(max(f))],'--g',...
             mean(vreal)*[1 1],[min(min(f)) max(max(f))],'--b',...
             sqrt(2)*mean(vreal)*[1 1],[min(min(f)) max(max(f))],'--g'); hold(hca,'on')
    plot(hca,v',f')    
    xlabel(hca,'v [km/s]')
    ylabel(hca,'psd')
    legend(hca,'<v_{ph}>','<v_{ph}>*sqrt(2)')
    legend(hca,'<v_{max\gamma}>','<v_{max\gamma}>*sqrt(2)')
    set(hca,'ylim',[1e-8 1e4],'xlim',[1e3 2e5],'xscale','log','yscale','log')
end
if 1 % dispersion relation
    hca=h(isub); isub=isub+1;
    plot(hca,z/zn,fr/fn,z/zn,fim/fn) % vph: ,z/zn,(fr/fn)./(z'/zn)/sqrt(2)    
    title(hca,[txtModel ' -> ' txtFile])
    ylabel(hca,'\omega/\omega_{pe}','fontsize',14);
    xlabel(hca,'k\lambda_{De}','fontsize',14);
    legend(hca,'\omega_r','\omega_i','v_{ph}/v_{te,bg}','location','northwest')
    %set(hca,'ylim',[0 1.2*max([fr/fn;fim/fn])],'xlim',[0 1.1*z([end])/zn])
end
if 1 % dispersion relation
    hca=h(isub); isub=isub+1;
    plot(hca,z/zn,(fr/fn)./(z'/zn)/sqrt(2))    
    %title(hca,[txtModel ' -> ' txtFile])
    ylabel(hca,'v_{ph}/v_{te,bg}','fontsize',14);
    xlabel(hca,'k\lambda_{De}','fontsize',14);    
    %set(hca,'ylim',[0 1.2*max([fr/fn;fim/fn])],'xlim',[0 1.1*z([end])/zn])
end
%xlabel(gca,'log_{10}(p)','fontsize',14);
%ylabel(gca,'log_{10}(z)','fontsize',14);
%zlabel(gca,'f/f_{cp}','fontsize',14);
%cax=colorbar;
%ylabel(cax,'f_i/f_{pp}')
%set(gcf,'paperpositionmode','auto');
%% Compare fsolve and WHAMP.
% Run solver_Buneman first, then this script
%ax_handles = findobj(gcf,'type','axes')
plot(h(1),z/zn,fr/fn,'k')
plot(h(2),z/zn,fim/fn,'k')
if 0 % check phase velocity
    %%
    vphs_inp=ginput(2);
    v_ph_1 = vphs_inp(1,2)/vphs_inp(1,1)/sqrt(2);
    v_ph_2 = vphs_inp(2,2)/vphs_inp(2,1)/sqrt(2);
    disp(['vph1 = ' num2str(v_ph_1) '*vtebg,  vph2 = ' num2str(v_ph_2) '*vtebg'])
end

%% Plot all surfaces which are within a certain folder
folder = 'Sx_R20/';
folderPath  = '/Users/Cecilia/Research/EH2/whamp/results/';
fileList = dir([folderPath folder '*.txt']);
nFiles = numel(fileList);

Ss = zeros(1,nFiles);
colors = jet(nFiles);get(0,'DefaultAxesColorOrder');
units=irf_units;

nPlots=3;
for kk = 1:nPlots; h(kk)=subplot(1,nPlots,kk); hold(h(kk),'on'); end

ifile=0;

for kk=4:nFiles;
    isub=1;    
    S = str2double(fileList(kk).name(7:10))/1000;
    
    file1   = load([folderPath folder fileList(kk).name]); ifile=ifile+1;
    [p,z,fr,fim]=m2xyz(file1); 
    vp=(fr/fn)./(z'/zn)/sqrt(2);
    vt = cn_eV2v(t(1)*1e3,'eV');
    vreal = vp*xvd*vt;
    v_mean = mean(vreal);
    v_fimax = mean(vreal(fim==max(fim)));    
    v_bun = (units.me/16/units.mp)^(1/3)*S; % normalized to vtbg
    tind = find(abs(z/zn-1)==min(abs(z/zn-1)));
    if 1 % dispersion relation
        hca=h(isub); isub=isub+1;
        plot(hca,z/zn,fr/fn,'color',colors(kk,:)*0.9)         
        ylabel(hca,'\omega_r/\omega_{pe}','fontsize',14);
        xlabel(hca,'k\lambda_{De}','fontsize',14);
        str1 = ['\leftarrow S=' num2str(S,'%.1f')];
        str2 = ['S=' num2str(S,'%.1f') ' \rightarrow'];
        text(z(tind)/zn,fr(tind)/fn,str1,'Parent',hca,'horizontalalignment','left')
        %set(hca,'ylim',[0 1.2*max([fr/fn;fim/fn])],'xlim',[0 1.1*z([end])/zn])
    end
    if 1 % dispersion relation
        hca=h(isub); isub=isub+1;
        plot(hca,z/zn,fim/fn,'color',colors(kk,:)*0.9)         
        ylabel(hca,'\omega_i/\omega_{pe}','fontsize',14);
        xlabel(hca,'k\lambda_{De}','fontsize',14);        
        %title(hca,'\omega_i')
        %set(hca,'ylim',[0 1.2*max([fr/fn;fim/fn])],'xlim',[0 1.1*z([end])/zn])
    end    
    if 1 % dispersion relation
        hca=h(isub); isub=isub+1;
        plot(hca,z/zn,(fr/fn)./(z'/zn)/sqrt(2),'color',colors(kk,:))
        plot(hca,z/zn,v_bun*z./z,'color',colors(kk,:),'linestyle','--')
        plot(hca,z/zn,(fr/fn)./(z'/zn)/sqrt(2),'color',colors(kk,:))
        %title(hca,[txtModel ' -> ' txtFile])
        ylabel(hca,'v_{ph}/v_{te,bg}','fontsize',14);
        xlabel(hca,'k\lambda_{De}','fontsize',14);  
        %set(hca,'yscale','log')
        %set(hca,'ylim',[0 1.2*max([fr/fn;fim/fn])],'xlim',[0 1.1*z([end])/zn])
    end
    
    %plot(h(1),z/zn,fr/fn,'color',colors(kk,:))
    %plot(h(2),z/zn,fim/fn,'color',colors(kk,:))            
end

for kk=1:nPlots
    set(h(kk),'xlim',[0,1.5])
end

plot(h(1),get(h(1),'xlim'),fce*1000/fpe*[1 1],'k')
text(0.1,fce*1000/fpe,'\omega_{ce}','Parent',h(1),'verticalalignment','bottom')

set(gcf,'position',[50 450 1070 500])



