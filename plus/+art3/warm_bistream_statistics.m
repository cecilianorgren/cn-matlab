%% Define filePath and fileName and fileIndex
fileDir = '/Users/Cecilia/Research/EH2/comp2stats/';
fileNames = {'warm_bistream.txt',...
            'modified_buneman.txt',...
            'modified_buneman2.txt',...
            'modified_buneman3.txt',...
            'electron_beam1.txt',...
            'test.txt'};
fileIndex = 5;        
filePath = [fileDir fileNames{fileIndex}];  

%% Model, generate some random values
N = 10; % number of distributions
instability = fileNames{fileIndex}(1:end-4);
newFile = 0;

switch instability
    case 'warm_bistream'
        % densities
        mR = 0.5; stdR = 0.1; % mean 0.5, standard deviation 0.1
        R = mR + stdR.*randn(N,1); 
        R(R>0.9) = mR;
        R(R<0.1) = mR;
        ni = 1*ones(N,1);
        ne1 = ni.*(1-R);
        ne2 = ni.*R;

        % temperatures
        mTi = 1000; stdTi = 140; Ti  = mTi  + stdTi*randn(N,1);
        mTe1 = 500; stdTe1 = 70; Te1 = mTe1 + stdTe1*randn(N,1);
        mTe2 = 500; stdTe2 = 70; Te2 = mTe2 + stdTe2*randn(N,1);

        % drift velocities
        mVdi = 0; stdVdi = 100;
        Vdi = mVdi + stdVdi*randn(N,1);
        mVde1 = cn_eV2v(400,'eV'); stdVde1 = cn_eV2v(70,'eV');
        Vde1 = mVde1 + stdVde1*randn(N,1); % km/s
        mVde2 = -cn_eV2v(400,'eV'); stdVde2 = cn_eV2v(70,'eV');
        Vde2 = mVde2 + stdVde2*randn(N,1); % km/s
    case 'modified_buneman'
        % densities
        mR = 0.5; stdR = 0.1; % mean 0.5, standard deviation 0.1
        R = mR + stdR.*randn(N,1); 
        R(R>0.9) = mR;
        R(R<0.1) = mR;
        ni = 1*ones(N,1);
        ne1 = ni.*(1-R);
        ne2 = ni.*R;

        % temperatures
        mTe1 = 1000; stdTe1 = 100;
        Te1 = mTe1 + stdTe1*randn(N,1);
        mTe2 = 50; stdTe2 = 10;
        Te2 = mTe2 + stdTe2*randn(N,1);
        mTi = 1000; stdTi = 140;
        Ti = mTi + stdTi*randn(N,1);

        % drift velocities
        mVdi = 0; stdVdi = 100;
        Vdi = mVdi + stdVdi*randn(N,1);
        mVde1 = 0; stdVde1 = 100;
        Vde1 = mVde1 + stdVde1*randn(N,1); % km/s
        mVde2 = cn_eV2v(mTe1,'eV'); stdVde2 = mVde2;
        Vde2 = mVde2 + stdVde2*randn(N,1); % km/s
    case 'modified_buneman2'
        % densities
        mR = 0.3; stdR = 0.1; % mean 0.5, standard deviation 0.1
        R = mR + stdR.*randn(N,1); 
        R(R>0.9) = mR;
        R(R<0) = 0.1;
        ni = 1*ones(N,1);
        ne1 = ni.*(1-R);
        ne2 = ni.*R;

        % temperatures
        mTe1 = 1000; stdTe1 = 100;
        Te1 = mTe1 + stdTe1*randn(N,1);
        mTe2 = 50; stdTe2 = 20;
        Te2 = mTe2 + stdTe2*randn(N,1);
        mTi = 1000; stdTi = 140;
        Ti = mTi + stdTi*randn(N,1);

        % drift velocities
        mVdi = 0; stdVdi = 100;
        Vdi = mVdi + stdVdi*randn(N,1);
        mVde1 = 0; stdVde1 = 100;
        Vde1 = mVde1 + stdVde1*randn(N,1); % km/s
        mVde2 = cn_eV2v(mTe1,'eV')*0.5; stdVde2 = mVde2*0.5;
        Vde2 = mVde2 + stdVde2*randn(N,1); % km/s 
    case 'modified_buneman3'
        % densities
        mR = 0.3; stdR = 0.1; % mean 0.5, standard deviation 0.1
        R = mR + stdR.*randn(N,1); 
        R(R>0.9) = mR;
        R(R<0) = 0.02;
        ni = 1*ones(N,1);
        ne1 = ni.*(1-R);
        ne2 = ni.*R;

        % temperatures
        mTe1 = 1000; stdTe1 = 100;
        Te1 = mTe1 + stdTe1*randn(N,1);
        mTe2 = 50; stdTe2 = 30;
        Te2 = mTe2 + stdTe2*randn(N,1);
        mTi = 1000; stdTi = 150;
        Ti = mTi + stdTi*randn(N,1);

        % drift velocities
        mVdi = 0; stdVdi = 100;
        Vdi = mVdi + stdVdi*randn(N,1);
        mVde1 = 0; stdVde1 = 100;
        Vde1 = mVde1 + stdVde1*randn(N,1); % km/s
        mVde2 = cn_eV2v(mTe1,'eV')*0.6; stdVde2 = mVde2*0.6;
        Vde2 = mVde2 + stdVde2*randn(N,1); % km/s    
    case 'electron_beam1'
        % densities
        mR = 0.1; stdR = 0.05; % mean 0.5, standard deviation 0.1
        R = mR + stdR.*randn(N,1); 
        R(R>0.9) = mR;
        R(R<0) = 0.005;
        ni = 1*ones(N,1);
        ne1 = ni.*(1-R);
        ne2 = ni.*R;

        % temperatures
        mTe1 = 1000; stdTe1 = 100;
        Te1 = mTe1 + stdTe1*randn(N,1);
        mTe2 = 30; stdTe2 = 20;
        Te2 = mTe2 + stdTe2*randn(N,1);
        mTi = 1000; stdTi = 150;
        Ti = mTi + stdTi*randn(N,1);

        % drift velocities
        mVdi = 0; stdVdi = 100;
        Vdi = mVdi + stdVdi*randn(N,1);
        mVde1 = 0; stdVde1 = 100;
        Vde1 = mVde1 + stdVde1*randn(N,1); % km/s
        mVde2 = cn_eV2v(mTe1,'eV'); stdVde2 = mVde2*0.5;
        Vde2 = mVde2 + stdVde2*randn(N,1); % km/s        
    case 'test'
        % densities
        mR = 0.3; stdR = 0.1; % mean 0.5, standard deviation 0.1
        R = mR + stdR.*randn(N,1); 
        R(R>0.9) = mR;
        R(R<0) = 0.02;
        ni = 1*ones(N,1);
        ne1 = ni.*(1-R);
        ne2 = ni.*R;

        % temperatures
        mTe1 = 1000; stdTe1 = 100;
        Te1 = mTe1 + stdTe1*randn(N,1);
        mTe2 = 50; stdTe2 = 30;
        Te2 = mTe2 + stdTe2*randn(N,1);
        mTi = 1000; stdTi = 150;
        Ti = mTi + stdTi*randn(N,1);

        % drift velocities
        mVdi = 0; stdVdi = 100;
        Vdi = mVdi + stdVdi*randn(N,1);
        mVde1 = 0; stdVde1 = 100;
        Vde1 = mVde1 + stdVde1*randn(N,1); % km/s
        mVde2 = cn_eV2v(mTe1,'eV')*0.6; stdVde2 = mVde2*0.6;
        Vde2 = mVde2 + stdVde2*randn(N,1); % km/s        
end

if newFile
    %%
    fid = fopen(filePath,'w');
    fileInputHeader = ['mR    stdR    mTi    stdTi    mTe1    stdTe1    mTe2    stdTe2     mVdi       stdVdi       mVde1       stdVde1      mVde2    stdVde2'];
    formatInputHeader = '\n%s %s %s %s %s %s %s %s %s %s %s %s %s %s';
    fprintf(fid,formatInputHeader,fileInputHeader);
    fclose(fid)
    fid = fopen(filePath,'a');
    fprintf(fid,'\n') 
    formatInputData = '%3.2f %5.2f  %5.0f  %5.0f    %5.0f   %5.0f     %5.0f   %5.0f     %8.0f   %8.0f     %8.0f    %8.0f      %8.0f    %8.0f\n\n';
    fprintf(fid,formatInputData,[mR    stdR    mTi    stdTi    mTe1    stdTe1    mTe2    stdTe2    mVdi    stdVdi    mVde1    stdVde1    mVde2    stdVde2]);

    fileHeader = 'ni    ne1   ne2   Ti    Te1  Te2  vdi    vde1    vde2    k     wr        wi       vph        vT1      vT2    vph_choose1  vph_choose2 ';
    formatHeader = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n';
    fprintf(fid,formatHeader,fileHeader);
    fprintf(fid,'\n')
    fclose(fid);
    type(filePath)
end
%% Display inputs
if 1 % show input
    %%
    nPl = 6;
    subplot(2,3,1); hist(R,20); set(gca,'xlim',[0 1]); xlabel('R');
    subplot(2,3,2); hist(Te1,20); set(gca,'xlim',[0 2*mTe1]); xlabel('Te1 [eV]');
    subplot(2,3,3); hist(Te2,20); set(gca,'xlim',[0 2*mTe2]); xlabel('Te2 [eV]');
    subplot(2,3,4); hist(Ti,20); set(gca,'xlim',[0 2000]); xlabel('Ti [eV]');
    subplot(2,3,5); [nb,xb] = hist(Vde1/1000,20); bh = bar(xb,nb); set(gca,'xlim',[-40 40]); xlabel('v_{de1}, v_{de2} [10^3 km/s]'); set(bh,'FaceColor',[0 1 0.5],'edgecolor',[0 1 0.5]); hold(gca,'on');
                    [nb,xb] = hist(Vde2/1000,20); bh = bar(xb,nb); set(gca,'xlim',[-40 40]); xlabel('v_{de1}, v_{de2} [10^3 km/s]'); set(bh,'FaceColor',[1 0.5 0,],'edgecolor',[1 0.5 0,]); hold(gca,'off');
    subplot(2,3,6); hist(Vdi,20); set(gca,'xlim',[-500 500]); xlabel('v_{di} [km/s]');
end

%% Initialize data to save
wr = zeros(N,1);
wi = zeros(N,1);
k = zeros(N,1);
vph = zeros(N,1);
vph_choose = nan(N,2);

0%% Loop over all distributions.
printStartNumber = 0; % (actually this + 1)
printName = 'electron_beam_1';
for ppp = 1:N


%% This is the script to produce figure and ask for ginput
no=1;
B=20;
% Ions
omega_pi = irf_plasma_calc(B,ni(ppp),no,Te1(ppp),Ti(ppp),'Fpp')*2*pi; % rad/s
vthi = irf_plasma_calc(B,ni(ppp),no,Te1(ppp),Ti(ppp),'Vtp'); % m/s 
vdi = Vdi(ppp)*1000; % m/s

% Background electrons
omega_pe1 = irf_plasma_calc(B,ne1(ppp),no,Te1(ppp),Ti(ppp),'Fpe')*2*pi; % rad/s
vthe1 = irf_plasma_calc(B,ne1(ppp),no,Te1(ppp),Ti(ppp),'Vte'); % m/s 
vde1 = Vde1(ppp)*1000; % m/s

% Beam electrons
omega_pe2 = irf_plasma_calc(B,ne2(ppp),no,Te2(ppp),Ti(ppp),'Fpe')*2*pi; % rad/s
vthe2 = irf_plasma_calc(B,ne2(ppp),no,Te2(ppp),Ti(ppp),'Vte'); % m/s 
vde2 = Vde2(ppp)*1000; % m/s

% Other
lamD = irf_plasma_calc(B,ni(ppp),no,Te1(ppp),Ti(ppp),'Ld'); % m
omega_bune = omega_pe1^(1/3)*omega_pi^(2/3)*16^(-1/3);

normal_frequency = omega_pe1;
normal_length = lamD;
normal_frequency_string = '\omega_{pe}';
normal_length_string = '\lambda_D';

fv = [omega_pi,omega_pe1,omega_pe2,vthi,vthe1,vthe2,vdi,vde1,vde2,lamD];

% set up figure for drawing results
close
figure(3)    
for nPlot = 1:4   
    h(nPlot) = subplot(1,4,nPlot); hold(h(nPlot),'on');        
    xlabel(h(nPlot),['k' normal_length_string],'Fontsize',14);
    ylabel(h(nPlot),['\omega/'  normal_frequency_string],'Fontsize',14);    
    set(h(nPlot),'XGrid','on','YGrid','on','Fontsize',14);
end
ylabel(h(3),'v_{ph} [km/s]')
title(h(1),'\Re (\omega)'); 
title(h(2),'\Im (\omega)')
title(h(3),'v_{ph}')
xlabel(h(4),'v [10^3 km/s]')
title(h(4),'Input distribution')
v_psd = -60000:100:60000;
e1_psd = cn.maxwellian(v_psd,Te1(ppp),ne1(ppp),Vde1(ppp),'e');
e2_psd = cn.maxwellian(v_psd,Te2(ppp),ne2(ppp),Vde2(ppp),'e');
i_psd = cn.maxwellian(v_psd,Ti(ppp),ni(ppp),Vdi(ppp),'i');
e_psd = e1_psd + e2_psd;
psd_lim = [0 max(e_psd)*1.2];
plot(h(4),v_psd,e_psd,v_psd,i_psd)
set(h(4),'ylim',psd_lim,'xlim',0.5e5*[-1 1])
axes(h(4))
distrStr = {['R=' num2str(R(ppp),'%.2f')],...
            ['Ti=' num2str(Ti(ppp),'%.0f')],...
            ['Te1=' num2str(Te1(ppp),'%.0f')],...
            ['Te2=' num2str(Te2(ppp),'%.0f')],...
            ['vdi=' num2str(Vdi(ppp),'%.0f')],...
            ['vde1=' num2str(Vde1(ppp),'%.0f')],...
            ['vde2=' num2str(Vde2(ppp),'%.0f')],...
            }; 
text(min(get(gca,'xlim')),max(get(gca,'ylim')),distrStr,'horizontalalignment','left','verticalalignment','top')

set(gcf,'position',[195 389 1455 376]);

%

ks = 1:0.05:1.8;
ks = 0.2:0.05:1.8;
ks = 0.01:0.05:1;
nk = numel(ks);
rights = zeros(1,numel(ks));
nTries = 10;

all_wr = zeros(nk,nTries);
all_wi = zeros(nk,nTries);
all_vph = zeros(nk,nTries);
all_k = zeros(nk,nTries);

for na=1:numel(ks);
    a = ks(na);
    
    % For each k*lD, make 20 random guesses within specified ranges
    for iii=1:nTries
        if iii<fix(nTries/2+1) %             
            x_real=0.1*randn*normal_frequency;             
            x_imag=0.05*randn*normal_frequency; % x_imag=2*rand*omega_pe;  
        else % around normal_frequency              
            x_real=(1.5*randn*normal_frequency);            
            x_imag=0.2*randn*normal_frequency; % x_imag=2*rand*omega_pe;
        end
        %af = @(temp) fff_two_stream(temp,a);
        %af = @(temp) fff_two_stream_test(temp,a);
        af = @(temp) myfun_20070831_EH(temp,a,fv);
        [x,FVAL,EXITFLAG] = fsolve(af,[x_real x_imag],optimoptions('fsolve','MaxIter',1000,'TolFun',1e-8,'TolX',1e-13,'MaxFunEvals',5000,'display','off'));
        [x_real, x_imag];
        all_wr(na,iii) = x(1); 
        all_wi(na,iii) = x(2); 
        all_k(na,iii) = ks(na);
        all_vph(na,iii) = x(1)/real(a)*normal_length/1000; % km/s
        
        if(EXITFLAG==1) % Found something, now plot            
            plot(h(1),real(a),x(1)/normal_frequency,'k.','Markersize',3);             
            plot(h(1),real(a),x(2)/normal_frequency,'k.','Markersize',3); hold on;
            plot(h(3),real(a),x(1)/real(a)*normal_length/1000,'k.','Markersize',3); hold on; 
            % only plot for reasonably omega_i > -0.6
            if(x(2)/normal_frequency>-0.6)
                if ( abs(x(1)/normal_frequency)<0.8)                    
                    plot(h(1),real(a),x(1)/normal_frequency,'r+','Markersize',6); hold on; 
                    plot(h(2),real(a),x(2)/normal_frequency,'r+','Markersize',6); hold on; 
                    plot(h(3),real(a),x(1)/real(a)*normal_length/1000,'r+','Markersize',6); hold on; 
                elseif 0;( abs(x(1)/normal_frequency)>0.8) % before it was 0.1                    
                    plot(h(1),real(a),x(1)/normal_frequency,'bd','Markersize',6); hold on;                    
                    plot(h(2),real(a),x(2)/normal_frequency,'bd','Markersize',6); hold on;    
                end
                if (x(2)/normal_frequency>0.0001)  % positive growth rate ?                  
                    rights(na)=rights(na)+1;
                    plot(h(1),real(a),x(1)/normal_frequency,'ch','Markersize',6); hold on;                    
                    plot(h(2),real(a),x(2)/normal_frequency,'ch','Markersize',6); hold on; 
                    plot(h(3),real(a),x(1)/real(a)*normal_length/1000,'ch','Markersize',6); hold on; 
                    disp(['x_real=' num2str(x_real,'%.4f') ', x_imag=' num2str(x_imag,'%.4f') ' => x_real=' num2str(x(1),'%.4f') ', x_imag=' num2str(x(2),'%.4f')])
                end
            end                   
            drawnow               
        end
    end
    rights;
    % Some analytical solution?
    %plot(h(1),real(a),sqrt(omega_pe*omega_pe + 3*(a/lamD)*(a/lamD)*vthe*vthe/2)/omega_pe  ,'r.','Markersize',3); hold on;
    %plot(h(1),real(a),(vSound*real(a)/2/lamD)/omega_pe,'g.','Markersize',3); hold on;   
end

if 0
    set(h(1),'ylim',[0 0.6],'xlim',[0 max(ks)]);
    set(h(2),'ylim',[-0.01 0.04],'xlim',[0 max(ks)])
end

if 1 % fix up figure and print
    %%
    set(h(1),'ylim',[-1 2])
    set(h(2),'ylim',[-0.2 max(get(h(2),'ylim'))])
    set(h(3),'ylim',[-0.2e5 0.2e5])
    
end

if 1
    %%
    nvs = input('How many vph would you like to add?');
    hold(h(4),'on')
    %set(h(3),'ylim',1e3*[-30 30])
    for kk = 1:nvs
        vph_choose_tmp = ginput(1);
        vph_choose(ppp,kk) = vph_choose_tmp(2);       
        plot(h(4),vph_choose(ppp,kk)*[1 1],[0 max(e_psd)],'r--')
    end
    %set(h(4),'ylim',[0 max(e_psd)])
    hold(h(4),'off')
end
if 0 % add phase velocities to distribution plot     
    nvs = input('How many vph would you like to add?')
    hold(h(4),'on')
    set(h(3),'ylim',1e3*[-30 30])
    for kk = 1:nvs
        vph_toplot = ginput(1);
        plot(h(4),vph_toplot(2)*[1 1],[0 max(e_psd)],'r--')
    end
    set(h(4),'ylim',[0 max(e_psd)])
    hold(h(4),'off')
else % automatic
    [max_value, max_index] = max(all_wi(:));
    wr(ppp) = all_wr(max_index);
    wi(ppp) = all_wi(max_index);
    vph(ppp) = all_vph(max_index);
    k(ppp) = all_k(max_index);
end
% print image
cn.print([printName '_' num2str(ppp+printStartNumber)])
end








%% Extract the trapping velocity.
% The closest difference between the phase velocity and the two beams
vT1 = abs(vph-Vde1*1e3); 
vT2 = abs(vph-Vde2*1e3); 
%% Write to file
fid = fopen(filePath,'a');
formatData = '%3.2f %5.2f %5.2f %5.0f %4.0f %4.0f %5.0f %7.0f %8.0f %5.2f %8.0f %8.0f %9.0f %9.0f %9.0f %9.0f %9.0f \n';
%sprintf(formatData,[ni ne1 ne2 Ti Te1 Te2 Vdi Vde1 Vde2 k wr wi vph vT1 vT2 vph_choose]')
fprintf(fid,formatData,[ni ne1 ne2 Ti Te1 Te2 Vdi Vde1 Vde2 k wr wi vph vT1 vT2 vph_choose]');
fclose(fid);    
type(filePath)

%% Read the collected data
fid = fopen(filePath,'r');
formatReadData = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n';
C = textscan(fid,formatReadData,'HeaderLines',1);

if 0 % Check trapping velocity
    %%
    R = C{3};
    v_eh = C{13};
    v_de1 = C{8};
    v_de2 = C{9};
    v_T1 = C{14};
    v_T2 = C{15};
    %%
    plot3(v_de1(v_eh<0),v_de2(v_eh<0),v_eh(v_eh<0),'o',v_de1(v_eh>0),v_de2(v_eh>0),v_eh(v_eh>0),'s')
    xlabel('vde1'); ylabel('vde2')
    set(gca,'xlim',[0 3]*1e4,'ylim',[-3 0]*1e4)
    %%
    plot(mean([v_de1 v_de2],2),v_eh,'o')
    xlabel('mean([vde1 vde2])'); ylabel('veh')
    %%
    plot(R,vde1,'ro',R,vde2,'go',R,v_eh,'bo',R,v_T,'mo')
    %legend('vde1','vde2','veh','vT')
end
    
% Read Daniels data
load('/Users/Cecilia/Research/EH2/Daniel/wavespeeds.mat')

% Compare theory to observations
plot(abs(C{13})*1e-2,C{14}*1e-2,'o',wavespeeds.veholes,wavespeeds.vtreholes,'gs'); 
xlabel('v_{ph}'); ylabel('v_{T}')


%% Compare warm bistream to modified buneman/electron beam
filePath = '/Users/Cecilia/Research/EH2/comp2stats/warm_bistream.txt';
fid = fopen(filePath,'r');
formatReadData = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n';
wbC = textscan(fid,formatReadData,'HeaderLines',1);

filePath = '/Users/Cecilia/Research/EH2/comp2stats/modified_buneman3.txt';
fid = fopen(filePath,'r');
formatReadData = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n';
mbC = textscan(fid,formatReadData,'HeaderLines',1);

mb_veh = mbC{13};
mb_vT1 = abs(mbC{16}-mbC{8}*1e3);
mb_vT2 = abs(mbC{16}-mbC{9}*1e3);
wb_veh = wbC{13};
wb_vT1 = abs(wbC{13}-wbC{8}*1e3);
wb_vT2 = abs(wbC{13}-wbC{9}*1e3);

loglog(abs(mb_veh),abs(mb_vT1),'ro',abs(mb_veh),abs(mb_vT2),'rs',...
     abs(wb_veh),abs(wb_vT1),'go',abs(wb_veh),abs(wb_vT2),'gs',...
     wavespeeds.veholes*1e3,wavespeeds.vtreholes*1e3,'bx',...
     [1e1 2e7],[1e1 2e7])
 xlabel('v_{ph}')
 ylabel('v_{T}')