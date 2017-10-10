%% THOR scatter plot (fp/fc,beta) - SW,MSH 1) get data
%Using Jan's database
% get MSH data
cd('/Users/Cecilia/Google Drive/THOR/Figures/matlab/')
load('../../Plasma parameters/databases/msht_stats_cl1_nov06jun08_cis_ok.mat')
B=pars.bmag;
beta_msh=pars.beta_perp;
n = pars.dens;
fp_kHz = 9*sqrt(n);
fc_kHz = 0.028*B;
fp2fc_msh = fp_kHz./fc_kHz;

% get solar wind data
tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 12 31 23 59 0])];
ff= irf_get_data_omni(tint,'B,n,beta');
load ff
fp_kHz = 9*sqrt(ff(:,3));
fc_kHz = 0.028*(ff(:,2));
fp2fc_sw = fp_kHz./fc_kHz;
beta_sw = ff(:,4);

%% Read MMS data
cd('/Users/andris/Dropbox (IRFU)/mms/SDP/db_1sec');
load_1sec_db

%% Prepare MMS data
db = dball;
n=db.data.ne;
T=db.data.te;
B=db.data.b;
B=sqrt(sum(B.^2,2));
ivok = sum(dball.data.vi.^2,2)<100^2; % v < 100km/s
inok = dball.data.ni<2;               % n < 2 cc
iok = ivok & inok;

beta = 0.4 * n(iok) .* T(iok) ./ B(iok).^2 ;
fp2fc = 9*n./(.028*B(iok));

plot(db.data.ni,db.data.ti,'.')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot(beta,fp2fc,'.');
set(gca,'YScale','log');set(gca,'XScale','log');
irf_figmenu

%% Read in MMS data

%% THOR scatter plot (fp/fc,beta) - SW,MSH 2) plot
%plot
figure;
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);

color_msh  = [1 0.5 0.5];
color_sw   = [0.5 0.5 0.9];
color_msph = [0.7 0.7 0.7];

scatter(beta_msh,fp2fc_msh,2,color_msh);
hold on;
scatter(beta_sw,fp2fc_sw,2,color_sw);
set(gca,'yscale','log')
set(gca,'xscale','log')

ylabel('f_{pe}/f_{ce}')
xlabel('\beta')

txt = struct('beta',[],'fp2fc',[],'region',[]);
txt(    1) = struct('beta',0.1  ,'fp2fc',10,  'region','Corona');
txt(end+1) = struct('beta',0.008,'fp2fc',1,   'region',{{'Lower','corona'}});
txt(end+1) = struct('beta',1    ,'fp2fc',220, 'region',{{'Outer','corona'}});
txt(end+1) = struct('beta',0.06 ,'fp2fc',0.6, 'region','Tokamaks');
txt(end+1) = struct('beta',0.3  ,'fp2fc',45,  'region','MRX(lab)');
txt(end+1) = struct('beta',0.2  ,'fp2fc',6,   'region','RFX(lab)');
txt(end+1) = struct('beta',3  ,  'fp2fc',26,  'region','laser plasma');
txt(end+1) = struct('beta',3  ,  'fp2fc',26,  'region','laser plasma');
txt(end+1) = struct('beta',10  , 'fp2fc',100, 'region',{{'Supernova','remnants'}});
txt(end+1) = struct('beta',20  , 'fp2fc',2e3, 'region',{{'Interstellar','medium'}});

for itxt = 1:numel(txt)
	lbl = txt(itxt);
	ht= text(lbl.beta,lbl.fp2fc,lbl.region);
	set(ht,'horizontalalignment','center');
end

text(1e-2,1.5e3,'THOR solar wind expected','color',color_sw,'fontsize',16);
text(1e-2,1e3,'THOR magnetosheath expected','color',color_msh,'fontsize',16);
text(1e-2,1e3/1.5,'THOR magnetosphere expected','color',color_msph,'fontsize',16);

ylim([3e-1 3e3])
xlim([3e-3 1e2])

grid on;
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
