% Distribution structures are named peaDEFlux?, peaDPFlux? and peaPSD?
% so do load peaDEFlux, peaDPFlux, peaPSD

% Plot 8 distributions for pitch angles 0, 90 and 180. Each distribution
% spans 4 seconds and are separated by 2 seconds, so they overlap with 2
% seconds.

if 0
 caa_download(tint,'C3_CP_PEA_3DXPH_cnts')
 caa_download(tint,'C4_CP_PEA_3DXPH_cnts')
 caa_download(tint,'C3_CP_PEA_3DXPL_cnts')
 caa_download(tint,'C4_CP_PEA_3DXPL_cnts')
end

figure(50)
nrow=4; ncol=5;
nplot=nrow*ncol;
for ii=1:nplot
    h(ii)=subplot(nrow,ncol,ii);
end

%tstart = toepoch([2007 08 31 10 17 39.3]); % C4
tstart = toepoch([2007 08 31 10 17 37.7]); % C3
tstart = toepoch([2007 08 31 10 17 36]); % C3
%tstop = toepoch([2007 08 31 10 17 50]); %tstart+20;peaPSD3.t(end);
tt = 1;
dt = 0.125;

tint3=tstart;
tint4=tstart;

%tint = tstart + [0 tt];
%%
for jj=1:10
if 1
for ii=1:nplot
    c_caa_plot_distribution_function(h(ii),'tint',tint3,'cross-section',peaPSD3_full);
    set(h(ii),'xlim',[4e1 3e4],'ylim',[1e-5 1e2],'xtick',[1e1 1e2 1e3 1e4 1e5])
    tint3 = tint3 + dt;
end
c_eval('cn.print(''C3_full_?'',''many_psd_full'')',jj)
end
if 0
for ii=1:nplot
    c_caa_plot_distribution_function(h(ii),'tint',tint4,'cross-section',peaPSD4);
    set(h(ii),'xlim',[4e1 3e4],'ylim',[1e-5 1e2],'xtick',[1e1 1e2 1e3 1e4 1e5])
    tint4 = tint4 + dt;
end
c_eval('cn.print(''C4_f_?'',''many_psd'')',jj)
end
end
%% Plot psd above and counts below
%%
tint=tstart;
for jj=1:10
if 1 % psd
for ii=1:(nplot/2)
    c_caa_plot_distribution_function(h(ii),'tint',tint,'cross-section',peaPSD4);
    set(h(ii),'xlim',[4e1 3e4],'ylim',[1e-5 1e2],'xtick',[1e1 1e2 1e3 1e4 1e5])
    
    c_caa_plot_distribution_function(h(ii+ncol*2),'tint',tint,'cross-section',peaCnts4);
    set(h(ii+ncol*2),'xlim',[4e1 3e4],'ylim',[1e-1 1e2],'xtick',[1e1 1e2 1e3 1e4 1e5])    

    tint = tint + dt;
end
c_eval('cn.print(''C4_psd_count?'',''psd&counts'')',jj)
end
end

%% Try just adding data for C3 and C4.
% combined data
ind4=find(peaCnts4.t>toepoch([2007 08 31 10 17 39.6]),1,'first')
ind3=find(peaCnts3.t>toepoch([2007 08 31 10 17 39.6]),1,'first')
%n?nting med intersekt kanske;

diff([peaPSD3.t(ind3) peaPSD4.t(ind4)])
%% Check spacecraft orientation
c_pl_sc_orient(3,toepoch([2007 08 31 10 17 38.2]))
c_pl_sc_orient(4,toepoch([2007 08 31 10 17 37.7]))

%% Compare distributions to see if they have the sam density, or count in 
% the PA180 direction
if 1 % the beam model
n  = [0.060 0.0025]*1e6; % density in m^-3
m  = [0 0]; % electrons
t  = [1.500 0.060]; % keV
vd = [0 -4];
d  = [1 1]; % no loss cone
a1 = [1 1]; % no anisotropy
a2 = [0 0];
pitchangles = [0 90 180];
plotoption = 1;
tbeam = toepoch([2007 08 31 10 17 38.2]);
end

if 0 % the beam model 2
n  = [0.050 0.004]*1e6; % density in m^-3
m  = [0 0]; % electrons
t  = [1.100 0.060]; % keV
vd = [0 -4];
d  = [1 1]; % no loss cone
a1 = [1.5 1]; % no anisotropy
a2 = [0 0];
pitchangles = [0 90 180];
plotoption = 1;
tbeam = toepoch([2007 08 31 10 17 38.2]);
end

if 0 % the broken beam model
n  = [0.055 0.0075]*1e6; % density in m^-3
m  = [0 0]; % electrons
t  = [1.500 0.12]; % keV
vd = [0 -1.3];
d  = [1 1]; % no loss cone
a1 = [1 1]; % no anisotropy
a2 = [0 0];
pitchangles = [0 90 180];
plotoption = 1;
tbeam = toepoch([2007 08 31 10 17 46.6]);
end

if 0 % the broken beam model, C3, another intermediate time interval
n  = [0.055 0.0075]*1e6; % density in m^-3
m  = [0 0]; % electrons
t  = [1.500 0.12]; % keV
vd = [0 -1.3];
d  = [1 1]; % no loss cone
a1 = [1 1]; % no anisotropy
a2 = [0 0];
pitchangles = [0 90 180];
plotoption = 1;
tbeam = toepoch([2007 08 31 10 17 42.37]);
end

if 0 % the broken beam model, C3, another intermediate time interval
n  = [0.055 0.0075]*1e6; % density in m^-3
m  = [0 0]; % electrons
t  = [1.500 0.12]; % keV
vd = [0 -1.3];
d  = [1 1]; % no loss cone
a1 = [1 1]; % no anisotropy
a2 = [0 0];
pitchangles = [0 90 180];
plotoption = 1;
tbeam = toepoch([2007 08 31 10 17 44.000]);
end

h_ax=axes; hold(h_ax,'on')
[h,f,Etot]=whamp.plot_f(n,m,t,vd,d,a1,a2,pitchangles,plotoption);
%c_caa_plot_distribution_function(h_ax,'tint',tbeam,'cross-section',peaPSD3);

set(h_ax,'xlim',[1e1 3e5],'ylim',[1e-5 1e2],...
         'yscale','log','xscale','log',...
         'xtick',[1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5]) 
     set(h(1),'marker','o')
     set(h(2),'marker','o')
     set(h(3),'marker','o')

text(2e1, 1e-4,distr.whamp_string(n,t,vd,a1))
hold(h_ax,'off')
%%

% Plot the same thing but with only one pitch angle and all in the same
% plot, in order to better see the development. Must I change
% c_caa_plot_distribution_dunction.m in order to do that? I think I must
% add the possibility to add arguments...

figure(51)
h=subplot(1,1,1);

tstart = toepoch([2007 08 31 10 17 34]);
tt = 4;
dt = 4;

tint = tstart + [0 4];

for ii=1:6    
    c_caa_plot_distribution_function(h,'tint',tint,'cross-section',peaCnts,...
        'pitchangle',180);
    hold(h,'on')
    tint = tint + dt;
end