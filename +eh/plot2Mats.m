for oo = 1:numel(axesObjs); try delete(axesObjs(oo)); end; end
clear axesObjs dataObjs
cd /Users/Cecilia/Data/BM/20070831/ % go to folder where data is saved
t1 = toepoch([2007 08 31 10 17 38.140]); % beam
t2 = toepoch([2007 08 31 10 17 39.689]); % beam
t3 = toepoch([2007 08 31 10 17 42.419]); % smearedbeam
t4 = toepoch([2007 08 31 10 17 43.840]); % smearedbeam
if ~exist('peaPSD3','var'); load peaPSD.mat; end % load data
pitchangles = [165];
% make figure
h = subplot(1,1,1); hold(h,'all'); set(h,'xscale','log','yscale','log')
% beam C3
c_caa_plot_distribution_function(h,'tint',t1,'cross-section',...
    'pitchangle',pitchangles,peaPSD3);
% beam C4
c_caa_plot_distribution_function(h,'tint',t2,'cross-section',...
    'pitchangle',pitchangles,peaPSD4);
% smeared C3
c_caa_plot_distribution_function(h,'tint',t3,'cross-section',...
    'pitchangle',pitchangles,peaPSD3);
% smeared beam C4
c_caa_plot_distribution_function(h,'tint',t4,'cross-section',...
    'pitchangle',pitchangles,peaPSD4);
hold(h,'off');

% Remake a little to do labeling better etc
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes  

% Delete all but one background line
delete(axesObjs(3))
delete(axesObjs(5))
delete(axesObjs(7))

c_eval('strTime? = irf_time(t?,''epoch2iso''); strTime?=strTime?(12:23);',1:4)
legend(['C3 - ' strTime1],['C4 - ' strTime2],['C3 - ' strTime3],['C4 - ' strTime4],'background','location','southwest')

title(h,'Particle distribution in direction antiparallel to B')
box(h,'on')
set(h,'ylim',[1e-5 1e2])
%%
%%
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
zdata = get(dataObjs, 'ZData');





