function gui_changeplotparameter(disprel,k,ope,oce,omref)
close all

h.fig = figure('position',[100 100 500 400]);
h.slider1 = uicontrol('style','slider','position',[20 20 450 20]);
h.slider2 = uicontrol('style','slider','position',[20 50 450 20]);                
h.ax=axes('position',[0.12    0.30    0.77    0.60]);

h.k=k;
h.omref=omref*ones(size(k));
h.ope=mean(ope);
h.oce=mean(oce);
h.disprel=disprel;
%h.ylim=[0 max([disprel(k,ope(1),oce(1)) disprel(k,ope(1),oce(2)),...
%               disprel(k,ope(2),oce(1)) disprel(k,ope(2),oce(2))])];

% initial plot
plotall(h);

% set up sliders
set(h.slider1,'callback',{@changeOpe,h});
set(h.slider2,'callback',{@changeOce,h});

range=ope; stepsize=[0.1 2];
set(h.slider1,'min',range(1),'max',range(2),...
    'sliderstep',stepsize/diff(range),'value',h.ope)
set(h.slider2,'min',range(1),'max',range(2),...
    'sliderstep',stepsize/diff(range),'value',h.oce)

function h = changeOpe(hObject,eventdata,h)
h.ope=get(h.slider1,'value');
h.oce=get(h.slider2,'value'); % why do i need to add these?
plotall(h);

function h = changeOce(hObject,eventdata,h)
h.ope=get(h.slider1,'value'); % why do i need to add these?
h.oce=get(h.slider2,'value');
plotall(h);

function plotall(h)
h.c=circlepos(h);
plot(h.ax,h.k,h.disprel(h.k,h.ope,h.oce),h.k,h.omref,h.c(1),h.c(2),'ro')
title(h.ax,'Localization of wave (given frequency) on dispersion curve as background parameter changes')
xlabel(h.ax,'k')
ylabel(h.ax,'\omega')
%h.ylim
%set(h.ax,'ylim',h.ylim)

function circlePosition = circlepos(h)
% find intersection or closest value in order to mark by circle
diff=abs(h.disprel(h.k,h.ope,h.oce)-h.omref);
minIndex = find(diff==min(diff));
if minIndex==1; minIndex=2; end
circlePosition = [h.k(minIndex) h.disprel(h.k(minIndex),h.ope,h.oce)];