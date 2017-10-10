%% Read data from list of points, ehpoints.txt
% First column is core, second is satellite
fid=fopen('/Users/Cecilia/TRED46/DATA/ehpoints.txt')
data=textscan(fid,'%f %f %*[^\n]');
uniquecores=unique(data{1});
ncore=length(uniquecores);
%% Load all satellites of the cores
format = repmat('%f',1,14*27); % 14*27 numbers
tic
for p=11:15%ncore
    p
eval(['fid=fopen(''/Users/Cecilia/TRED46/VirtualSatellite/VirtualSatelliteTraces',num2str(uniquecores(p)),'.txt'');'])
Pos=textscan(fid, '%f%f%f', 27, 'delimiter', '\n', 'headerlines', 0);
eval(['Position(',num2str(p),',:,:)=[Pos{1} Pos{2} Pos{3}];'])
Data=textscan(fid, format, 'delimiter', '\n', 'headerlines', 0);
% every point has 14 own columns for each time step
ind=reshape(1:length(Data),length(Data)/27,27);
for k=1:27
    % Structure: Quantity(core,sc (1-27),time,xyz)
    eval(['B(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(1,k)),'} Data{1,',num2str(ind(2,k)),'} Data{1,',num2str(ind(3,k)),'}];']) 
    eval(['E(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(4,k)),'} Data{1,',num2str(ind(5,k)),'} Data{1,',num2str(ind(6,k)),'}];']) 
    eval(['Je(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(7,k)),'} Data{1,',num2str(ind(8,k)),'} Data{1,',num2str(ind(9,k)),'}];']) 
    eval(['Ji(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(10,k)),'} Data{1,',num2str(ind(11,k)),'} Data{1,',num2str(ind(12,k)),'}];']) 
    eval(['Ne(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(13,k)),'}];']) 
    eval(['Ni(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(14,k)),'}];']) 
end
end
toc
%% Calculate Bnorm, Epar, Jipar, Jepar
for k=1:ncore
    for l=1:27
        Epar(k,l,:,1)=irf_dot(squeeze(E(k,l,:,:)),irf_norm(squeeze(B(k,l,:,:))));
        Jipar(k,l,:,1)=irf_dot(squeeze(Ji(k,l,:,:)),irf_norm(squeeze(B(k,l,:,:))));
        Jepar(k,l,:,1)=irf_dot(squeeze(Je(k,l,:,:)),irf_norm(squeeze(B(k,l,:,:))));        
        Bnorm(k,l,:,:)=irf_norm(squeeze(B(k,l,:,:)));
    end
end
%% Plot all Epar and 'comp' within cores + 3dplot with locations
% One panel for each core
opi=50;
t=tred_norm(tocolumn(1:size(squeeze(Epar(1,1,:,:)))),'t',opi);

figure(99);h=irf_plot(ncore);
for k=1:ncore
    hca=h(k);
    uniquecores(k);
    isc=find(data{1}==uniquecores(k));
    nsc=length(isc);
    for p=1:nsc        
        %data{2}(isc(p))
        irf_plot(hca,[t tred_norm(squeeze(Epar(k,data{2}(isc(p)),:,:)),'E',opi)]); hold(hca,'on');
    end
    ylabel(hca,[' ',num2str(uniquecores(k)),' '])
end
%irf_zoom(h,'x',[40 t(end)]);
%eval(['print -dpng /Users/Cecilia/TRED46/Satellitbilder/epar/ov_epar.png']);

%%
% Locations
figure(98); 
for k=1:ncore
    isc=find(data{1}==uniquecores(k));
    nsc=length(isc);
    for p=1:nsc
        sat=data{2}(isc(p)); % satellite
        cor=data{1}(isc(p)); % satellite
        
        quiver3(Position(k,sat,1),Position(k,sat,2),Position(k,sat,3),...
                Bnorm(k,sat,16000,1)*0.3,Bnorm(k,sat,16000,2)*0.3,Bnorm(k,sat,16000,3)*0.3); hold on;        
        text(Position(k,sat,1),Position(k,sat,2),Position(k,sat,3),[num2str(cor),'-',num2str(sat)]);
    end
end
xlabel('x: xgsm');ylabel('y: zgms');zlabel('z: -ygsm')
%% Limit the data a little bit and plot
% One panel for each core
t=tocolumn(1:size(squeeze(Epar(1,1,:,:))));
cores=[437 532 533];
ncores=length(cores);
figure;h=irf_plot(ncores);
for k=1:ncores
    hca=h(k);
    uniquecores(k);
    isc=find(data{1}==cores(k));
    nsc=length(isc);
    coren=find(uniquecores==cores(k))
    
    %colors1=[ones(1,ceil(nsc/2)), ones(1,floor(nsc/2))*0.1];
    %colors2=rand(1,nsc);%[ones(1,ceil(nsc/2))*0.1, ones(1,floor(nsc/2))];
    %colors3=rand(1,nsc);%[linspace(0,1,ceil(nsc/2)) linspace(0,1,floor(nsc/2))];
    line=[cellstr(repmat('-',ceil(nsc/3),1)); cellstr(repmat('-.',ceil(nsc/3),1));...
          cellstr(repmat('--',ceil(nsc/3),1))];
    colors=[1 0 0; 0 1 0; 0 0 1; 0.7 0 0.9; 0.7 0.7 0; 0 0.7 0.7;...
            0 0 0; 1 0 0.5; 0.5 1 0; 0.5 0 1; 0.8 0.5 0.2; 0.5 0.4 0.9;...
            0.8 0.8 0.5; 0.3 0.7 0.4; 0.8 0.3 0.3; 1 0.4 0.7; 0.4 0.7 1];
        colors=repmat(colors,2,1);
    
    for p=1:nsc        
        %set(hca,'colorOrder',[colors1(p) colors2(p) colors3(p)])
        set(hca,'colorOrder',colors(p,:))
        irf_plot(hca,[t squeeze(Epar(coren,data{2}(isc(p)),:,:))],line{p}); hold(hca,'on');
    end
    ylabel(hca,[' ',num2str(cores(k)),' '])
    legendCell=cellstr(num2str(data{2}(isc)));
    legend(hca,legendCell)
end

% Locations
figure; 
for k=1:ncores
    isc=find(data{1}==cores(k));
    nsc=length(isc);
    for p=1:nsc
        sat=data{2}(isc(p)); % satellite
        cor=data{1}(isc(p)); % satellite
        
        quiver3(Position(k,sat,1),Position(k,sat,2),Position(k,sat,3),...
                Bnorm(k,sat,16600,1)*0.3,Bnorm(k,sat,16600,2)*0.3,Bnorm(k,sat,16600,3)*0.3); hold on;        
        text(Position(k,sat,1),Position(k,sat,2),Position(k,sat,3),[num2str(cor),'-',num2str(sat)]);
    end
end
xlabel('x: xgsm');ylabel('y: zgms');zlabel('z: -ygsm')
%% Plot locations of satellites that has observed an eh, only one time.
figure; 
for k=1:ncore
    isc=find(data{1}==uniquecores(k));
    nsc=length(isc);
    for p=1:nsc
        sat=data{2}(isc(p)); % satellite
        cor=data{1}(isc(p)); % satellite       
        quiver3(Position(k,sat,1),Position(k,sat,2),Position(k,sat,3),...
                Bnorm(k,sat,16000,1)*0.3,Bnorm(k,sat,16000,2)*0.3,Bnorm(k,sat,16000,3)*0.3); hold on;        
        text(Position(k,sat,1),Position(k,sat,2),Position(k,sat,3),[num2str(cor),'-',num2str(sat)]);
    end
end
xlabel('x: xgsm');ylabel('y: zgms');zlabel('z: -ygsm')
%% Use two satellites to calculate distance along B-field
% Give core, two satellites and time
core=532; sc=[14 15]; time=16450;
ic=find(uniquecores==core);
t=tocolumn(1:size(E,3));
% Calculations
dist=torow(squeeze(Position(ic,sc(2),:)-Position(ic,sc(1),:)));
Bdir=mean([torow(squeeze(Bnorm(ic,sc(2),time,:)));torow(squeeze(Bnorm(ic,sc(1),time,:)))]);
distB=dot(Bdir,dist);
% Find timelag
[xcorrs, lags]=xcorr(Epar(ic,sc(2),:),Epar(ic,sc(1),:),70);
[maxcorr, lagind]=max(squeeze(xcorrs));
tlag=lags(lagind);
subplot(2,1,1);plot(lags,squeeze(xcorrs))
subplot(2,1,2);
irf_plot({[t squeeze(Epar(ic,sc(1),:,:))],[t squeeze(Epar(ic,sc(2),:,:))],...
        [t squeeze(Epar(ic,sc(2),:,:))]},'comp','dt', [0 0 tlag])
legend([' ',num2str(core),'-',num2str(sc(1)),' '],...
       [' ',num2str(core),'-',num2str(sc(2)),' '],...
       [' ',num2str(core),'-',num2str(sc(2)),' {shift}'])
title([num2str(sc(1)),': [',num2str(Position(ic,sc(1),1)),' ',num2str(Position(ic,sc(1),2)),' ',num2str(Position(ic,sc(1),3)),']  ',...
       num2str(sc(2)),': [',num2str(Position(ic,sc(2),1)),' ',num2str(Position(ic,sc(2),2)),' ',num2str(Position(ic,sc(2),3)),']']);
    % Find velocity, distance/timelag
irf_units
velocity=distB/tlag;
vkms=velocity*Units.c*1e-3; %km/s
% Plot match

%% Use four stellites and c_4_v
% For this we need sc separated in all three dimensions, plane is not
% enough.
% Need core, four satellites and four times when the same eh was observed
cores=[532 532 532 533];%[437 532 532 533]; 
sc=[7 8 15 13];%[7 15 13 21];%[25 7 8 13];
cores=[532 532 533 533];sc=[6 15 20 21];
cores=[532 532 533 533];sc=[14 15 22 21];
for k=1:4; ic(k)=find(uniquecores==cores(k)); end
%[~,ic,~]=intersect([data{1}(:) data{2}(:)],[cores' sc'],'rows');
% Locations, can not be in a plane
figure; 
for k=1:4
    isc=find(uniquecores==cores(k));
    nsc=length(isc);
    for p=1:nsc
        sat=sc(k); % satellite
        cor=isc; % core
        
        quiver3(Position(cor,sat,1),Position(cor,sat,2),Position(cor,sat,3),...
                Bnorm(cor,sat,16000,1)*0.3,Bnorm(cor,sat,16000,2)*0.3,Bnorm(cor,sat,16000,3)*0.3); hold on;        
        text(Position(cor,sat,1),Position(cor,sat,2),Position(cor,sat,3),[num2str(cor),'-',num2str(sat)]);
    end
end
xlabel('x: xgsm');ylabel('y: zgms');zlabel('z: -ygsm');axis equal

t=tocolumn(1:size(squeeze(Epar(1,1,:,:))));
figure;h=irf_plot({[t squeeze(Epar(ic(1),sc(1),:,:))],[t squeeze(Epar(ic(2),sc(2),:,:))],[t squeeze(Epar(ic(3),sc(3),:,:))],[t squeeze(Epar(ic(4),sc(4),:,:))]},'comp')
legend(h,'1','2','3','4')
%% Then zoom in and click on times
[times,~]=ginput(4);
scale=1;
velocity=c_4_v([times(1) scale*squeeze(Position(ic(1),sc(1),:))'],...
               [times(2) scale*squeeze(Position(ic(2),sc(2),:))'],...
               [times(3) scale*squeeze(Position(ic(3),sc(3),:))'],...
               [times(4) scale*squeeze(Position(ic(4),sc(4),:))']);
%% Define normalizations
%% Create overview plot 
for p=23%:length(data{1})
sc(1)=data{2}(p);
core=data{1}(p);    
sc=7;core=532;
figure(10);
set(gcf,'PaperPositionMode','auto');
h=irf_plot(13);
set(gcf,'Position',[1 1 838 955])
isub=1;
t=tred_norm(tocolumn(squeeze(1:size(E,3))),'t',opi);

cr=find(uniquecores==core);
opi=50; % Hz
if 1 % B
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(B(cr,sc(1),:,:)),'B',opi)]);
    ylabel(hca,'B [nT]')
end
if 1 % E
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(E(cr,sc(1),:,:)),'E',opi)]);
    ylabel(hca,'E [mV/m]')    
end
if 1 % Epar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Epar(cr,sc(1),:)),'E',opi)]);
    ylabel(hca,'E_{||} [mV/m]')    
end
if 1 % Ji
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Ji(cr,sc(1),:,:)),'hej',opi)]);
    ylabel(hca,'J_i')
end
if 1 % Je
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Je(cr,sc(1),:,:)),'hej',opi)]);    
    ylabel(hca,'J_e')
end
if 1 % Ni
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Ni(cr,sc(1),:)),'Ni',opi)]);
    ylabel(hca,'N_i [cc]')
end
if 1 % Ne
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Ne(cr,sc(1),:)),'Ni',opi)]);
    ylabel(hca,'N_e [cc]')
end
if 1 % Ntot
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Ne(cr,sc(1),:))+squeeze(Ni(cr,sc(1),:)),'Ne',opi)]);
    ylabel(hca,'N_{tot} [cc]')
end
if 1 % Ve
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Je(cr,sc(1),:,:))./repmat(squeeze(Ne(cr,sc(1),:)),1,3),'Je',opi)]);
    ylabel(hca,'V_e [km/s]')
end
if 1 % Vi
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Ji(cr,sc(1),:,:))./repmat(squeeze(Ni(cr,sc(1),:)),1,3),'Ji',opi)]);
    ylabel(hca,'V_i [km/s]')
end    
if 1 % Vepar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Jepar(cr,sc(1),:,:))./squeeze(Ne(cr,sc(1),:)),'Je',opi)]);
    ylabel(hca,'V_{e,||} [km/s]')
end
if 1 % Vipar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Jipar(cr,sc(1),:,:))./squeeze(Ni(cr,sc(1),:)),'Ji',opi)]);
    ylabel(hca,'V_{i,||} [km/s]')
end  
if 1 % deltaVpar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Jipar(cr,sc(1),:,:))./squeeze(Ni(cr,sc(1),:)),'Ji',opi)+tred_norm(squeeze(Jepar(cr,sc(1),:,:))./squeeze(Ne(cr,sc(1),:)),'Ji',opi)]);
    ylabel(hca,'\Delta V_{i,e,||} [km/s]')
end  

title(h(1),[' ',num2str(core),'-',num2str(sc(1)),'    Position: ',...
            num2str(Position(cr,sc,1)),' ',num2str(Position(cr,sc,2)),...
            ' ',num2str(Position(cr,sc,3))])
irf_zoom(h,'x',[0 t(end)])
irf_zoom(h,'y')
eval(['print -dpng /Users/Cecilia/TRED46/Satellitbilder/ov_',num2str(core),'-',num2str(sc),'.png']);

irf_zoom(h,'x',[26 t(end)])
irf_zoom(h,'y')
eval(['print -dpng /Users/Cecilia/TRED46/Satellitbilder/ov_',num2str(core),'-',num2str(sc),'_zoom.png']);
end