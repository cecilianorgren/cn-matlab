function tred_load_n_plot_ov(core,sc)
format = repmat('%f',1,14*27); % 14*27 numbers
ncore=length(core);

for p=1:ncore
eval(['fid=fopen(''/Users/Cecilia/TRED46/VirtualSatellite/VirtualSatelliteTraces',num2str(core(p)),'.txt'');'])
Pos(p,:)=textscan(fid, '%f%f%f', 27, 'delimiter', '\n', 'headerlines', 0);
eval(['Position(',num2str(p),',:,:)=[Pos{1} Pos{2} Pos{3}];'])
Data=textscan(fid, format, 'delimiter', '\n', 'headerlines', 0);
% every point has 14 own columns for each time step
ind=reshape(1:length(Data),length(Data)/27,27);
    % Structure: Quantity(core,sc (1-27),time,xyz)
    eval(['B(',num2str(p),',:,:) = [Data{1,',num2str(ind(1,sc(p))),'} Data{1,',num2str(ind(2,sc(p))),'} Data{1,',num2str(ind(3,sc(p))),'}];']) 
    eval(['E(',num2str(p),',:,:) = [Data{1,',num2str(ind(4,sc(p))),'} Data{1,',num2str(ind(5,sc(p))),'} Data{1,',num2str(ind(6,sc(p))),'}];']) 
    eval(['Je(',num2str(p),',:,:) = [Data{1,',num2str(ind(7,sc(p))),'} Data{1,',num2str(ind(8,sc(p))),'} Data{1,',num2str(ind(9,sc(p))),'}];']) 
    eval(['Ji(',num2str(p),',:,:) = [Data{1,',num2str(ind(10,sc(p))),'} Data{1,',num2str(ind(11,sc(p))),'} Data{1,',num2str(ind(12,sc(p))),'}];']) 
    eval(['Ne(',num2str(p),',:,:) = [Data{1,',num2str(ind(13,sc(p))),'}];']) 
    eval(['Ni(',num2str(p),',:,:) = [Data{1,',num2str(ind(14,sc(p))),'}];']) 
end

%% Calculate Bnorm, Epar, Jipar, Jepar
for k=1:ncore
        Epar(k,:,1)=irf_dot(squeeze(E(k,:,:)),irf_norm(squeeze(B(k,:,:))));
        Jipar(k,:,1)=irf_dot(squeeze(Ji(k,:,:)),irf_norm(squeeze(B(k,:,:))));
        Jepar(k,:,1)=irf_dot(squeeze(Je(k,:,:)),irf_norm(squeeze(B(k,:,:))));        
        Bnorm(k,:,:)=irf_norm(squeeze(B(k,:,:)));

end

%% Create overview plot 
for p=1:ncore
figure;
set(gcf,'PaperPositionMode','auto');
h=irf_plot(11);
set(gcf,'Position',[1 1 838 955])
isub=1;

opi=50; % Hz
t=tred_norm(tocolumn(squeeze(1:size(E,2))),'t',opi);
if 1 % B
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(B(p,:,:)),'B',opi)]);
    ylabel(hca,'B [nT]')
end
if 1 % E
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(E(p,:,:)),'E',opi)]);
    ylabel(hca,'E [mV/m]')    
end
if 1 % Epar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tocolumn(tred_norm(squeeze(Epar(p,:)),'E',opi))]);
    ylabel(hca,'E_{||} [mV/m]')    
end
if 0 % Ji
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Ji(p,:,:)),'Ji',opi)]);
    ylabel(hca,'J_i')
end
if 0 % Je
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Je(p,:,:)),'Je',opi)]);    
    ylabel(hca,'J_e')
end
if 1 % Ni
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tocolumn(tred_norm(squeeze(Ni(p,:)),'Ni',opi))]);
    ylabel(hca,'N_i [cc]')
end
if 1 % Ne
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tocolumn(tred_norm(squeeze(Ne(p,:)),'Ni',opi))]);
    ylabel(hca,'N_e [cc]')
end
if 1 % Ntot
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tocolumn(tred_norm(squeeze(Ne(p,:))+squeeze(Ni(p,:)),'Ne',opi))]);
    ylabel(hca,'N_{tot} [cc]')
end
if 1 % Ve
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Je(p,:,:))./repmat(squeeze(tocolumn(Ne(p,:))),1,3),'Je',opi)]);
    ylabel(hca,'V_e [km/s]')
end
if 1 % Vi
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tred_norm(squeeze(Ji(p,:,:))./repmat(squeeze(tocolumn(Ni(p,:))),1,3),'Ji',opi)]);
    ylabel(hca,'V_i [km/s]')
end    
if 1 % Vepar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tocolumn(tred_norm(squeeze(Jepar(p,:,:))./squeeze(Ne(p,:)),'Je',opi))]);
    ylabel(hca,'V_{e,||} [km/s]')
end
if 1 % Vipar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tocolumn(tred_norm(squeeze(Jipar(p,:,:))./squeeze(Ni(p,:)),'Ji',opi))]);
    ylabel(hca,'V_{i,||} [km/s]')
end  
if 1 % deltaVpar
    hca=h(isub);isub=isub+1;
    irf_plot(hca,[t tocolumn(tred_norm(squeeze(Jipar(p,:,:))./squeeze(Ni(p,:)),'Ji',opi)+tred_norm(squeeze(Jepar(p,:,:))./squeeze(Ne(p,:)),'Ji',opi))]);
    ylabel(hca,'\Delta V_{i,e,||} [km/s]')
end  

title(h(1),[' ',num2str(core(p)),'-',num2str(sc(p)),'    Position: ',...
            num2str(Pos{1}(sc(p))),' ',num2str(Pos{2}(sc(p))),...
            ' ',num2str(Pos{3}(sc(p)))])
        set(gcf,'PaperPositionMode','auto');
end
    