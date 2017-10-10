function return_j=cn_check_j(tint,gsmB3,gsmB4,gsePos3,gsePos4,cosys,method)
% %%%%% check j %%%%%%
% inputs: tint, gsmB3, gsmB4, gsePos3, gsePos4, coordinate system, method
% method = 'av' - fixed coordinate system during the whole time interval
% method = 'mov' - the coordinate system moves with b for each time step

%% Clear variables
%clear b3 b4 b gsmdb z zav sagsm3 sagsm4 sagsmav spgsmav facb3 facb4 facdb
%clear gsmdist facdist facdistav facj gsmj
%%
%t1=toepoch([2007 08 31 10 13 10]); t2=toepoch([2007 08 31 10 13 50]);
%t1=toepoch([2007 08 31 10 19 00]); t2=toepoch([2007 08 31 10 19 10]);
%% 

t1=tint(1); t2=tint(2);
%%
%te1=toepoch(t1);
%te2=toepoch(t2);
%%

% b-field
gsmb3=irf_tlim(gsmB3(:,1:5),t1,t2);   % b3
gsmb4=irf_tlim(gsmB4(:,1:5),t1,t2);   % b4
gsmdb=irf_add(1,gsmb3(:,1:4),-1,gsmb4(:,1:4));  % db
gsmb=irf_add(0.5,gsmb3,0.5,gsmb4);   % bav
z=irf_norm(gsmb(:,1:4));       % xav for each time step

%% positions
c_eval('gsmPos?=c_coord_trans(''gse'',''gsm'',cn_toepoch(t1,t2,gsePos?),''cl_id'',?);',3:4)
c_eval('gsmPos?=irf_resamp(gsmPos?,gsmB?);',3:4);
gsmdist=irf_add(1,gsmPos3,-1,gsmPos4);


%% Field aligned coordinate system
if strcmp(lower(cosys),'fac')
    c_eval('sagsm?=c_coord_trans(''dsi'',''gsm'',[z(:,1) repmat([0 0 1],size(gsmb?,1),1)],''cl_id'',?);',3:4)
                            % spin axis for each sc and each time step
    sagsm=irf_add(0.5,sagsm3,0.5,sagsm4); % spin axis average between sc

    if strcmp(lower(method),'av') % fixed coordinate system during the whole time interval
        z=irf_norm(mean(z(:,2:4),1));           % xav during time interval
        sagsm=irf_norm(mean(sagsm(:,2:4),1));   % spin axis av during time interval
    elseif strcmp(lower(method),'mov') % the coordinate system moves with b for each time step
        z=z(:,1:4);
        sagsm=irf_norm(sagsm(:,1:4));
    end
            
    y=irf_norm(irf_cross(sagsm,z)); % in spin plane 
    y=irf_cross(z,irf_cross([-0.09 -0.52 0.85],z));    
    x=irf_cross(y,z);
    facb33=irf_lmn(gsmb3,x,y,z);
    facb3=irf_lmn(gsmb3,x,y,z);
    facb4=irf_lmn(gsmb4,x,y,z);
    facdb=irf_lmn(gsmdb,x,y,z); % fac delta b based on bav and spav
    if size(facb3,2)<5
        facb3=irf_abs(facb3);
        facb4=irf_abs(facb4);
    end

    facdist=irf_lmn(gsmdist,x,y,z);
    
    %% current
    facj=irf_abs(cn_j(facdb,facdist)); % A
    facj(:,2:5)=facj(:,2:5)*1e9; % nA
    return_j=facj;
else
    gsmj=irf_abs(cn_j(gsmdb,gsmdist)); % A
    gsmj(:,2:5)=gsmj(:,2:5)*1e9; % nA 
    method='';
    yy=gsmb;
    yy(:,2:3)=yy(:,2:3)*0;
    v3=[0.8154    0.4595    0.3523];
    v3=[-0.8060    0.5386    0.2454];
    gsmb
    yy=irf_norm([gsmb(:,1) repmat(v3,size(gsmb,1),1)]);
    return_j=irf_lmn(gsmj,gsmb,yy);
    %return_j=gsmj;
end 



%% what to plot
if strcmp(lower(cosys),'gsm')
    figure;h=irf_plot(8);isub=1;
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{gsmb3(:,[1 2]),gsmb4(:,[1 2])},'comp');
    ylabel(hca,'B_x  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{gsmb3(:,[1 3]),gsmb4(:,[1 3])},'comp');
    ylabel(hca,'B_y  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{gsmb3(:,[1 4]),gsmb4(:,[1 4])},'comp');
    ylabel(hca,'B_z  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{gsmb3(:,[1 5]),gsmb4(:,[1 5])},'comp');
    ylabel(hca,'|B|  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,gsmdb(:,1:4));ylabel(hca,'\Delta B  [nT]');irf_legend(hca,{'x','y','z'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,gsmdist);ylabel(hca,'s/c separation  [km]');irf_legend(hca,{'x','y','z'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,gsmj(:,1:4));ylabel(hca,'j  [nA]');irf_legend(hca,{'x','y','z','|j|'},[0.02 0.04])
    irf_zoom(h,'x',tint);
    hca=h(isub);isub=isub+1;
    irf_plot(hca,return_j(:,1:4));ylabel(hca,'j  [nA]');irf_legend(hca,{'x_{perp}','y_{perp}','z (B)','|j|'},[0.02 0.04])
    irf_zoom(h,'x',tint);
    title(h(1),upper(cosys));
    
elseif strcmp(lower(cosys),'fac')
    figure;h=irf_plot(7);isub=1;
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{facb3(:,[1 2]),facb4(:,[1 2])},'comp');
    ylabel(hca,'B_x  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{facb3(:,[1 3]),facb4(:,[1 3])},'comp');
    ylabel(hca,'B_y  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{facb3(:,[1 4]),facb4(:,[1 4])},'comp');
    ylabel(hca,'B_z  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,{facb3(:,[1 5]),facb4(:,[1 5])},'comp');
    ylabel(hca,'|B|  [nT]');irf_legend(hca,{'C3','C4'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,facdb(:,1:4));ylabel(hca,'\Delta B  [nT]');irf_legend(hca,{'x','y (spin plane)','z (B)'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,facdist);ylabel(hca,'s/c separation  [km]');irf_legend(hca,{'x','y (spin plane)','z (B)'},[0.02 0.04])
    hca=h(isub);isub=isub+1;
    irf_plot(hca,facj(:,1:4));ylabel(hca,'j  [nA]');irf_legend(hca,{'x','y (spin plane)','z (B)','|j|'},[0.02 0.04])
    irf_zoom(h,'x',tint);
    title(h(1),upper([cosys,'   ',method]));
    
end
%% plot
%% print
if 1 % make date/time strings
    t1str_p=datestr(epoch2date(t1),'HHMMSSFFF');
    t2str_p=datestr(epoch2date(t2),'MMSSFFF');
end
set(gcf,'PaperPositionMode','auto');
eval(['print -depsc2 /Users/Cecilia/Dropbox/Cecilia/EH/current/',t1str_p,'-',t2str_p,'_j_',cosys,'_',method,'.eps']);
