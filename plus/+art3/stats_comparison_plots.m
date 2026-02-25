% stats_comparison_plots
% load daniels observation data
load('/Users/Cecilia/Research/EH2/Daniel/wavespeeds.mat')

dg.vphEH = wavespeeds.veholes;
dg.vTEH = wavespeeds.vtreholes;
dg.vteEH = wavespeeds.vtheholes;
dg.vphES = wavespeeds.vESwaves;
dg.vTES = wavespeeds.vtrESwaves;
dg.vteES = wavespeeds.vthESwaves;

% choose
filePaths = {'/Users/Cecilia/Research/EH2/comp2stats/warm_bistream.txt',...
             '/Users/Cecilia/Research/EH2/comp2stats/modified_buneman.txt',...
             '/Users/Cecilia/Research/EH2/comp2stats/modified_buneman2.txt',...
             '/Users/Cecilia/Research/EH2/comp2stats/modified_buneman3.txt',...
             '/Users/Cecilia/Research/EH2/comp2stats/electron_beam1.txt'};
formatReadData = {'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',...
                  '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f NaN NaN\n',...
                  '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f NaN NaN\n',...
                  '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',...
                  '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n'};         
fileIndex = 1;
%HeaderLines = [0 0 0 0 3];         

filePath = filePaths{fileIndex};
nHeader = [17 15 15 17 17];
startData = [3 0 0 0 3];
fid = fopen(filePath,'r');
headerC = textscan(fid,'%s',nHeader(fileIndex),'HeaderLines',startData(fileIndex)); headerC{:};
formatReadData = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n';
C = textscan(fid,formatReadData,'HeaderLines',1);


% Construct data structure
clear data
data.file = filePaths{fileIndex};
data.number = 1: numel(C{1});
data.vde1 = C{8};
data.vde2 = C{9};
% if vph is not in between vde1 and vde2, put those values to NaN
vvv = C{13}; 
vvv(vvv<min([C{8} C{9}],[],2)) = NaN;
vvv(vvv>max([C{8} C{9}],[],2)) = NaN;    
data.vph_auto = vvv;
data.vph_manual = C{16};
data.vT1_auto = abs(data.vph_auto-C{8});
data.vT2_auto = abs(data.vph_auto-C{9});
data.vT1_manual = abs(data.vph_manual-C{8});
data.vT2_manual = abs(data.vph_manual-C{9});
data.vTmin_auto = min([data.vT1_auto data.vT2_auto],[],2);
data.vTmax_auto = min([data.vT1_auto data.vT2_auto],[],2);
data.vTmin_manual = min([data.vT1_manual data.vT2_manual],[],2);
data.vTmax_manual = min([data.vT1_manual data.vT2_manual],[],2);


switch fileIndex
    case 1; data.vT = data.vTmin_manual; wb = data;
    case 2; data.vT = data.vT2_manual; mb = data;
    case 3; data.vT = data.vT2_manual; mb = data;
    case 4; data.vT = data.vT2_manual; mb = data;
    case 5; data.vT = data.vT2_manual; eb = data;
end
        
%% new, based on data structure

nPl = 4;
nrows = 2;
ncols = ceil(nPl/nrows);
for kk=1:nPl; h(kk) = subplot(nrows,ncols,kk); end
isub=1;

if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plot(hca,data.number,data.vde1,...
             data.number,data.vde2,...
             data.number,data.vph_auto,'o-',...
             data.number,data.vph_manual,'x--',...
             data.number,data.vT1_manual,'s',...
             data.number,data.vT2_manual,'^')
    legend(hca,'v_{de,1}','v_{de,2}','v_{ph,auto}','v_{ph,manual}','v_{T,1}','v_{T,2}','location','bestoutside')
    xlabel(hca,'Event #')
    ylabel(hca,'v [km/s]')
    title(hca,filePaths{fileIndex})
end
if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plot(hca,data.number,data.vde1,...
             data.number,data.vde2,...             
             data.number,data.vph_manual,'x--')%,...
             %data.number'*[1 1],data.vph_manual*[1 1] + data.vTmin_manual*[-1 1],'c-')
    hold(hca,'on')
    errorbar(hca,data.number,data.vph_manual,data.vTmin_manual,'c')     
    hold(hca,'off')
    legend(hca,'v_{de,1}','v_{de,2}','v_{ph,manual}','v_{T,min,range}','location','bestoutside')
    xlabel(hca,'Event #')
    ylabel(hca,'v [km/s]')
    title(hca,'Appropriate v_T for 2 counter-streaming electron populations')
end
if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plot(hca,data.number,data.vde1,...
             data.number,data.vde2,...             
             data.number,data.vph_manual,'x--')%,...
             %data.number'*[1 1],data.vph_manual*[1 1] + data.vTmin_manual*[-1 1],'c-')
    hold(hca,'on')
    errorbar(hca,data.number,data.vph_manual,data.vT2_manual,'c')     
    hold(hca,'off')
    legend(hca,'v_{de,1}','v_{de,2}','v_{ph,manual}','v_{T,2,range}','location','bestoutside')
    xlabel(hca,'Event #')
    ylabel(hca,'v [km/s]')
    title(hca,'Appropriate v_T for only one drifting electron population')
end
if 1 % vph,auto and vph,manual vs vT2 (the one relevant for mod Buneman)
    hca = h(isub); isub=isub+1;   
    plot(hca,data.vph_auto,data.vT2_manual,'o',data.vph_manual,data.vT2_manual,'x')
    xlabel(hca,'v_{ph} [km/s]'); ylabel(hca,'v_{T} [km/s]')
    legend(hca,'v_{ph,auto}','v_{ph,maunal}','v_{ph}=v_T')
    maxv = max([get(hca,'ylim') get(hca,'xlim')]);
    set(hca,'xlim',[0 maxv],'ylim',[0 maxv])
end

%% Compare 3 files

data1 = wb;
data2 = mb;
data3 = eb;

nPl = 2;
nrows = 1;
ncols = ceil(nPl/nrows);
for kk=1:nPl; h(kk) = subplot(nrows,ncols,kk); end
isub=1;


if 1 % vph,auto and vph,manual vs vT2 (the one relevant for mod Buneman)
    hca = h(isub); isub=isub+1;   
    multx = 1e-1;
    multy = 1e-1;
    plot(hca,data1.vph_manual*multx,data1.vT*multy,'o',...
             data2.vph_manual*multx,data2.vT*multy,'s',...
             data3.vph_manual*multx,data3.vT*multy,'^',...
             dg.vphEH,dg.vTEH,'c*',...
             dg.vphES,dg.vTES,'g*')
    xlabel(hca,'v_{ph} [km/s]'); ylabel(hca,'v_{T} [km/s]')    
    maxv = max([get(hca,'ylim') get(hca,'xlim')]);
    hold(hca,'on')
    plot(hca,[0 maxv],[0 maxv])
    hold(hca,'off')
    legend(hca,'electron bi-stream','slow electron beam','fast electron beam','v_{ph}=v_T','location','best')
    set(hca,'xlim',[0 maxv],'ylim',[0 maxv])
end
if 1 % vph,auto and vph,manual vs vT2 (the one relevant for mod Buneman)
    hca = h(isub); isub=isub+1;   
    multx = 15e-2;
    multy = 3e-2;
    plot(hca,dg.vphEH./dg.vteEH,dg.vTEH./dg.vteEH,'c*',...
             dg.vphES./dg.vteES,dg.vTES./dg.vteES,'g*')
    xlabel(hca,'v_{ph}/v_{te,bg}'); ylabel(hca,'v_{T}/v_{te,bg}')
    legend(hca,'EH','ES')
    maxv = max([get(hca,'ylim') get(hca,'xlim')]);
    set(hca,'xlim',[0 maxv],'ylim',[0 maxv])
end

%%
if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plx = 1: numel(C{1});
    ply1 = C{8};
    ply2 = C{9};
    ply3 = C{13};
    plot(hca,plx,ply1,plx,ply2,plx,ply3)
    legend(hca,'v_{de1}','v_{de2}','v_{ph}')
end
if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plx = C{13}; % vph,auto [km/s]
    ply = C{16};  % vph,manual [km/s]
    plot(hca,[1e2 2e4],[1e2 2e4],'r',plx,ply,'o')
    xlabel(hca,'v_{ph,auto} [km/s]')
    ylabel(hca,'v_{ph,manual} [km/s]')
    set(hca,'yscale','log','xscale','log','ylim',get(hca,'xlim'))
end
if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plx = C{13}; % vph,auto [km/s]
    ply = C{16};  % vph,manual [km/s]
    plot(hca,[1e2 2e4],[1e2 2e4],'r',plx,ply,'o')
    xlabel(hca,'v_{ph,auto} [km/s]')
    ylabel(hca,'v_{ph,manual} [km/s]')
    set(hca,'yscale','log','xscale','log','ylim',get(hca,'xlim'))
end
if 1 % vph,auto vs vde,2
    %%
    hca = h(isub); isub=isub+1;
    plx = C{13}; % vph,auto [km/s]
    ply1 = C{9};  % vde2 [km/s]
    ply2 = C{8};  % vde2 [km/s]
    plot(hca,plx*1e-3,ply1*1e-3,'bo',plx*1e-3,ply2*1e-3,'ro')
    xlabel(hca,'v_{ph,auto} [10^3 km/s]')
    ylabel(hca,'v_{de,1(red),2(blue)} [10^3 km/s]')
    set(hca,'yscale','lin','xscale','lin')
end
if 1
    hca = h(isub); isub=isub+1;
    plx = C{16}; % vph,manual [km/s]
    ply = abs(C{9}-C{16});  % vde2-vph,manual [km/s]
    plot(hca,[1e2 1.3e4],[1e2 1.3e4],'r',plx,ply,'o',wavespeeds.veholes*1e0,wavespeeds.vtreholes*1e0,'gx')
    xlabel(hca,'v_{ph,manual}')
    ylabel(hca,'v_T = |v_{de,2}-v_{ph,manual}|')
    set(hca,'yscale','lin','xscale','lin','ylim',get(hca,'xlim'))
end
if 0
    hca = h(isub); isub=isub+1;
    plx = C{16}/1000; % vph [km/s]
    ply = abs(C{9}-C{16}/1000);  % vT [km/s]
    plot(hca,plx,ply,'o')
end
if 1
    hca = h(isub); isub=isub+1;
    fac=1e-2;
    plx = abs(C{13})*fac; % vph,auto [km/s]    
    ply1  = abs(C{13}-C{8})*fac;  % vph,auto-vde1 [km/s]
    ply2 = abs(C{13}-C{9})*fac;  % vph,auto-vde2 [km/s]
    plot(hca,plx,ply1,'ro',plx,ply2,'bo',wavespeeds.veholes*1e0,wavespeeds.vtreholes*1e0,'gx')
    xlabel(hca,'|v_{ph,auto}|')
    ylabel(hca,'|v_T| = |v_{de,1,2}-v_{ph,auto}|')
    set(hca,'yscale','lin','xscale','lin','ylim',get(hca,'xlim'))
end
if 1
    hca = h(isub); isub=isub+1;
    
    Units=irf_units; % read in standard units
    Me=Units.me;
    Mp=Units.mp;
    c=Units.c;
    e=Units.e;
    epso=Units.eps0;
    mu0=Units.mu0;
    Mp_Me = Mp/Me; % ratio of proton and electron mass;
    
    fac=1e-2;
    plx = C{3}./C{1}; % R
    Vte1 = c*sqrt(1-1./(C{5}.*e./(Me*c^2)+1).^2);
    ply = abs(C{9}./Vte1*1e3);
    plot(hca,plx,ply,'kx')
    xlabel(hca,'R')
    ylabel(hca,'S = v_{de,2}/v_{te,1}')
    set(hca,'yscale','lin','xscale','lin','ylim',get(hca,'xlim'))
end    
if 1
    hca = h(isub); isub=isub+1;
    
    v_psd = -60000:100:60000;
    n_psd = numel(C{1});
    e_psd = zeros(n_psd,numel(v_psd));
    e_psd_turned = zeros(n_psd,numel(v_psd));
    for kk=1:n_psd        
        e1_psd = cn.maxwellian(v_psd,C{5}(kk),C{2}(kk),C{8}(kk),'e');
        e2_psd = cn.maxwellian(v_psd,C{6}(kk),C{3}(kk),C{9}(kk),'e'); 
        e_psd_temp = e1_psd + e2_psd;                
        e_psd(kk,:) = e_psd_temp;
        if v_psd(e_psd_temp==max(e_psd_temp))>0
            
            e_psd_temp = e_psd_temp(end:-1:1);
        end
        e_psd_turned(kk,:) = e_psd_temp;
    end
    plot(hca,v_psd*1e-3,e_psd')
    psd_lim = [0 max(e_psd(:))*1.2];
    set(hca,'xlim',50*[-1 1],'ylim',psd_lim,'yscale','lin','xscale','lin')    
    title(hca,'Input electron distributions')
    xlabel(hca,'v [10^3 km/s]')
    
    hca = h(isub); isub=isub+1;
    plot(hca,v_psd*1e-3,e_psd_turned')
    psd_lim = [0 max(e_psd(:))*1.2];
    set(hca,'xlim',50*[-1 1],'ylim',psd_lim,'yscale','lin','xscale','lin') 
    title(hca,'Input electron distributions, flipped')
    xlabel(hca,'v [10^3 km/s]')
end






%% old

if 0 
nPl = 7;
for kk=1:nPl; h(kk) = subplot(3,3,kk); end
isub=1;

if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plx = 1: numel(C{1});
    ply1 = C{8};
    ply2 = C{9};
    ply3 = C{13};
    plot(hca,plx,ply1,plx,ply2,plx,ply3)
    legend(hca,'v_{de1}','v_{de2}','v_{ph}')
end
if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plx = C{13}; % vph,auto [km/s]
    ply = C{16};  % vph,manual [km/s]
    plot(hca,[1e2 2e4],[1e2 2e4],'r',plx,ply,'o')
    xlabel(hca,'v_{ph,auto} [km/s]')
    ylabel(hca,'v_{ph,manual} [km/s]')
    set(hca,'yscale','log','xscale','log','ylim',get(hca,'xlim'))
end
if 1 % vph,auto vs vph,manual
    hca = h(isub); isub=isub+1;
    plx = C{13}; % vph,auto [km/s]
    ply = C{16};  % vph,manual [km/s]
    plot(hca,[1e2 2e4],[1e2 2e4],'r',plx,ply,'o')
    xlabel(hca,'v_{ph,auto} [km/s]')
    ylabel(hca,'v_{ph,manual} [km/s]')
    set(hca,'yscale','log','xscale','log','ylim',get(hca,'xlim'))
end
if 1 % vph,auto vs vde,2
    %%
    hca = h(isub); isub=isub+1;
    plx = C{13}; % vph,auto [km/s]
    ply1 = C{9};  % vde2 [km/s]
    ply2 = C{8};  % vde2 [km/s]
    plot(hca,plx*1e-3,ply1*1e-3,'bo',plx*1e-3,ply2*1e-3,'ro')
    xlabel(hca,'v_{ph,auto} [10^3 km/s]')
    ylabel(hca,'v_{de,1(red),2(blue)} [10^3 km/s]')
    set(hca,'yscale','lin','xscale','lin')
end
if 1
    hca = h(isub); isub=isub+1;
    plx = C{16}; % vph,manual [km/s]
    ply = abs(C{9}-C{16});  % vde2-vph,manual [km/s]
    plot(hca,[1e2 1.3e4],[1e2 1.3e4],'r',plx,ply,'o',wavespeeds.veholes*1e0,wavespeeds.vtreholes*1e0,'gx')
    xlabel(hca,'v_{ph,manual}')
    ylabel(hca,'v_T = |v_{de,2}-v_{ph,manual}|')
    set(hca,'yscale','lin','xscale','lin','ylim',get(hca,'xlim'))
end
if 0
    hca = h(isub); isub=isub+1;
    plx = C{16}/1000; % vph [km/s]
    ply = abs(C{9}-C{16}/1000);  % vT [km/s]
    plot(hca,plx,ply,'o')
end
if 1
    hca = h(isub); isub=isub+1;
    fac=1e-2;
    plx = abs(C{13})*fac; % vph,auto [km/s]    
    ply1  = abs(C{13}-C{8})*fac;  % vph,auto-vde1 [km/s]
    ply2 = abs(C{13}-C{9})*fac;  % vph,auto-vde2 [km/s]
    plot(hca,plx,ply1,'ro',plx,ply2,'bo',wavespeeds.veholes*1e0,wavespeeds.vtreholes*1e0,'gx')
    xlabel(hca,'|v_{ph,auto}|')
    ylabel(hca,'|v_T| = |v_{de,1,2}-v_{ph,auto}|')
    set(hca,'yscale','lin','xscale','lin','ylim',get(hca,'xlim'))
end
if 1
    hca = h(isub); isub=isub+1;
    
    Units=irf_units; % read in standard units
    Me=Units.me;
    Mp=Units.mp;
    c=Units.c;
    e=Units.e;
    epso=Units.eps0;
    mu0=Units.mu0;
    Mp_Me = Mp/Me; % ratio of proton and electron mass;
    
    fac=1e-2;
    plx = C{3}./C{1}; % R
    Vte1 = c*sqrt(1-1./(C{5}.*e./(Me*c^2)+1).^2);
    ply = abs(C{9}./Vte1*1e3);
    plot(hca,plx,ply,'kx')
    xlabel(hca,'R')
    ylabel(hca,'S = v_{de,2}/v_{te,1}')
    set(hca,'yscale','lin','xscale','lin','ylim',get(hca,'xlim'))
end    
if 1
    hca = h(isub); isub=isub+1;
    
    v_psd = -60000:100:60000;
    n_psd = numel(C{1});
    e_psd = zeros(n_psd,numel(v_psd));
    e_psd_turned = zeros(n_psd,numel(v_psd));
    for kk=1:n_psd        
        e1_psd = cn.maxwellian(v_psd,C{5}(kk),C{2}(kk),C{8}(kk),'e');
        e2_psd = cn.maxwellian(v_psd,C{6}(kk),C{3}(kk),C{9}(kk),'e'); 
        e_psd_temp = e1_psd + e2_psd;                
        e_psd(kk,:) = e_psd_temp;
        if v_psd(e_psd_temp==max(e_psd_temp))>0
            
            e_psd_temp = e_psd_temp(end:-1:1);
        end
        e_psd_turned(kk,:) = e_psd_temp;
    end
    plot(hca,v_psd*1e-3,e_psd')
    psd_lim = [0 max(e_psd(:))*1.2];
    set(hca,'xlim',50*[-1 1],'ylim',psd_lim,'yscale','lin','xscale','lin')    
    title(hca,'Input electron distributions')
    xlabel(hca,'v [10^3 km/s]')
    
    hca = h(isub); isub=isub+1;
    plot(hca,v_psd*1e-3,e_psd_turned')
    psd_lim = [0 max(e_psd(:))*1.2];
    set(hca,'xlim',50*[-1 1],'ylim',psd_lim,'yscale','lin','xscale','lin') 
    title(hca,'Input electron distributions, flipped')
    xlabel(hca,'v [10^3 km/s]')
end


end



