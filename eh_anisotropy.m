if 0 % % % % % % make anisotropy plot
% load
varname=irf_ssub('Data__C?_CP_PEA_PITCH_3DXH_DEFlux',4);

res=c_caa_construct_subspin_res_data(varname); 

% divide into parallel antiparallel and perpendicular
% 1:2 is 0-30
% 3:10 is 30-150
% 11:12 is 150-180
res_size=size(res.data);
flux=zeros(res_size(1),3,res_size(3));
flux(:,1,:)=nanmean(res.data(:,1:2,:),2);
flux(:,2,:)=nanmean(res.data(:,3:10,:),2);
flux(:,3,:)=nanmean(res.data(:,11:12,:),2);
%
for k=1:3
specr{k}.t=res.tt;
specr{k}.f=res.en;
specr{k}.p=squeeze(flux(:,k,:));
end
figure;

h=irf_plot(3);
for k=1:3;
irf_spectrogram(h(k),specr{k});
end
end

%% MUST MAKE SPIN DATA, NTO SUBSPIN
if ~exist('diE3','var'); load matlabE; end
if ~exist('gsmB3','var'); load matlabB; end
sc=3;
%varname=irf_ssub('Data__C?_CP_PEA_PITCH_3DXH_DEFlux',sc);
varname=irf_ssub('C?_CP_PEA_PITCH_3DXH_DEFlux',sc);
xtra_vn=''; % could be _HEEA or _LEEA when using pitch full
[var,dobj,varmat,varunits]=c_caa_var_get(['Data',xtra_vn,'__',varname]);
varname2=irf_ssub(['Sweep_Energy',xtra_vn,'__',varname],sc);
%varname3=irf_ssub('Sweep_PitchAngle__C?_CP_PEA_PITCH_3DXH_DEFlux',sc);
energies=c_caa_var_get(varname2,'mat');
energy_units=c_caa_var_get(varname2,'units'); % These give eV, should be keV
%piichangle=c_caa_var_get(varname3,'mat');
%
or_var=var;
% Taking away all fluxes that are below one, these are something else,
% what? FILLVALS
%ind_nodata=find(var.data(:,:,:)<0);
%var.data(ind_nodata)=NaN;

% Taking away all fluxes that are below one, these are something else,
% what?
var=or_var;
var.data(or_var.data == or_var.FILLVAL)=NaN;

% taking mean of asimuthal angles
deflux=squeeze(nanmean(var.data(:,:,:,:),2));

%
% Making the basic flux structure
b_flux.t=energies(:,1);
b_flux.f=energies(1,2:end); nan_en=~isnan(b_flux.f);
b_flux.f=b_flux.f(nan_en);
b_flux.p=deflux(:,:,nan_en);
%lowe_ind=find(b_flux.f>0.009e4);b_flux.f=b_flux.f(lowe_ind); b_flux.p=b_flux.p(:,:,lowe_ind);
b_flux.f_units='keV';
%%
if 1 % combine c3 and c4
    b_flux.f=b_flux3.f;
    b_flux.t=b_flux3.t;
    b_flux.p=zeros(size(b_flux3));
    
end
%%
% Taking away all fluxes that are below one, these are something else,
% what?
%ind_nodata=find(b_flux.p(:,:,:)<1);
%b_flux.p(ind_nodata)=NaN;


% This is the division I make 7.5/22.5/37.5 etc. (this is the center pitch
% angle.
clear ind;
ind{1}=1:2;     % 0-30
ind{2}=3:10;    % 30-150
ind{3}=11:12;   % 150-180

for k=1:3
    flux{k}.t=b_flux.t;
    flux{k}.f=b_flux.f;
    flux{k}.p=squeeze(nanmean(b_flux.p(:,ind{k},:),2));
    flux{k}.p_units=varunits;
    flux{k}.f_units='keV'; % energy_units; % gives eV which is wrong
end
%
% now to the anisotropies
aind=0;
if 0  % (par+apar)/2per-1    
    aind=aind+1;
    an{aind}.p=(flux{1}.p+flux{3}.p)./2./flux{2}.p-1;
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='(par+apar)/ \newline 2per-1';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 0 % par/apar-1
    aind=aind+1;
    an{aind}.p=(flux{1}.p./flux{3}.p)-1;
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='par/apar-1';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 0 % log(par/apar)
    aind=aind+1;
    an{aind}.p=log((flux{1}.p./flux{3}.p));
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='log(par/apar)';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-2 2];
end
if 0 % par/per-1
    aind=aind+1;
    an{aind}.p=(flux{1}.p./flux{2}.p)-1;
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='par/per-1';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 0 % apar/per-1
    aind=aind+1;
    an{aind}.p=(flux{3}.p./flux{2}.p)-1;
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='apar/per-1';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 0 % (par)/(tot)
    aind=aind+1;
    an{aind}.p=flux{1}.p./(flux{1}.p+flux{2}.p+flux{3}.p);
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='par/tot';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 0 % (per)/(tot)
    aind=aind+1;
    an{aind}.p=flux{2}.p./(flux{1}.p+flux{2}.p+flux{3}.p);
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='per/tot';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 0 % (par)/(tot)
    aind=aind+1;
    an{aind}.p=flux{3}.p./(flux{1}.p+flux{2}.p+flux{3}.p);
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='apar/tot';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 0 % (per-par)/(par+par) , recommended by huishan
    aind=aind+1;
    an{aind}.p=(flux{2}.p-0.5*(flux{1}.p+flux{3}.p))./(flux{2}.p+0.5*(flux{1}.p+flux{3}.p));
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='(per-par)/\newline (par+par)';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 0 % (par-apar)/(par+apar)
    aind=aind+1;
    an{aind}.p=(flux{1}.p-flux{3}.p)./(flux{1}.p+flux{3}.p);
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='(par-apar)/\newline(par+apar)';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end
if 1 % log(per/(par+apar)/2) , recommended by huishan
    aind=aind+1;
    an{aind}.p=log(flux{2}.p./(0.5*(flux{1}.p+flux{3}.p)));
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='log(2*per/\newline (par+apar)';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-2 2];
end
if 1 % log(par/apar) , recommended by huishan
    aind=aind+1;
    an{aind}.p=log(flux{1}.p./flux{3}.p);
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='log(par/apar)';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-2 2];
end
if 0 % alfa=fpar/fper-1 or fper/fpar-1 , the one huishan used, andris idea
    flux_par_apar_mean=(flux{1}.p+flux{3}.p)/2;
    ind_alfa_par=find(flux_par_apar_mean>=flux{2}.p);
    ind_alfa_per=find(flux{2}.p>flux_par_apar_mean);
    aniso=zeros(size(squeeze(flux{2}.p)));
    aniso(ind_alfa_par)=1-flux{2}.p(ind_alfa_par)./(flux{1}.p(ind_alfa_par)+flux{3}.p(ind_alfa_par))./2;
    aniso(ind_alfa_per)=(flux{1}.p(ind_alfa_per)+flux{3}.p(ind_alfa_per))./flux{2}.p(ind_alfa_per)./2./-1;
    aind=aind+1;
    an{aind}.p=aniso;
    %log(flux{2}.p(ind_alfa_par)./(0.5*(flux{1}.p(ind_alfa_par)+flux{3}.p(ind_alfa_par))));
    an{aind}.t=b_flux.t;
    an{aind}.f=b_flux.f;
    an{aind}.p_units='\alpha= \newline par/per-1 \newline or =per/par-1';
    an{aind}.f_units=b_flux.f_units;
    an{aind}.c_axis=[-1 1];
end

% make the plot
% time interval
tint=bm_tint;
tint=tint(9,:);

nf=3; % the 3 fluxes
nani=aind; % the anisotropies
ne=2; % 2 E-fields + deltaB
figure(77);set(77,'position',[1 -21 838 955]);
h=irf_plot(nf+nani+ne);

for k=1:nf; % basic fluxes
    irf_plot(h(k),flux{k});
    set(h(k),'yscale','log');
    set(h(k),'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    hcmap(k)=colorbar('peer',h(k));
    ylabel(h(k),['Energy \newline [',flux{k}.f_units,']'],'fontsize',10)
    ylabel(hcmap(k),flux{k}.p_units,'fontsize',10)
    caxis(h(k),[3 7])
end
for k=1:nani % anisotropies
    irf_plot(h(nf+k),an{k});
    set(h(nf+k),'yscale','log');
    set(h(nf+k),'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    hcmap(nf+k)=colorbar('peer',h(nf+k));
    ylabel(h(nf+k),['Energy \newline [',an{k}.f_units,']'],'fontsize',10)
    ylabel(hcmap(nf+k),an{k}.p_units,'fontsize',10)
    caxis(h(nf+k),an{k}.c_axis)
end
if 1 % E-fields + deltaB
    isub=nf+nani+1;
    hca=h(isub);isub=isub+1; 
    irf_plot(hca,{irf_tlim(diE3(:,[1 2]),tint),irf_tlim(diE4(:,[1 2]),tint)},'comp');
    irf_legend(hca,{'C3','C4'},[0.02 0.9]); ylabel(hca,'E_{X,ISR2}');grid(hca,'off');
    hca=h(isub);isub=isub+1; 
    irf_plot(hca,{irf_tlim(diE3(:,[1 3]),tint),irf_tlim(diE4(:,[1 3]),tint)},'comp');
    irf_legend(hca,{'C3','C4'},[0.02 0.9]); ylabel(hca,'E_{Y,ISR2} ');grid(hca,'off');
    if 0
    hca=h(isub);isub=isub+1; 
    irf_plot(hca,irf_add(1,irf_tlim(gsmB3,tint),-1,irf_tlim(gsmB4,tint)));
    irf_legend(hca,{'x','y','z'},[0.02 0.9]); ylabel(hca,'dB_{GSM} ');grid(hca,'off');
    irf_zoom(hca,'y',[-0.6 0.6])
    end
end
if 1 % special corrections to limits etc
    for k=1:3; clims(k,:)=caxis(h(k)); end        
    %climits=[min(clims(:,1)) max(clims(:,2))];
    climits=[4 8];
    for k=1:3; caxis(h(k),climits); end
    %caxis(h(9),[-1 1])
    %caxis(h(12),[-1 1])
end
if 1 % custom colormap
    it=0:.02:1;it=it';
    xcm=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 ...
        flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]];
    clear it;
    colormap(xcm)
end
if 1 % mark eh intervals
    [tint_eh nscobs comments] = eh_tint;
    scncobs=nscobs;
    % mark all eh's, and especially the presented one.
    % prio 1/2/3/4/5/rubbish/dl/presented interval
    color={'green','green','green','yellow','yellow','white','blue','red'};
    for p = 1:size(tint_eh,1)        
        cin=scncobs(p);        
        irf_pl_mark(h,tint_eh{p},color{cin})               
    end    
end

title(h(1),[varname,'  ',xtra_vn])
linkaxes(h,'x')
irf_zoom(h,'x',[gsmB3(1,1) gsmB3(end,1)])
cn_plot_axis_align(h);
set(gcf,'PaperPositionMode','auto');
if 0 % printing three different zooms
    eval(['print -dpng /Users/Cecilia/EH/Pics/ANI_',varname,xtra_vn,'_1.png'])
    
    tint2=[toepoch([2007 08 31 10 14 00]) toepoch([2007 08 31 10 20 00])];
    irf_zoom(h,'x',tint2);
    eval(['print -dpng /Users/Cecilia/EH/Pics/ANI_',varname,xtra_vn,'_2.png'])
    
    tint3=[toepoch([2007 08 31 10 16 30]) toepoch([2007 08 31 10 19 30])];
    irf_zoom(h,'x',tint3)
    eval(['print -dpng /Users/Cecilia/EH/Pics/ANI_',varname,xtra_vn,'_3.png'])
end