%% Make energy flux from particle flux
varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',3);
[var,dobj,varmat,varunits]=c_caa_var_get(varname);
energies=varmat.dep_x{2}.data; % 2166x44 energy
phi=varmat.dep_x{1}.data; % 2166x12 pitch angles
flux=varmat.data; % particle flux
time=varmat.t;
flux(isnan(flux))=0.000;
fluxsum0=sum(flux,2);
fluxsum=zeros(size(fluxsum0,1),size(fluxsum0,3));

for k=1:44;
fluxsum(:,k)=fluxsum0(:,1,k);
end
fluxsum(:,26:end)=[];
energies(:,26:end)=[];
energies(isnan(energies))=0;
eng=log10(fluxsum)+log10(energies);
eng2=10.^eng;
eng2(isinf(eng2))=9;
%eng(isinf(eng))=0;

figure;%h=irf_plot(2);
pc=pcolor(time(:,1),energies(1,:),eng'); shading flat; 
%hcb = colorbar('peer',pc);
%ylabel(hcb,{'log_{10} dEF';'keV/cm^2s sr keV'},'fontsize',12)
     
irf_colormap(gca,'space');
irf_timeaxis
set(gca,'yscale','log','TickDir','Out')

   
%%  Plot subtimescale electron data
varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',3);
[var,dobj,varmat,varunits]=c_caa_var_get(varname);
energies=double(varmat.dep_x{2}.data); % 2166x44 energy
phi=double(varmat.dep_x{1}.data); % 2166x12 pitch angles
flux=double(varmat.data); % particle flux
time=double(varmat.t);

%%
tint3=[cn_toepoch(t1) cn_toepoch(t2)];
tint3=tint;
res=c_caa_construct_subspin_res_data(irf_ssub('Data__C3_CP_PEA_3DXPH_PSD',3));
[delmett,ind]=irf_tlim(res.tt,tint3);
specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
  
