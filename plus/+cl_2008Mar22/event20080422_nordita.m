% run_single_event
% Specify data folder
%cd /Users/Cecilia/Data/Cluster/20070417/
%savePath = '/Users/Cecilia/Research/LH2/20070417/singles/';
savePath = '/Users/Cecilia/Research/LH2/20080422/nordita/';

% Input parameters
sc = 1;
%tint = toepoch([2007 04 17 15 34 22.02; 2007 04 17 15 34 22.10])';
%tint = toepoch([2007 04 17 15 32 50.70; 2007 04 17 15 32 51.00])';
tint_inout = toepoch([2008 04 22 17 50 00;2008 04 22 18 15 00])';
tint_out = toepoch([2008 04 22 17 35 30;2008 04 22 17 37 30])';
em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])';
em_tint = toepoch([2008 04 22 17 37 13.0;2008 04 22 17 37 14.0])';
em_tint = toepoch([2008 04 22 17 37 14.1;2008 04 22 17 37 14.3])';
em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])';
%em_tint = toepoch([2008 04 22 18 06 05;2008 04 22 18 06 07])';
csys = 'gsm';
flim = 0.1;%.1 + [-.1 .1];
v=linspace(20,800,30);
tint = em_tint;
%save_diE1 = diE1;
%diE1(:,4)=0;
% Obtain correlation parameters
tool.single_event
%diE1 = save_diE1;
doPrint = 0;
if 1
    %%
    tint_zoom = tint_out;
    h=irf_plot(6,'newfigure');
    hca = irf_panel('B'); c_eval('irf_plot(hca,irf_tlim(irf_abs(gsmB?),tint_zoom+[-4 4]));',sc); hca.YLabel.String = 'B [nT]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    hca = irf_panel('Vi'); c_eval('irf_plot(hca,irf_tlim(gsmVi?,tint_zoom+[-4 4]));',sc); hca.YLabel.String = 'V_i [km/s]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    %irf_legend(hca,{['X','Y']},[0.02 0.02])
    %hold(hca,'on')
    %co = [0 0 0; 0 0 1; 1 0 0;     
    % h = ax;        
    % set(h,'ColorOrder',co)
    %c_eval('irf_plot(hca,irf_tlim(gsmVi?,tint_zoom));',sc); hca.YLabel.String = 'V_i [km/s]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    hca = irf_panel('E'); c_eval('irf_plot(hca,irf_tlim(gsmE?,tint_zoom+[-4 4]));',sc); hca.YLabel.String = 'E [mv/m]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    irf_zoom(h(1:3),'x',tint_zoom)
    irf_pl_mark(h(1:3),tint)
    delete(h(4))
    hca = axes('Position',h(5).Position); delete(h(5)); h(5) = hca;
    irf_match_phibe_vis(hca,'best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir,flim,flh_loc);
    
    irf_legend(h(5),['E_n = ' num2str(avEn,'%.f') ' mV/m   E_k = ' num2str(avEk,'%.f') ' mV/m   B_0 = ' num2str(B0,'%.f') ' nT   n = ' num2str(n_loc,'%.1f') ' cc' ],[0.02 -.5])
    irf_legend(h(5),['v_{ExB} = [' num2str(irf_norm(vExB),'%.2f') ']\times' num2str(irf_abs(vExB,1),'%.f') 'km/s' ],[0.02 -.75])
    irf_legend(h(5),['v-v_{ExB} = [' num2str(irf_norm(direction*velocity-vExB),'%.2f') ']\times' num2str(irf_abs(direction*velocity-vExB,1),'%.f') ' km/s' ],[0.02 -1.0])
    irf_legend(h(5),['\theta_{v>v_{ExB}} = ' num2str(acosd(direction*irf_norm(vExB)'),'%.f') ' ^{\circ}'],[0.02 -1.25])
    %text(h(5).XLim(1),h(5).YLim(2)*1,'a')
    delete(h(6))
    irf_plot_zoomin_lines_between_panels(h(3),h(5))
    title(h(1),['C' num2str(sc) '  ' upper(csys)])
    irf_zoom(h(1:3),'x',tint_zoom)
    if doPrint; cn.print([ 'C' num2str(sc) '_' csys '_inclov_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_best'],'path',savePath); end
else
    %%
    h=irf_plot(4,'newfigure'); set(gcf,'position',[10   239   565   466])
    delete(h(1)); delete(h(3)); delete(h(4));
    irf_match_phibe_vis(h(2),'best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir,flim,flh_loc);
    
    irf_legend(h(2),['E_n = ' num2str(avEn,'%.f') ' mV/m   E_k = ' num2str(avEk,'%.f') ' mV/m   B_0 = ' num2str(B0,'%.f') ' nT   n = ' num2str(n_loc,'%.1f') ' cc' ],[0.02 -.5])
    irf_legend(h(2),['v_{ExB} = [' num2str(irf_norm(vExB),'%.2f') ']\times' num2str(irf_abs(vExB,1),'%.f') 'km/s' ],[0.02 -.75])
    irf_legend(h(2),['v-v_{ExB} = [' num2str(irf_norm(direction*velocity-irf_norm(vExB)),'%.2f') ']\times' num2str(irf_abs(direction*velocity-vExB,1),'%.f') ' km/s' ],[0.02 -1.0])
    irf_legend(h(2),['\theta_{v>v_{ExB}} = ' num2str(acosd(direction*irf_norm(vExB)'),'%.f') ' ^{\circ}'],[0.02 -1.25])
    
    if doPrint; cn.print([ 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_best'],'path',savePath); end
end
%% Visualize
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass,B0);      
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B,v,n_loc);

%% Print figures
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_dir_aaaaaaa.gif'],'DelayTime',0.01,'LoopCount',inf);
imwrite(gif_stuff_v.im,gif_stuff_v.map,[savePath 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_v_.gif'],'DelayTime',0.01,'LoopCount',inf);


%% Add the event and save data to TimeTable
comment = 'nordita';
toSave = 'single';
tool.add_event

%% 
em_tint = toepoch([2008 04 22 17 37 12.2;2008 04 22 17 37 12.8])';
h = irf_plot({irf_abs(gsmB1),irf_abs(gsmVi1),gsmE1,[gsmE1(:,1) ang1],diE1,irf_abs(diB1),diE1par});
h(1).YLabel.String = 'B GSM [nT)';
h(2).YLabel.String = 'Vi GSM [km/s)';
h(3).YLabel.String = 'E GSM [mV/m)';
h(4).YLabel.String = '\theta_{B>SP} [deg]';
h(5).YLabel.String = 'E ISR2 [mV/m)';
h(6).YLabel.String = 'B ISR2 [nT)';
h(7).YLabel.String = 'E_{ExB=0} ISR2 [mV/m)';
for kk = 1:6; irf_legend(h(kk),{'x','y','z'},[0.99 0.95]); end
irf_zoom(h,'x',em_tint)
irf_zoom(h,'y')

%% Check projection of B on spin plane
em_tint = toepoch([2008 04 22 17 37 12.2;2008 04 22 17 37 12.8])';
tt = mean(em_tint);
Bii = irf_resamp(diB1,tt,'nesrest'); Bii = Bii(2:4)
para = irf_norm(Bii);
perp = irf_norm(cross([0 0 1],para))
proj = cross(perp,[0 0 1])
acosd(perp*proj')
acosd(perp*para')
acosd(proj*para')

%% Plot magnetic field directions before and after crossing part.
t_msp = toepoch([2008 04 22 17 35 30])';
t_sh = toepoch([2008 04 22 17 37 30])';
t_fr = toepoch([2008 04 22 17 37 12.60])';
B_msp = irf_resamp(gsmB1,t_msp,'nearest'); B_msp = B_msp(2:4)
B_sh = irf_resamp(gsmB1,t_sh,'nearest'); B_sh = B_sh(2:4)
B_fr = irf_resamp(gsmB1,t_fr,'nearest'); B_fr = B_fr(2:4)
B = [B_msp; B_fr; B_sh];
%quiver3(-B(:,1),-B(:,2),-B(:,2),B(:,1),B(:,2),B(:,3))
quiver3(-B(:,1)*0,-B(:,2)*0,-B(:,2)*0,B(:,1),B(:,2),B(:,3))
view([1 0 0]); xlabel('x');ylabel('y');zlabel('z')
text(B(1,1),B(1,2),B(1,3),'MSP')
text(B(2,1),B(2,2),B(2,3),'FL')
text(B(3,1),B(3,2),B(3,3),'SH')
%%
quiver3(-B(:,1)*0,-B(:,2)*0,-B(:,2)*0,B(:,1),B(:,2),B(:,3))
quiver3(-B(:,1)*0,-B(:,2)*0,-B(:,2)*0,B(:,1),B(:,2),B(:,3))
quiver3(-B(:,1)*0,-B(:,2)*0,-B(:,2)*0,B(:,1),B(:,2),B(:,3))
%% Plot with STAFF
z=irf_norm(B_fr); % B/z direction, tries*3
y=irf_norm(irf_cross(irf_cross(z,direction),z)); % perp1
x=irf_norm(irf_cross(y,z)); % perp2

sysB3 = [gsmB3(:,1) gsmB3(:,2:4)*[x' y' z']];
%sysVi3 = [gsmVi3(:,1) gsmVi3(:,2:4)*[x' y' z']];
sysE3 = [gsmE3(:,1) gsmE3(:,2:4)*[x' y' z']];
fl2_tint = toepoch([2008 04 22 17 37 12.34;2008 04 22 17 37 12.54])';
fl1_tint = toepoch([2008 04 22 17 37 12.54;2008 04 22 17 37 12.75])';


h=irf_plot(3,'newfigure');
irf_plot({irf_tlim(irf_filt(sysB3,0,1,5),fl2_tint),irf_tlim(phi_B,fl2_tint)})
h=irf_plot(3,'newfigure');
irf_plot({irf_tlim(irf_filt(sysB3,0,1,5),fl1_tint),irf_tlim(phi_B,fl1_tint)})

%% Estimate current
deltaB = 10e-9; % T
deltaxk = 15e3; % m
mu0 = 4*pi*1e-7;
Je = deltaB/deltaxk/mu0; % A
Je_muA = Je*1e6
Je_nA = Je*1e9
nJ = 33*1e6; % m^-3
e = 1.6e-19;

veJ = -Je/(nJ*e)

%%
c_eval('P?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',sclist);
%%
em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])';
[mvaB1td,l,v]=irf_minvar(irf_tlim(gsmB1,em_tint),'td');
[mvaB1n0,l,v]=irf_minvar(irf_tlim(gsmB1,em_tint),'<Bn>=0');
[mvaB1uc,l,v]=irf_minvar(irf_tlim(gsmB1,em_tint));

%%
irf_plot({mvaB1td,...
          mvaB1n0,...
          mvaB1uc,...
          irf_tlim(gsmB1,mvaB1td([1 end],1)),...
          irf_filt(irf_tlim(gsmB1,mvaB1td([1 end],1)),0.1,0,450,5),...
          irf_filt(irf_tlim(gsmB1,mvaB1td([1 end],1)),1,0,450,5),...
          irf_filt(irf_tlim(gsmB1,mvaB1td([1 end],1)),2,0,450,5),...
          irf_filt(irf_tlim(gsmB1,mvaB1td([1 end],1)),3,0,450,5),...
          irf_tlim(irf_filt(gsmB1,3,0,450,5),mvaB1td([1 end],1))})

%%
vBz = [Bz(:,1) Bz(:,2)*z(1,:)];
bdiff = irf_add(1,irf_tlim(gsmB1,vBz([1 end],1)),-1,vBz);
h=irf_plot({irf_abs(gsmB1),Bz,irf_abs(bdiff)});
em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])';
irf_zoom(h,'x',em_tint)
irf_zoom(h,'y')
%%
h=irf_plot({gsmB1,bdiff},'comp')
irf_legend(h(1),{'B1_{GSM}','B1_{GSM}-\delta B_z'},[0.02 0.95])
irf_legend(h(2),{'B1_{GSM}','B1_{GSM}-\delta B_z'},[0.02 0.95])
irf_legend(h(3),{'B1_{GSM}','B1_{GSM}-\delta B_z'},[0.02 0.95])
ylabel(h(1),'B_x');ylabel(h(2),'B_y');ylabel(h(3),'B_z')
irf_zoom(h,'x',em_tint)
irf_zoom(h,'y')
%%
h = irf_plot({irf_abs(gsmB1),gsmVi1,gsmE1,flh1,re1,vte1,peaNe1});
ylabel(h(1),'B [nT]');
ylabel(h(2),'V_i [km/s]');
ylabel(h(3),'E [mV/m]');
ylabel(h(4),'f_{LH} [Hz]');
ylabel(h(5),'\rho_{e} [km]');
ylabel(h(6),'v_{te} [km/s]');
ylabel(h(7),'n_{e} [cm^{-3}]');
em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])';
%irf_zoom(h,'x',em_tint)
title(h(1),'C1')
irf_zoom(h,'y')
%%
h = irf_plot(7);
hca = irf_panel('B'); irf_plot(hca,gsmB1); ylabel(hca,'B [nT]');
hca = irf_panel('Vi'); irf_plot(hca,gsmVi1); ylabel(hca,'V_i [km/s]');
hca = irf_panel('Te'); irf_plot(hca,parTe1); ylabel(hca,'T_e [eV]');
hca = irf_panel('ne'); irf_plot(hca,peaNe1); ylabel(hca,'n_e [cm^{-3}]');
hca = irf_panel('partPr'); irf_plot(hca,[partPr1(:,1) partPr1(:,2)*1e9]); ylabel(hca,'P_{e} [nPa]');
hca = irf_panel('magPr'); irf_plot(hca,[magPr1(:,1) magPr1(:,2)*1e9]); ylabel(hca,'P_{B} [nPa]');
hca = irf_panel('beta'); irf_plot(hca,beta1); ylabel(hca,'\beta_e'); hca.YScale = 'log'; hca.YTick = 10.^[-3:1:2]; hca.YLim = [2e-3 2e1];

em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])';
title(h(1),'C1')
irf_zoom(h,'y')
irf_zoom(h,'x',toepoch([2008 04 22 17 10 00 ;2008 04 22 18 30 00])')
%% Intgrate Faradays law to compare the amplitudes of B
% divE=-dB/dt
% units: mV/m*km/s*s
%farB = irf_integrate([dEn(:,1) dEn(:,1+i_dir)*velocity]);
farB = [dEn(:,1) dEn(:,1+i_dir)/velocity];
% Plot
h=irf_plot(5);

hca=irf_panel('E_n'); 
irf_plot(hca,dEn(:,[1 i_dir+1]))
ylabel(hca,'E_n [mV/m]')

hca=irf_panel('E_k'); 
irf_plot(hca,dEk(:,[1 i_dir+1]))
ylabel(hca,'E_k [mV/m]')

hca=irf_panel('B_z'); 
irf_plot(hca,Bz)
ylabel(hca,'\delta B_z [nT]')

hca=irf_panel('\phi'); 
irf_plot(hca,{phi_E(:,[1 1+i_v]),phi_B},'comp')
ylabel(hca,'\phi [V]')

hca=irf_panel('Faraday B'); 
irf_plot(hca,farB)
ylabel(hca,'B')

irf_zoom(h,'x',em_tint)

%% Momentum balance
% after running tool.run_single_event
% mn*dv/dt = T1 + T2 + T3 (T2 = -grad(P))
T1 = units.e*n_loc*1e6*y(i_dir,:)*avEn*1e-3;
T3 = cross(Vi_loc*1e3,B0*z(i_dir,:)*1e-9)*units.e*n_loc*1e6;

T2 = -T1 -T3

T1norm = irf_norm(T1)
T2norm = irf_norm(T2)
T3norm = irf_norm(T3)
T1abs = irf_abs(T1,1)
T2abs = irf_abs(T2,1)
T3abs = irf_abs(T3,1)

%% Ohm's law terms

c_eval('gsmVixB?= irf_tappl(irf_cross(gsmVi?,gsmB?),''*-1e-3'');',[1 3])

h=irf_plot(8);

hca=irf_panel('B'); 
irf_plot(hca,gsmB1)
ylabel(hca,'B [nT]')

hca=irf_panel('Vi'); 
irf_plot(hca,gsmVi1)
ylabel(hca,'V_i [km/s]')

VidotB = irf_dot(irf_norm(gsmB1),irf_norm(gsmVi1));
angViB = [VidotB(:,1) acosd(VidotB(:,2))];
hca=irf_panel('angle Vi B'); 
irf_plot(hca,angViB)
ylabel(hca,'\theta_{V_i>B} [\circ]')

hca=irf_panel('E'); 
irf_plot(hca,gsmE1)
ylabel(hca,'E [mV/m]')

hca=irf_panel('VixB'); 
irf_plot(hca,gsmVixB1)
ylabel(hca,'-V_i\times B [mV/m]')

hca=irf_panel('VixB and E'); 
irf_plot(hca,irf_filt(gsmE1,0,1,450,5)); hold(hca,'on');
irf_plot(hca,gsmVixB1); hold(hca,'off');
irf_legend(hca,{'E_x','E_y','E_z','(-V_i\times B)_x','(-V_i\times B)_y','(-V_i\times B)_z'},[0.02 0.95])
ylabel(hca,'E [mV/m]')

hca=irf_panel('n'); 
irf_plot(hca,{irf_tappl(scpNe1,'*0.25'),peaNe1,whiNe1},'comp')
irf_legend(hca,{'EFW','PEA','WHI'},[0.02 0.95])
ylabel(hca,'N_e [cm^{-3}]')

hca=irf_panel('VixB + E'); 
LHS = irf_abs(irf_add(1,irf_filt(gsmE1,0,1,450,5),-1,gsmVixB1));
irf_plot(hca,LHS)
irf_legend(hca,{'x','y','z'},[0.02 0.95])
ylabel(hca,'E+V_i\times B [mV/m]')


irf_zoom(h,'y')
irf_zoom(h,'x',toepoch([2008 04 22 17 10 00 ;2008 04 22 18 30 00])')
%%
hca=irf_panel('E_n'); 
irf_plot(hca,dEn(:,[1 i_dir+1]))
ylabel(hca,'E_n [mV/m]')

hca=irf_panel('E_k'); 
irf_plot(hca,dEk(:,[1 i_dir+1]))
ylabel(hca,'E_k [mV/m]')


%%
h=irf_plot(4);
irf_plot(h(1),{irf_tappl(scpNe1,'*0.25'),peaNe1,whiNe1},'comp')
irf_plot(h(2),{irf_tappl(scpNe2,'*0.25'),peaNe2,whiNe2},'comp')
irf_plot(h(3),{irf_tappl(scpNe3,'*0.25'),peaNe3,whiNe3},'comp')
irf_plot(h(4),{irf_tappl(scpNe4,'*0.25'),peaNe4,whiNe4},'comp')
ylabel(h(1),'C1')
ylabel(h(2),'C2')
ylabel(h(3),'C3')
ylabel(h(4),'C4')
irf_legend(h(1),{'EFW','PEA','WHI'},[0.02 0.95])
irf_legend(h(2),{'EFW','PEA','WHI'},[0.02 0.95])
irf_legend(h(3),{'EFW','PEA','WHI'},[0.02 0.95])
irf_legend(h(4),{'EFW','PEA','WHI'},[0.02 0.95])
irf_zoom(h,'x',tint)

for ii=1:4;
    set(h(ii),'yscale','log','ylim',[0.5 150],'ytick',[1 10 100])
end

%% Calculate j and grad
c_eval('gradR? = irf_resamp(gsmR?,scpNe?);',sclist);
c_eval('gradR? = [gradR?(:,1:4) scpNe?(:,2)];',sclist);
gradPhi = c_4_gradphi(gradR1,gradR2,gradR3,gradR4);

C3norm = '*0.25';
C4norm = '*0.34';

scpNediff34 = irf_add(1,irf_tappl(scpNe4,C4norm),-1,irf_tappl(scpNe3,C3norm));
peaNediff34 = irf_add(1,peaNe4,-1,peaNe3);
whiNediff34 = irf_add(1,whiNe4,-1,whiNe3);

h=irf_plot(4);
irf_plot(h(1),{irf_tappl(scpNe3,C3norm),peaNe3,whiNe3},'comp'); 
irf_plot(h(2),{irf_tappl(scpNe4,C4norm),peaNe4,whiNe4},'comp')
irf_plot(h(3),{irf_tappl(scpNe3,C3norm),irf_tappl(scpNe4,C4norm)},'comp')
irf_plot(h(4),{scpNediff34,peaNediff34,whiNediff34},'comp')

ylabel(h(1),'N_e [cm^{-3}] (C3)')
ylabel(h(2),'N_e [cm^{-3}] (C4)')
ylabel(h(3),'N_e [cm^{-3}] (EFW)')
ylabel(h(4),'\Delta N_e [cm^{-3}]')

irf_legend(h(1),{'EFW','PEA','WHI'},[0.02 0.95])
irf_legend(h(2),{'EFW','PEA','WHI'},[0.02 0.95])
irf_legend(h(3),{'C3','C4'},[0.02 0.95])
irf_legend(h(4),{'EFW','PEA','WHI'},[0.02 0.95])

%irf_zoom(h,'x',tint_zoom)
irf_zoom(h,'x',tint)
irf_zoom(h,'y')

%%
tint = toepoch([2008 04 22 17 16 45;2008 04 22 17 17 15])';
[minvarE,l,v] = irf_minvar(irf_tlim(irf_filt(irf_tlim(gsmE1,tint+[-10 10]),0,0.5,450,5),tint));

figure;
irf_plot(minvarE)
title(['x=[' num2str(v(1,:),'%.2f') '],  y=[' num2str(v(1,:),'%.2f') '],  z=[' num2str(v(1,:),'%.2f') ']'])
irf_legend({'x','y','z'},[0.02 0.95])
ylabel('E [mV/m]')

