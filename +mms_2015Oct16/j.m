% Load B
tint_J = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');
%c_eval('[Bxyz?,dobj_Bxyz?] = cn_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint_J(1));',sc)
%c_eval('Bxyz? = Bxyz?.resample(Bxyz1);',1:4);

R  = mms.get_data('R_gse',tint_J);
c_eval('Rxyz? = TSeries(R.time,R.gseR?,''to'',1);',1:4);
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',1:4);
R0 = (Rxyz1.data+Rxyz2.data+Rxyz3.data+Rxyz4.data)/4;
c_eval('relR? = Rxyz?-R0;',sc)
% Assuming GSE and DMPA are the same coordinate system.
[j,divB,B,jxB,divTshear,divPb] = c_4_j('Rxyz?','Bxyz?');

if 1
    %%
    % Assuming GSE and DMPA are the same coordinate system.
    [j,divB,B,jxB,divTshear,divPb] = c_4_j('gseR?','dfg?');
    
    Bav = (dfg1.data+dfg2.data+dfg3.data+dfg4.data)/4;
    Bav = TSeries(dfg1.time,Bav,'to',1);

    divovercurl = divB;
    divovercurl.data = abs(divovercurl.data)./j.abs.data;

    % Transform current density into field-aligned coordinates
    SCpos = [0 1 0];

    Bmag = dfg1.abs.data;
    Rpar = dfg1.data./[Bmag Bmag Bmag];
    Rperpy = irf_cross(Rpar,SCpos);
    Rmag   = irf_abs(Rperpy,1);
    Rperpy = Rperpy./[Rmag Rmag Rmag];
    Rperpx = irf_cross(Rperpy, Rpar);
    Rmag   = irf_abs(Rperpx,1);
    Rperpx = Rperpx./[Rmag Rmag Rmag];

    jpar = dot(Rpar,j.data,2);
    jperp = dot(Rperpx,j.data,2);
    jperp2 = dot(Rperpy,j.data,2);

    jfac = TSeries(j.time,[jperp jperp2 jpar],'to',1);

end


%% Make plot
h = irf_plot(5,'newfigure');

hca = irf_panel('B1');
irf_plot(hca,Bxyz1);
ylabel(hca,{'B_{DMPA}','[nT]'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
irf_legend(hca,'(a)',[0.99 0.98],'color','k')

hca = irf_panel('R1');
irf_plot(hca,gseR1/units.RE*1e3);
ylabel(hca,{'R_{GSE}','[RE]'},'Interpreter','tex');
irf_legend(hca,{'X','Y','X'},[0.88 0.10])
irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca = irf_panel('J');
%j.data = j.data*1e9;
irf_plot(hca,j*1e9);
ylabel(hca,{'J_{DMPA}','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
irf_legend(hca,'(c)',[0.99 0.98],'color','k')

hca = irf_panel('Jfac');
%jfac.data = jfac.data*1e9;
irf_plot(hca,jfac*1e9);
ylabel(hca,{'J_{FAC}','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'J_{\perp 1}','J_{\perp 2}','J_{||}'},[0.88 0.10])
irf_legend(hca,'(d)',[0.99 0.98],'color','k')

hca = irf_panel('divovercurl');
irf_plot(hca,divovercurl);
ylabel(hca,{'|\nabla . B|','|\nabla \times B|'},'Interpreter','tex');
irf_legend(hca,'(e)',[0.99 0.98],'color','k')
set(hca,'yscale','lin');

irf_zoom(h,'x',tint_J)
irf_zoom(h,'y')
%%
hca = irf_panel('jxB');
jxB.data = jxB.data./[ne.data ne.data ne.data]; 
jxB.data = jxB.data/1.6e-19/1000; %Convert to (mV/m)
jxB.data(find(abs(jxB.data) > 100)) = NaN; % Remove some questionable fields
irf_plot(hca,jxB);
ylabel(hca,{'J \times B/n_{e} q_{e}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,'(f)',[0.99 0.98],'color','k')

%%
irf_zoom(h,'x',tint_J)
irf_zoom(h,'y')

% Roy Torberts requested intervals
tint1 = irf.tint('2015-10-16T13:00:00',8*60);
tint2 = irf.tint('2015-10-16T11:25:00',60*7);
irf_pl_mark(h,tint1.epochUnix')
irf_pl_mark(h,tint2.epochUnix')
