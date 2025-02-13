% A routine that identifies strong magnetic nulls (potential EDR) from brst l2pre data.
% Written by Elin Eriksson
%
%Important: This uses brst data so tetrahedron quality is available and is
%used, thus this is only compatible with data on server /data/mms currently.
%
% Loads data and finds nulls that has a BoxLim of 70km and a Current
% Lim of 500E-9 [nA/m^2 if B is in nT and R is in km].
%Plots B, J (curlometer method), E, JxB electric field, and J.E as an average
%of all spacecraft values and a seperate plot with the same variables where only the E and
%B are from a specific spacecraft given by the number of ic.

%To-Do: Fix data gap so that the function looks for more data after the
%first gap.
%% Set time interval to look in, the specific spacecraft of interest and set limits for Null method.
% taken from
ic=1; %Gives number of spacecraft where density is taken for Hall field calculations.
Tint  = irf.tint('2015-10-16T10:33:40Z/2015-10-16T10:33:55Z');
boxLim=70;
currentLim=500E-9;
smallInterval=false; %Option to make the intervals just around the nulls
avgPlot=false; %Makes the plot with average E- and B-fields
%% Load magnetic field and spacecraft positional data
%mms.db_init('local_file_db','/data/mms'); %Point towards the folder with the mms data

% Magnetic Field
disp('Loading Magnetic fields');
c_eval('B?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_gse'',Tint);',1:4);
%c_eval('B?=removerepeatpnts(B?)',[1:4]); %Removes data points with too small a time difference between them which buggs the resample functions
%Resamples according to the time line where all spacecraft has data
timeStart=max([B1.time.start.epochUnix B2.time.start.epochUnix B3.time.start.epochUnix B4.time.start.epochUnix],[],2);
timeStart=EpochUnix(timeStart);
timeStop=min([B1.time.stop.epochUnix B2.time.stop.epochUnix B3.time.stop.epochUnix B4.time.stop.epochUnix],[],2);
timeStop=EpochUnix(timeStop);
timeLog=B1.time >= timeStart & B1.time <=timeStop;
newTime=B1.time(timeLog,:);
c_eval('B? = B?.resample(newTime);',1:4);

% Spacecraft Position
disp('Loading Spacecraft Position');
R  = mms.get_data('R_gse',Tint);%Cailbrated position
if length(R.gseR1(1,:))==4 || length(R.gseR2(1,:))==4 || length(R.gseR3(1,:))==4 || length(R.gseR4(1,:))==4
    c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?(:,1:3));',1:4);
else
c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?);',1:4);
end
clear R
% Checks if there is any position data missing in that case go to predicted
% spacecraft position
if isempty(R1) || isempty(R2) || isempty(R3) || isempty(R4)
    disp('Loading predicted spacecraft position');
    c_eval('R?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_ql_pos_gse'',Tint);',1:4) %Predicted s/c position
end
%If the spacecraft data is still missing for a spacecraft give error
if isempty(R1) || isempty(R2) || isempty(R3) || isempty(R4)
    error('Missing spacecraft position from at least one spacecraft');
else
    c_eval('R? = R?.resample(B1);',1:4);
end
%% Looks for Nulls and if they are found makes average calculations and Loads electric fields, density and calculates J and JxB with curlometer method
%Quality data comes 2 days late
% Load quality of tetrahedron
quality=mms.db_get_variable('mms_ancillary_defq','quality',Tint); %Every 30s
if isempty(quality.quality)
        error('No tetrahedron quality available');
        
    else
        quality=irf.ts_scalar(EpochTT(quality.time),quality.quality);
      
        quality=quality.resample(B1);
        tetrahedronBad= quality.data < 0.7;
        % Removes all time steps with bad tetrahedron quality
        c_eval('R?_null = [R?.time.epochUnix double(R?.data)];',1:4);
        c_eval('B?_null = [B?.time.epochUnix double(B?.data)];',1:4);
        
        c_eval('R?_null(tetrahedronBad,:)=[];',1:4);
        c_eval('B?_null(tetrahedronBad,:)=[];',1:4);
        
        if isempty(B1_null) || isempty(B2_null) || isempty(B3_null) || isempty(B4_null)
            error('Tetrahedron quality is not good enough to search for Nulls');
        else
            disp('Looking for Nulls');
            Nulls=c_4_null(R1_null,R2_null,R3_null,R4_null,B1_null,B2_null,B3_null,B4_null,'boxLim',boxLim,'strong',currentLim);
        end
end
if isempty(Nulls.t)
    time_int=irf_time(Tint,'tint>utc');
    errormessage=['No Nulls were found with strong currents between ',time_int];
    error(errormessage);
else
    %% Average calculations. Loads electric fields, density and calculates J and JxB with curlometer method
    
    disp('Calculates average Bfield')
    Bav = (B1.data+B2.data+B3.data+B4.data)/4;
    Bav=[B1.time.epochUnix double(Bav)];
    disp('Loads Electric fields')
    % Electric field
    c_eval('E?_orig=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',Tint);',1:4);
    timeStartE=max([E1_orig.time.start.epochUnix E2_orig.time.start.epochUnix E3_orig.time.start.epochUnix E4_orig.time.start.epochUnix],[],2);
    timeStartE=EpochUnix(timeStartE);
    timeStopE=min([E1_orig.time.stop.epochUnix E2_orig.time.stop.epochUnix E3_orig.time.stop.epochUnix E4_orig.time.stop.epochUnix],[],2);
    timeStopE=EpochUnix(timeStopE);
    timeLogE=E1_orig.time >= timeStartE & E1_orig.time <=timeStopE;
    newTimeE=E1_orig.time(timeLogE,:);
    c_eval('E? = E?_orig.resample(newTimeE);',1:4);
    %Removes 1.5 mV/m offset from Ex for all spacecraft
    c_eval('E?.data(:,1) = E?.data(:,1)-1.5;',1:4);
    c_eval('E? = E?.resample(B1);',1:4); %Resampling to B-field due to the index searching later to make smaller filesizes.
    
    % Average E-field for EdotJ
    disp('Calculates average E-fields')
    Eav = (E1.data+E2.data+E3.data+E4.data)/4;
    Eav=[E1.time.epochUnix double(Eav)];
    
    %Special satellite
    c_eval('Efield=E?;',ic);
    Efield=[Efield.time.epochUnix double(Efield.data)];
    c_eval('Bfield=B?;',ic);
    Bfield=[Bfield.time.epochUnix double(Bfield.data)];
    
    disp(['Loads density from spacecraft',num2str(ic)])
    %Density from spacecraft ic for Hall field calculation
    c_eval('ni=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_numberDensity'',Tint);',ic);
    ni = TSeries(ni.time,ni.data,'to',1);
    ni = ni.resample(B1);
    ni=[ni.time.epochUnix double(ni.data)];
  
    disp('Calculates current')
    % Assuming GSE and DMPA are the same coordinate system calculates j and jXB.
    [j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
    divovercurl = divB;
    divovercurl.data = abs(divovercurl.data)./j.abs.data;
    divovercurl=[B1.time.epochUnix double(divovercurl.data)];
    
    j=[j.time.epochUnix double(j.data)];
    j(:,2:4) = j(:,2:4).*1e9;
    jxB=[jxB.time.epochUnix double(jxB.data)];
    jxB(:,2:4) = jxB(:,2:4)./[ni(:,2) ni(:,2) ni(:,2)];
    jxB(:,2:4) = jxB(:,2:4)./1.6e-19./1000; %Convert to (mV/m)
    jxB((irf_abs(jxB(:,2:4),1) > 100),1:4) = NaN; % Remove some questionable fields
    disp('Calculates EdotJ')
    % Calculates EdotJ
    EdotJ = dot(Eav(:,2:4),j(:,2:4),2)./1000; %J (nA/m^2), E (mV/m), E.J (nW/m^3)
    EdotJ = [Eav(:,1) EdotJ];
end

%% Only picks out a smaller time interval around the nulls (to keep filesize small)
%Possibility to pick out interesting time interval to keep the file size small
if smallInterval
    if isempty(Nulls.t)
        error('No Nulls are found')
    else
        index=B1.time(:,1).epochUnix>(Nulls.t(1,1)-1) & B1.time(:,1).epochUnix<Nulls.t(end,1)+1;
        Efield=Efield(index,:);
        Bfield=Bfield(index,:);
        Bav=Bav(index,:);
        Eav=Eav(index,:);
        ni=ni(index,:);
        j=j(index,:);
        jxB=jxB(index,:);
        EdotJ=EdotJ(index,:);
        divovercurl=divovercurl(index,:);
    end
end
%% Plot data should only be used if nulls are found. Plots are only focused around the intervals of nulls found because of above segment.
if avPlot
if isempty(Nulls.t)
    error('Could not plot because no nulls have been found')
else
    h = irf_plot(6,'newfigure');
    
    hca = irf_panel('BMMSav');
    irf_plot(hca,Bav);
    ylabel(hca,{'B_{DMPA}','(nT)'},'Interpreter','tex');
    irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
    %irf_legend(hca,'(a)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('J');
    irf_plot(hca,j);
    ylabel(hca,{'J_{DMPA}','(nA m^{-2})'},'Interpreter','tex');
    irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
    %irf_legend(hca,'(c)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('divovercurl');
    irf_plot(hca,divovercurl);
    ylabel(hca,{'$\frac{|\nabla \cdot \mathbf{B}|} {|\nabla \times \mathbf{B}|}$'},'Interpreter','latex','Fontsize',18);
    %irf_legend(hca,'(e)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('EMMSav');
    irf_plot(hca,Eav);
    ylabel(hca,{'E_{DSL}','(mV m^{-1})'},'Interpreter','tex');
    irf_legend(hca,{'E_{x}','E_{y}','E_{z}'},[0.88 0.10])
    %irf_legend(hca,'(b)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('jxB');
    irf_plot(hca,jxB);
    ylabel(hca,{'J \times B/n_{e} q_{e}','(mV m^{-1})'},'Interpreter','tex');
    %irf_legend(hca,'(f)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('jdotE');
    irf_plot(hca,EdotJ);
    ylabel(hca,{'E . J','(nW m^{-3})'},'Interpreter','tex');
    %irf_legend(hca,'(g)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    title(h(1),strcat('MMS averaged - Current density and fields'));
    tmarks=Nulls.t; %Makes tmarks for all null data points.
    irf_plot_axis_align(1,h);
    irf_pl_mark(h,tmarks,[0.8 0.8 0.8]) %Setting the color lines to grey
    irf_pl_number_subplots(h,[0.99, 0.95]);
    tint_zoom=[Bav(1,1) Bav(end,1)];
    irf_zoom(h,'x',tint_zoom);
end
end
%% Plot data for specific spacecraft given by ic.
if isempty(Nulls.t)
    error('Could not plot because no nulls have been found')
else
    h = irf_plot(6,'newfigure');
    
    hca = irf_panel('BMMS');
    irf_plot(hca,Bfield);
    ylabel(hca,{'B_{DMPA}','(nT)'},'Interpreter','tex');
    irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
    %irf_legend(hca,'(a)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('J');
    irf_plot(hca,j);
    ylabel(hca,{'J_{DMPA}','(nA m^{-2})'},'Interpreter','tex');
    irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
    %irf_legend(hca,'(c)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('divovercurl');
    irf_plot(hca,divovercurl);
    ylabel(hca,{'$\frac{|\nabla \cdot \mathbf{B}|} {|\nabla \times \mathbf{B}|}$'},'Interpreter','latex','Fontsize',18);
    %irf_legend(hca,'(e)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('EMMS');
    irf_plot(hca,Efield);
    ylabel(hca,{'E_{DSL}','(mV m^{-1})'},'Interpreter','tex');
    irf_legend(hca,{'E_{x}','E_{y}','E_{z}'},[0.88 0.10])
    %irf_legend(hca,'(b)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('jxB');
    irf_plot(hca,jxB);
    ylabel(hca,{'J \times B/n_{e} q_{e}','(mV m^{-1})'},'Interpreter','tex');
    %irf_legend(hca,'(f)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    hca = irf_panel('jdotE');
    irf_plot(hca,EdotJ);
    ylabel(hca,{'E . J','(nW m^{-3})'},'Interpreter','tex');
    %irf_legend(hca,'(g)',[0.99 0.98],'color','k')
    grid(hca,'off');
    
    title(h(1),strcat(['MMS', num2str(ic)]));
    tmarks=Nulls.t; %Makes tmarks for all null data points.
    irf_plot_axis_align(1,h);
    irf_pl_mark(h,tmarks,[0.8 0.8 0.8]) %Setting the color lines to grey
    irf_pl_number_subplots(h,[0.99, 0.95]);
    tint_zoom=[Bfield(1,1) Bfield(end,1)];
    irf_zoom(h,'x',tint_zoom);
end

% 2) from terminal convert to eps file without white margins
% > epstool --copy --bbox Jan11pic10.eps Jan11pic10_crop.eps
% 3) convert eps file to pdf, result is in Current_crop.pdf
% > ps2pdf -dEPSFitPage -dEPSCrop -dAutoRotatePages=/None Jan11pic10_crop.eps

