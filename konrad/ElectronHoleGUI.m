function [] = customGUIv3(varargin)
%% NEW CHANGES:
% - Now works for 2-spacecraft manually and automatically.
%% TO BE ADDED:
% - Make it so one can load a large Tint and then select smaller intervals
% with zooming on the main plot, e.g. E = Ewide.tlim(whatever goes
% here), and a button to return to E.
% - Check box for plasma parameters.
% - Ion velocity to enable ion frame.
% - Check box for plot L_perp.
% - Output spacecraft separations in all dimensions.
% - Calculate L_perp.
% - Make it 3 spacecraft. (the difficulty lies in fixing peak_filter and
%   matching_EH codes...)

%% Main Window  %     'position',[1000 1000 1400 1000],...
S.fh = figure('units','pixels',...
    'position',[1 1 1280 720],...%get(groot,'ScreenSize'),...
    'menubar','none',...
    'numbertitle','off',...
    'name','customGUI',...
    'resize','on');

% Graphical items...

%% Filter options
S.filter_panel = uipanel('Title','Filter options','FontSize',12,...
    'units','normalized','pos',[0.675 0.7 0.15 0.3]);

% Filter input boxes:
delta = 0.05/(1.5);
alpha = 0.035/(1.5);
w = 0.2/(1.5);   
S.filter_ed(1) = uicontrol(S.filter_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.5 1-w-alpha 0.4 w],...
    'string','30');
S.filter_ed(2) = uicontrol(S.filter_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.5 1-2*w-alpha-delta 0.4 w],...
    'string','0');
S.filter_ed(3) = uicontrol(S.filter_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.5 1-3*w-alpha-2*delta 0.4 w],...
    'string','');
S.filter_ed(4) = uicontrol(S.filter_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.5 1-4*w-alpha-3*delta 0.4 w],...
    'string','5');
%Filter input boxes text:

S.filter_tx(1) = uicontrol(S.filter_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-w-alpha 0.4 w],...
                    'string',sprintf('fmin'),...
                    'fontsize',12);   
S.filter_tx(2) = uicontrol(S.filter_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-2*w-delta-alpha 0.4 w],...
                    'string',sprintf('fmax'),...
                    'fontsize',12);   
S.filter_tx(3) = uicontrol(S.filter_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-3*w-2*delta-alpha 0.4 w],...
                    'string',sprintf('Sample freq.'),...
                    'fontsize',12);   
S.filter_tx(4) = uicontrol(S.filter_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-4*w-3*delta-alpha 0.4 w],...
                    'string',sprintf('Order'),...
                    'fontsize',12);   

% Apply filter button:

S.apply_filter = uicontrol(S.filter_panel,'style','push',...
                  'unit','normalized',...
                  'position',[0.25,0.05,0.5,0.2],...
                  'string','Apply filter',...
                  'call',{@Apply_filter,S});
%% Result options:
S.result_panel = uipanel('Title','Result options','FontSize',12,...
    'units','normalized','pos',[0.675 0.4 0.15 0.3]);              

S.plot_check = uicontrol(S.result_panel,'style','check',...
                 'unit','normalized',...
                 'position',[0.1 0.8 0.6 0.1],...
                 'string','Plot correlation',...
                 'fontsize',12);
             
S.result_txt(1) = uicontrol(S.result_panel,'style','text',...
                    'unit','normalized',...
                    'position',[-0.04 0.3 0.7 0.1],...
                    'string',sprintf('Output variable name:'),...
                    'fontsize',12);   

S.result_txt(2) = uicontrol(S.result_panel,'style','text',...
                    'unit','normalized',...
                    'position',[-0.04 0.6 0.7 0.1],...
                    'string',sprintf('\x3D5 int. range (\x0B1tpp)'),...
                    'fontsize',12);   
                
S.result_option(1) = uicontrol(S.result_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.65 0.275 0.35 0.15],...
    'string','EHProps');

S.result_option(2) = uicontrol(S.result_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.65 0.6 0.35 0.15],...
    'string','1.5');

%% Peak Finding Method
S.bg1 = uibuttongroup('Title','Peak finding method','units','normalized',...
    'pos',[0.825 0.9 0.15 0.1]);
S.rd(1) = uicontrol(S.bg1,...
    'style','rad',...
    'unit','pix',...
    'position',[20 35 70 20],...
    'string','Automatic');
S.rd(2) = uicontrol(S.bg1,...
    'style','rad',...
    'unit','pix',...
    'position',[20 10 70 20],...
    'string','Manual');

%% Manual Peak Finding Options
S.manual_input_panel = uipanel('Title','Manual peak finding options','visible','off','FontSize',12,...
    'units','normalized','pos',[0.825 0.7 0.15 0.2]);

S.manual_tx(1) = uicontrol(S.manual_input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 0.5 0.6 0.2],...
                    'string',sprintf('Window size (index)'),...
                    'fontsize',12);
S.manual_ed(1) = uicontrol(S.manual_input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 0.5 0.3 0.2],...
    'string','10');                
                
%% Automatic Peak Finding Options
S.input_panel = uipanel('Title','Automatic peak finding options','FontSize',12,...
    'units','normalized','pos',[0.825 0.2 0.15 0.7]);

delta = 0.015;
alpha = 0.025;
w = 0.09;     
S.tx(1) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-w-alpha 0.55 w],...
                    'string',sprintf('\x3C3-crit [\x3C3(E)]'),...
                    'fontsize',12);   
S.tx(2) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-2*w-delta-alpha 0.55 w],...
                    'string',sprintf('min(E\x22A5 / E\x2223\x2223)'),...
                    'fontsize',12);   
S.tx(3) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-3*w-2*delta-alpha 0.55 w],...
                    'string',sprintf('min(|E|)'),...
                    'fontsize',12);   
S.tx(4) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-4*w-3*delta-alpha 0.55 w],...
                    'string',sprintf('Spacecraft delay (ms)'),...
                    'fontsize',12);   
S.tx(5) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-5*w-4*delta-alpha 0.55 w],...
                    'string',sprintf('Minimum criteria duration (index)'),...
                    'fontsize',12);                   
S.tx(6) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-6*w-5*delta-alpha 0.55 w],...
                    'string',sprintf('min(|E\x22A5|) (mV/m)'),...
                    'fontsize',12); 
S.tx(7) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-7*w-6*delta-alpha 0.55 w],...
                    'string',sprintf('max(tpp) (s)'),...
                    'fontsize',12);
S.tx(8) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-8*w-7*delta-alpha 0.55 w],...
                    'string',sprintf('First peak sign (+/-1)'),...
                    'fontsize',12);
S.tx(9) = uicontrol(S.input_panel,'style','text',...
                    'unit','normalized',...
                    'position',[0 1-9*w-8*delta-alpha 0.55 w],...
                    'string',sprintf('Spacecraft (e.g. 13, 1234)'),...
                    'fontsize',12);
delta = 0.05;
alpha = 0.025;
w = (1-9*delta)/10;    
S.ed(1) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-w-alpha 0.3 w],...
    'string','1.5');
S.ed(2) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-2*w-alpha-delta 0.3 w],...
    'string','0');
S.ed(3) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-3*w-alpha-2*delta 0.3 w],...
    'string','0.3');
S.ed(4) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-4*w-alpha-3*delta 0.3 w],...
    'string','15');
S.ed(5) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-5*w-alpha-4*delta 0.3 w],...
    'string','2');
S.ed(6) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-6*w-alpha-5*delta 0.3 w],...
    'string','0');
S.ed(7) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-7*w-alpha-6*delta 0.3 w],...
    'string','0.005');
S.ed(8) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-8*w-alpha-7*delta 0.3 w],...
    'string','1');
S.ed(9) = uicontrol(S.input_panel,'style','edit',...
    'unit','normalized',...
    'position',[0.6 1-9*w-alpha-8*delta 0.3 w],...
    'string','1234');

%% Find peaks button
S.select_PEAK = uicontrol('style','push',...
                  'unit','normalized',...
                  'position',[0.2,0.05,0.2,0.1],...
                  'string','Select Peaks',...
                  'call',{@select_PEAKS_call,S});

%% Analyse EHs
S.calc_props = uicontrol('style','push',...
                  'unit','normalized',...
                  'position',[0.4,0.05,0.2,0.1],...
                  'string','Analyse holes',...
                  'call',{@calc_EH_props,S});

%% Plot
ysize = 0.16;
y0pos = 0.36;
x0pos = 0.05;
S.ax1 = axes('units','normalized',...
    'position',[x0pos y0pos+3*ysize 0.6, ysize],...
    'fontsize',8,...
    'nextplot','replacechildren');
S.ax2 = axes('units','normalized',...
    'position',[x0pos y0pos+2*ysize 0.6, ysize],...
    'fontsize',8,...
    'nextplot','replacechildren');
S.ax3 = axes('units','normalized',...
    'position',[x0pos y0pos+ysize 0.6, ysize],...
    'fontsize',8,...
    'nextplot','replacechildren');
S.ax4 = axes('units','normalized',...
    'position',[x0pos y0pos 0.6, ysize],...
    'fontsize',8,...
    'nextplot','replacechildren');

%% GUI Callbacks and other functions..

set(S.fh,'toolbar','figure');
set(S.fh,'menubar','figure');
set(S.rd(1),'callback',{@automatic_peak_check,S})  % Set callback.
set(S.rd(2),'callback',{@manual_peak_check,S})  % Set callback.
% set(S.plot_check,'callback',{@plot_check,S})  % Set callback.

%% Startup stuff

% Setting up default options.
global run_data
ic=1:4;
if isa(varargin{1},'char') && ~isempty(strfind(varargin{1},'?'))
    varString = varargin{1};
    c_eval('run_data.Exyz?=evalin(''base'',irf_ssub(varString,?));');
    temp_tint = [run_data.Exyz1.time(1),run_data.Exyz1.time(end)];
    run_data.Tint=temp_tint;
    if length(varargin) == 2
        if isa(varargin{2},'char') && ~isempty(strfind(varargin{2},'?'))
            varString = varargin{2};
            c_eval('R?=evalin(''base'',irf_ssub(varString,?));');
            c_eval('Bgse?=mms.get_data(''B_gse_fgm_brst_l2'',temp_tint+[-1,1],?);',ic);
            
            c_eval('Rr?=R?.resample(Bgse1);',ic);
            c_eval('Rfac?=irf_convert_fac(Rr?,Bgse1,[1 0 0]);',ic);
            
            
            c_eval('dRfac?!=Rfac?-Rfac!;',1,2:4);
            c_eval('dRfac?!=Rfac?-Rfac!;',2,3:4);
            c_eval('dRfac?!=Rfac?-Rfac!;',3,4);
            
            r_1=TSeries(dRfac12.time,zeros(length(dRfac12.time),3));
            c_eval('r_?=TSeries(dRfac!?.time,[-dRfac!?.data(:,1),-dRfac!?.data(:,2),-dRfac!?.data(:,3)]);',[2,3,4],1);
            c_eval('r_resamp?=TSeries(r_?.time,r_?.data(:,3));',1:4);
            c_eval('r_resamp?=r_resamp?.resample(run_data.Exyz?);',1:4);
            c_eval('r_resamp_all?=TSeries(r_?.time,r_?.data);',1:4);
            c_eval('r_resamp_all?=r_resamp_all?.resample(Efac?);',1:4);
            c_eval('run_data.sc_pos?=r_resamp_all?;',1:4);
        else
            disp('Incorrect input. \n Should either be: \n')
            disp('gui(''E?'',''R?'') where E? and R? are TSeries, or: \n')
            disp('gui(''E?'')');
        end
    else %Download position
        temp_tint = [run_data.Exyz1.time(1),run_data.Exyz1.time(end)];
        run_data.Tint=temp_tint;
        Tintl = temp_tint+[-60 60];
        R = mms.get_data('R_gse',Tintl);
        c_eval('R? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);
        c_eval('Bgse?=mms.get_data(''B_gse_fgm_brst_l2'',temp_tint+[-1,1],?);',ic);
                
        c_eval('Rr?=R?.resample(Bgse1);',ic);
        c_eval('Rfac?=irf_convert_fac(Rr?,Bgse1,[1 0 0]);',ic);
        
        
        c_eval('dRfac?!=Rfac?-Rfac!;',1,2:4);
        c_eval('dRfac?!=Rfac?-Rfac!;',2,3:4);
        c_eval('dRfac?!=Rfac?-Rfac!;',3,4);
        
        r_1=TSeries(dRfac12.time,zeros(length(dRfac12.time),3));
        c_eval('r_?=TSeries(dRfac!?.time,[-dRfac!?.data(:,1),-dRfac!?.data(:,2),-dRfac!?.data(:,3)]);',[2,3,4],1);
        c_eval('r_resamp?=TSeries(r_?.time,r_?.data(:,3));',1:4);
        c_eval('r_resamp?=r_resamp?.resample(run_data.Exyz?);',1:4);
        c_eval('r_resamp_all?=TSeries(r_?.time,r_?.data);',1:4);
        c_eval('r_resamp_all?=r_resamp_all?.resample(run_data.Exyz?);',1:4);
        c_eval('run_data.sc_pos?=r_resamp_all?;',1:4);
        
    end
    fmin=30;
    fmax=0;
    order=5;
    c_eval('Bxyz?=mms.get_data(''B_dmpa_fgm_brst_l2'',run_data.Tint,?);',ic);
    c_eval('run_data.Bxyz?=Bxyz?;',ic);
    c_eval('run_data.Exyz_filt?=run_data.Exyz?.filt(fmin,fmax,[],order);',ic);
    c_eval('run_data.Exyz_nofilt? = irf_convert_fac(run_data.Exyz?,Bxyz?,[1 0 0]);',ic);
    c_eval('run_data.Exyz_filt? = irf_convert_fac(run_data.Exyz_filt?,Bxyz?,[1 0 0]);',ic);

elseif isa(varargin{1},'EpochTT')
    Tint = varargin{1};
    run_data.Tint=Tint;
    c_eval('Bxyz?=mms.get_data(''B_dmpa_fgm_brst_l2'',Tint,?);',ic);
    c_eval('Exyz?=mms.get_data(''E_dsl_edp_brst_l2'',Tint,?);',ic);
    fmin=30;
    fmax=0;
    order=5;
    c_eval('Exyz_filt?=Exyz?.filt(fmin,fmax,[],order);',ic);
    c_eval('run_data.Bxyz? = Bxyz?;',ic);
    c_eval('run_data.Exyz_nofilt? = irf_convert_fac(Exyz?,Bxyz?,[1 0 0]);',ic);
    c_eval('run_data.Exyz_filt? = irf_convert_fac(Exyz_filt?,Bxyz?,[1 0 0]);',ic);
   
    
    temp_tint = [Exyz1.time(1),Exyz1.time(end)];
    Tintl = temp_tint+[-60 60];
    R = mms.get_data('R_gse',Tintl);
    c_eval('R? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);
    c_eval('Bgse?=mms.get_data(''B_gse_fgm_brst_l2'',temp_tint+[-1,1],?);',ic);
    
    c_eval('Rr?=R?.resample(Bgse1);',ic);
    c_eval('Rfac?=irf_convert_fac(Rr?,Bgse1,[1 0 0]);',ic);
    
    c_eval('dRfac?!=Rfac?-Rfac!;',1,2:4);
    c_eval('dRfac?!=Rfac?-Rfac!;',2,3:4);
    c_eval('dRfac?!=Rfac?-Rfac!;',3,4);
    
    r_1=TSeries(dRfac12.time,zeros(length(dRfac12.time),3));
    c_eval('r_?=TSeries(dRfac!?.time,[-dRfac!?.data(:,1),-dRfac!?.data(:,2),-dRfac!?.data(:,3)]);',[2,3,4],1);
    c_eval('r_resamp?=TSeries(r_?.time,r_?.data(:,3));',1:4);
    c_eval('r_resamp?=r_resamp?.resample(Exyz?);',1:4);
    c_eval('r_resamp_all?=TSeries(r_?.time,r_?.data);',1:4);
    c_eval('r_resamp_all?=r_resamp_all?.resample(Exyz?);',1:4);
    c_eval('run_data.sc_pos?=r_resamp_all?;',1:4);
    
    
end
%Initial properites:
run_data.peak_method = 'Automatic';
run_data.corr_plot=0;

S.combined_axes = [S.ax1, S.ax2, S.ax3, S.ax4];
c_eval('irf_plot(S.combined_axes(?),run_data.Exyz_nofilt?.z,''k'');',ic);
c_eval('ylabel(S.combined_axes(?),''E_{||} (mV/m)'',''interpreter'',''tex'');',ic);
c_eval('irf_zoom(S.combined_axes(?),''x'',run_data.Tint);',ic);
c_eval('irf_legend(S.combined_axes(?),''mms?'',[0.98 0.12],''fontsize'',12);',ic);
c_eval('box(S.combined_axes(?),''on'');',ic);




function [] = load_auto_options_call(varargin)
% Callback for pushbutton.
global run_data
run_data.peak_method = 'Automatic';
S=varargin{3};

function [] = automatic_peak_check(varargin)
global run_data
% Callback for pushbutton.
S=varargin{3};
run_data.peak_method = 'Automatic';
set(S.input_panel,'visible','on')
set(S.manual_input_panel,'visible','off');

function [] = manual_peak_check(varargin)
global run_data
% Callback for pushbutton.
S=varargin{3};
run_data.peak_method = 'Manual';
set(S.input_panel,'visible','off')
set(S.manual_input_panel,'visible','on')

function [] = Apply_filter(varargin)
S=varargin{3};
global run_data;
ic=1:4;
plot_start_nr = length(S.fh.Children);
filter_options = str2double(get(S.filter_ed(1:4),'string'));  % The numbers to operate on.
fmin = filter_options(1);
fmax = filter_options(2);
sample_freq = filter_options(3);
if isnan(sample_freq)
    sample_freq=[];
end
order=filter_options(4);
if fmin==0 && fmax == 0
    c_eval('Exyz_filt?=run_data.Exyz_nofilt?;',ic);
else
c_eval('run_data.Exyz_filt?=run_data.Exyz_nofilt?.filt(fmin,fmax,sample_freq,order);',ic);
end
c_eval('cla(S.fh.Children(plot_start_nr+1-?));',ic);
c_eval('irf_plot(S.fh.Children(plot_start_nr+1-?),run_data.Exyz_filt?.z,''k'');',ic);
c_eval('irf_zoom(S.fh.Children(plot_start_nr+1-?),''y'');',ic);
c_eval('ylabel(S.fh.Children(plot_start_nr+1-?),''E_{||} (mV/m)'',''interpreter'',''tex'');',ic);
c_eval('irf_legend(S.fh.Children(plot_start_nr+1-?),''mms?'',[0.98 0.12],''fontsize'',12);',ic);

function [] = select_PEAKS_call(varargin)
S=varargin{3};
global run_data;
plot_start_nr = length(S.fh.Children);
switch run_data.peak_method
    case 'Manual'
        manual_options = str2double(get(S.manual_ed(1),'string'));
        run_data.manual_sens = manual_options(1);
        manual_peaks;
    case 'Automatic'
        auto_options = str2double(get(S.ed(1:8),'string'));  % The numbers to operate on.
        spacecraft_string = get(S.ed(9),'string');
        jSC=1;
        for iSC=1:4
            if contains(spacecraft_string,num2str(iSC))
               selected_spacecraft(jSC)=iSC;
               jSC=jSC+1;
            end
        end
        run_data.use_spacecraft = selected_spacecraft;
        run_data.sigma_crit = auto_options(1);
        run_data.E_perp_ratio = auto_options(2);
        run_data.E_crit = auto_options(3);
        run_data.spacecraft_delay = auto_options(4);
        run_data.minimum_duration = auto_options(5);
        run_data.E_perp_crit = auto_options(6);
        run_data.tpp_crit = auto_options(7);
        if auto_options(8) >=0
            run_data.peak_order = [1,-1];
        else
            run_data.peak_order = [-1,1];
        end
        automatic_peaks;
end

% Plotta peaks här
for iSC = 1:4
    c_eval('Edata = run_data.Exyz_filt?.data(:,3);',iSC);
    c_eval('peaks{iSC} = run_data.peaks?.data;',iSC);
    if ~isempty(peaks{iSC})
        for iPeak = 1:length(peaks{iSC})
            c_eval('indx?(iPeak) = find(Edata == peaks{?}(iPeak));',iSC);
        end
    else
        c_eval('indx?=[];');
    end
end
colours=[1,0,0;1,0,0;0,1,0;0,1,0;0,0,1;0,0,1;0,1,1;0,1,1;1,0,1;1,0,1;0.4,0.3,0;0.4,0.3,0];
c_eval('hold(S.fh.Children(plot_start_nr+1-?),''on'');',1:4);
for iSC=1:4
    colournr=1;
    c_eval('peaks = run_data.peaks?;',iSC);
    if ~isempty(peaks)
        for ii=1:1:length(peaks)
            irf_plot(S.fh.Children(plot_start_nr+1-iSC),peaks(ii),'o','color',colours(colournr,:));
            colournr=colournr+1;
            if colournr==13
                colournr=1;
            end
        end
    end
end
c_eval('ylabel(S.fh.Children(plot_start_nr+1-?),''E_{||} (mV/m)'',''interpreter'',''tex'');',1:4);
% Gör allt jox så att det funkar med 2,3, eller 4 sc... fyfan.. Manual
% borde jag inte behöva göra just nått med, och jag har redan automatic för
% 2, så behöver bara göra en för 3.

function manual_peaks()
global run_data;
% selection_sens = run_data.
ic=1:4;
c_eval('Efac?=run_data.Exyz_filt?;',1:4);
manual_sens = run_data.manual_sens;
dt = median(Efac1.time(2:end)-Efac1.time(1:end-1));
c_eval('PEAK_TIMES?=int64([]);',ic);
c_eval('PEAK_DATA?=[];',ic);
currkey=0;
colours=[1,0,0;1,0,0;0,1,0;0,1,0;0,0,1;0,0,1;0,1,1;0,1,1;1,0,1;1,0,1;0.4,0.3,0;0.4,0.3,0];
colour_nr=1;
Tint=EpochTT([Efac1.time.epoch(1);Efac1.time.epoch(end)]);
h=irf_plot(4,'newfigure');
scrsz = get(groot,'ScreenSize');
set(h(1).Parent,'Position',scrsz);
c_eval('irf_plot(h(?),Efac?.z);',ic);
c_eval('ylabel(h(?),'' '')',ic);
irf_zoom(h(1:4),'x',Tint);
deltaT = 0.005; %Seconds forgiveness interval in peak selection.


for ii=1:4
%     clear PEAK_TIMES PEAK_DATA plot_points
    PEAK_TIMES = int64([]);
    PEAK_DATA = [];
    jj=1;
    currkey=0;
    colour_nr=1;
    hold(h(ii),'on');
    click_str = sprintf('Click near peaks in PANEL %d, SPACE to undo, press ENTER when finished.',ii);
    title(h(1),click_str)
    c_eval('use_E = Efac?;',ii);
    while currkey~=1
        [x1,y1,button]=ginput(1);
        
        if ~isempty(x1)
            if button==1
                
                guess_time = Efac1.time(1) + x1;
                if manual_sens==0
                    [~,closest_point] = min(abs(use_E.time-guess_time));
                    E_short = TSeries(use_E.time(closest_point),use_E.data(closest_point,:));
                else
                E_short = use_E.tlim(guess_time+[-manual_sens*dt, manual_sens*dt]);
                end
                guess_amp = y1;

                
                [~, closest_data]= min((E_short.data(:,3)-guess_amp).^2);
                
                
                data_sign = sign(E_short.data(closest_data,3));
                data_vector = closest_data+[(-20:1:20)];
                data_vector(data_vector>length(E_short.data(:,3)))=[]; %Constrain to inside interval
                data_vector(data_vector<1)=[];
                switch data_sign
                    case 1
                        [extremeval,extreme_index] = max(E_short.data(data_vector,3));
                    case -1
                        [extremeval,extreme_index] = min(E_short.data(data_vector,3));
                end
                PEAK_TIMES(jj)=E_short.time.epoch(data_vector(extreme_index));
                PEAK_DATA(jj) = E_short.data(data_vector(extreme_index),3);
                output_str = sprintf('Nr of selected peaks: %d',length(PEAK_DATA));
                title(h(1),{click_str,output_str});
                
                plot_point{jj} = TSeries(EpochTT(PEAK_TIMES(jj)),PEAK_DATA(jj));
                
                irf_plot(h(ii),plot_point{jj},'o','color',colours(colour_nr,:))
                colour_nr = (colour_nr+1)*(colour_nr~=12)+1*(colour_nr==12);
                forward = 1;
            else
                forward = 0;
                if jj>1
                    PEAK_TIMES(jj-1)=[];
                    PEAK_DATA(jj-1)=[];
                    output_str = sprintf('Nr of selected peaks: %d',length(PEAK_DATA));
                    title(h(1),{click_str,output_str});
                    irf_plot(h(ii),plot_point{jj-1},'o','color',[0.5,0.5,0.5])
                    colour_nr = (colour_nr~=1)*(colour_nr-1)+(colour_nr==1)*12;
                else
                    title(h(1),'All selections removed...')
                    colour_nr = 1;
                    jj=1;
                end
            end
            if forward
                jj=jj+1;
            else
                if jj>1
                    jj=jj-1;
                end
            end

        end
        if jj==2
            output_str = sprintf('Nr of selected peaks: %d',length(PEAK_DATA));
            title(h(1),{click_str,output_str})
        end
        
        currkey=get(gcf,'CurrentKey');
        if strcmp(currkey, 'return')
            currkey=1;
        else
            currkey=0;
        end
        output_str = sprintf('Nr of selected peaks: %d',length(PEAK_DATA));
%         disp(output_str);
    end
    c_eval('PEAK_TIMES?=PEAK_TIMES;',ii);
    c_eval('PEAK_DATA?=PEAK_DATA;',ii);
    if ii<4
        title(h(1),'Complete, press SPACE to move on to next panel...');
        pause; % wait for a keypress
    else
        title(h(1),'Complete!');
    end
end
[sorted_times1, order1]=sort(PEAK_TIMES1);
[sorted_times2, order2]=sort(PEAK_TIMES2);
[sorted_times3, order3]=sort(PEAK_TIMES3);
[sorted_times4, order4]=sort(PEAK_TIMES4);

sorted_data1 = PEAK_DATA1(order1);
sorted_data2 = PEAK_DATA2(order2);
sorted_data3 = PEAK_DATA3(order3);
sorted_data4 = PEAK_DATA4(order4);

manual_peaksTS1 = TSeries(EpochTT(sorted_times1'),sorted_data1');
manual_peaksTS2 = TSeries(EpochTT(sorted_times2'),sorted_data2');
manual_peaksTS3 = TSeries(EpochTT(sorted_times3'),sorted_data3');
manual_peaksTS4 = TSeries(EpochTT(sorted_times4'),sorted_data4');
logical_sc = [1*~isempty(manual_peaksTS1),2*~isempty(manual_peaksTS2),...
    3*~isempty(manual_peaksTS3),4*~isempty(manual_peaksTS4)];
logical_sc(logical_sc==0)=[];
run_data.use_spacecraft=logical_sc;
c_eval('run_data.peaks?=manual_peaksTS?;',ic);
close(h(1).Parent);


function automatic_peaks()
global run_data;
ic=run_data.use_spacecraft;
sigma_factor=run_data.sigma_crit;
E_crit = run_data.E_crit;
min_dur = run_data.minimum_duration; 
spacecraft_delay = run_data.spacecraft_delay*10^6;
perp_par_ratio = run_data.E_perp_ratio;
E_perp_crit = run_data.E_perp_crit; %Fixa som input...
dTcrit = run_data.tpp_crit; %Fixa som input...
Epeak_order=run_data.peak_order;
c_eval('Efac?=run_data.Exyz_filt?;',ic);

%%Standard deviation calculations...
std_range=200;
std_stepsize=50;

%Calculate approximate values of the standard deviation to be used later.
for iSC=1:length(ic)
    sc=ic(iSC);
    c_eval('vec_L=length(Efac?.data(:,3));',sc);
    c_eval('glob_std_crit?=0;',sc);%std(Efac?.data(:,3));',sc_num);
    c_eval('glob_stdTS?=TSeries([Efac?.time.start;Efac?.time.stop],[glob_std_crit?;glob_std_crit?]);',sc);
    for ii=1:std_stepsize:vec_L
        
        if ii<=std_range
            c_eval('std_crit?(ii,1)=std(Efac?.data(1:ii+std_range,3));',sc);
            
        elseif ii>std_range && ii+std_range<=vec_L
            c_eval('std_crit?(ii,1)=std(Efac?.data(ii-std_range:ii+std_range,3));',sc);
            
        elseif ii+std_range>vec_L
            c_eval('std_crit?(ii,1)=std(Efac?.data(ii-std_range:vec_L,3));',sc);
        end
    end
    c_eval('if length(std_crit?)<vec_L, std_crit?(end+1:vec_L,1)=0; end',sc);
    c_eval('interp_points=find(std_crit?==0);',sc);
    c_eval('data_points=find(std_crit?~=0);',sc);
    c_eval('data_values=std_crit?(data_points);',sc);
    c_eval('interp_values=interp1(data_points,data_values,interp_points);',sc);
    c_eval('all_std_points=[data_points;interp_points]; all_std_values=[data_values;interp_values];',sc);
    c_eval('[sorted_points,sorting]=sort(all_std_points);',sc);
    c_eval('sorted_values=all_std_values(sorting,:);',sc);
    
    c_eval('interp_std_crit?=sorted_values;',sc);
end

%First peak selection:
% Peak-finding algorithm.
c_eval('crit_matrix1{?}=((abs(Efac?.data(:,3)))<(sigma_factor*interp_std_crit?));',ic);
c_eval('crit_matrix2{?}=((abs(Efac?.data(:,3)))<E_crit);',ic);


%Check that peaks are actually detected:
std_flag=false(1,length(ic));
E_crit_flag=false(1,length(ic));
for ii=1:length(ic)
    if sum(crit_matrix1{ic(ii)}(20:end-20))==length(crit_matrix1{ic(ii)}(20:end-20))
       std_flag(ii)=true;
       str1=sprintf('No peaks fulfilling \x3C3-crit for MMS%d. Consider decreasing \x3C3-crit.',ic(ii));
       
    end
    if sum(crit_matrix2{ic(ii)})==length(crit_matrix2{ic(ii)})
       E_crit_flag(ii)=true; 
       str2=sprintf('No peaks larger than min(|E|) for MMS%d. Consider decreasing min(|E|).',ic(ii));
      
    end
    if std_flag(ii)==true
    
        if E_crit_flag(ii)==true
            error(sprintf([str1,'\n',str2]))
            
        end
        error(str1)
    elseif E_crit_flag(ii)==true
        error(sprintf(str2))        
    end
    
end



%We want to find 0s in the resulting matrix.
c_eval('criteria_matrix{?}=crit_matrix1{?}+crit_matrix2{?};',ic);
c_eval('criteria_matrix{?}(criteria_matrix{?}==2)=1;',ic);
empty_peaks_flag=0;
for JJ=1:length(ic)
    filtered_peaksT=int64([]);
    filtered_peaksD=[];
    c_eval('filtered_peaksTS?=TSeries(EpochTT(filtered_peaksT),filtered_peaksD);',ic(JJ));
    dsig = diff([1 (criteria_matrix{ic(JJ)})' 1]);
    startIndex = find(dsig < 0);
    endIndex = find(dsig > 0)-1;
    duration = endIndex-startIndex+1;
    stringIndex = (duration >= min_dur);
    startIndex = startIndex(stringIndex);
    endIndex = endIndex(stringIndex);
    if~isempty(startIndex)
        for ii=1:length(startIndex)
            temp_vector{ii}=startIndex(ii):endIndex(ii);
            c_eval('[peaks?(ii),peak_indices?(ii)]=max(Efac?.data(temp_vector{ii},3));',ic(JJ));
            c_eval('if peaks?(ii)<0, peaks?(ii)=nan; end',ic(JJ));
            c_eval('[minpeaks?(ii),minpeak_indices?(ii)]=min(Efac?.data(temp_vector{ii},3));',ic(JJ));
            c_eval('if minpeaks?(ii)>0, minpeaks?(ii)=nan; end',ic(JJ));
            
            c_eval('peak_times?(ii)=(Efac?.time.epoch(temp_vector{ii}(peak_indices?(ii))));',ic(JJ));
            c_eval('minpeak_times?(ii)=(Efac?.time.epoch(temp_vector{ii}(minpeak_indices?(ii))));',ic(JJ));
        end
        c_eval('combPeaks?=[peaks?'';minpeaks?''];',ic(JJ));
        c_eval('combTimes?=[peak_times?'';minpeak_times?''];',ic(JJ));
        c_eval('combTimes?(isnan(combPeaks?))=[];',ic(JJ));
        c_eval('combPeaks?(isnan(combPeaks?))=[];',ic(JJ));
        c_eval('[sortedTimes,TimeOrder]=sort(combTimes?);',ic(JJ));
        c_eval('sortedPeaks=combPeaks?(TimeOrder,:);',ic(JJ));
        c_eval('empty_peaks_flag=isempty(sortedPeaks);');
        if ~empty_peaks_flag
            c_eval('peaksTS?=TSeries(EpochTT(sortedTimes),sortedPeaks);',ic(JJ));
            length_loop=3;
            kk=1;
            ii=1;
            while ii<length_loop
                if ii==1
                    c_eval('length_loop=length(peaksTS?.data);',ic(JJ));
                end
                if length_loop>1
                    c_eval('dT=peaksTS?.time(ii+1)-peaksTS?.time(ii);',ic(JJ));
                    c_eval('sgn(1)=sign(peaksTS?.data(ii)); sgn(2)=sign(peaksTS?.data(ii+1));',ic(JJ));
                    c_eval('p2p=abs(peaksTS?.data(ii+1)-peaksTS?.data(ii));',ic(JJ));
                    c_eval('minpeak=min(abs([peaksTS?.data(ii+1),peaksTS?.data(ii)]));',ic(JJ));
                    if (dT<dTcrit*2*minpeak/p2p) && ((sgn(1)==Epeak_order(1) && sgn(2)==Epeak_order(2)))
                        % Peaks in a pair must have opposite signs, and
                        % pairs with similar field strength is allowed a
                        % larger peak to peak time. (This essentially
                        % removes noise-peak matches.
                        c_eval('filtered_peaksT(kk)=peaksTS?.time.epoch(ii);',ic(JJ));
                        c_eval('filtered_peaksD(kk)=peaksTS?.data(ii);',ic(JJ));
                        c_eval('filtered_peaksT(kk+1)=peaksTS?.time.epoch(ii+1);',ic(JJ));
                        c_eval('filtered_peaksD(kk+1)=peaksTS?.data(ii+1);',ic(JJ));
                        kk=kk+2;
                        ii=ii+2;
                    else
                        ii=ii+1;
                        
                    end
                end
            end
        end
        c_eval('filtered_peaksTS?=TSeries(EpochTT(filtered_peaksT)'',filtered_peaksD'');',ic(JJ));
        c_eval('if isempty(filtered_peaksTS?), empty_peaks_flag=1; end;',ic(JJ));
        
    else
        empty_peaks_flag=1;
    end
end

%Filters to select only the multi-spacecraft electron holes. This should be
%generalized to the two-spacecraft case. Could just be done using "if" and
%different functions, but yeah...

switch length(ic)
    case 2
        c_eval('filtered_peaks1 = filtered_peaksTS?;',ic(1));
        c_eval('filtered_peaks2 = filtered_peaksTS?;',ic(2));
        [filtered2_peaksTS1,filtered2_peaksTS2,~]=peak_filter_2sc(filtered_peaks1,filtered_peaks2,[0,0,0,0],spacecraft_delay);
        empty_peaks_flag=isempty(filtered2_peaksTS1); 
        c_eval('r_par1 = TSeries(run_data.sc_pos?.time(:,1),run_data.sc_pos?.data(:,3));',ic(1));
        c_eval('r_par2 = TSeries(run_data.sc_pos?.time(:,1),run_data.sc_pos?.data(:,3));',ic(2));
        
        [matchingESW1,matchingESW2]=label_EHs_2sc(filtered2_peaksTS1,filtered2_peaksTS2,r_par1,r_par2,spacecraft_delay); %#ok<ASGLU>
        c_eval('matchingESW_TS?=matchingESW1;',ic(1));
        c_eval('matchingESW_TS?=matchingESW2;',ic(2));
        
    case 4
        [filtered2_peaksTS1,filtered2_peaksTS2,filtered2_peaksTS3,filtered2_peaksTS4,~]=peak_filter_v2(filtered_peaksTS1,filtered_peaksTS2,filtered_peaksTS3,filtered_peaksTS4,[0,0,0,0],spacecraft_delay);
        empty_peaks_flag=isempty(filtered2_peaksTS1);
        
        c_eval('r_par? = TSeries(run_data.sc_pos?.time(:,1),run_data.sc_pos?.data(:,3));',ic);
        [matchingESW_TS1,matchingESW_TS2,matchingESW_TS3,matchingESW_TS4]=label_EHs(filtered2_peaksTS1,filtered2_peaksTS2,filtered2_peaksTS3,filtered2_peaksTS4,r_par1,r_par2,r_par3,r_par4,spacecraft_delay); %#ok<ASGLU>
        
end
% Calculate perpendicular field of EH. (Related to when multiple spacecraft
% fultill the distance criteria to be "forced" to be correct I think.
c_eval('EppIndx?=[]; for ii=1:length(matchingESW_TS?), [indices,~]=find(Efac?.time.epoch==matchingESW_TS?.time.epoch(ii)); EppIndx?(end+1:end+length(indices))=indices; end',ic);
c_eval('Eperp1?=Efac?.x;',ic);
c_eval('Eperp2?=Efac?.y;',ic);
c_eval('Efac_perp?=TSeries(Efac?.time,sqrt(Eperp1?.data.^2+Eperp2?.data.^2));',ic);

for jj=1:length(ic)
    c_eval('loop_length=length(matchingESW_TS?.time);',ic(jj));
    c_eval('used_peakTS=matchingESW_TS?;',ic(jj));
    if ~isempty(used_peakTS.time.epoch)
        
        for ii=1:loop_length/2
            c_eval('[~,E_perpidx?(ii)]=min(abs(Efac?.data(EppIndx?(2*ii-1):EppIndx?(2*ii),3)));',ic(jj))
            % Take the perpendicular field as close as possible to E_||=0;
            
            c_eval('Epp?(ii)=Efac?.data(EppIndx?(2*ii),3)-Efac?.data(EppIndx?(2*ii-1),3);',ic(jj));
            c_eval('E_perpidx_adj?(ii)=E_perpidx?(ii)+EppIndx?(2*ii-1)-1;',ic(jj));
            c_eval('E_p1?(ii)=Eperp1?.data(E_perpidx_adj?(ii));',ic(jj));
            c_eval('E_p2?(ii)=Eperp2?.data(E_perpidx_adj?(ii));',ic(jj));
            c_eval('E_pmag?(ii)=Efac_perp?.data(E_perpidx_adj?(ii));',ic(jj));
            c_eval('Tpp?(ii)=matchingESW_TS?.time(2*ii)-matchingESW_TS?.time(2*ii-1);',ic(jj));
            c_eval('Tpp?t(ii)=(matchingESW_TS?.time.epoch(2*ii)+matchingESW_TS?.time.epoch(2*ii-1))/2;',ic(jj));
            c_eval('Epp?(ii)=Efac?.data(EppIndx?(2*ii),3)-Efac?.data(EppIndx?(2*ii-1),3);',ic(jj));
            
            
        end
        c_eval('EperpTS?=TSeries(EpochTT(Tpp?t)'',[(E_p1?)'', (E_p2?)'',E_pmag?'']);',ic(jj));
        c_eval('TppTS?=TSeries(EpochTT(Tpp?t)'',Tpp?'');',ic(jj));
        c_eval('EppTS?=TSeries(EpochTT(Tpp?t)'',abs(Epp?)'');',ic(jj));
        c_eval('EperpTS?=TSeries(EpochTT(Tpp?t)'',[(E_p1?)'', (E_p2?)'',E_pmag?'']);',ic(jj));
    end
end

% Filter to find ESWs with good perp fields.

c_eval('p2pfield?=EppTS?.data;',ic);
c_eval('good_peaks?=zeros(size(EppTS!.data));',1:4,ic(1));
c_eval('perpfield?=EperpTS?.data(:,3);',ic);
c_eval('good_peaks?=((perpfield?./p2pfield?)>perp_par_ratio).*(perpfield?>E_perp_crit);',ic);
c_eval('good_peak_indx?=good_peaks?.*[1:length(good_peaks?)]'';',ic);
consistently_good=sum([good_peaks1,good_peaks2,good_peaks3,good_peaks4],2); %
selected_peaks=consistently_good==length(ic); 
% Above could perhaps be changed to for example 3. It is the number of
% spacecraft fulfilling the perpendicual field constraints. But depending
% on where the spacecraft pass the EH they might measure a very weak
% perpendicular field. So, maybe this should just be set to 2 or something
% like that. Maybe even 1?

c_eval('selected_peak_indx=sort([2*good_peak_indx?(selected_peaks)-1;2*good_peak_indx?(selected_peaks)]);',ic(1));
c_eval('filtered3_peaksTS?=TSeries(EpochTT(int64([])),[]);',1:4);
c_eval('filtered3_peaksTS?=TSeries(matchingESW_TS?.time(selected_peak_indx), matchingESW_TS?.data(selected_peak_indx));',ic);
c_eval('if isempty(filtered3_peaksTS?), empty_peaks_flag=1; end',ic);
%Should make notes at these with more information if they end up being
%empty. i.e. where they become empty...

run_data.peaks1 = filtered3_peaksTS1;
run_data.peaks2 = filtered3_peaksTS2;
run_data.peaks3 = filtered3_peaksTS3;
run_data.peaks4 = filtered3_peaksTS4;

function []=calc_EH_props(varargin)
global run_data;
S=varargin{3};
%Load selected spacecraft and create combinations.
ic = run_data.use_spacecraft;
all_sc_flag = length(ic)==4;
sc_combs=nchoosek(ic,2);
[nrOfCombs,~]=size(sc_combs);
c_eval('nrOfEHs = length(run_data.peaks?)/2;',ic(1));
 %Flag for plotting x-corr results.
switch get(S.plot_check,'Value')
    case 1
        run_data.corr_plot=1;
    case 0
        run_data.corr_plot=0;
end
Phi_int_range=str2double(get(S.result_option(2),'string'));

corr_plot = run_data.corr_plot;
% Allocate with NaNs
tpp = NaN(nrOfEHs,4);
Eperp1=NaN(nrOfEHs,4);
Eperp2=NaN(nrOfEHs,4);
Epp=NaN(nrOfEHs,4);
EH_times=zeros(nrOfEHs,4,'int64');

%Download particle data.
Tint=run_data.Tint;
c_eval('Te?=mms.get_data(''Te_dbcs_fpi_brst_l2'',Tint,?);',ic);
c_eval('Ti?=mms.get_data(''Ti_dbcs_fpi_brst_l2'',Tint,?);',ic);
c_eval('Ne?=mms.get_data(''Ne_fpi_brst_l2'',Tint,?);',ic);
c_eval('Vi_dbcs_fpi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',Tint,?);',ic);

switch length(ic)
    case 2
    c_eval('Tii1 = Ti?;',ic(1));
    c_eval('Tii2 = Ti?;',ic(2));
        if isempty(Tii1) || isempty(Tii2)
            c_eval('Ti?=mms.get_data(''Ti_dbcs_fpi_brst_l2'',Tint+[-1,1],?);',ic);
            c_eval('Vi_dbcs_fpi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',Tint+[-1,1],?);',ic);
        end        
    case 4
        if isempty(Ti1) || isempty(Ti2) || isempty(Ti3) || isempty(Ti4)
            c_eval('Ti?=mms.get_data(''Ti_dbcs_fpi_brst_l2'',Tint+[-1,1],?);',ic);
            c_eval('Vi_dbcs_fpi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',Tint+[-1,1],?);',ic);
        end
        
end
c_eval('Teavg?=Te?.trace.data/3;',ic);
c_eval('Teavg?=TSeries(Te?.time,Teavg?);',ic);
c_eval('Tiavg?=Ti?.trace.data/3;',ic);
c_eval('Tiavg?=TSeries(Ti?.time,Tiavg?);',ic);
c_eval('Vi? = irf_convert_fac(Vi_dbcs_fpi?,run_data.Bxyz?,[1 0 0]);',ic);

Me=9.10938356e-31; %Electron mass (kg)
Mp=1.672621898e-27; %Proton mass (kg)
c=299792458; %Speed of light (m/s)
qe=1.6021766208e-19; %Electron charge (C)
epso=1/((4*pi*10^-7)*c^2); %Dielectric constant (F/m)

c_eval('Ne?=Ne?*10^6;',ic) %m^-3
c_eval('wpe?=sqrt(Ne?*qe^2/Me/epso);',ic); %rad/s
c_eval('Vte? = c*sqrt(1-1/(Teavg?*qe/(Me*c^2)+1).^2);',ic); % m/s (relativ. correct), particle with Vte has energy e*Te
c_eval('Ld? = 10^-3*(Vte?/wpe?/sqrt(2));',ic); %km
c_eval('TempeV?=Teavg?;',ic);
c_eval('Ne?=Ne?*10^-6;',ic);
c_eval('Te_ion?=Teavg?.resample(Tiavg?);',ic);
c_eval('c_is? = sqrt((Te_ion?*qe+3*Tiavg?*qe)/Mp)*10^-3;',ic); %km/s ion sound speed.
c_eval('B?=run_data.Bxyz?.abs*10^-9;',ic) %T
c_eval('wce?=qe*B?/Me;',ic); %rad/s



v_ion = NaN(nrOfEHs,4);
T_ion=NaN(nrOfEHs,4);
T_electron=NaN(nrOfEHs,4);
w_cyclotron=NaN(nrOfEHs,4);
w_plasma=NaN(nrOfEHs,4);
T_ion=NaN(nrOfEHs,4);
Bangle=NaN(nrOfEHs,4);
n_electron=NaN(nrOfEHs,4);
c_ion=NaN(nrOfEHs,4);
Ld=NaN(nrOfEHs,4);
B_field=cell(nrOfEHs,4);

% Step through each EH.
for iEH=1:nrOfEHs
    
    % Peak-to-peak time in s:
    c_eval('tpp(iEH,?) = double(run_data.peaks?.time.epoch(2*iEH)-run_data.peaks?.time.epoch(2*iEH-1))*10^-9;',ic);
    % Peak-to-peak E-field.
    c_eval('Epp(iEH,?) = run_data.peaks?.data(2*iEH)-run_data.peaks?.data(2*iEH-1);',ic);
    
    c_eval('E_tint? = [run_data.peaks?.time(2*iEH-1);run_data.peaks?.time(2*iEH)];',ic); %Time interval between peaks.
    c_eval('E_short? = run_data.Exyz_filt?.tlim(E_tint?);',ic); %Selection of E-field data between peaks.
    c_eval('[~,EH_middle_time(?)] = min(abs(E_short?.data(:,3)));',ic);
    % Evaluate the perpendicular field at z=0 (centre of EH crossing).
    c_eval('Eperp1(iEH,?) = E_short?.data(EH_middle_time(?),1);',ic);
    c_eval('Eperp2(iEH,?) = E_short?.data(EH_middle_time(?),2);',ic);
    c_eval('EH_times(iEH,?) = E_short?.time.epoch(EH_middle_time(?));',ic);
    
    %Find spacecraft position wrt mms 1 as they observe the EHs
    c_eval('[a?,pos_indx?]=min(abs(run_data.sc_pos?.time.epoch-EH_times(iEH,?)));',ic);
    c_eval('sc_pos?(iEH,:) = run_data.sc_pos?(pos_indx?).data;',ic);
    
    %Calculate 3D-velocity
    if all_sc_flag
        Rmat=[sc_pos2(iEH,:);sc_pos3(iEH,:);sc_pos4(iEH,:)];
        tvec = double([EH_times(iEH,2)-EH_times(iEH,1);EH_times(iEH,3)-EH_times(iEH,1);...
            EH_times(iEH,4)-EH_times(iEH,1)])*10^-9;
        nu_vec=(Rmat\tvec)';
        v_vec = nu_vec/(norm(nu_vec)^2);
        v_3d(iEH,:)=v_vec;
    end
    %Quantities evaluated at centre of EH: Plasma params (L_D, wpe,
    %wce, Te, Ti, ne, ni), E_perp, B, ..?
    c_eval('[~,Ti_index(?)]=min(abs(Tiavg?.time.epoch-EH_times(iEH,?)));',ic);
    c_eval('[~,Te_index(?)]=min(abs(Teavg?.time.epoch-EH_times(iEH,?)));',ic);
    c_eval('[~,B_index(?)]=min(abs(run_data.Bxyz?.time.epoch-EH_times(iEH,?)));',ic);
    
    c_eval('v_ion(iEH,?)=Vi?.data(Ti_index(?),3);',ic);
    c_eval('T_ion(iEH,?)=Tiavg?.data(Ti_index(?));',ic);
    c_eval('T_electron(iEH,?)=Teavg?.data(Te_index(?));',ic);
    c_eval('w_cyclotron(iEH,?)=wce?.data(B_index(?));',ic);
    c_eval('w_plasma(iEH,?)=wpe?.data(Te_index(?));',ic);
    c_eval('B_field{iEH,?}=run_data.Bxyz?.data(B_index(?),:);',ic);
    c_eval('Bangle(iEH,?)=atand(B_field{iEH,?}(3)/sqrt(B_field{iEH,?}(1)^2+B_field{iEH,?}(2)^2));',ic);
    c_eval('n_electron(iEH,?)=Ne?.data(Te_index(?));',ic);
    c_eval('c_ion(iEH,?)=c_is?.data(Ti_index(?));',ic);
    c_eval('Ld(iEH,?)=Ld?.data(Te_index(?));',ic);
end

% Get FPI data, Te, Ti, ne, ni, ve, vi etc.
% Calculate plasma freq, cyclo freq, debye lenght.
% peak-2-peak-time as well.

% For every EH, do x-corr.
for iEH=1:nrOfEHs
    
    % Do cross-correlation for every unique spacecraft combination
    for ic_loop=1:nrOfCombs
        iicc=sc_combs(ic_loop,:);
        iiiccc{iEH,ic_loop}=iicc; %Save combinations for later.
        %Find the "number" of the sc combination, e.g. [1,2] is number 1,
        %[1,3] is 2, ... [3,4] is 6.
        R_pair=find_scpair(iicc); 

        legends_label1{iEH,ic_loop} = sprintf('mms %01d',iicc(1)); %For plotting purposes.
        legends_label2{iEH,ic_loop} = sprintf('mms %01d',iicc(2));
        
        % Select the two E-field vectors to be correlated (E_||).
        c_eval('E_first=TSeries(run_data.Exyz_filt?.time, run_data.Exyz_filt?.data(:,3));',iicc(1));
        c_eval('E_second=TSeries(run_data.Exyz_filt?.time, run_data.Exyz_filt?.data(:,3));',iicc(2));
        
        %Create time interval between observations for later use:
        Tint=irf.tint(min(E_first.time.epoch(1),E_second.time.epoch(1)),max(E_first.time.epoch(end),E_second.time.epoch(end)));
        
        %Resample fields
        E_second=E_second.resample(E_first);
        
        %Create time interval for cross-correlation:
        T_PP = [tpp(iEH,iicc(1)),tpp(iEH,iicc(2))]; %peak2peak times for the two spacecraft
        T1 = EpochTT(EH_times(iEH,iicc(1))); %Time of centre of EH crossing for sc1
        T2 = EpochTT(EH_times(iEH,iicc(2))); % sc2.
        %Save time interval for nice plots:
        [sort_times,time_order] = sort([T1,T2]);
        plot_tint{iEH,ic_loop} = sort_times + [-2*T_PP(time_order(1)),2*T_PP(time_order(2))];
        
        %Generate vectors to be cross-correlated.
        eh1_times = E_first.time;
        eh1_data = E_first.data;
        eh1_tint = T1+[-1.4*double(T_PP(1)),1.4*double(T_PP(1))]; %MAKE THIS PROPERTY AN OPTION!!!!!
        [~,eh1_start]=min(abs((eh1_times) - (eh1_tint(1))));
        [~,eh1_stop]=min(abs((eh1_times) - (eh1_tint(2))));
        eh1_TSeries = TSeries(eh1_times,[zeros(eh1_start-1,1);eh1_data(eh1_start:eh1_stop);zeros(length(eh1_times)-eh1_stop,1)]);
        E_first = eh1_TSeries; %Replace with shortened vector.
        
        %Same for second E-field.
        eh2_times = E_second.time;
        eh2_data = E_second.data;
        eh2_tint = T2+[-1.4*T_PP(2),1.4*T_PP(2)];
        [~,eh2_start]=min(abs((eh2_times) - (eh2_tint(1))));
        [~,eh2_stop]=min(abs((eh2_times) - (eh2_tint(2))));
        eh2_TSeries = TSeries(eh2_times,[zeros(eh2_start-1,1);eh2_data(eh2_start:eh2_stop);zeros(length(eh2_times)-eh2_stop,1)]);
        E_second = eh2_TSeries;
        dt=E_second.time(2)-E_second.time(1);
        
        %Normalize the vectors for cross-correlation :: MAKE AN OPTION
        normE_first=(E_first.data-median(E_first.data))./(max(E_first.data));
        normE_second=(E_second.data-median(E_second.data))./max(E_second.data);
        
        [acor,lag]=xcorr(normE_first,normE_second,'coeff');
        
        corr=acor; %Full cross correlation.
        [~,I] = max(acor);
        lagDiff=lag(I);
        maxcorr=max(corr); %Correlation coefficient.
        
        % Velocity calculations
        lagTS=TSeries(T1,(lagDiff*dt)');
        
        %This calculates the velocity of the ESW assuming it travels along R_||
        [Rpar,~,~] = mms.mms4_displacement(Tint+[-10,10]);

        Rpar_diff = Rpar.data(:,R_pair);

        delta_Rpar = TSeries(Rpar.time,Rpar_diff);
        delta_Rpar= delta_Rpar.resample(E_first);
        
        dR_times = delta_Rpar.time.epoch;
        dR_values = delta_Rpar.data;
        
        dR=TSeries(EpochTT(dR_times),dR_values);
        dR2=dR.resample(lagTS);
        % %
        vel = TSeries(lagTS.time,(dR2.data)./(lagTS.data)); %km/s;


        if corr_plot
            v1_plot{iEH,ic_loop} = TSeries(E_first.time,normE_first);
            v3_plot{iEH,ic_loop} = TSeries(E_second.time,normE_second);
            v3_shift{iEH,ic_loop} = TSeries(E_second.time+lagDiff*dt,normE_second);
            E_v1_plot{iEH,ic_loop}=E_first;
            E_v2_plot{iEH,ic_loop}=E_second;
            deltat_plot{iEH,ic_loop}=lagDiff*dt*10^3;
        end
        deltat(iEH,ic_loop)=lagDiff*dt*10^3;
        lag_val(iEH,ic_loop)=abs(lagDiff);
        vel_val(iEH,ic_loop)=vel.data;
        corr_val(iEH,ic_loop)=maxcorr;
        scSep_val(iEH,ic_loop)=dR2.data;
    end

end

%% Weigh velocities together:
%Use baselines with more than 9 deltaIndex, if meadian is below 9, use instead best 3 baselines.
for iEH=1:nrOfEHs
   
    median_lag(iEH) = median(lag_val(iEH,:));
    if median_lag(iEH)>=9
        
        lag_weight(iEH,:)=lag_val(iEH,:)>=9;
    else        
        lag_weight(iEH,:)=lag_val(iEH,:)>=median_lag(iEH);
    end
    
end
%To avoid NaN, set any infinite velocities to 0, they will not be used 
%anyway since they have a corresponding value of lag_val = 0.
vel_val(vel_val==inf)=0;
vel_val(vel_val==-inf)=0; 


dist_exp_w = 1; %Exponent to weight of spacecraft separation.
Phi_int_range=str2double(get(S.result_option(2),'string'));
% Number of peak-to-peak times from the centre of the EH we carry out integration
                              % to calculate potential. i.e. int_time=[T1-int_t_factor*tpp, T1+int_t_factor*tpp]
for iEH=1:nrOfEHs
    %Weigh velocities together based on spacecraft separation and corr coeff.
    vel_combined(iEH,1)=sum(lag_weight(iEH,:).*vel_val(iEH,:).*corr_val(iEH,:).*(abs(scSep_val(iEH,:)).^dist_exp_w))/(sum(lag_weight(iEH,:).*corr_val(iEH,:).*(abs(scSep_val(iEH,:)).^dist_exp_w)));
    lpp(iEH,:) = abs(vel_combined(iEH,1)).*tpp(iEH,:);
    
    %Calculate the integrated potential
    for iSC=1:length(ic)
        sc = ic(iSC);
        T_PP = tpp(iEH,sc);
        T = EH_times(iEH,sc);
        
        c_eval('eh_times = run_data.Exyz_filt?.time;',sc);
        c_eval('eh_data = run_data.Exyz_filt?.data(:,3);',sc);
        eh_tint = EpochTT(T)+[-Phi_int_range*T_PP,Phi_int_range*T_PP];
        [~,eh_start]=min(abs((eh_times) - (eh_tint(1))));
        [~,eh_stop]=min(abs((eh_times) - (eh_tint(2))));
        eh_TSeries = TSeries(eh_times,[zeros(eh_start-1,1);eh_data(eh_start:eh_stop);zeros(length(eh_times)-eh_stop,1)]);
        
        int_pot{iEH,sc} = -irf_integrate(eh_TSeries*(-vel_combined(iEH)));
        integrated_pot(iEH,sc)=max(int_pot{iEH,sc}.data);
    end
    
    
    
end
if corr_plot % This is not generalized for 2-3 spacecraft.. only 4...
    for iEH=1:nrOfEHs
        for ic_loop=1:nrOfCombs
            
            iicc=iiiccc{iEH,ic_loop};
            h=irf_plot(3,'newfigure');
            xSize=800; ySize=1000;
            set(gcf,'Position',[10 10 xSize ySize]);
            xwidth = 0.86;
            ywidth = 0.095;
            
            irf_plot(h(1),E_v1_plot{iEH,ic_loop});
            hold(h(1),'on');
            irf_plot(h(1),E_v2_plot{iEH,ic_loop});
            ylabel(h(1),'E_{||} (mV/m)','fontsize',18,'interpreter','tex');
            irf_legend(h(1),{legends_label1{iEH,ic_loop},legends_label2{iEH,ic_loop}},[0.98 0.12]);
            irf_legend(h(1),sprintf('EH # %i',iEH),[0.98 0.88]);
            
            irf_plot(h(2),v1_plot{iEH,ic_loop});
            hold(h(2),'on');
            irf_plot(h(2),v3_shift{iEH,ic_loop});
            ylabel(h(2),'E_{||} (mV/m)','fontsize',18,'interpreter','tex');
            irf_legend(h(2),{legends_label1{iEH,ic_loop},legends_label2{iEH,ic_loop}},[0.98 0.12]);
            irf_legend(h(2),sprintf('dt = %.2f ms', deltat_plot{iEH,ic_loop}),[0.28 0.88]);
            irf_legend(h(2),sprintf('v=%.2f km/s', vel_val(iEH,ic_loop)),[0.28,0.78]);
            
            irf_plot(h(3),int_pot{iEH,iicc(1)})
            hold(h(3),'on');
            irf_plot(h(3),int_pot{iEH,iicc(2)})
            ylabel(h(3),'\Phi (V)','interpreter','tex');
            
            irf_zoom(h(1:3),'x',plot_tint{iEH,ic_loop});
            irf_zoom(h(1:3),'y');
            % Add spacecraft separation and some other stuff too.
        end
    end
end

% Do x-correlation to obtain v, Phi, L_parallel.
outputvar = get(S.result_option(1),'string');
if isempty(outputvar)
    outputvar = 'EHProps';
end

        lag_val(iEH,ic_loop)=abs(lagDiff);
        vel_val(iEH,ic_loop)=vel.data;
        corr_val(iEH,ic_loop)=maxcorr;
        scSep_val(iEH,ic_loop)=dR2.data;
        lag_weight;
        
output.plasma.Ti=T_ion;
output.plasma.Te=T_electron;
output.plasma.Ld = Ld;
output.plasma.wpe=w_plasma;
output.plasma.wce=w_cyclotron;
output.plasma.cis=c_ion;
output.plasma.ne = n_electron;
output.plasma.B = B_field;
output.plasma.B_angle = Bangle;
output.plasma.vi = v_ion;

output.xcorr.dIndex = lag_val;
output.xcorr.dt = deltat;
output.xcorr.lag_weight = lag_weight;
output.xcorr.dR = scSep_val;
output.xcorr.v = vel_val;
output.xcorr.corr = corr_val;

output.v(:,1)=vel_combined;
if all_sc_flag
    output.v3d = v_3d;
end
output.potential = integrated_pot;
output.tpp = tpp;
output.Eperp1 = Eperp1;
output.Eperp2 = Eperp2;
output.Epp = Epp;
output.lpp = lpp;
output.obs_time=EH_times;
output.v_norm_cis = (vel_combined-mean(v_ion,2))./mean(c_ion,2);

assignin('base', outputvar, output);



%% Help functions:


%Finds the "number" of the spacecraft pair. Ex input: [1,2] gives pair_nr
%1, and pair_name 12. Or [3,4] gives pair_nr = 6, pair_name =34.
function [pair_nr,pair_name]=find_scpair(spacecraft)

if spacecraft(1)==1
    if spacecraft(2)==2
        pair_nr=1;
        pair_name=12;
    elseif spacecraft(2)==3
        pair_nr=2;
        pair_name=13;
    elseif spacecraft(2)==4
        pair_nr=3;
        pair_name=14;
    end
elseif spacecraft(1)==2
    
    if spacecraft(2)==3
        pair_nr=4;
        pair_name=23;
    elseif spacecraft(2)==4
        pair_nr=5;
        pair_name=24;
    end
elseif spacecraft(1)==3
    
    pair_nr=6;
    pair_name=34;
end

function [filtered_peaks1,filtered_peaks2,flag]=peak_filter_2sc(input_peaksTS1,input_peaksTS2,flag,dT_crit)

%Checking input arguments and defining what the filter should do.
emptySC=zeros(1,4);
emptyFlag=zeros(1,4);
c_eval('if isempty(input_peaksTS?), emptyFlag=1; emptySC(?)=?; end',[1:2]);

if emptyFlag
    c_eval('filtered_peaks?=[];',[1:2]);
else
    
    if flag(1)==1 && flag(2)==1
        c_eval('filtered_peaks?=input_peaksTS?;',[1:2]);
    else
        
        if ~(length(input_peaksTS1)==length(input_peaksTS2))
            
            if norm(flag)==0
                flag=[0,0];
                [maxPeaks,max_sc]=max([length(input_peaksTS1),length(input_peaksTS2)]);
            else
                max_sc=find(flag==0,1);
                c_eval('maxPeaks=length(input_peaksTS?);',max_sc);
            end
        else
            max_sc=find(flag==0,1);
            c_eval('maxPeaks=length(input_peaksTS?);',max_sc);
            flag(max_sc)=1;
        end
        sc=[1,2];
        c_sc=sc; c_sc(max_sc)=[];
        
        [~,emptyindex]=intersect(c_sc,emptySC);
        c_sc(emptyindex)=[];
        
        emptySC(emptySC==0)=[];
        nonemptySC=sc; nonemptySC(emptySC)=[];
        
        %Returning an empty TS for empty input TS...
        c_eval('filtered_peaks?=input_peaksTS?;',emptySC);
        
        %Splitting into positive and negative peaks:
        c_eval('pos_peaksIdx?=find(input_peaksTS?.data>0);',sc);
        c_eval('pos_peaksTS?=TSeries(input_peaksTS?.time(pos_peaksIdx?),input_peaksTS?.data(pos_peaksIdx?));',sc);
        c_eval('neg_peaksIdx?=find(input_peaksTS?.data<0);',sc);
        c_eval('neg_peaksTS?=TSeries(input_peaksTS?.time(neg_peaksIdx?),input_peaksTS?.data(neg_peaksIdx?));',sc);
        
        
        %Finding peaks that are not common to all signals.
        c_eval('pos_lonely_peak=zeros(1,length(pos_peaksIdx?));',max_sc);
        c_eval('neg_lonely_peak=zeros(1,length(neg_peaksIdx?));',max_sc);
        for ii=1:maxPeaks/2
            
            c_eval('dTpos(ii,!)=min(abs(pos_peaksTS?.time.epoch(ii)-pos_peaksTS!.time.epoch));',max_sc,c_sc);
            c_eval('if dTpos(ii,?)>dT_crit, pos_lonely_peak(ii)=1; end',c_sc);
            c_eval('dTneg(ii,!)=min(abs(neg_peaksTS?.time.epoch(ii)-neg_peaksTS!.time.epoch));',max_sc,c_sc);
            c_eval('if dTneg(ii,?)>dT_crit, neg_lonely_peak(ii)=1; end',c_sc);
            
        end
        lonely_peaks=neg_lonely_peak+pos_lonely_peak;
        lonely_peaks(lonely_peaks==2)=1;
        
        c_eval('pos_T_vec=pos_peaksTS?.time.epoch;',max_sc);
        c_eval('pos_D_vec=pos_peaksTS?.data;',max_sc);
        c_eval('neg_T_vec=neg_peaksTS?.time.epoch;',max_sc);
        c_eval('neg_D_vec=neg_peaksTS?.data;',max_sc);
        if max(lonely_peaks)>0
            pos_T_vec(lonely_peaks==1)=[];
            pos_D_vec(lonely_peaks==1)=[];
            neg_T_vec(lonely_peaks==1)=[];
            neg_D_vec(lonely_peaks==1)=[];
            
            [T_vec,order]=sort([pos_T_vec;neg_T_vec]);
            temp_vec=[pos_D_vec;neg_D_vec];
            D_vec=temp_vec(order,:);
            
            c_eval('input_peaksTS?=TSeries(EpochTT(T_vec),D_vec);',max_sc);
            flag=[0,0];
            c_eval('filtered_peaks?=input_peaksTS?;',nonemptySC);
            
            [filtered_peaks1,filtered_peaks2,flag]=peak_filter_2sc(input_peaksTS1,input_peaksTS2,flag,dT_crit);
            
        else
            flag(max_sc)=1;
            c_eval('filtered_peaks?=input_peaksTS?;',1:2);
            [filtered_peaks1,filtered_peaks2,flag]=peak_filter_2sc(input_peaksTS1,input_peaksTS2,flag,dT_crit);
            
        end
    end
end

function [matchingESW_TS1,matchingESW_TS2]=label_EHs_2sc(filtered2_peaksTS1,filtered2_peaksTS2,r_resamp1,r_resamp2,spacecraft_delay)
%This function first selects the spacecraft with the fewest solitary waves.
%The code then goes through the input peaks of the selected spacecraft one
%by one, and matches it to the solitary waves in the other spacecrafts'
%time series. The matching works the following way: First, for each
%solitary wave in the "minimum" spaceraft time series, a interval of
%+/-spacecraft_delay is created. All electron holes within this interval
%are selected for the other spacecraft. Then all possible combinations of
%the possible matches are formed. If the supposed match has the same (or
%opposite) sequence in time as the spacecraft have in space (along B), they
%are "approved". If only one such "approved" combination occurs, then they
%are said to be matching, and labeld as one electron hole appearing on all
%four spacecraft. However, if there are multiple "approved" combinations,
%the one time delays most fitting to the spatial separation (assuming
%constant velocity) is chosen as the best fit. 

ic=[1,2];
eswNr=1; numb=1;
c_eval('matching_ESW?t=int64([]);',ic);
c_eval('matching_ESW?d=[];',ic);

%Making a note of the spacecraft with fewest solitary waves.
[min_peak_nr,min_peaks_sc]=min([length(filtered2_peaksTS1),length(filtered2_peaksTS2)]);
for ii=1:2:min_peak_nr-1
    
    approved_peak=[0,0]; 
    temp_dist_check=[];
    temp_time_check=[];
    complete_flag=0;
    c_eval('tempESW_time=filtered2_peaksTS?.time.epoch(ii);',min_peaks_sc);
    c_eval('temp?times=filtered2_peaksTS?.time.epoch;',ic);
    c_eval('dt?=abs(temp?times-tempESW_time);',ic);
    c_eval('inside_int?=find(dt?<=spacecraft_delay);',ic);
    
    c_eval('if ~mod(inside_int?(1),2), temp_vec=inside_int?; inside_int?=[inside_int?(1)-1;temp_vec]; end',ic); %Make sure it start with odd and ends with even..
    c_eval('if mod(inside_int?(end),2), inside_int?(end+1)=inside_int?(end)+1; end',ic);
    
    
    c_eval('possible_ESW.mms?{ii}=TSeries(filtered2_peaksTS?.time(inside_int?),filtered2_peaksTS?.data(inside_int?));',ic);
    c_eval('possible_ESW.mms?{ii}=TSeries(filtered2_peaksTS?.time([ii;ii+1]),filtered2_peaksTS?.data([ii;ii+1]));',min_peaks_sc);
    
    [min_possible_nr,min_possible_nr_sc]=min([length(possible_ESW.mms1{ii}),length(possible_ESW.mms2{ii})]);
    for sc_loop=1:2
        c_eval('min_possible_nr=length(possible_ESW.mms?{ii});',sc_loop);
        for possible_peak_nr=1:min_possible_nr/2
            c_eval('temp_peakT=possible_ESW.mms?{ii}.time.epoch(possible_peak_nr*2-1:possible_peak_nr*2);',sc_loop);
            c_eval('temp_peakD=possible_ESW.mms?{ii}.data(possible_peak_nr*2-1:possible_peak_nr*2);',sc_loop);
            c_eval('isolated_peaks?{ii}{possible_peak_nr}=temp_peakT;',sc_loop);
            c_eval('isolated_peaksD?{ii}{possible_peak_nr}=temp_peakD;',sc_loop);
        end
    end
    
    elements = {1:length(isolated_peaks1{ii}), 1:length(isolated_peaks2{ii})}; %cell array with N vectors to combine
    combinations = cell(1, numel(elements)); %set up the varargout result
    [combinations{:}] = ndgrid(elements{:});
    combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
    result_combs = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.
    [comb_nr,~]=size(result_combs);
    dTdR_ratio=cell(1,comb_nr);
    dRdT_ratio=cell(1,comb_nr);
    sgn_flag=zeros(1,comb_nr);
    for comb_loop=1:comb_nr
        approved_peak=[0,0];
        temp_dist_check=[]; %%%%%%%%% NEW CHANGE
        temp_time_check=[];
        
        c_eval('[r_times?,r_indices?] = intersect(r_resamp?.time.epoch,isolated_peaks?{ii}{result_combs(comb_loop,?)});',ic);
        R_values=mean([r_resamp1.data(r_indices1),r_resamp2.data(r_indices2)]);
        [sorted_R_values,mms_R_order]=sort(R_values);
        [ESW_times{comb_loop},ESW_time_orders{comb_loop}]=sort([r_times1,r_times2],2);
        [~,sort_back]=sort(ESW_time_orders{comb_loop},2);
        
        
        for esw_Peak=1:2
            peak_sign=[1,-1];
            c_eval('if ~(sign(isolated_peaksD?{ii}{result_combs(comb_loop,?)}(esw_Peak))==sign(isolated_peaksD!{ii}{result_combs(comb_loop,!)}(esw_Peak))), sgn_flag(comb_loop)=1; end',ic, min_peaks_sc);
           
            c_eval('sign_thing = sign(isolated_peaksD?{ii}{result_combs(comb_loop,?)}(esw_Peak));',min_peaks_sc);
%             approved_peak(esw_Peak)=1;
            temp_ESW_time_orders=ESW_time_orders{comb_loop}(esw_Peak,:);
            
            if mean(1*(mms_R_order==temp_ESW_time_orders))==1*(sign_thing==1*peak_sign(esw_Peak))
                
                approved_peak(esw_Peak)=1;
            
            elseif mean(1*(mms_R_order==flip(temp_ESW_time_orders)))==1*(sign_thing==-1*peak_sign(esw_Peak))
                
                approved_peak(esw_Peak)=1;
            end
%             if mean(1*(mms_R_order==temp_ESW_time_orders))==1 || mean(1*(mms_R_order==flip(temp_ESW_time_orders)))==1
%                 approved_peak(esw_Peak)=1;
%                 
%             end
            
            c_eval('dT_?!=isolated_peaks?{ii}{result_combs(comb_loop,?)}(esw_Peak)-isolated_peaks!{ii}{result_combs(comb_loop,!)}(esw_Peak);',1,2);
            c_eval('if abs(double(dT_?!))*10^-6<1.8, temp_time_check(end+1,1:2)=[?,!]; end',1,2);
            c_eval('dR_?!=R_values(?)-R_values(!);',1,2);
            c_eval('if abs(dR_?!)<0.1, temp_dist_check(end+1,1:2)=[?,!]; end',1,2);

            
            if ~isempty(temp_time_check)
                [rows,cols]=size(temp_time_check);
                if rows>=3
                    approved_peak(esw_Peak)=1;
                end
                for jj=1:rows
                    
                    [~,swapping_index]=intersect(ESW_time_orders{comb_loop}(esw_Peak,:),temp_time_check(jj,:));
                    temp_mms_time_order=ESW_time_orders{comb_loop}(esw_Peak,:);
                    temp_mms_time_order([swapping_index(1),swapping_index(2)])=temp_mms_time_order([swapping_index(2),swapping_index(1)]);
                    
                    if (mean(1*(temp_mms_time_order==temp_ESW_time_orders))==1)*(sign_thing==1*peak_sign(esw_Peak)) || (mean(1*(temp_mms_time_order==flip(temp_ESW_time_orders)))==1)*(sign_thing==-1*peak_sign(esw_Peak))
                        approved_peak(esw_Peak)=1;
                    end
                    
                end
            end
            if ~isempty(temp_dist_check)
                [rows,cols]=size(temp_dist_check);
                for jj=1:rows
                    
                    [~,swapping_index]=intersect(mms_R_order,temp_dist_check(jj,:));
                    temp_mms_R_order=mms_R_order;
                    temp_mms_R_order([swapping_index(1),swapping_index(2)])=temp_mms_R_order([swapping_index(2),swapping_index(1)]);
                    
                    if (mean(1*(temp_mms_R_order==temp_ESW_time_orders))==1)*(sign_thing==1*peak_sign(esw_Peak)) || (mean(1*(temp_mms_R_order==flip(temp_ESW_time_orders)))==1)*(sign_thing==-1*peak_sign(esw_Peak))
                        approved_peak(esw_Peak)=1;
                    end
                    
                end
            end
        end
        
        if mean(approved_peak)==1 && sgn_flag(comb_loop)==0
            
            complete_flag(comb_loop)=comb_loop;
            temp_ESW_times1=ESW_times{comb_loop}(1,:);
            temp_ESW_times2=ESW_times{comb_loop}(2,:);
            temp_ESW_times1_sorted=temp_ESW_times1(:,sort_back(1,:));
            temp_ESW_times2_sorted=temp_ESW_times2(:,sort_back(2,:));
            temp_ESW_times_sorted{comb_loop}=[temp_ESW_times1_sorted;temp_ESW_times2_sorted];
            c_eval('comb_dT_?!{comb_loop}=abs(temp_ESW_times_sorted{comb_loop}(:,?)-temp_ESW_times_sorted{comb_loop}(:,!))'';',1,2);
            c_eval('comb_dR_?{comb_loop}=abs(dR_?);',12);
            c_eval('dRdT_ratio{comb_loop}(end+1,1:2)=double(comb_dT_?{comb_loop})./double(comb_dR_?{comb_loop});',12);
            
        end
    end
    
    complete_flag(complete_flag==0)=[];
    if length(complete_flag)==1
        %     Select the best match here::,
        best_combinations=result_combs(complete_flag,:);
        c_eval('matching_ESW?t(eswNr:eswNr+1)=isolated_peaks?{ii}{best_combinations(?)};',ic);
        c_eval('matching_ESW?d(eswNr:eswNr+1)=isolated_peaksD?{ii}{best_combinations(?)};',ic);
        
        best_dRdT{numb}=dRdT_ratio{complete_flag}(:,1);
        numb=numb+1;
        eswNr=eswNr+2;
    elseif length(complete_flag)>1
        asd1=[dRdT_ratio{complete_flag}];
        asd2=asd1(:,2:2:end);
        asd1(:,2:2:end)=[];
        korv1=std(asd1);
        korv2=std(asd2);
        [~,inx]=min(korv1);
        best_combination_nr=complete_flag(inx);
        best_combinations=result_combs(best_combination_nr,:);
        c_eval('matching_ESW?t(eswNr:eswNr+1)=isolated_peaks?{ii}{best_combinations(?)};',ic);
        c_eval('matching_ESW?d(eswNr:eswNr+1)=isolated_peaksD?{ii}{best_combinations(?)};',ic);
        best_dRdT{numb}=dRdT_ratio{best_combination_nr}(:,1);
        numb=numb+1;
        eswNr=eswNr+2;
        
    end
    %
    %     [length(matching_ESW1t),length(matching_ESW2t),length(matching_ESW3t),length(matching_ESW4t)]
     
end

%Find and remove duplicates based on best dRdT thingie. This can't
%be the optimal way of doing this...

c_eval('temp_ESWd?=matching_ESW?d;',1:2);
c_eval('temp_ESWt?=matching_ESW?t;',1:2);
if ~isempty(temp_ESWd1)
temp_best_dRdT=best_dRdT;
for loop_variable=1:2
    c_eval('[unique_times?, unique_indices?]=unique(temp_ESWt?);',loop_variable);
    c_eval('dups?=1:length(temp_ESWt?);',loop_variable);
    c_eval('uniqs?=intersect(unique_indices?,dups?);',loop_variable);
    c_eval('dups?(uniqs?)=[];',loop_variable);
    c_eval('dup_times?=temp_ESWt?(dups?(1:2:end-1));',loop_variable);
    c_eval('nrOfDups=length(dup_times?);',loop_variable);
    c_eval('dup_indx=[];',loop_variable);
    c_eval('corresponding_ESW?=[];',loop_variable);
    for duplicates=1:nrOfDups
        c_eval('[~,dup_indx]=find(dup_times?(duplicates)==temp_ESWt?);',loop_variable);
        c_eval('corresponding_ESW?=(dup_indx+1)/2;',loop_variable);
        
        c_eval('[min_std,min_indx]=min(std([temp_best_dRdT{corresponding_ESW?}]));',loop_variable);
        remove_indx=dup_indx;
        remove_indx(min_indx)=[];
        c_eval('temp_ESWd?([remove_indx,remove_indx+1])=[];',ic)
        c_eval('temp_ESWt?([remove_indx,remove_indx+1])=[];',ic)
        temp_best_dRdT((remove_indx+1)/2)=[];
    end
end
end
c_eval('matchingESW_TS?=TSeries(EpochTT(temp_ESWt?)'',temp_ESWd?'');',ic);

function [filtered_peaks1,filtered_peaks2,filtered_peaks3,filtered_peaks4,flag]=peak_filter_v2(input_peaksTS1,input_peaksTS2,input_peaksTS3,input_peaksTS4,flag,dT_crit)
%Removes peaks that do not appear in all four signals. Only works for
%well-spaced peaks since it's using dT_crit as the criteria, could probably
%be generalized by using the time shifts from cross-correlation but it
%could still be problematic. I apologize in advance for the messy code.
%Checking input arguments and defining what the filter should do.
emptySC=zeros(1,4);
emptyFlag=zeros(1,4);
c_eval('if isempty(input_peaksTS?), emptyFlag=1; emptySC(?)=?; end',[1:4]);

if emptyFlag
    c_eval('filtered_peaks?=[];',[1:4]);
else
    
    if flag(1)==1 && flag(2)==1 && flag(3)==1 && flag(4)==1
        c_eval('filtered_peaks?=input_peaksTS?;',[1:4]);
    else
        
        if ~(length(input_peaksTS1)==length(input_peaksTS2)==length(input_peaksTS3)==length(input_peaksTS4))
            
            if norm(flag)==0
                flag=[0,0,0,0];
                [maxPeaks,max_sc]=max([length(input_peaksTS1),length(input_peaksTS2),length(input_peaksTS3),length(input_peaksTS4)]);
            else
                max_sc=find(flag==0,1);
                c_eval('maxPeaks=length(input_peaksTS?);',max_sc);
            end
        else
            max_sc=find(flag==0,1);
            c_eval('maxPeaks=length(input_peaksTS?);',max_sc);
            flag(max_sc)=1;
        end
        sc=[1,2,3,4];
        c_sc=sc; c_sc(max_sc)=[];
        
        [~,emptyindex]=intersect(c_sc,emptySC);
        c_sc(emptyindex)=[];
        
        emptySC(emptySC==0)=[];
        nonemptySC=sc; nonemptySC(emptySC)=[];
        
        %Returning an empty TS for empty input TS...
        c_eval('filtered_peaks?=input_peaksTS?;',emptySC);
        
        %Splitting into positive and negative peaks:
        c_eval('pos_peaksIdx?=find(input_peaksTS?.data>0);',sc);
        c_eval('pos_peaksTS?=TSeries(input_peaksTS?.time(pos_peaksIdx?),input_peaksTS?.data(pos_peaksIdx?));',sc);
        c_eval('neg_peaksIdx?=find(input_peaksTS?.data<0);',sc);
        c_eval('neg_peaksTS?=TSeries(input_peaksTS?.time(neg_peaksIdx?),input_peaksTS?.data(neg_peaksIdx?));',sc);
        
        
        %Finding peaks that are not common to all signals.
        c_eval('pos_lonely_peak=zeros(1,length(pos_peaksIdx?));',max_sc);
        c_eval('neg_lonely_peak=zeros(1,length(neg_peaksIdx?));',max_sc);
        for ii=1:maxPeaks/2
            
            c_eval('dTpos(ii,!)=min(abs(pos_peaksTS?.time.epoch(ii)-pos_peaksTS!.time.epoch));',max_sc,c_sc);
            c_eval('if dTpos(ii,?)>dT_crit, pos_lonely_peak(ii)=1; end',c_sc);
            c_eval('dTneg(ii,!)=min(abs(neg_peaksTS?.time.epoch(ii)-neg_peaksTS!.time.epoch));',max_sc,c_sc);
            c_eval('if dTneg(ii,?)>dT_crit, neg_lonely_peak(ii)=1; end',c_sc);
            
        end
        lonely_peaks=neg_lonely_peak+pos_lonely_peak;
        lonely_peaks(lonely_peaks==2)=1;
        
        c_eval('pos_T_vec=pos_peaksTS?.time.epoch;',max_sc);
        c_eval('pos_D_vec=pos_peaksTS?.data;',max_sc);
        c_eval('neg_T_vec=neg_peaksTS?.time.epoch;',max_sc);
        c_eval('neg_D_vec=neg_peaksTS?.data;',max_sc);
        if max(lonely_peaks)>0
            pos_T_vec(lonely_peaks==1)=[];
            pos_D_vec(lonely_peaks==1)=[];
            neg_T_vec(lonely_peaks==1)=[];
            neg_D_vec(lonely_peaks==1)=[];
            
            [T_vec,order]=sort([pos_T_vec;neg_T_vec]);
            temp_vec=[pos_D_vec;neg_D_vec];
            D_vec=temp_vec(order,:);
            
            c_eval('input_peaksTS?=TSeries(EpochTT(T_vec),D_vec);',max_sc);
            flag=[0,0,0,0];
            c_eval('filtered_peaks?=input_peaksTS?;',nonemptySC);
            
            [filtered_peaks1,filtered_peaks2,filtered_peaks3,filtered_peaks4,flag]=peak_filter_v2(input_peaksTS1,input_peaksTS2,input_peaksTS3,input_peaksTS4,flag,dT_crit);
            
        else
            flag(max_sc)=1;
            c_eval('filtered_peaks?=input_peaksTS?;',[1:4]);
            [filtered_peaks1,filtered_peaks2,filtered_peaks3,filtered_peaks4,flag]=peak_filter_v2(input_peaksTS1,input_peaksTS2,input_peaksTS3,input_peaksTS4,flag,dT_crit);
            
        end
    end
end

function [matchingESW_TS1,matchingESW_TS2,matchingESW_TS3,matchingESW_TS4]=label_EHs(filtered2_peaksTS1,filtered2_peaksTS2,filtered2_peaksTS3,filtered2_peaksTS4,r_resamp1,r_resamp2,r_resamp3,r_resamp4,spacecraft_delay)
%This function first selects the spacecraft with the fewest solitary waves.
%The code then goes through the input peaks of the selected spacecraft one
%by one, and matches it to the solitary waves in the other spacecrafts'
%time series. The matching works the following way: First, for each
%solitary wave in the "minimum" spaceraft time series, a interval of
%+/-spacecraft_delay is created. All electron holes within this interval
%are selected for the other spacecraft. Then all possible combinations of
%the possible matches are formed. If the supposed match has the same (or
%opposite) sequence in time as the spacecraft have in space (along B), they
%are "approved". If only one such "approved" combination occurs, then they
%are said to be matching, and labeld as one electron hole appearing on all
%four spacecraft. However, if there are multiple "approved" combinations,
%the one time delays most fitting to the spatial separation (assuming
%constant velocity) is chosen as the best fit. This code is a real mess
%unfortunately, but it works well enough. So just sweep it under the rug...


eswNr=1; numb=1; ic=1:4;
c_eval('matching_ESW?t=int64([]);',ic);
c_eval('matching_ESW?d=[];',ic);

%Making a note of the spacecraft with fewest solitary waves.
[min_peak_nr,min_peaks_sc]=min([length(filtered2_peaksTS1),length(filtered2_peaksTS2),length(filtered2_peaksTS3),length(filtered2_peaksTS4)]);
for ii=1:2:min_peak_nr-1
    
    approved_peak=[0,0]; 
    temp_dist_check=[];
    temp_time_check=[];
    complete_flag=0;
    c_eval('starting_sign = sign(filtered2_peaksTS?.data(ii));',min_peaks_sc);
    c_eval('tempESW_time=filtered2_peaksTS?.time.epoch(ii);',min_peaks_sc);
    c_eval('temp?times=filtered2_peaksTS?.time.epoch;',ic);
    c_eval('dt?=abs(temp?times-tempESW_time);',ic);
    c_eval('inside_int?=find(dt?<=spacecraft_delay);',ic);
    c_eval('inside_int?=[ii;ii+1];',min_peaks_sc);
    
    c_eval('if ~mod(inside_int?(1),2), temp_vec=inside_int?; inside_int?=[inside_int?(1)-1;temp_vec]; end',ic); %Make sure it start with odd and ends with even..
    c_eval('if mod(inside_int?(end),2), inside_int?(end+1)=inside_int?(end)+1; end',ic);
    
    
    c_eval('possible_ESW.mms?{ii}=TSeries(filtered2_peaksTS?.time(inside_int?),filtered2_peaksTS?.data(inside_int?));',ic);
    c_eval('possible_ESW.mms?{ii}=TSeries(filtered2_peaksTS?.time([ii;ii+1]),filtered2_peaksTS?.data([ii;ii+1]));',min_peaks_sc);
    
    [min_possible_nr,min_possible_nr_sc]=min([length(possible_ESW.mms1{ii}),length(possible_ESW.mms2{ii}),length(possible_ESW.mms3{ii}),length(possible_ESW.mms4{ii})]);
    for sc_loop=1:4
        c_eval('min_possible_nr=length(possible_ESW.mms?{ii});',sc_loop);
        for possible_peak_nr=1:min_possible_nr/2
            c_eval('temp_peakT=possible_ESW.mms?{ii}.time.epoch(possible_peak_nr*2-1:possible_peak_nr*2);',sc_loop);
            c_eval('temp_peakD=possible_ESW.mms?{ii}.data(possible_peak_nr*2-1:possible_peak_nr*2);',sc_loop);
            c_eval('isolated_peaks?{ii}{possible_peak_nr}=temp_peakT;',sc_loop);
            c_eval('isolated_peaksD?{ii}{possible_peak_nr}=temp_peakD;',sc_loop);
        end
    end
    
    elements = {1:length(isolated_peaks1{ii}), 1:length(isolated_peaks2{ii}), 1:length(isolated_peaks3{ii}), 1:length(isolated_peaks4{ii})}; %cell array with N vectors to combine
    combinations = cell(1, numel(elements)); %set up the varargout result
    [combinations{:}] = ndgrid(elements{:});
    combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
    result_combs = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.
    [comb_nr,~]=size(result_combs);
    dTdR_ratio=cell(1,comb_nr);
    dRdT_ratio=cell(1,comb_nr);
    sgn_flag=zeros(1,comb_nr);
    for comb_loop=1:comb_nr
        approved_peak=[0,0];
        temp_dist_check=[]; 
        temp_time_check=[];
        
        c_eval('[r_times?,r_indices?] = intersect(r_resamp?.time.epoch,isolated_peaks?{ii}{result_combs(comb_loop,?)});',ic);
        R_values=mean([r_resamp1.data(r_indices1),r_resamp2.data(r_indices2),r_resamp3.data(r_indices3),r_resamp4.data(r_indices4)]);
        [sorted_R_values,mms_R_order]=sort(R_values);
        [ESW_times{comb_loop},ESW_time_orders{comb_loop}]=sort([r_times1,r_times2,r_times3,r_times4],2);
        [~,sort_back]=sort(ESW_time_orders{comb_loop},2);
        dist(1,:) = abs([R_values(1)-R_values(2),R_values(1)-R_values(3),R_values(1)-R_values(4)]);
        dist(2,:) = abs([R_values(1)-R_values(2),R_values(2)-R_values(3),R_values(2)-R_values(4)]);
        dist(3,:) = abs([R_values(1)-R_values(3),R_values(2)-R_values(3),R_values(3)-R_values(4)]);
        dist(4,:) = abs([R_values(1)-R_values(4),R_values(2)-R_values(4),R_values(3)-R_values(4)]);
        
        dist_crit = 3.5; %km
        total_dist = sum(dist,2);
        forced_order_check = dist>dist_crit; %Check if SC distances is larger than criteria.
        forced_order_fulfilled = sum(forced_order_check,2)==3; %If one spacecraft is further than criteria away from all the others it must be in the right place.
        [~,ok_sc] = find(forced_order_fulfilled'==1);
        %         if numel(ok_sc)>1
        %             [~,best_sc]=max(total_dist);
        %             ok_sc=best_sc;
        %         end
        for iEHpos = 1:length(ok_sc)
            [~,forced_appearance_pos(iEHpos)] = find(mms_R_order==ok_sc(iEHpos));
        end
        for esw_Peak=1:2
            peak_sign=[1,-1];
            temp_time_check=[];
            c_eval('if ~(sign(isolated_peaksD?{ii}{result_combs(comb_loop,?)}(esw_Peak))==sign(isolated_peaksD!{ii}{result_combs(comb_loop,!)}(esw_Peak))), sgn_flag(comb_loop)=1; end',ic, min_peaks_sc);
           
            c_eval('sign_thing = sign(isolated_peaksD?{ii}{result_combs(comb_loop,?)}(esw_Peak));',min_peaks_sc);
%             approved_peak(esw_Peak)=1;
            temp_ESW_time_orders=ESW_time_orders{comb_loop}(esw_Peak,:);
            
            if mean(1*(mms_R_order==temp_ESW_time_orders))==1*(sign_thing==1*peak_sign(esw_Peak))
                
                approved_peak(esw_Peak)=1;
            
            elseif mean(1*(mms_R_order==flip(temp_ESW_time_orders)))==1*(sign_thing==-1*peak_sign(esw_Peak))
                
                approved_peak(esw_Peak)=1;
            end
%             if mean(1*(mms_R_order==temp_ESW_time_orders))==1 || mean(1*(mms_R_order==flip(temp_ESW_time_orders)))==1
%                 approved_peak(esw_Peak)=1;
%                 
%             end
            
            c_eval('dT_?!=isolated_peaks?{ii}{result_combs(comb_loop,?)}(esw_Peak)-isolated_peaks!{ii}{result_combs(comb_loop,!)}(esw_Peak);',1,2:4);
            c_eval('dT_?!=isolated_peaks?{ii}{result_combs(comb_loop,?)}(esw_Peak)-isolated_peaks!{ii}{result_combs(comb_loop,!)}(esw_Peak);',2,3:4);
            c_eval('dT_?!=isolated_peaks?{ii}{result_combs(comb_loop,?)}(esw_Peak)-isolated_peaks!{ii}{result_combs(comb_loop,!)}(esw_Peak);',3,4);
            c_eval('if abs(double(dT_?!))*10^-6<1.8, temp_time_check(end+1,1:2)=[?,!]; end',1,2:4);
            c_eval('if abs(double(dT_?!))*10^-6<1.8, temp_time_check(end+1,1:2)=[?,!]; end',2,3:4);
            c_eval('if abs(double(dT_?!))*10^-6<1.8, temp_time_check(end+1,1:2)=[?,!]; end',3,4);
            
            c_eval('dR_?!=R_values(?)-R_values(!);',1,2:4);
            c_eval('dR_?!=R_values(?)-R_values(!);',2,3:4);
            c_eval('dR_?!=R_values(?)-R_values(!);',3,4);
            c_eval('if abs(dR_?!)<0.1, temp_dist_check(end+1,1:2)=[?,!]; end',1,2:4);
            c_eval('if abs(dR_?!)<0.1, temp_dist_check(end+1,1:2)=[?,!]; end',2,3:4);
            c_eval('if abs(dR_?!)<0.1, temp_dist_check(end+1,1:2)=[?,!]; end',3,4);
            
            if ~isempty(temp_time_check)
                [rows,cols]=size(temp_time_check);
                if rows>=3
                    
                  
                    goodpossc=[];
                    goodpos=[];
                    forced_posflag=zeros(1,length(ok_sc));
                    if ~isempty(ok_sc)
                        for jPOS=1:length(ok_sc)
                            for kSC = 1:4
                                if kSC==ok_sc
                                    goodpos=forced_appearance_pos; % The selected sc appears as number "goodpos(kk)"
                                    goodpossc=kSC; % The good spacecraft appearing at "goodpos(kk)" is "goodpossc(kk)"
                                    
                                end
                            end
                            if ~isempty(goodpossc)
                                switch goodpos
                                    
                                    case 1
                                        
                                        if (temp_ESW_time_orders(1)==goodpossc)*starting_sign==1
                                            %                                         approved_peak(esw_Peak)=1;
                                            forced_posflag(jPOS)=1;
                                            
                                        elseif (temp_ESW_time_orders(4)==goodpossc)*starting_sign==-1
                                            
                                            forced_posflag(jPOS)=1;
                                        end
                                        
                                    case 2
                                        if (temp_ESW_time_orders(2)==goodpossc)*starting_sign==1
                                            forced_posflag(jPOS)=1;
                                        elseif (temp_ESW_time_orders(3)==goodpossc)*starting_sign==-1
                                            forced_posflag(jPOS)=1;
                                            
                                        end
                                        
                                        
                                    case 3
                                        if (temp_ESW_time_orders(3)==goodpossc)*starting_sign==1
                                            
                                            forced_posflag(jPOS)=1;
                                        elseif (temp_ESW_time_orders(2)==goodpossc)*starting_sign==-1
                                            forced_posflag(jPOS)=1;
                                        end
                                        
                                    case 4
                                        if (temp_ESW_time_orders(4)==goodpossc)*starting_sign==1
                                            forced_posflag(jPOS)=1;
                                            
                                        elseif (temp_ESW_time_orders(1)==goodpossc)*starting_sign==-1
                                            
                                            forced_posflag(jPOS)=1;
                                        end
                                end
                                
                            else
                                approved_peak(esw_Peak)=1;
                            end
                        end
                        if mean(forced_posflag)==1
                            approved_peak(esw_Peak)=1;
                        end
                    else
                        approved_peak(esw_Peak)=1;
                    end
                end
                for jj=1:rows
                    
                    [~,swapping_index]=intersect(ESW_time_orders{comb_loop}(esw_Peak,:),temp_time_check(jj,:));
                    temp_mms_time_order=ESW_time_orders{comb_loop}(esw_Peak,:);
                    temp_mms_time_order([swapping_index(1),swapping_index(2)])=temp_mms_time_order([swapping_index(2),swapping_index(1)]);
                    
                    if (mean(1*(temp_mms_time_order==temp_ESW_time_orders))==1)*(sign_thing==1*peak_sign(esw_Peak)) || (mean(1*(temp_mms_time_order==flip(temp_ESW_time_orders)))==1)*(sign_thing==-1*peak_sign(esw_Peak))
                        approved_peak(esw_Peak)=1;
                    end
                    
                end
            end
            if ~isempty(temp_dist_check)
                [rows,cols]=size(temp_dist_check);
                for jj=1:rows
                    
                    [~,swapping_index]=intersect(mms_R_order,temp_dist_check(jj,:));
                    temp_mms_R_order=mms_R_order;
                    temp_mms_R_order([swapping_index(1),swapping_index(2)])=temp_mms_R_order([swapping_index(2),swapping_index(1)]);
                    
                    if (mean(1*(temp_mms_R_order==temp_ESW_time_orders))==1)*(sign_thing==1*peak_sign(esw_Peak)) || (mean(1*(temp_mms_R_order==flip(temp_ESW_time_orders)))==1)*(sign_thing==-1*peak_sign(esw_Peak))
                        approved_peak(esw_Peak)=1;
                    end
                    
                end
            end
        end
        
        if mean(approved_peak)==1 && sgn_flag(comb_loop)==0
            
            complete_flag(comb_loop)=comb_loop;
            temp_ESW_times1=ESW_times{comb_loop}(1,:);
            temp_ESW_times2=ESW_times{comb_loop}(2,:);
            temp_ESW_times1_sorted=temp_ESW_times1(:,sort_back(1,:));
            temp_ESW_times2_sorted=temp_ESW_times2(:,sort_back(2,:));
            temp_ESW_times_sorted{comb_loop}=[temp_ESW_times1_sorted;temp_ESW_times2_sorted];
            c_eval('comb_dT_?!{comb_loop}=abs(temp_ESW_times_sorted{comb_loop}(:,?)-temp_ESW_times_sorted{comb_loop}(:,!))'';',1,2:4);
            c_eval('comb_dT_?!{comb_loop}=abs(temp_ESW_times_sorted{comb_loop}(:,?)-temp_ESW_times_sorted{comb_loop}(:,!))'';',2,3:4);
            c_eval('comb_dT_?!{comb_loop}=abs(temp_ESW_times_sorted{comb_loop}(:,?)-temp_ESW_times_sorted{comb_loop}(:,!))'';',3,4);
            c_eval('comb_dR_?{comb_loop}=abs(dR_?);',[12,13,14,23,24,34]);
            c_eval('dRdT_ratio{comb_loop}(end+1,1:2)=double(comb_dT_?{comb_loop})./double(comb_dR_?{comb_loop});',[12,13,14,23,24,34]);
            
        end
    end
    
    complete_flag(complete_flag==0)=[];
    if length(complete_flag)==1
        %     Välj bästa matchningen här::,
        best_combinations=result_combs(complete_flag,:);
        c_eval('matching_ESW?t(eswNr:eswNr+1)=isolated_peaks?{ii}{best_combinations(?)};',ic);
        c_eval('matching_ESW?d(eswNr:eswNr+1)=isolated_peaksD?{ii}{best_combinations(?)};',ic);
        
        best_dRdT{numb}=dRdT_ratio{complete_flag}(:,1);
        numb=numb+1;
        eswNr=eswNr+2;
    elseif length(complete_flag)>1
        asd1=[dRdT_ratio{complete_flag}];
        asd2=asd1(:,2:2:end);
        asd1(:,2:2:end)=[];
        korv1=std(asd1);
        korv2=std(asd2);
        [~,inx]=min(korv1);
        best_combination_nr=complete_flag(inx);
        best_combinations=result_combs(best_combination_nr,:);
        c_eval('matching_ESW?t(eswNr:eswNr+1)=isolated_peaks?{ii}{best_combinations(?)};',ic);
        c_eval('matching_ESW?d(eswNr:eswNr+1)=isolated_peaksD?{ii}{best_combinations(?)};',ic);
        best_dRdT{numb}=dRdT_ratio{best_combination_nr}(:,1);
        numb=numb+1;
        eswNr=eswNr+2;
        
    end
    %
    %     [length(matching_ESW1t),length(matching_ESW2t),length(matching_ESW3t),length(matching_ESW4t)]
     
end

%Find and remove duplicates based on best dRdT thingie. This can't
%be the optimal way of doing this...

c_eval('temp_ESWd?=matching_ESW?d;',1:4);
c_eval('temp_ESWt?=matching_ESW?t;',1:4);
if ~isempty(temp_ESWd1)
temp_best_dRdT=best_dRdT;
for loop_variable=1:4
    c_eval('[unique_times?, unique_indices?]=unique(temp_ESWt?);',loop_variable);
    c_eval('dups?=1:length(temp_ESWt?);',loop_variable);
    c_eval('uniqs?=intersect(unique_indices?,dups?);',loop_variable);
    c_eval('dups?(uniqs?)=[];',loop_variable);
    c_eval('dup_times?=temp_ESWt?(dups?(1:2:end-1));',loop_variable);
    c_eval('nrOfDups=length(dup_times?);',loop_variable);
    c_eval('dup_indx=[];',loop_variable);
    c_eval('corresponding_ESW?=[];',loop_variable);
    for duplicates=1:nrOfDups
        c_eval('[~,dup_indx]=find(dup_times?(duplicates)==temp_ESWt?);',loop_variable);
        c_eval('corresponding_ESW?=(dup_indx+1)/2;',loop_variable);
        
        c_eval('[min_std,min_indx]=min(std([temp_best_dRdT{corresponding_ESW?}]));',loop_variable);
        remove_indx=dup_indx;
        remove_indx(min_indx)=[];
        c_eval('temp_ESWd?([remove_indx,remove_indx+1])=[];',ic)
        c_eval('temp_ESWt?([remove_indx,remove_indx+1])=[];',ic)
        temp_best_dRdT((remove_indx+1)/2)=[];
    end
end
end
c_eval('matchingESW_TS?=TSeries(EpochTT(temp_ESWt?)'',temp_ESWd?'');',ic);






