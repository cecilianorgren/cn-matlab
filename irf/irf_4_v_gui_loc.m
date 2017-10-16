function out = irf_4_v_gui_loc(varargin)
%IRF_4_V_GUI interactive discontinuity analyzer for Cluster and MMS.
%
%   IRF_4_V_GUI('B?',mission) use interactive discontinuity analyzer on
%   variables B1,B2,B3 and B4 in the Workspace. The function loads position
%   data depending on the variable mission. mission is either 'Cluster' or
%   'MMS'.
%   
%   IRF_4_V_GUI(B1,B2,B3,B4,mission) uses variables B1,B2,B3 and B4.
%   
%   IRF_4_V_GUI('B?','R?') uses variables B1,B2,B3,B4 and R1,R2,R3,R4 in
%   Workspace. R? is a Nx4 matrix which contains time data (Cluster
%   format). R? can also be TSeries objects.
%
%   IRF_4_V_GUI(B1,B2,B3,B4,R1,R2,R3,R4) Also possible.
%
%   IRF_4_V_GUI(...,mission,column) Allows to specify mission and column
%   used for plotting.
%
%   The function first tries to load position data from local disk. For
%   Cluster, data might be downloaded.
%
%   Examples:
%       irf_4_v_gui_loc('B?','mms',3) % Loads position data of MMS from disk.
%
%       % First read predicted position data from dfg cdf file.
%       c_eval('R?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_ql_pos_gse'',tint);',1:4)
%       irf_4_v_gui_loc('B?','R?') % Input both field and position
%
%   See also: IRF_4_V

%persistent ud

%% Input-------------
if nargin == 0
    % Display help
    help irf_4_v_gui
    return;
elseif nargin == 1,
    error;
elseif nargin==2 && ~ischar(varargin{1}) &&  all(ishandle(varargin{1}))
    % Only action input
    hfig=varargin{1};
    action = varargin{2};
    ud=get(hfig,'userdata');
    ud.flag_first_call = 0;
else
    ud = [];
    ud.flag_first_call = 1;
    ud.hfig=figure;
    if(nargin>=2 && nargin<=4)
        if ischar(varargin{1}) && ~isempty(strfind(varargin{1},'?'))
            ud.variable_str = varargin{1};
            c_eval('ud.var?=evalin(''base'',irf_ssub(ud.variable_str,?));');
        end
        if ischar(varargin{2}) && ~isempty(strfind(varargin{2},'?'))
            ud.pos_str = varargin{2};
            c_eval('R?=evalin(''base'',irf_ssub(ud.pos_str,?));');
            if isa(R1,'TSeries')
                ud = r_ts2mat(ud,R1,R2,R3,R4);
            else
                c_eval('ud.pos?=R?;',1:4)
            end
            
        end
        % Set the rest of parameters if applicable
        ud = set_col_and_sc(ud,varargin{2:end});
        
    elseif(nargin>=5 && nargin<=7)
        ud.var1 = varargin{1};
        ud.var2 = varargin{2};
        ud.var3 = varargin{3};
        ud.var4 = varargin{4};
        ud.variable_str = [inputname(1) '..' inputname(4)];
        ud = set_col_and_sc(ud,varargin{5:end});
        
    elseif(nargin>=8 && nargin<=10)
        ud.var1 = varargin{1};
        ud.var2 = varargin{2};
        ud.var3 = varargin{3};
        ud.var4 = varargin{4};
        ud.variable_str = [inputname(1) '..' inputname(4)];
        
        R1 = varargin{5};
        R2 = varargin{6};
        R3 = varargin{7};
        R4 = varargin{8};
        
        if isa(R1,'TSeries')
            ud = r_ts2mat(ud,R1,R2,R3,R4);
        else
            c_eval('ud.pos?=R?;',1:4)
        end
        
        ud = set_col_and_sc(ud,varargin{9:end});
    
    else
        error('Unknown input type.')
    end
    set(ud.hfig,'userdata',ud)
end
%----------------------

%% Actions
if ud.flag_first_call
    % Initialize gui and read position data if not inputted.
    ud=irf_4_v_gui_loc(ud.hfig,'init');
    
    if ~is_pos_ok(ud)
        if ~isfield(ud,'sc')
            error('Wrong format of position data')
        end
        if strcmpi(ud.sc,'mms')
            ud = get_mms_pos(ud);
        elseif strcmpi(ud.sc,'cluster')
            ud = get_c_pos(ud);
        end
    end
    
else
    switch action
        case 'init'
            ud = init_gui(ud);
        case {'c1','c2','c3','c4','c5','c6'}
            ud.var_col=str2double(action(2:end));
            set(gcf,'userdata',ud);
            ud=irf_4_v_gui_loc(ud.hfig,'update_var_col');
        case 'dt'
            ud = v_from_dt(ud);
        case 'v'
            ud = dt_from_v(ud);
        case 'autoY'
            ud = auto_y(ud);
        case 'click_times'
            ud = click_times(ud);
        case 'distance'
            ud = time2distance(ud);
        case 'new_var_enter'
            ud = new_var_enter(ud);
        case 'new_var'
            ud = new_var(ud);
        case 'update_var_col'
            ud = update_var_col(ud);
        case 'update_boundary'
            ud = update_boundary(ud);
        otherwise % Mostly for debugging
            error(['Not implemented action: ', action])
    end
end
set(ud.hfig,'userdata',ud);
if nargout == 1
    out = ud;
end
end

%% Action functions
function ud = init_gui(ud)
%   Initiates the GUI and creates all menu items and input boxes.

irf_figmenu;
set(gcf,'position',get(gcf,'position').*[1 1 1.4 1])
set(gcf,'color','white'); % white background for figures (default is grey)
set(gcf,'userdata',ud); % because irf_pl_tx can also changed userdata)

% Subplot for field data
h(1)=subplot(3,4,[1 2 3]);
%h(1)=subplot(3,1,1);
irf_pl_tx(h(1),ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col);zoom(h(1),'on');
ylabel(h(1),var_label(ud.variable_str,ud.var_col));

%Subplot for time-shifted field data
h(2)=subplot(3,4,[5 6 7]);
%h(2)=subplot(3,1,2);
irf_pl_tx(h(2),ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col);
ylabel(h(2),var_label(ud.variable_str,ud.var_col));

% Subplot spacecraft position
h(3)=subplot(3,4,[4]);
%h(3)=subplot(6,4,23);
hold(h(3),'on')
plot(h(3),mean(ud.relpos1(:,2)),mean(ud.relpos1(:,3)),'ks','markersize',10)
plot(h(3),mean(ud.relpos2(:,2)),mean(ud.relpos2(:,3)),'rd','markersize',10)
plot(h(3),mean(ud.relpos3(:,2)),mean(ud.relpos3(:,3)),'go','markersize',10)
plot(h(3),mean(ud.relpos4(:,2)),mean(ud.relpos4(:,3)),'bv','markersize',10)
xlabel(h(3),'x')
ylabel(h(3),'y')
box(h(3),'on')
axis(h(3),'equal')
set(h(3),'XLim',max(abs(get(h(3),'XLim')))*[-1 1]*1.1)
set(h(3),'YLim',max(abs(get(h(3),'YLim')))*[-1 1]*1.1)
hold(h(3),'off')

% Subplot spacecraft position
h(4)=subplot(3,4,[8]);
%h(4)=subplot(6,4,24);
hold(h(4),'on')
plot(h(4),mean(ud.relpos1(:,2)),mean(ud.relpos1(:,4)),'ks','markersize',10)
plot(h(4),mean(ud.relpos2(:,2)),mean(ud.relpos2(:,4)),'rd','markersize',10)
plot(h(4),mean(ud.relpos3(:,2)),mean(ud.relpos3(:,4)),'go','markersize',10)
plot(h(4),mean(ud.relpos4(:,2)),mean(ud.relpos4(:,4)),'bv','markersize',10)
xlabel(h(4),'x')
ylabel(h(4),'z')
box(h(4),'on')
axis(h(4),'equal')
set(h(4),'XLim',max(abs(get(h(4),'XLim')))*[-1 1]*1.1)
set(h(4),'YLim',max(abs(get(h(4),'YLim')))*[-1 1]*1.1)
hold(h(4),'off')

% Subplot spacecraft position
h(5)=subplot(3,4,[12]);
hca=h(5);
%h(5)=subplot(6,4,24);
hold(hca,'on')
plot(hca,mean(ud.relpos1(:,3)),mean(ud.relpos1(:,4)),'ks','markersize',10)
plot(hca,mean(ud.relpos2(:,3)),mean(ud.relpos2(:,4)),'rd','markersize',10)
plot(hca,mean(ud.relpos3(:,3)),mean(ud.relpos3(:,4)),'go','markersize',10)
plot(hca,mean(ud.relpos4(:,3)),mean(ud.relpos4(:,4)),'bv','markersize',10)
xlabel(hca,'y')
ylabel(hca,'z')
box(hca,'on')
axis(hca,'equal')
set(hca,'XLim',max(abs(get(hca,'XLim')))*[-1 1]*1.1)
set(hca,'YLim',max(abs(get(hca,'YLim')))*[-1 1]*1.1)
hold(hca,'off')

% add information to the plot
irf_legend(0,['irf\_4\_v\_gui ' datestr(now)],[0.01 0.99],'fontsize',7); 
irf_legend(h(1),sc_label(ud),[1, 1.1],'color','cluster');

hh=h(1,1);  % use the first subplot to estimate available time interval
xl=get(hh,'XLim');
hc=get(hh,'Children');
xd=get(hc(end),'XData');
avail=[min([xl xd]) max([xl xd])];
%presel=xl;

dt = 0.02*diff(avail);
xlim = [avail(1)-dt avail(2)+dt];
ttics = timeaxis(xlim);

ud.tlim = avail;
ud.h=h;

% dt text input
xp=0.05;yp=0.25;
uicontrol('style', 'text', 'string', '[dt1 dt2 dt3 dt4] =','units','normalized','position', [xp yp 0.15 0.03]);
ud.dt_input = uicontrol('style', 'edit', ...
    'string', '[0 0 0 0]', ...
    'callback', 'irf_4_v_gui_loc(gcf,''dt'')', ...
    'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.2 0.05]);
ud.dt=[0 0 0 0]; % default values

% Velocity text input
xp=0.05;yp=0.2;
uicontrol('style', 'text', 'string', '[vx vy vz] km/s =','units','normalized','position', [xp yp 0.15 0.03]);
ud.v = uicontrol('style', 'edit', ...
    'string', '0*[0 0 0]', ...
    'callback', 'irf_4_v_gui_loc(gcf,''v'')', ...
    'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.2 0.05]);

% Low pass filter text input.
xp=0.05;yp=0.15;
uicontrol('style', 'text', 'string', 'Low pass filter f/Fs = ','units','normalized','position', [xp yp 0.15 0.03]);
ud.filter = uicontrol('style', 'edit', ...
    'string', '1', ...
    'callback', 'irf_4_v_gui_loc(gcf,''dt'')', ...
    'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.1 0.05]);

% GSM checkbox
xp=0.05;yp=0.10;
ud.coord_sys = uicontrol('style', 'checkbox', ...
    'string', 'velocity in GSM', ...
    'callback', 'irf_4_v_gui_loc(gcf,''dt'')', ...
    'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.3 0.05]);

% Reference satellite text input.
xp=0.05;yp=0.05;
uicontrol('style', 'text', 'string', 'Reference satellite ','units','normalized','position', [xp yp 0.15 0.03]);
ud.ref_satellite = uicontrol('style', 'edit', ...
    'string', '1', ...
    'callback', 'irf_4_v_gui_loc(gcf,''v'')', ...
    'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.1 0.05]);

uimenu('label','Auto &YLim','accelerator','y','callback','irf_4_v_gui_loc(gcf,''autoY'')');
uimenu('label','&Distance','accelerator','d','callback','irf_4_v_gui_loc(gcf,''distance'')');
uimenu('label','Click&Times','accelerator','t','callback','irf_4_v_gui_loc(gcf,''click_times'')');
uimenu('label','New&Variable','accelerator','v','callback','irf_4_v_gui_loc(gcf,''new_var_enter'')');
ud.columns=uimenu('label','&Columns','accelerator','c');
%ud.t_start_epoch = h(1).Parent.UserData.t_start_epoch;
 figud = get(get(h(1),'Parent'),'UserData');
 ud.t_start_epoch = figud.t_start_epoch;
%uimenu(ud.column

if isa(ud.var1,'TSeries'), nCol = size(ud.var1.data,2);
else nCol = size(ud.var1,2);
end

for j_col=1:nCol
    eval_str=['ud.hcol(j_col)=uimenu(ud.columns,''label'',''' num2str(j_col) ''',''callback'',''irf_4_v_gui_loc(gcf,''''c' num2str(j_col) ''''')'');'];
    eval(eval_str);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(ud.hfig,'userdata',ud);
end


function ud = v_from_dt(ud)
%   Calculates velocity from time difference in text box.

hh=ud.h(1);  % use the first subplot to estimate available time interval
xl=get(hh,'XLim')+ud.t_start_epoch;
yl=get(hh,'YLim');
%hc=get(hh,'Children');
ud.dt=eval(['[' get(ud.dt_input,'string') ']']);
tstr=['[' num2str(ud.dt,'%9.2f') '] s'];
t=0.5*(xl(1)+xl(2))+ud.dt;
if max(abs(ud.dt))==0
    vstr='0 * [0 0 0]';
    v = NaN;
else
    v=irf_4_v(ud.pos1,ud.pos2,ud.pos3,ud.pos4,t);
    if strcmpi(coord_sys(ud),'gsm')
        v = irf_gse2gsm([mean(t),v]);
        v = v(2:4);
    end
    vstr=[num2str(norm(v),3) ' * [' num2str(v./norm(v),'%6.2f') ']'];
end
%set(ud.v,'string',vstr);
ud.vVector=v;
ud.nVector=v./norm(v);
if eval(get(ud.filter,'string'))<1
    x1=ud.var1;Fs=1/(x1(2,1)-x1(1,1));flim=Fs*eval(get(ud.filter,'string'));
    c_eval('x?=irf_tlim(ud.var?,xl+[-20/Fs 20/Fs]);x?=irf_filt(x?,0,flim,Fs,5);');
    irf_pl_tx(ud.h(2),'x?',ud.var_col,ud.dt);
else
    irf_pl_tx(ud.h(2),'ud.var?',ud.var_col,ud.dt);
end
irf_zoom(ud.h(2),'x',xl);
irf_zoom(ud.h(2),'y',yl);
irf_timeaxis(ud.h(2));
text(.5,-.6,['t_{2nd panedel} = t_{1st panel} - dt\newline  dt = ' tstr '\newline V_{discontinuity}=' vstr ' km/s ' coord_sys(ud) ],'units','normalized','verticalalignment','top','paren',ud.h(2));
set(ud.hfig,'userdata',ud);
ud = update_boundary(ud);
end


function ud = dt_from_v(ud)
%   Calculates time difference from velocity in text box.

hh = ud.h(1);  % use the first subplot to estimate available time interval
xl = get(hh,'XLim'); yl = get(hh,'YLim');
hc = get(hh,'Children');
v = eval(['[' get(ud.v,'string') ']']);
t = ud.t_start_epoch +0.5*(xl(1) +xl(2));
if max(abs(v))==0
    dt=[0 0 0 0];
else
    if strcmpi(coord_sys(ud),'gsm')
        v = irf_gse2gsm([t(1) v], -1);
        v = v(2:4);
    end
    
    % The actual calculation
    dt=irf_4_v(ud.pos1,ud.pos2,ud.pos3,ud.pos4,[t v]); 
    
    ref_satellite_string=get(ud.ref_satellite,'string');
    ref_satellite=str2double(ref_satellite_string);
    if ref_satellite<1 || ref_satellite>4, ref_satellite=1;end
    dt=dt-dt(ref_satellite);
end
tstr=['[' num2str(dt,'%9.2f') ']'];
if norm(v) > 0,
    vstr=[num2str(norm(v),3) ' * [' num2str(v./norm(v),'%6.2f') ']'];
else
    vstr='0*[0 0 0]';
end
set(ud.dt_input,'string',tstr);
%set(ud.v,'string',vstr); % print the normalized string
ud.vVector=v; % update userdata velocity 
ud.nVector=v./norm(v); % and direction

if eval(get(ud.filter,'string'))<1,
    x1=ud.var1;Fs=1/(x1(2,1)-x1(1,1));
    flim=Fs*eval(get(ud.filter,'string'));
    c_eval('x?=irf_tlim(var?,xl+[-20/Fs 20/Fs]);x?=irf_filt(x?,0,flim,Fs,5);');
    irf_pl_tx(ud.h(2),x1,x2,x3,x4,ud.var_col,dt);
else
    irf_pl_tx(ud.h(2),ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col,dt);
end
axis(ud.h(2),[xl yl]);
irf_timeaxis(ud.h(2));
text(.5,-.6,['t_{2nd panedel} = t_{1st panel} - dt\newline dt = ' tstr '\newline V_{discontinuity}=' vstr ' km/s ' coord_sys(ud)],'units','normalized','verticalalignment','top','paren',ud.h(2));
set(ud.hfig,'userdata',ud);
ud = update_boundary(ud);
end


function ud = click_times(ud) 
%   Starts prompt to click to specify times for analysis.

zoom(ud.h(1),'off');
if (~isfield(ud,'ic') || isempty(ud.ic)), ud.ic=0;ud.dtv=[];end
if ud.ic==0,
    set(ud.hfig,'windowbuttondownfcn', 'irf_4_v_gui_loc(gcf,''click_times'')');
    ud.ic=1;
else
    p = get(ud.h(1), 'currentpoint');
    ud.dtv(ud.ic)=p(1);
    ud.ic=ud.ic+1;
end
title(['click on s/c ' num2str(ud.ic)]);
if ud.ic==5,
    set(ud.hfig,'windowbuttondownfcn', '');
    title('');
    ud.ic=0;
    ud.dt=ud.dtv-ud.dtv(1);
    tstr=['[' num2str(ud.dt,'%10.2f') ']'];
    set(ud.dt_input,'string',tstr);
    zoom(ud.h(1),'on');
    set(ud.hfig,'userdata',ud);
    ud=irf_4_v_gui_loc(ud.hfig,'dt');
end
end


function ud = auto_y(ud) 
%   Updates the limits of the y-axis for both panels.

for h=ud.h(1:2),
    set(h,'YLimMode','auto');
end
end


function ud = time2distance(ud)
%   Converts time to distance given velocity for the time-shifted panel.

hh=ud.h(1);  % use the first subplot to estimate available time interval
xl=get(hh,'XLim');
v=eval(['[' get(ud.v,'string') ']']);
tcenter = mean(xl);distance=norm(v)*diff(xl)/2;logd=log10(distance);
if logd>round(logd), dx=10^(round(logd))/2;
else dx=10^(round(logd))/5;
end
xticks=(-30:30)*dx/norm(v)/5+tcenter;
xticklabels=cell(size(-30:30));
for j=-30:30, xticklabels{j+31}=' ';end
for j=-6:6, xticklabels{j*5+31}=num2str(j*dx);end
set(ud.h(2),'xtick',xticks,'xticklabel',xticklabels);
xlabel(ud.h(2),'km');
end


function ud = new_var_enter(ud)
%   Input for new variable to be read from Workspace

xx=inputdlg('Enter new variable mask. Examples: B? or R? or P?p1','**',1,{'B?'});
ud.variable_str=xx{1};
set(ud.hfig,'userdata',ud);
ud = new_var(ud);
end

%% Auxiliary functions
function ud = new_var(ud)
%   Reads in new variable from Workspace.

evalin('base',['if ~exist(''' irf_ssub(ud.variable_str,1) '''), c_load(''' ud.variable_str ''');end' ]);
c_eval('ud.var?=evalin(''base'',irf_ssub(ud.variable_str,?));');
if ud.var_col > size(ud.var1,2), ud.var_col=2;end % in case new variable has less columns
if ud.flag_first_call,
    set(ud.hfig,'userdata',ud);
    irf_4_v_gui_loc(ud.hfig,'init');
else
    if ishandle(ud.h(1)),
        for j_col=2:size(ud.var1,2)
            if j_col<=length(ud.hcol),
                if ishandle(ud.hcol(j_col)),
                    set(ud.hcol(j_col),'enable','on')
                else
                    eval_str=['ud.hcol(j_col)=uimenu(ud.columns,''label'',''' num2str(j_col) ''',''callback'',''irf_4_v_gui_loc(gcf,''''c' num2str(j_col) ''''')'');'];
                    eval(eval_str);
                end
            else
                eval_str=['ud.hcol(j_col)=uimenu(ud.columns,''label'',''' num2str(j_col) ''',''callback'',''irf_4_v_gui_loc(gcf,''''c' num2str(j_col) ''''')'');'];
                eval(eval_str);
            end
            for jj_col=(size(ud.var1,2)+1):length(ud.hcol)
                set(ud.hcol(jj_col),'enable','off')
            end
        end
        set(ud.hfig,'userdata',ud);
        irf_4_v_gui_loc(ud.hfig,'update_var_col');
    else
        set(ud.hfig,'userdata',ud);
        irf_4_v_gui_loc(ud.hfig,'init');
    end
end
end


function ud = update_var_col(ud)
%   Updates which column is used for plotting.

hca=ud.h(1);
xl=get(ud.h(1),'XLim');
irf_pl_tx(hca,ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col);zoom(hca,'on');
ylabel(ud.h(1),var_label(ud.variable_str,ud.var_col));
axis(hca,[xl(1) xl(2) 0 1]);
irf_legend(hca,sc_label(ud),[1, 1.1],'color','cluster');
irf_zoom(hca,'y');irf_timeaxis(hca)
irf_pl_tx(ud.h(2),ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col,ud.dt);
ylabel(ud.h(2),var_label(ud.variable_str,ud.var_col));
end


function ud = update_boundary(ud)
  % plot or update the boundary orientation and velocity direction
  
  % Delete old boundary and velocity quiver
  delete(findobj(ud.h(3).Children,'Tag','boundary'));
  delete(findobj(ud.h(3).Children,'Tag','quiver'));
  delete(findobj(ud.h(4).Children,'Tag','boundary'));
  delete(findobj(ud.h(4).Children,'Tag','quiver'));
  delete(findobj(ud.h(5).Children,'Tag','boundary'));
  delete(findobj(ud.h(5).Children,'Tag','quiver'));
  
  % Equations for plane
  funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
  funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
  funZ = @(x,y,n) (-n(1)*x-n(2)*y)/n(3);
  
  % Plot boundary plane and velocity quiver
  % xy-plane
  hca = ud.h(3); 
  xl = get(hca,'XLim');
  yl = get(hca,'YLim');
  if ud.nVector(1) == 0 && ud.nVector(2) == 0
    plotX = [NaN NaN];
    plotY = [NaN NaN];
  elseif ud.nVector(2) == 0    
    plotX = [0 0];
    plotY = yl;
  else
    plotX = xl;
    plotY = funY(xl,[0 0],ud.nVector);    
  end
  
  hold(hca,'on')
  plot(hca,plotX,plotY,'color',0.1*[0 0 0],'tag','boundary')
  quiver(hca,0,0,ud.nVector(1),ud.nVector(2),4,'tag','quiver','color',[0 0 0])
  hold(hca,'off')
  set(hca,'XLim',xl,'YLim',yl)
  
  % xz-plane
  hca = ud.h(4); 
  xl = get(hca,'XLim');
  yl = get(hca,'YLim');
  if ud.nVector(1) == 0 && ud.nVector(3) == 0
    plotX = [NaN NaN];
    plotZ = [NaN NaN];
  elseif ud.nVector(3) == 0    
    plotX = [0 0];
    plotZ = yl;
  else
    plotX = xl;
    plotZ = funZ(xl,[0 0],ud.nVector);
  end  

  hold(hca,'on')
  plot(hca,plotX,plotZ,'color',0.1*[0 0 0],'tag','boundary')
  quiver(hca,0,0,ud.nVector(1),ud.nVector(3),4,'tag','quiver','color',[0 0 0])
  hold(hca,'off')
  set(hca,'XLim',xl,'YLim',yl)  
  
  % yz-plane
  hca = ud.h(5); 
  xl = get(hca,'XLim');
  yl = get(hca,'YLim');
  if ud.nVector(2) == 0 && ud.nVector(3) == 0
    plotY = [NaN NaN];
    plotZ = [NaN NaN];
  elseif ud.nVector(3) == 0    
    plotY = [0 0];
    plotZ = yl;
  else
    plotY = xl;
    plotZ = funZ([0 0],yl,ud.nVector);
  end  

  hold(hca,'on')
  plot(hca,plotY,plotZ,'color',0.1*[0 0 0],'tag','boundary')
  quiver(hca,0,0,ud.nVector(2),ud.nVector(3),4,'tag','quiver','color',[0 0 0])
  hold(hca,'off')
  set(hca,'XLim',xl,'YLim',yl) 
  1;
end


function label = var_label(var_str,var_col)
%   Returns variable labels.

iVecComponent = var_col; % number of vector component
dd=c_desc(irf_ssub(var_str,1));
if isempty(dd)
    label=[var_str '[' num2str(iVecComponent) ']'];
else
    if numel(dd.units)==1,
        labUnit = dd.units{1};
    else
        labUnit = dd.units{iVecComponent};
    end
    if numel(dd.labels)==1,
        labVar = dd.labels{1};
    else
        labVar = dd.labels{iVecComponent};
    end
    if isfield(dd,'col_labels')
        colLabels=dd.col_labels{1};
        labVar = [labVar colLabels{iVecComponent}];
    end
    label=[labVar '[' labUnit ']'];
end
end


function label = sc_label(ud)
%   Returns labels for the spacecraft.

if ~isfield(ud,'sc')
    label = {'SC1','SC2','SC3','SC4'};
else
    if strcmpi(ud.sc,'Cluster')
        label = {'C1','C2','C3','C4'};
    elseif strcmpi(ud.sc,'MMS')
        label = {'MMS1','MMS2','MMS3','MMS4'};
    end
end
end


function coord_sys_label = coord_sys(ud)
%   Returns coordinate system string 'GSE' or 'GSM'

coord_sys_flag=get(ud.coord_sys,'value');
if coord_sys_flag == 1
    coord_sys_label='GSM';
else
    coord_sys_label='GSE';
end
end


function answer=is_pos_ok(ud)
%   Checks if position data in ud is ok. Returns true or false.

if(isfield(ud,'pos1') && isfield(ud,'pos2') && isfield(ud,'pos3') && isfield(ud,'pos4'))
    s1 = size(ud.pos1,2);
    s2 = size(ud.pos2,2);
    s3 = size(ud.pos3,2);
    s4 = size(ud.pos4,2);
    if (s1==4 && s2==4 && s3==4 && s4==4)
        answer = true;
    else
        answer = false;
    end
else
    answer = false;
end
end


function read_RV_from_caa_stream(tint)
%   Downloads data from CAA and reads it.

currentDir = pwd;
tempDir = tempname;
mkdir(tempDir);
cd(tempDir);
caa_download(tint,'CL_SP_AUX','stream');
cd('CAA/CL_SP_AUX');
d=dir('*.cef.gz');
cefFile = d.name;
pos = c_caa_cef_var_get('sc_r_xyz_gse',cefFile);
for sc='1234'
    tempR = c_caa_cef_var_get(['sc_dr' sc '_xyz_gse'],cefFile);
    ud.(['pos' sc])=pos+[zeros(size(pos,1),1) tempR(:,2:end)];
end

cd(currentDir);
rmdir(tempDir,'s');
end


function ud = set_col_and_sc(ud,x1,x2)
%   sets the values ud.var_col and ud.sc appropriately, the order does not matter.

% Supported spacecraft
psc = {'mms','cluster'};

if nargin >= 2
    if ischar(x1)
        if ismember(x1,psc)
            ud.sc = x1;
        end
        ud.var_col = 1;
    else
        ud.var_col = x1;
    end
end

if nargin == 3
    if ischar(x2)
        if ismember(x2,psc)
            ud.sc = x2;
        end
    else
        ud.var_col = x2;
    end
end

if nargin == 1
    ud.var_col = 1;
end

end

% Position TSeries to Matricies
function ud = r_ts2mat(ud,R1,R2,R3,R4)
c_eval('ud.pos?=[R?.time.epochUnix,double(R?.data(:,1:3))];',1:4)
% add relative positions, for plot
R0 = (R1.resample(R1.time)+R2.resample(R1.time)+R3.resample(R1.time)+R4.resample(R1.time))/4;
c_eval('RR? = R?.resample(R0)-R0;',1:4)
c_eval('ud.relpos?=[RR?.time.epochUnix,double(RR?.data(:,1:3))];',1:4)
end


%% Get position functions
function ud = get_mms_pos(ud) 
%   Reads position data for MMS and stores it in ud.

% Gets position data for all data +-2 min for good measure
tstart = ud.var1.time(1);
tstop = ud.var1.time(end);

Tint = irf.tint(tstart,tstop);
R  = mms.get_data('R_gse',Tint);

T = irf_time(R.time.epoch,'ttns>epoch');

ud.pos1 = [T,R.gseR1];
ud.pos2 = [T,R.gseR2];
ud.pos3 = [T,R.gseR3];
ud.pos4 = [T,R.gseR4];

if ~is_pos_ok(ud)
    error('Unable to read MMS position data.')
end

end

function ud = get_c_pos(ud)
%   Reads position data for Cluster and stores it in ud.

irf.log('warning','Trying to read Cluster position')
tint = ud.tlim + ud.t_start_epoch + [-120, 120];
if ~is_pos_ok(ud) && exist('./mR.mat','file'),
    irf.log('warning','Trying to read position from mR.mat')
    load mR R1 R2 R3 R4;
    c_eval('R.R?=R?;');
end
if ~is_pos_ok(ud)
    irf.log('warning','Trying to read CAA files C?_CP_AUX_POSGSE_1M...')
    var = {'sc_r_xyz_gse__C1_CP_AUX_POSGSE_1M','sc_r_xyz_gse__C2_CP_AUX_POSGSE_1M','sc_r_xyz_gse__C3_CP_AUX_POSGSE_1M','sc_r_xyz_gse__C4_CP_AUX_POSGSE_1M'};
    ttt=c_caa_var_get(var,'mat','tint',tint);
    ud.pos1 = ttt{1}; ud.pos2 = ttt{2}; ud.pos3 = ttt{3}; ud.pos4 = ttt{4};
end
if ~is_pos_ok(ud) && exist('CAA/CL_SP_AUX','dir')==7,
    irf.log('warning','Trying to read CAA files CL_CP_AUX ...')
    R.R=irf_get_data('sc_r_xyz_gse__CL_SP_AUX','caa','mat');
    if ~isempty(R.R)
        c_eval('ud.pos?=irf_get_data(''sc_dr?_xyz_gse__CL_SP_AUX'',''caa'',''mat'');');
    end
end
if ~is_pos_ok(ud)
    irf.log('warning','Streaming s/c position from CAA');
    read_RV_from_caa_stream(tint);
end
if ~is_pos_ok(ud)
    disp('loading position from isdat. Connecting ...');
    tmin=min(tint);tmax=max(tint);
    for ic='1234'
        [tt,temp] = irf_isdat_get(['Cluster/' ic '/ephemeris/position'], tmin-30, tmax-tmin+30);
        ud.(['pos' ic])=[tt temp];
        fprintf('%s',ud.(['pos' ic]));
    end
    disp('');
end
if ~is_pos_ok(ud)
    irf.log('warning','!!! Could not obtain Cluster position data !!!');
    return
end
end