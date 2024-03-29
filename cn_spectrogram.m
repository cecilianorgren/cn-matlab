function [ax,cax] = cn_spectrogram(varargin)
% Slight modfification to irf_spectrogram so that the higher energy channel
% is seen. irf_spectrogram uses pcolor, which skipes the highest one. By
% using surf, one can instead specify the grid_edges.

[h,args,nargs] = axescheck(varargin{:});
if isempty(h)
  fig = get(groot,'CurrentFigure'); h = findobj(gcf,'type','axes'); h = h(end:-1:1);
end

%% Defaults
flagLog = true;            % want log10(data) dy default
f_multiplier = 1;          % default value using Hz units when units not specified, can be overwritten later if kHz makes labels more reasonable
fitColorbarLabel = true;   % fit font size of colorbar label to fit into axes size
showColorbar = true;

specrec = args{1};

%specrec = args{1};
  if ~isfield(specrec,'dt'), specrec.dt=[];end
  if ~isfield(specrec,'df'), specrec.df=[];end
  for iArgs = 2:numel(args)
    flagValue = args{iArgs};
    if ischar(flagValue)
      switch lower(flagValue)
        case 'donotfitcolorbarlabel'
          fitColorbarLabel = false;
        case 'log'
          flagLog = true;
        case 'lin'
          flagLog = false;
        case 'donotshowcolorbar'
          showColorbar = false;
        otherwise
          errStr= ['irf_spectrogram(), unknown flag:' flagValue];
          irf.log('critical',errStr);
          error('irf_spectrogram:unknown_flag',errStr);
      end
    end
  end


specrec.t = double(specrec.t);
specrec.f = double(specrec.f);

% specrec.p can be cell array
if iscell(specrec.p)
  ncomp=length(specrec.p);
elseif isnumeric(specrec.p)
  ncomp=1;
  specrec.p={specrec.p};
else
  disp('WARNING: cannot interpret input parameters in irf_spectrogram, returning.')
  return
end

ndata = length(specrec.t);
if ndata<1, if nargout>0, hout=h; end, return, end


% Initiate figure if handles not given
if isempty(h)
  h=irf_plot(ncomp,'newfigure');
end



% If H is specified, but is shorter than NCOMP, we plot just first
% length(H) spectra
for comp=1:min(length(h),ncomp)
  
  specrec.p{comp}(isnan(specrec.p{comp})) = NaN; % WHY is this done? NaN = NaN already.
  
  % get start time of plot
  ud = get(gcf,'userdata');
  ii = find(~isnan(specrec.t));
  if isfield(ud,'t_start_epoch')
    t_start_epoch = double(ud.t_start_epoch);
  elseif specrec.t(ii(1))> 1e8
    % Set start_epoch if time is in isdat epoch
    % Warn about changing t_start_epoch
    t_start_epoch = double(specrec.t(ii(1)));
    ud.t_start_epoch = t_start_epoch; set(gcf,'userdata',ud);
    irf.log('notice',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
  else
    t_start_epoch = double(0);
  end
  
  % Special case when we have only one spectrum
  % We duplicate it
  if ndata==1
    specrec.dt = double(.5/specrec.f(2));
    %		specrec.t = [specrec.t-dt; specrec.t+dt];
    %		specrec.p(comp) = {[specrec.p{comp}; specrec.p{comp}]};
  end
  if ~isfield(specrec,'f_unit') && ~isfield(specrec,'f_label') % if not specified assume units are Hz
    if max(specrec.f) > 2000 % check whether to use kHz
      specrec.f=specrec.f*double(1e-3);
      f_multiplier=1e-3;
      specrec.f_unit='kHz';
    else
      specrec.f_unit='Hz';
    end
    if ~isfield(specrec,'f_label')
      specrec.f_label=['f [' specrec.f_unit ']'];
    end
  end
  
  if min(size(specrec.f))==1, ff=double(specrec.f(:))';
  else
    ff=double(specrec.f);
  end % if f vector make it row vector
  tt=double(specrec.t_edges(:));
  pp=specrec.p{comp};
  ff=[specrec.f_edges; specrec.f_edges(end,:)];
 
  
  tag=get(h(comp),'tag'); % keep tag during plotting
  ud=get(h(comp),'userdata'); % keep tag during plotting
 
  tt = double(tt-t_start_epoch);
  tt = repmat(tt,1,size(ff,2));
  if ~flagLog || any(min(pp)<0) % spectra include negative values linear spectrogram
    surf(h(comp),tt,ff,ff*0,double(pp'))
  else
    surf(h(comp),tt,ff,ff*0,log10(double(pp)))
  end
  view(h(comp),[0 0 1]);
 
  set(h(comp),'tag',tag);
  set(h(comp),'userdata',ud);
  zoom_in_if_necessary(h(comp)); %
  
  shading(h(comp),'flat')
  set(h(comp),'TickDir','out')
  %check ylabel
  if ~isfield(specrec,'f_label')
    if ~isfield(specrec,'f_unit')
      specrec.f_unit='a.u.';
    end
    specrec.f_label=['[' specrec.f_unit ']'];
  end
  ylabel(h(comp),specrec.f_label)
  
  if showColorbar
    if isfield(specrec,'p_label')
      if isa(h(comp),'handle'), hcb = colorbar(h(comp)); % HG2
      else, hcb = colorbar('peer',h(comp));
      end
      drawnow
      posCb = get(hcb,'Position');
      posAx = get(h(comp),'Position');
      drawnow
      set(hcb,'TickDir','out','Position',...
        [posCb(1) posCb(2)+posCb(4)*0.05 posCb(3)*.75 posCb(4)*0.9])
      set(h(comp),'Position',[posAx(1) posAx(2) (posCb(1)-posAx(1))*0.97 posAx(4)])
      ylabel(hcb,specrec.p_label);
      if fitColorbarLabel
        irf_colorbar_fit_label_height(hcb);
      end
    end
  end
  if comp==min(length(h),ncomp), irf_timeaxis;
  else, set(h(comp),'XTicklabel','')
  end
end

if nargout>0, hout=h; end

function zoom_in_if_necessary(h)
ud=get(h,'userdata');
if isfield(ud,'zoom_x')
  disp('zooming in the updated plot')
  irf_zoom(h,'x',ud.zoom_x);
  if ud.zoom_x(1) > 1e8 && ud.zoom_x(1) < 1e10 % isdat epoch
    irf_timeaxis(h,'nolabel');
  end
end
