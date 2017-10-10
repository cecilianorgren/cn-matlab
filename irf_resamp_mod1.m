function out = irf_resamp(x,y,varargin)
%IRF_RESAMP   Resample X to the time line of Y
%
% if sampling of X is more than two times higher than Y, we average X,
% otherwise we interpolate X.
%
% out = irf_resamp(X,Y,[METHOD],['fsample',FSAMPLE],['window',WIN],
%                      ['thresh',THRESH],['median'],['max'])
% method - method of interpolation 'spline', 'linear' etc. (default 'linear')
%          if method is given then interpolate independant of sampling
% thresh - points above STD*THRESH are disregarded for averaging
% fsample - sampling frequency of the Y signal, 1/window
% window - length of the averaging window, 1/fsample
% median - use median instead of mean when averaging
% max    - return max within each averaging window, rather than mean
%
% See also INTERP1

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,8)

have_options = 0;
args = varargin; 
if nargin > 2, have_options = 1; end

% Default values that can be override by options
sfy = [];
thresh = 0;
method = '';
flag_do='check'; % if no method check if interpolate or average
median_flag=0;
max_flag=0;
%flag_output = 'tseries'; % default

while have_options
	l = 1;
	switch(lower(args{1}))
		case {'nearest','linear','spline','pchip','cubic','v5cubic'}
			method = args{1};
			flag_do='interpolation'; % if method is given do interpolation
		case 'method'
			if length(args)>1
				if ischar(args{2})
					method = args{2};
					l = 2;
					flag_do='interpolation'; % if method is given do interpolation
				else irf.log('critical','wrongArgType : METHOD must be numeric')
				end
			else irf.log('critical','wrongArgType : METHOD value is missing')
			end
		case {'fs','fsample'}
			if length(args)>1
				if ~isempty(sfy)
          msgS = 'FSAMPLE/WINDOW already specified';
          irf.log('critical',msgS), error(msgS)
				end
				if isnumeric(args{2})
					sfy = args{2};
					l = 2;
				else irf.log('critical','wrongArgType : FSAMPLE must be numeric')
				end
			else irf.log('critical','wrongArgType : FSAMPLE value is missing')
			end
		case {'win','window'}
			if length(args)>1
				if ~isempty(sfy)
					msgS = 'FSAMPLE/WINDOW already specified';
          irf.log('critical',msgS), error(msgS)
				end
				if isnumeric(args{2})
					sfy = 1/args{2};
					l = 2;
				else irf.log('critical','wrongArgType : WINDOW must be numeric')
				end
			else irf.log('critical','wrongArgType : WINDOW value is missing')
			end
		case {'thresh','threshold'}
            if length(args)>1
                if isnumeric(args{2})
                    thresh = args{2};
                    l = 2;
                else irf.log('critical','wrongArgType : THRESHOLD must be numeric')
                end
            else irf.log('critical','wrongArgType : THRESHOLD value is missing')
            end
        case 'median'
            median_flag=1;
        case 'max'
            max_flag=1;
		otherwise
			irf.log('warning',['Skipping parameter ''' args{1} ''''])
			args = args(2:end);
	end
	args = args(l+1:end);
	if isempty(args), break, end
end

% return in case inputs are empty
% changed numel(x) == 0 to isempty(x), since it also works on TSeries
if isempty(x) || isempty(y)
    out=[];
    irf.log('warining','Some of input is empty, returning empty ouput');
    return;
end

% Check format of input
% x, old time and data
if isa(x,'TSeries')
  flag_output = 'tseries';
  oldData = x.data;
  oldTime = x.time.epoch;
elseif isa(x,'numeric') % old format [time data]
  flag_output = 'mat';
  oldData = x(:,2:end);
  oldTime = x(:,1);
else
  msgS = 'Cannot recognize input x (old data), it is neither a TSeries nor a numeric array';
  irf.log('critical',msgS), error(msgS)
end

% y, new time
if isa(y,'struct'),
  if isfield(y,'t'), t=y.t; t=t(:);
  else
    msgS = 'Input is structure without time field';
    irf.log('critical',msgS), error(msgS)
  end
elseif isa(y,'TSeries')
  type_epoch = class(y);
  t = y.time.epoch;
elseif isa(y,'GenericTimeArray')
  type_epoch = class(y);
  t = y.epoch;  
elseif isa(y,'numeric')
  if size(y,2)==1, t = y(:);  % y is only time
  else t = y(:,1); t = t(:);  % first column of y is time
  end
else 
  msgS = 'Cannotrecognize input y (new time), it is neither a TSeries nor a numeric array';
  irf.log('critical',msgS), error(msgS)
end

% Same timeline - no need to do anything
% old data already divided into time and data
if length(oldTime)==length(t) && all(oldTime==t)
  irf.log('notice','New and old timelines are identical - no resampling needed')
  out = x; 
  return
end

% if X (oldData) has only one point, this is a trivial case and we 
% return directly there is only one number
% the single data point is duplicated into the new timeline
if numel(oldTime) == 1,       
  out = construct_output(t,repmat(oldData,numel(t),1,1));
  return
end

ndata = length(t);
if strcmp(flag_do,'check'), % Check if interpolation or average
	if ndata>1 
		% If more than one output time check sampling frequencies
		% to decide interpolation/average
		
		% Guess samplings frequency for Y
    % This becomes problematic when having time in nanoseconds, the
    % division with dt becomes = 0. The format int64 does not allow for
    % floating numbers. I make it a double
		if isempty(sfy)
			sfy1 = (1/double(t(2) - t(1)));
			if ndata==2, sfy = sfy1;
			else
				not_found = 1; cur = 3; MAXTRY = 10;
				while (not_found && cur<=ndata && cur-3<MAXTRY)
					sfy = 1/(t(cur) - t(cur-1));
                    if abs(sfy-sfy1)<sfy*0.001
                        not_found = 0;
                        sfy = (sfy+sfy1)/2;
                        break
                    end
                    sfy1=sfy;
					cur = cur + 1;
				end
				if not_found
					sfy = sfy1;
					irf.log('warning',	sprintf(...
						'Cannot guess sampling frequency. Tried %d times',MAXTRY));
				end
			end
			clear sfy1
		end
		
		if length(oldTime(:,1))/(oldTime(end,1) - oldTime(1,1)) > 2*sfy
			flag_do='average';
			irf.log('warning','Using averages in irf_resamp.');
		else
			flag_do='interpolation';
		end
	else
		flag_do='interpolation';  % If one output time then do interpolation
	end
end

if strcmp(flag_do,'average')
  % This part is not made suitable for tensor 2 data
    dt2 = .5/sfy; % Half interval
    if median_flag || max_flag || (exist('irf_average_mx','file')~=3)
        if (~median_flag && ~max_flag), irf.log('warning','cannot find mex file, defaulting to Matlab code.')
        end
        newData = zeros(ndata,size(oldData,2),size(oldData,3));
        %out(:,1) = t;
        for j=1:ndata
            ii = find(t(:,1) <=  t(j) + dt2 & t(:,1) >  t(j) - dt2);
            if isempty(ii), newData(j,:) = NaN;
            else
                if thresh % Throw away points above THERESH*STD()
                    sdev = std(oldData(ii,:));
                    mm = mean(oldData(ii,:));
                    if any(~isnan(sdev))
                        for k=1:length(sdev)
                            if ~isnan(sdev(k))
                                kk = find( abs( oldData(ii,k) -mm(k) ) <= thresh*sdev(k));
                                %disp(sprintf(...
                                %	'interval(%d) : disregarding %d 0f %d points',...
                                %	j, length(ii)-length(kk),length(ii)));
                                if ~isempty(kk)
                                    if median_flag, newData(j,k) = median(oldData(ii(kk),k));
                                    elseif max_flag, newData(j,k) = max(oldData(ii(kk),k));
                                    else newData(j,k) = mean(oldData(ii(kk),k));
                                    end
                                end
                            else newdata(j,k) = NaN;
                            end
                        end
                    else newData(j,:) = NaN;
                    end
                else
                    if median_flag, newData(j,:,:) = median(oldData(ii,:,:));
                    elseif max_flag, newData(j,:,:) = max(oldData(ii,:,:));
                    else newData(j,:,:) = mean(oldData(ii,:,:));
                    end
                end
            end
        end
        out = construct_output(t,newData);
    else
        if ( (t(1,1) > t(end) + dt2) || (t(end,1) <= t(1) - dt2) )
            irf.log('warning','Interval mismatch - empty return')           
            out = construct_output([],[]);
        else          
            % input to irf_average_mx needs to be double
            if strcmp(flag_output,'tseries') && x.tensorOrder == 2 
              % must divide tensor into vectors
              tmpData = nan(size(oldData));
              for ii = 1:size(oldData,2)
                %irf_resamp([tData squeeze(dataTmp(:,ii,:))], newTimeTmp, varargin{:});
                % Maybe insert warning here that data is converted to double, if it is 
                % not already double
                newTmpData = irf_average_mx([double(oldTime-oldTime(1)) squeeze(oldData(:,ii,:))],double(t-oldTime(1)),double(dt2),thresh);              
                tmpData(:,ii,:) = newTmpData(:,2:end);
              end   
              newData = tmpData;
              %outTmp = irf_average_mx(oldData,double(t-t(1)),double(dt2),thresh);              
            else
              % Maybe insert warning here that data is converted to double, if it is 
              % not already double
              time1 = double(oldTime-oldTime(1));
              time2 = double(t-oldTime(1));
              tmpData = irf_average_mx([time1 double(oldData)],time2,double(dt2),thresh); 
              %tmpData = irf_average_mx([double(oldTime-oldTime(1)) double(oldData)],double(t-oldTime(1)),double(dt2),thresh); 
              newData = tmpData(:,2:end);
            end
            out = construct_output(t,newData);
        end
    end
elseif strcmp(flag_do,'interpolation'),
  if nargin < 3 || isempty(method), method = 'linear'; end

  % If time series agree, no interpolation is necessary.
  if size(oldData,1)==size(t,1), if oldTime==t(:,1), out = construct_output(oldTime,oldData); return, end, end

  %out = [t interp1(x(:,1),x(:,2:end),t,method,'extrap')]; 
  
  % Maybe insert warning here that data is converted to double, if it is 
  % not already double
  out = construct_output(t,interp1(double(oldTime-oldTime(1)),double(oldData),double(t-oldTime(1)),method,'extrap'));
end

function out = construct_output(t,newData)
  % Instead of checking each time is the input was a TSeries or not, a flag
  % is defined once, and then this functions called and constructs the
  % appropriate output.
  if isempty(t) || isempty(newData)
    out = [];
  end
  switch flag_output
    case 'tseries'
      switch lower(type_epoch)
        case 'epochtt', genericTime = EpochTT(t);
        case 'epochunix',  genericTime = EpochUnix(t);
      end
      % TODO: pass on all the metadata
      % Did not remember the good way to do this   
      out = x.clone(genericTime,newData);
    case 'mat'
      out = [t newData];
  end

end % end of construct_output()

end % end of irf_resamp()
