function varargout = fit_gaussian(tsE,tsL,varargin)
% Returns the 

if isnumeric(tsL) && numel(tsL == 1) % tsL is a speed
  L = (tsE.time - tsE.time(1))*tsL;
elseif isa(tsL,'TSeries')
  tsL = tsL.resample(tsE);
  L = tsL.data;
end
E = tsE.data;

if isempty(varargin) % initial guess of half width
  imax = find(E == max(E));
  imin = find(E == min(E));
  guessL = 0.5*abs(L(imax)-L(imin));
else 
  guessL = varargin{1};
end

guessL;
varargout = {guessL};