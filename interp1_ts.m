function out = interp1_ts(TS,varargin)

dimData = ndims(TS.data)-1; % remove time dim
sizeData = size(TS.data); 
sizeData = sizeData(2:end); % remove time size

data = TS.data;
new_data = data;

for iDim = 1:dimData
  for iSize = 1:sizeData(iDim)
    new_data(:,iSize) = smooth(squeeze(data(:,iSize)),varargin{:});
  end
end

out = TS;
out.data = new_data;