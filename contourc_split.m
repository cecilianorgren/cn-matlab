function [x,y] = contourc_split(varargin)

cline = contourc(varargin{:});

iLine = 0;
while not(isempty(cline))
  iLine = iLine + 1;
  C(iLine) = cline(1,1);
  nSeg = cline(2,1);
  x_tmp = cline(1,1+(1:nSeg));  
  y_tmp = cline(2,1+(1:nSeg));  
  cline(:,1:(nSeg+1)) = [];
  
  x{iLine} = x_tmp;
  y{iLine} = y_tmp;
end