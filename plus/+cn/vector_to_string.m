function out = vector_to_string(vector,varargin)
% CN.VECTOR_TO_STRING Transforms vector to string.
%   CN.VECTOR_TO_STRING(vector,format)
%       vector - 1x3 vector
%       format - specify format to be used in num2str
%                OR is number of decimals, 1 is default
%   Example:
%       title(gca,['x=' cn.vector_to_string(x) '  y=' cn.vector_to_string(y) '  z=' cn.vector_to_string(y)])

if nargin>1
    if ischar(varargin{1})
        format=varargin{1};
    elseif isnumeric(varargin{1})
        format=['%.',num2str(varargin{1},'%.0f'),'f'];
    else
        error('Format in unknown format ;). Try again.')
        return;
    end
else
    format='%.1f';
end
    
    
out = ['[',num2str(vector(1),format),' ',...
           num2str(vector(2),format),' ',...
           num2str(vector(3),format),']'];

           