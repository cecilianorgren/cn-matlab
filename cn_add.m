function sum=cn_add(a,at,b,bt,varargin)
varargin
if size(varargin,2)==0
    sum=a*at+b*bt;
elseif size(varargin,2)==2
    sum=cn_add(1,cn_add(a,at,b,bt),varargin{1},varargin{2});
else
    sum=cn_add(1,cn_add(a,at,b,bt),varargin{1},varargin{2},varargin{3:end});
end


if 0
if 0% mod(nargin,2)
    disp('Uneven number of input arguments.');
    return;
end
%varargin
iscell(varargin)
sum=varargin{1}+varargin{2}
varargout={sum, varargin{3:end}}
if nargin<2
    return;
end
if size(varargout,2)==2
    return
end
cn_add(varargout);
%cn_add(varargin{3:end})
%for k=1:nargin/2
    %sum=cn_add(varargin
end