function varargout = max(A,varargin)
% find outliers
if nargin
    N = varargin{1};
    orgA = A;
    A=abs(A);
    maxes = [];
    rows=[];
    cols=[];
    for k=1:N
        [row col] = find(A==max(max(A)));
        maxes = [maxes orgA(row,col)];
        rows=[rows row];
        cols=[cols col];
        A(A==max(max(A)))=NaN;
    end
    out=maxes;
    varargout={maxes,rows,cols};
else    
    maxInd = find(abs(A)==max(abs(A)));
    out = A(maxInd);
    out = max(out); % choose the positive if they have exact same value
    varargout{1} = out;
end