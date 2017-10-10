function out = abs(in,flag)
% always assume first column is time
% this can be overridden by giving a flag = 1

if nargin > 1 % flag is given
    if flag == 1 % first column is not time series
        in = in;
    end
else
    in = in(:,[2:end]); % take away time column
end

out = sqrt(sum(in.^2,2));
        
        
