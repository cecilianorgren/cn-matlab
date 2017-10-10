function out = R(varargin)
% Calculates R or n1/n2.
%   if input is R, n/n1 is output
%   if input is n1 and n2, output is R

switch numel(varargin)
    case 1 % input is R
        R = varargin{1};
        n1 = R;
        n2 = 1-R;
        out = n1/n2;        
    case 2 % input is n1,n2
        n1 = varargin{1};
        n2 = varargin{2};
        R = n1/(n1+n2);
        out = R;        
end
        
disp(['R = ' num2str(R) ', n1 = ' num2str(n1) ', n2 = ' num2str(n2) ', ntot = ' num2str(n1+n2) ', n1/n2 = ' num2str(n1/n2)])        