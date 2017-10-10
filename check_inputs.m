function check_inputs(A, B, varargin)
minargs=2;  maxargs=5;

% Number of inputs must be >=minargs and <=maxargs.
nargchk(nargin,minargs, maxargs)

fprintf('Received 2 required, %d optional inputs.\n\n', ...
    size(varargin, 2))