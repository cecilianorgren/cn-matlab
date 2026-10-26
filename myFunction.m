function myFunction(varargin)
    % 1. Create the inputParser object
    p = inputParser();
    
    % 2. Define positional arguments (value only, no name string needed)
    addOptional(p, 'firstValue', 1, @isnumeric);
    addOptional(p, 'secondValue', 10, @isnumeric);
    
    % 3. Parse the input cell array (varargin)
    parse(p, varargin{:});
    
    % 4. Access the parsed results
    val1 = p.Results.firstValue;
    val2 = p.Results.secondValue;
    
    fprintf('Value 1: %d, Value 2: %d\n', val1, val2);
end
