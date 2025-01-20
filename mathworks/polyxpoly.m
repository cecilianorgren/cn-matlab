%% [xi,yi] = my_polyxpoly(x1,y1,x2,y2,acc)
%
%%% --- INPUTS --- %%%
%   acc ---> accuracy of vectors, to what decimal place, default 1e-4
%   x1, y1 ---> mx1 or 1xm vectors
%   x2, y2 ---> nx1 or 1xn vectors
%
%%% --- OUTPUTS --- %%%
%   xi, yi ---> are a kx1 vectors repesenting all the intercepts
%
% m does not have to equal n
% m and n both have to be greater than or equal to 2
%
% acc is desired accuracy, will defualt to 1e-4, i.e. rounding to nearest
% 0.0001
%
% creating my own version of ...
%
% "[xi,yi] = polyxpoly(x1,y1,x2,y2) returns the intersection points of two 
% polylines in a planar, Cartesian system. x1 and y1 are vectors containing 
% the x- and y-coordinates of the vertices in the first polyline, and x2 
% and y2 contain the vertices in the second polyline. The output variables, 
% xi and yi, are column vectors containing the x- and y-coordinates of each 
% point at which a segment of the first polyline intersects a segment of 
% the second."
%
% From: http://www.mathworks.com/help/map/ref/polyxpoly.html
%
% x1, y1 must be vectors of the same size and x2, y2 as well, but both sets
% don't need to be.  They must at least be a vector with 2 cells (i.e. a
% line)
%
%%% --- EXAMPLE --- %%%
% >> x1 = [0,1:11,10]; 
% >> y1 = [0 1 1 5 -1 5 0 1 2 10 0 0 -1]; 
% >> x2 = [0,1:2:9,10,11,10]; 
% >> y2 = [0 1 10 0 7 4 0 0 -1];
% >> plot(x1,y1); hold on; plot(x2,y2,'r');
%
% See also length, mean, floor, meshgrid, reshape, find, sort, nargin,
% nargout
function [xi,yi] = polyxpoly(x1,y1,x2,y2,acc)
    % default inputs    
    if(nargin < 5)
        % do nothing this is correct
        acc = 1e-4; % accuracy ammount
    end
    % demmo mode
    if(nargout == 0)
        MODE = 'DEMO';
        if(nargin == 0)
            x1 = [-1,0,1,0,-1]; 
            y1 = [0,1,0,-1,0]; 
            x2 = -1:0.001:1;    % [sec] time
            f = 1;              % [Hz] frequency
            y2 = sin((2*pi*f).*x2);
        end
        x1_old = x1;
        y1_old = y1;
        x2_old = x2;
        y2_old = y2;
    else
        MODE = 'NORM';
    end
    % making sure x and y are same size vector
    test_vecs(x1,y1);
    test_vecs(x2,y2);
    % making all vector inputs in column vectors
    x1 = col_vec(x1,acc);
    y1 = col_vec(y1,acc);
    x2 = col_vec(x2,acc);
    y2 = col_vec(y2,acc);
    
    % finding slopes
    m1 = slope_vec(x1,y1);          % slopes of first line
    m2 = slope_vec(x2,y2);          % slopes of second line
    
    % finding intercepts
    b1 = intercept_vec(x1,y1,m1);   % intercepts of first line
    b2 = intercept_vec(x2,y2,m2);   % intercepts of second line
    
    % creating meshes
    [m1,m2] = meshgrid(m1,m2);
    [b1,b2] = meshgrid(b1,b2);
    [x1,x2] = meshgrid(x1,x2);
    [y1,y2] = meshgrid(y1,y2);
    
    % --- same points --- %
    x_test = x1 == x2;
    y_test = y1 == y2;
    same_points = x_test & y_test;
    IDX = find(same_points == 1);
    x_same = x1(IDX);
    y_same = y1(IDX);
    
    % --- intercepts possibilities --- %
    % x intercepts
    xi = (b2 - b1)./(m1 - m2);
    yi = (m2.*b1-m1.*b2)./(m2 - m1);
    
    % x intercepts
    [xi_1] = intercept_matrix(xi,x1);
    [xi_2] = intercept_matrix(xi,x2);
    xi_valid = xi_1 & xi_2; % logical locations of valid x intercepts
    
    % y intercepts
    [yi_1] = intercept_matrix(yi,y1);
    [yi_2] = intercept_matrix(yi,y2);
    yi_valid = yi_1 & yi_2; % logical locations of valid y intercepts
     
    % all valid possibilities, will have to check for x = 0 and y = 0 if
    % they are a possibility as well later
    xi = xi .* (xi_valid);
    [r,c] = find(isnan(xi));
    xi(r,c) = 0;
    yi = yi .* (yi_valid);
    [r,c] = find(isnan(yi));
    yi(r,c) = 0;
    
    % cleaning up xi and yi
    IDX = find(xi_valid | yi_valid == 1);
    xi = [x_same;xi(IDX)];
    yi = [y_same;yi(IDX)];
    
    if(strcmp(MODE,'DEMO'))
        figure();   
        plot(x1_old,y1_old);    hold on;   
        plot(x2_old,y2_old,'r');
        plot(xi,yi,'mo');       hold off;   
        xlabel('x');
        ylabel('y');
        title([num2str(length(xi)),' number intercept points found']);
        [xi yi] %#ok<NOPRT>
    elseif(strcmp(MODE,'NORM'))
        % sorting based on x column
        [xi,IDX] = sort(xi);
        yi = yi(IDX);
    else
        help('my_polyxpoly');
        error('something has gone wrong...');
    end
end
% finding intercepts matrix
function [xi_1] = intercept_matrix(xi,x1)
    test0 = (xi>x1(1:end-1,1:end-1)) & (xi<x1(2:end,2:end));  % assending
    test1 = (xi<x1(1:end-1,1:end-1)) & (xi>x1(2:end,2:end));  % decending
    xi_1 = test0 | test1;
end
% finding intercepts
function [b] = intercept_vec(x,y,m)
    % remember there is 1 less point in the slope vector because this is
    % the slope between two points and thus we'll have the same numbers of
    % intercepts, e.g. if x and y are a nx1 vector, m and b will be mx1
    % where m = n - 1
    b = y(1:end-1) - m.*x(1:end-1);
end
% finding slope vector
function [m] = slope_vec(x,y)
    m = diff(y) ./ diff(x); % slope vector
end
% creating all column vectors for work
function [x_new] = col_vec(x_old,acc)
    [r,c] = size(x_old);
    
    if((r == 1) && (c > 1))
        % got a row vector
        x_new = x_old'; % transpose to get column
    elseif((r > 1) && (c == 1))
        % got a column vector
        x_new = x_old;  % doing nothing
    else
        error(['Data passed in is not a vector [r,c] = [', ...
            num2str(r),',',num2str(c),'] or is 1x1']);
    end
    x_new = round(x_new ./ acc) .* acc;
end
% are x and y the same size?
function test_vecs(x,y)
    test = size(x) == size(y);
    test = floor(mean(test));
    if(test == 1)
        % they are the same time, we are all good
    else
        error('a pair of x and y vectors are not the same size');
    end
end