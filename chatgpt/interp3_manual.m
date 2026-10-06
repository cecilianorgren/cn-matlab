function data_new = interp3_manual(D1,D2,D3,data,newD1,newD2,newD3,method)
%INTERP3_MANUAL Sequential 3-D interpolation using interp1.
%
%   data_new = interp3_manual(D1,D2,D3,data,newD1,newD2,newD3,method)
%
%   Interpolates data defined on the grid (D1,D2,D3) onto
%   (newD1,newD2,newD3), applying interp1 sequentially along
%   each dimension.
%
%   method can be e.g. 'linear', 'pchip', or 'spline'.

if nargin < 8
  method = 'pchip';
end

%% Interpolate along dimension 1
tmp = interp1(D1, data, newD1, method);

%% Interpolate along dimension 2
tmp = permute(tmp,[2 1 3]);
tmp = interp1(D2, tmp, newD2, method);
tmp = permute(tmp,[2 1 3]);

%% Interpolate along dimension 3
tmp = permute(tmp,[3 1 2]);
tmp = interp1(D3, tmp, newD3, method);
data_new = permute(tmp,[2 3 1]);

end
