function varargout = find_peaks_iteratively_mp(MP,dv,N,w)
% Find peaks iteratively in 3D matrix. Macroparticle structure as input.
%   1. Make F.
%   2. Find max.
%   3. Remove w surrounding point: F(i-w:1:i+w,...,...) = 0. This
%      corresponds to a velocity interval of +- w*dv-
%   4. Repeat 2-3 N times.

if numel(dv) == 1; dv = repmat(dv,N,1); end
V = [MP.vx, MP.vy, MP.vz];

vabs = vecnorm(V,2,2);
vmax = max(vabs);
vlim = ceil(vmax/max(dv))*max(dv); % round up to nerest integer of dv

vec{1} = -vlim:dv(1):vlim;
vec{2} = -vlim:dv(2):vlim;
vec{3} = -vlim:dv(3):vlim;

[X,Y,Z] = ndgrid(vec{:});
XYZ = [X(:) Y(:) Z(:)];

dn = MP.df.*MP.dv;
[n,edges,mid] = histcn(V,vec{:},'AccumData',dn);
F = n/prod(dv);

%f_smooth = imgaussfilt3(F, 1);
%f_smooth = F;
%isosurface(f_smooth, 0.3*max(f_smooth(:)))

F = imgaussfilt3(F, 1);

% Find peaks
i = zeros(N,1);
j = zeros(N,1);
k = zeros(N,1);
Fiter = F;
for ipeak = 1:N
  [maxval, idx] = max(Fiter(:));
  [i(ipeak),j(ipeak),k(ipeak)] = ind2sub(size(Fiter), idx);
  ix = i(ipeak)+(-w:1:w);
  iy = j(ipeak)+(-w:1:w);
  iz = k(ipeak)+(-w:1:w);
  ix(ix<1) = [];
  iy(iy<1) = [];
  iz(iz<1) = [];
  ix(ix>size(F,1)) = [];
  iy(iy>size(F,2)) = [];
  iz(iz>size(F,3)) = [];
  Fiter(ix,iy,iz) = 0;
end

if nargout == 3
  varargout{1} = i;
  varargout{2} = j;
  varargout{3} = k;
elseif nargout == 4
  varargout{1} = i;
  varargout{2} = j;
  varargout{3} = k;
  out.F = F;
  out.vec = edges;
  varargout{4} = out;
end