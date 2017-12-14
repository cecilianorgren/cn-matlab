function out = cmap(map)
% cmap = cn.cmap(choice)
%   Choice can be 'islands', 'pinkpurple', 'bluered', 'bluered2',
%   'bluered3','fall' or whatever (then the default is returned).

%path = '/Users/Cecilia/MATLAB/colormaps/';
path = datastore('colormaps','path');
switch map
    case 'white_jet'
        out = load([path 'white_jet.mat']);
        out = out.cmap;
    case 'white_parula'
        out = load([path 'white_parula.mat']);
        out = out.cmap;
    case 'discrete_jet'
        out = load([path 'discrete_jet.mat']);
        out = out.cmap;
    case 'islands'
        out = load([path 'cmap_hc.mat']);
        out = out.cmap;
    case 'bluepink'
        out = load([path 'cmap_bluepink.mat']);
        out = out.cmap;
    case 'purplepink'
        out = load([path 'cmap_purplepink.mat']);
        out = out.cmap;        
    case 'pinkpurple'
        out = load([path 'cmap_pp.mat']);
        out = out.cmap;
    case 'poynting'
        out = load([path 'irf_poynting.mat']);
        out = out.cmap;
    case 'bluered'
        out = load([path 'cmap_bluered.mat']);
        out = out.cmap;
    case 'bluered2'
        out = load([path 'cmap_redblue2.mat']);
        out = out.cmap;
    case 'bluered3'
        out = load([path 'cmap_redblue3.mat']);
        out = out.cmap;
    case 'whitered'
        out = load([path 'cmap_whitered.mat']);
        out = out.cmap;
    case 'fall'
        out = load([path 'cmap_fall.mat']);
        out = out.cmap;
    otherwise
        disp('Didn''t find the colormap. Returning default.');
        out = load([path 'cmap_hc']);
end
        
        
        