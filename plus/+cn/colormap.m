function cmap1=irf_colormap(varargin)
% IRF_COLORMAP return colormap by name or apply and freeze the colormap
%
% IRF_COLORMAP(colormap_name)
%  Colormap_names:
%       'standard'  - (default), same as 'space','cmap' (commonly used showing space data)
%       'poynting'  - white in center and blue/green for negative and red/black for positive values
%       'poynting_gray'  - gray in center and blue/green for negative and red/black for positive values
%
% IRF_COLORMAP(AX,colormap_name) - apply colormap to axis AX
%

% $Id$

[ax,args,nargs] = axescheck(varargin{:});
% double check for "double handle" in first argument
if isempty(ax)
    possible_handles = varargin{1};
    for ii = 1:numel(varargin{1})
        if ishghandle(possible_handles,'axes')
            ax(ii) = possible_handles(ii);
        end
    end
    [~,args,nargs] = axescheck(varargin{2:end});
end


if nargs == 0, % show only help
    help irf_colormap;
    return
end

% check which axis to apply
if isempty(ax), 
    axes(gca);
else
    axes(ax(1));
end

colormap_name=args{1};

load caa/cmap.mat % default map
if nargs > 0,
    switch lower(colormap_name)
        case 'poynting'
            it=0:.02:1;it=it(:);
            cmap=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
        case {'poynting_grey','poynting_gray'}
            it=0:.02:1;it=it(:);
            cmap=[ [0*it flipud(it) it];...
				[it*.8 it*.8 it*0+1];...
				[it*0+1 flipud(it*.8) flipud(it*.8)];...
				[flipud(it) 0*it 0*it]]; 
			clear it;
        case 'solo'
            it=0:.02:1;it=it(:);
            cmap=[ [it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
    end
end

if nargout == 0, % freeze others and apply the colormap to given axes 
    all_axes = findall(gcf,'type','axes');
    all_plot = findall(gcf,'type','axes','tag','');
    all_cbar = findobj(gcf,'tag','Colorbar');       
    nonactive_plot = all_plot;
    nonactive_cbar = all_cbar;
    for ii = 1:numel(ax)
        active_plot = ax(ii);
        active_cbar = findobj(gcf, 'Type', 'axes', 'Tag', 'Colorbar','Axes', active_plot);    
        nonactive_plot(find(nonactive_plot==active_plot)) = [];
        nonactive_cbar(find(nonactive_cbar==active_cbar)) = [];
    end   
    
    for ii = 1:numel(nonactive_plot); freezeColors(nonactive_plot(ii)); end      
    for ii = 1:numel(nonactive_cbar), % workaround cbfreeze bug that cbfreeze removes cblabel
        hcb = nonactive_cbar(ii);
        hy=get(hcb,'ylabel');
        ylabel_string=get(hy,'string');
        ylabel_fontsize=get(hy,'fontsize');
        cbar_position = get(hcb,'position');
        new_hcb = cbfreeze(hcb);
        new_hy=get(new_hcb,'ylabel');        
        set(new_hy,'string',ylabel_string,'fontsize',ylabel_fontsize);
        set(new_hcb,'position',cbar_position);
    end
    colormap(active_plot,cmap);   
elseif nargout == 1, % only return colormap
    cmap1=cmap;
end

