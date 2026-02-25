function ax = cn_plot3d(varargin)
% CN_PLOT3D Plot 3D arrows.
%
% Example cn_plot3d(ax,origin,vector,str,linespec)
% 
%   ax - axis   
%   origin - Kx3 matrix where each row is speciefiec the origin, or 1x3
%            vector where, in that case same origin for all vectors
%   vector - Lx3 matric where each row specifies the vector
%   str - Mx1 cell or row vector where each cell contains vector 
%         description
%   linespec - Nx1 cell or row vector where each cell contains line
%              specifications, for example 'b'. If the 
%
%   L>=K,M,N ! 
%           If L>K, the rest of the entries is filled up with the last
%           entry for the origin. If isempty(origin), then origin=[0 0 0].
%           If L>M, the rest of the entries is filled up with ' '. If
%           isempty(str), then str=' '.
%           If L>N, the rest of the entries is filled up with the last entry
%           for the line specification. If isempty(str), then linespec=' '.
%
%   The data can at anytime be accessed through get(gcf,'userdata').

[ax,args,nargs] = axescheck(varargin{:});
original_args=args;

% Origins
origin=args{1};
no=size(origin,1);

% Directions
vector=args{2};
nv=size(vector,1);

% Fill up origins if less than number of vectors are given
if no == 0 % empty entry
    origin=zeros(nv,3);
elseif no < nv % fill up rest with last entry
    origin=[origin;repmat(origin(end,:),nv-no,1)];
end

% Scalings
scale=args{3};
ns=size(scale,1);

% Fill up scalings if less than number of vectors are given
if ns == 0 % empty entry, make scale=1
    scale=ones(nv,3);
elseif ns < nv % fill up rest with last entry
    scale=[scale;repmat(scale(end,:),nv-no,1)];
end

% Labels
str=args{4};
ns=max(size(str,1));

if ns == 0 % empty input
    clear str;
    [str{1:nv}]=deal(' ');
else % otherwise
    if ischar(str) % convert to cell array
        str=cellstr(str);
    end
    [str{end+1:nv}]=deal(' ');
end

% Linespec
linespec=args{5};
nl=max(size(linespec));

if nl == 0 % empty input
    clear linespec;
    [linespec{1:nv}]=deal(' ');
else % otherwise
    if ischar(linespec) % convert to cell array
        linespec=cellstr(linespec);
    end
    [linespec{end+1:nv}]=deal(linespec{end});
end

if nargs>4
    flag_remove=1;
    remove=args{end}; % cell array with
    if ischar(remove)
        remove=cellstr(remove);
    end
    nr=max(size(remove));
else
    flag_remove=0;
end

% initialize figure, remove data and replot
if isempty(ax) 
    figure;
    plot(1);
    ax=gca;
    ud=get(gcf,'userdata');
    ud.origin=origin;
    ud.vector=vector;
    ud.scale=scale;
    ud.str=str;
    ud.linespec=linespec; 
else
    ud=get(gcf,'userdata');
    ud.origin=[ud.origin;origin];
    ud.vector=[ud.vector;vector];
    ud.scale=[ud.scale;scale];
    ud.str=[ud.str,str];
    ud.linespec=[ud.linespec,linespec];
end
nud=max(size(ud.vector));
nplot=nud+1-nv:nud;
if flag_remove
    ri=[];
    for k=1:nr
        ri=[ri,find(strcmp(ud.str,remove{k}))];        
    end
    nri=length(ri);
    ud.origin(ri,:)=[];
    ud.vector(ri,:)=[];
    ud.scale(ri,:)=[];
    ud.str(ri)=[];
    ud.linespec(ri)=[];
    nud=max(size(ud.str)); 
    nplot=1:nud;
end

for k=nplot
    quiver3(ax,ud.origin(k,1),ud.origin(k,2),ud.origin(k,3),...
               ud.vector(k,1),ud.vector(k,2),ud.vector(k,3),...
               ud.scale(k),...
               ud.linespec{k})
    text(ud.origin(k,1)+ud.vector(k,1)*ud.scale(k),...
         ud.origin(k,2)+ud.vector(k,2)*ud.scale(k),...
         ud.origin(k,3)+ud.vector(k,3)*ud.scale(k),...
         ud.str{k})
    hold(gca,'on')
end
axis equal;
set(gcf,'userdata',ud)