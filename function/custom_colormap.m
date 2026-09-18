function cmap = custom_colormap(colors,varargin)
%CUSTOM_COLORMAP Interpolate a colormap from user-defined RGB nodes.
%   CMAP=CUSTOM_COLORMAP(COLORS) returns 256 colors and applies the map to
%   the current axes. COLORS is an N-by-3 array in [0,1] or [0,255].
%
%   Name-value options:
%     'n'         Number of output colors (default 256)
%     'apply'     'on' or 'off' (default 'on')
%     'positions' Strictly increasing node positions spanning [0,1]

p=inputParser;
addRequired(p,'colors',@(x) isnumeric(x) && size(x,2)==3 && size(x,1)>=3);
addParameter(p,'n',256,@(x) isscalar(x) && x>0 && fix(x)==x);
addParameter(p,'apply','on',@(x) (ischar(x) || isstring(x)) && ...
    any(strcmpi(x,{'on','off'})));
addParameter(p,'positions',[],@(x) isnumeric(x) && isvector(x) && ...
    numel(x)==size(colors,1));
parse(p,colors,varargin{:});

n=p.Results.n;
applyFlag=p.Results.apply;
pos=p.Results.positions;

colors=double(colors);
if max(colors(:))>1
    colors=colors/255;
end
colors=max(0,min(1,colors));

k=size(colors,1);
if isempty(pos)
    pos=linspace(0,1,k);
else
    pos=pos(:)';
    if any(diff(pos)<=0)
        error('custom_colormap:InvalidPositions', ...
            'Positions must be strictly increasing.');
    end
    if pos(1)~=0 || pos(end)~=1
        warning('custom_colormap:RescaledPositions', ...
            'Positions were rescaled to span [0,1].');
        pos=(pos-pos(1))/(pos(end)-pos(1));
    end
end

xi=linspace(0,1,n);
r=interp1(pos,colors(:,1),xi,'linear');
g=interp1(pos,colors(:,2),xi,'linear');
b=interp1(pos,colors(:,3),xi,'linear');
cmap=[r(:),g(:),b(:)];

if strcmpi(applyFlag,'on')
    colormap(cmap);
end
end
