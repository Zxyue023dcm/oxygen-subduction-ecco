function [DIV,gX,gY] = divdffs(datax,distx,datay,disty,method)
%DIVDFFS Calculate horizontal divergence on a two-dimensional grid.
%   Longitude is the first dimension and latitude is the second. METHOD is
%   either 'backward' or 'gradient'.

lon_l=distx;
lat_l=disty;

if strcmp(method,'backward')
    delX=diff(datax,1,1);
    delx=diff(lon_l,1,1);
    gX=delX(:,2:end)./delx(:,2:end);

    delY=diff(datay,1,2);
    dely=diff(lat_l,1,2);
    gY=delY(2:end,:)./dely(2:end,:);

elseif strcmp(method,'gradient')
    guHmlx1=nan(size(datax));
    guHmly1=nan(size(datay));
    for p=1:size(lat_l,2)
        hx=lon_l(:,p);
        hy=lat_l(1,:);
        [~,px]=gradient(datax,hy,hx);
        [qy,~]=gradient(datay,hy,hx);
        guHmlx1(:,p)=px(:,p);
        guHmly1(:,p)=qy(:,p);
    end
    gX=guHmlx1(2:end-1,:);
    gY=guHmly1(2:end-1,:);
else
    error('divdffs:UnknownMethod', ...
        'Method must be ''backward'' or ''gradient''.');
end

DIV=gX+gY;
end
