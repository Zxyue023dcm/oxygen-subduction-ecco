function [CURL,gXy,gYx] = curldffs_region(datax,distx,datay,disty)
%CURLDFFS_REGION Calculate curl terms on a two-dimensional regional grid.
%   Longitude is the first dimension and latitude is the second.

lon_l=distx;
lat_l=disty;
gXy1=nan(size(datax));
gYx1=nan(size(datay));

for p=1:size(lat_l,2)
    hx=lon_l(:,p);
    hy=lat_l(1,:);
    [py,~]=gradient(datax,hy,hx);
    [~,qx]=gradient(datay,hy,hx);
    gXy1(:,p)=py(:,p);
    gYx1(:,p)=qx(:,p);
end

gXy=gXy1;
gYx=gYx1;
CURL=gYx-gXy;
end

