function [lon_l,lat_l] = distcoordinate(lon,lat)
%DISTCOORDINATE Build metric coordinates from longitude and latitude.
%   LON and LAT are column vectors. LON_L and LAT_L are in metres.

nlon=size(lon,1);nlat=size(lat,1);
distlat=zeros(nlon,nlat-1);
distlon=zeros(nlon-1,nlat);


for i=1:nlat-1
    distlat(:,i)=sw_dist([lat(i,1) lat(i+1,1)], ...
        [lon(1,1) lon(1,1)],'km');
end
for i=1:nlat
    distlon(:,i)=sw_dist([lat(i,1) lat(i,1)], ...
        [lon(1,1) lon(2,1)],'km');
end
lat_l=cumsum([zeros(nlon,1) distlat],2)*1000;
lon_l=cumsum([zeros(1,nlat);distlon],1)*1000;

end

