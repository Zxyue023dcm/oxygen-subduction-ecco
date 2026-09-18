function [field,lon,lat,dep] = read_ecco_monthly(filePath,variableName)
%READ_ECCO_MONTHLY Read one monthly ECCO field on the 0.5-degree grid.
%   FIELD is returned as longitude-by-latitude-by-depth-by-time.

longitude=double(ncread(filePath,'longitude'));
latitude=double(ncread(filePath,'latitude'));
dep=double(ncread(filePath,'Z'));
field=double(ncread(filePath,variableName));

fieldSize=size(field);
lonDim=find(fieldSize==numel(longitude),1);
latDim=find(fieldSize==numel(latitude),1);
depDim=find(fieldSize==numel(dep),1);
assert(~isempty(lonDim) && ~isempty(latDim) && ~isempty(depDim) && ...
    numel(unique([lonDim latDim depDim]))==3, ...
    'Cannot identify longitude, latitude, and depth dimensions in %s.',filePath);

order=[lonDim latDim depDim setdiff(1:ndims(field), ...
    [lonDim latDim depDim],'stable')];
field=permute(field,order);
field=reshape(field,[numel(longitude) numel(latitude) numel(dep) ...
    numel(field)/(numel(longitude)*numel(latitude)*numel(dep))]);
assert(size(field,4)==1,'Expected one monthly time record in %s.',filePath);

[lat,lon]=meshgrid(latitude,longitude);
dep=dep(:);
end
