%% Prepare GOBAI-O2 fields on a 0-360 degree longitude grid
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
inputDir=fullfile(projectRoot,'input_data');
functionDir=fullfile(projectRoot,'function');
addpath(inputDir);
addpath(genpath(fullfile(functionDir,'gsw_matlab_v3_06_16_1')));

sourceFile=fullfile(inputDir,'GOBAI-O2-v2.2.nc');
lon=ncread(sourceFile,'lon');
lat=ncread(sourceFile,'lat');
pres=ncread(sourceFile,'pres');
oxy=ncread(sourceFile,'oxy');

[pres,lat]=meshgrid(pres,lat);
depth=-gsw_z_from_p(pres,lat);
lon=cat(1,lon(341:360),lon(1:340));
oxy=cat(1,oxy(341:360,:,:,:),oxy(1:340,:,:,:));
[LAT,LON]=meshgrid(lat(:,1),0.5:1:359.5);
