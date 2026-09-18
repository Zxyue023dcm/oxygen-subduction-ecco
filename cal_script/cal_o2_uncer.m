%% Extract GOBAI-O2 means and uncertainties at Hmax and Hmax+10 m
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
inputDir=fullfile(projectRoot,'input_data');
functionDir=fullfile(projectRoot,'function');
resultDir=fullfile(projectRoot,'result_data');
addpath(inputDir);
addpath(genpath(functionDir));

sourceFile=fullfile(inputDir,'GOBAI-O2-v2.2.nc');
lon=ncread(sourceFile,'lon');
lat=ncread(sourceFile,'lat');
pres=ncread(sourceFile,'pres');
oxy=ncread(sourceFile,'oxy');
uncer=ncread(sourceFile,'uncer');

[pres,lat]=meshgrid(pres,lat);
depth=-gsw_z_from_p(pres,lat);

% Rotate longitude from -180:180 to 0:360.
lon=cat(1,lon(341:360),lon(1:340));
oxy=cat(1,oxy(341:360,:,:,:),oxy(1:340,:,:,:));
uncer=cat(1,uncer(341:360,:,:,:),uncer(1:340,:,:,:));
[LAT,LON]=meshgrid(lat(:,1),0.5:1:359.5);

load(fullfile(inputDir,'PhscData_monthly.mat'),'Hmax','oceanlon','oceanlat');
mMLD=griddata(oceanlon,oceanlat,Hmax,double(LON),double(LAT),'linear');

outputSize=[size(oxy,1) size(oxy,2) size(oxy,4)];
oxy_in_mml=nan(outputSize,'like',oxy);
oxy_in_mml10=nan(outputSize,'like',oxy);
unc_in_mml=nan(outputSize,'like',uncer);
unc_in_mml10=nan(outputSize,'like',uncer);

for t=1:size(oxy,4)
    disp(t);tic
    for i=1:size(oxy,1)
        for j=1:size(oxy,2)
            ml=mMLD(i,j);
            if isnan(ml)
                continue
            end
            [~,d]=min(abs(depth(j,:)-ml));
            [~,d1]=min(abs(depth(j,:)-ml-10));
            oxy_in_mml(i,j,t)=weighted_average_depth(depth(j,1:d), ...
                squeeze(oxy(i,j,1:d,t)));
            oxy_in_mml10(i,j,t)=weighted_average_depth(depth(j,d:d1), ...
                squeeze(oxy(i,j,d:d1,t)));
            unc_in_mml(i,j,t)=weighted_average_depth(depth(j,1:d), ...
                squeeze(uncer(i,j,1:d,t)));
            unc_in_mml10(i,j,t)=weighted_average_depth(depth(j,d:d1), ...
                squeeze(uncer(i,j,d:d1,t)));
        end
    end
    toc
end

save(fullfile(resultDir,'UncO2_Hmax.mat'),'oxy_in_mml','oxy_in_mml10', ...
    'unc_in_mml','unc_in_mml10','LAT','LON','-v7.3');
