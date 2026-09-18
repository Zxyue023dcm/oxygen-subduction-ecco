%% Calculate Hml, Hmax, and seasonal thermocline thickness
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
sourceDir=fullfile(projectRoot,'input_data','raw');
functionDir=fullfile(projectRoot,'function');
outputFile=fullfile(projectRoot,'input_data','PhscData_monthly.mat');
addpath(genpath(functionDir));
assert(isfolder(sourceDir),'ECCO raw-data folder not found: %s',sourceDir);
assert(isfile(outputFile),'Prepared physical-data file not found: %s',outputFile);

years=2004:2017;
mxl=[];

for yr=1:numel(years)
    for mon=1:12
        fprintf('Calculating Hml for %04d-%02d\n',years(yr),mon);tic
        thetaFile=fullfile(sourceDir,sprintf('THETA_%04d_%02d.nc',years(yr),mon));
        saltFile=fullfile(sourceDir,sprintf('SALT_%04d_%02d.nc',years(yr),mon));
        assert(isfile(thetaFile),'ECCO file not found: %s',thetaFile);
        assert(isfile(saltFile),'ECCO file not found: %s',saltFile);

        [theta,lon,lat,dep]=read_ecco_monthly(thetaFile,'THETA');
        [salt,lonSalt,latSalt,depSalt]=read_ecco_monthly(saltFile,'SALT');
        assert(isequaln(lon,lonSalt) && isequaln(lat,latSalt) && ...
            isequaln(dep,depSalt),'THETA and SALT grids do not match for %04d-%02d.', ...
            years(yr),mon);

        if isempty(mxl)
            mxl=nan(size(theta,1),size(theta,2),12,numel(years));
        end
        mxl(:,:,mon,yr)=calculate_mxl(theta,salt,dep,lon,lat);
        clear theta salt lonSalt latSalt depSalt
        toc
    end
end

mxl=reshape(mxl,[size(mxl,1) size(mxl,2) 12*numel(years)]);
halfLon=size(mxl,1)/2;
assert(mod(size(mxl,1),2)==0,'The longitude dimension must have an even length.');
Hml=cat(1,mxl(halfLon+1:end,:,:),mxl(1:halfLon,:,:));
save(outputFile,'Hml','-append');

Hmax=max(Hml,[],3);
save(outputFile,'Hmax','-append');

Hst=Hmax-Hml;
save(outputFile,'Hst','-append');
