%% Calculate monthly oxygen-subduction flux and all decomposition terms
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
inputDir=fullfile(projectRoot,'input_data');
functionDir=fullfile(projectRoot,'function');
resultDir=fullfile(projectRoot,'result_data');
physicalFile=fullfile(inputDir,'PhscData_monthly.mat');
oxygenFile=fullfile(inputDir,'oxydata_monthly.mat');
outputFile=fullfile(resultDir,'OSresults_monthly.mat');

assert(isfolder(functionDir),'Function folder not found: %s',functionDir);
assert(isfile(physicalFile),'Physical-data file not found: %s',physicalFile);
assert(isfile(oxygenFile),'Oxygen-data file not found: %s',oxygenFile);
addpath(genpath(functionDir));

physicalVariables={'Hml','Hst','evel1','evel2','nvel1','nvel2', ...
    'wvelb2','estar1','estar2','nstar1','nstar2','prho','prho10'};
oxygenVariables={'oxy','oxy10','oceanlon','oceanlat'};

fprintf('Loading physical and oxygen data...\n');
physical=load(physicalFile,physicalVariables{:});
oxygen=load(oxygenFile,oxygenVariables{:});
missingPhysical=setdiff(physicalVariables,fieldnames(physical));
missingOxygen=setdiff(oxygenVariables,fieldnames(oxygen));
assert(isempty(missingPhysical),'Missing physical variable(s): %s', ...
    strjoin(missingPhysical,', '));
assert(isempty(missingOxygen),'Missing oxygen variable(s): %s', ...
    strjoin(missingOxygen,', '));

fprintf('Calculating monthly oxygen flux...\n');
[OS,OSw,OSeddy,OSlavoX,OSlavoY]=oxygen_flux( ...
    oxygen.oceanlon,oxygen.oceanlat,physical.Hml,physical.Hst, ...
    physical.evel1,physical.evel2,physical.nvel1,physical.nvel2, ...
    physical.wvelb2,physical.estar1,physical.estar2, ...
    physical.nstar1,physical.nstar2,physical.prho,physical.prho10, ...
    oxygen.oxy,oxygen.oxy10);

oceanlon=oxygen.oceanlon;
oceanlat=oxygen.oceanlat;
save(outputFile,'OS','OSw','OSeddy','OSlavoX','OSlavoY', ...
    'oceanlon','oceanlat','-v7.3');
fprintf('Saved monthly oxygen-flux results to:\n%s\n',outputFile);
