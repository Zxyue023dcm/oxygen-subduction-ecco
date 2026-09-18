%% Calculate seasonal oxygen-subduction centroid latitudes
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
functionDir=fullfile(projectRoot,'function');
resultFile=fullfile(projectRoot,'result_data','OSresults_monthly.mat');
addpath(genpath(functionDir));

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');
BASIN=BASIN'==2;

m=matfile(resultFile);
oceanlon=m.oceanlon;
oceanlat=m.oceanlat;
Sox=m.OS;

% Shift to December-February seasons and average all 14 years.
Sox=cat(3,Sox(:,:,end),Sox(:,:,1:end-1));
Sox=squeeze(mean(reshape(Sox,[size(Sox,1) size(Sox,2) 3 4 14]), ...
    [3 5],'omitnan'));
Sox=permute(Sox(:,:,[1 3]),[2 1 3]);
Sox(repmat(~BASIN,[1 1 2]))=nan;

latBands=[-66 -10;10 63];
bandNames={'SH','NH'};
pacificLon=[122 292];
power=1;
qInfluence=1;
negFactor=[0 0.5 1];
central_LAT=nan(4,numel(negFactor));
coherence=nan(4,numel(negFactor));

for i=1:numel(negFactor)
    S=positiveCentroidBands(Sox,oceanlat(1,:),oceanlon(:,1), ...
        pacificLon,power,negFactor(i),BASIN,latBands,bandNames,qInfluence);
    central_LAT(:,i)=S.latSph(:);
    coherence(:,i)=S.coherence(:);
end

disp(central_LAT);
disp(coherence);
