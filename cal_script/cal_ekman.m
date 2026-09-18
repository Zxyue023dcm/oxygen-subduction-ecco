%% Calculate monthly and climatological Ekman pumping velocity
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
functionDir=fullfile(projectRoot,'function');
inputFile=fullfile(projectRoot,'input_data','PhscData_monthly.mat');
outputFile=fullfile(projectRoot,'result_data','Wek.mat');
addpath(genpath(functionDir));

m=matfile(inputFile);
oceanlon=m.oceanlon;oceanlat=m.oceanlat;
taue=m.taue;taun=m.taun;

Wek=cal_Wek(oceanlon,oceanlat,taue,taun);
Wek_0=Wek;
Wek(:,170:191,:)=nan; % Mask unrealistically large values near f=0.
save(outputFile,'Wek_0','Wek','-v7.3');

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');
WEKm=mean(reshape(Wek,[720 360 12 14]),4,'omitnan');
WEKm(repmat(BASIN,[1 1 12])~=2)=nan;
WEKm=reshape(WEKm,[720*4 90 12]);
WEKm=squeeze(mean(WEKm,1,'omitnan'));
WEKm(WEKm==0)=nan;
save(outputFile,'WEKm','-append');
