%% Identify the component most correlated with total oxygen subduction
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
functionDir=fullfile(projectRoot,'function');
resultDir=fullfile(projectRoot,'result_data');
addpath(genpath(functionDir));

m=matfile(fullfile(resultDir,'OSresults_monthly.mat'));
oceanlon=m.oceanlon;oceanlat=m.oceanlat;
OS=cat(4,m.OS,m.OSlavoX,m.OSlavoY,m.OSw,m.OSeddy);
OS=reshape(OS,[size(OS,1) size(OS,2) 12*14 5]);

for cpn=1:size(OS,4)
    for i=1:size(OS,1)
        for j=1:size(OS,2)
            OS(i,j,:,cpn)=detrend(squeeze(OS(i,j,:,cpn)));
        end
    end
end

R=nan(size(OS,1),size(OS,2),4);
pValue=nan(size(R));
Main=nan(size(OS,1),size(OS,2));
MainPV=false(size(Main));

for i=1:size(OS,1)
    for j=1:size(OS,2)
        x=squeeze(OS(i,j,:,1));
        aver=nan(4,1);
        for cpn=2:5
            y=squeeze(OS(i,j,:,cpn));
            [a,b]=corrcoef(x',y');
            R(i,j,cpn-1)=a(1,2);
            pValue(i,j,cpn-1)=b(1,2);
            aver(cpn-1)=median(y);
        end

        localR=squeeze(R(i,j,:));
        nn=find(localR==max(localR));
        if numel(nn)==1
            Main(i,j)=nn;
        elseif numel(nn)>1
            [~,largestTerm]=max(abs(aver(nn)));
            Main(i,j)=nn(largestTerm);
        end
        MainPV(i,j)=all(squeeze(pValue(i,j,:))>0.05);
    end
end

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');
Main(MainPV)=nan;
Main(BASIN~=2)=nan;
[lat,lon]=meshgrid(-89:2:89,1:2:359);
maxR=griddata(oceanlon,oceanlat,Main,lon,lat,'nearest');

outputFile=fullfile(resultDir,'seas_corr_maxR.mat');
save(outputFile,'oceanlon','oceanlat','Main','maxR','lat','lon','-v7.3');
save(outputFile,'R','pValue','-append');
