%% Calculate oxygen-flux statistics in five eddy-dominated regions
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
functionDir=fullfile(projectRoot,'function');
resultDir=fullfile(projectRoot,'result_data');
addpath(genpath(functionDir));

load(fullfile(resultDir,'seas_corr_maxR.mat'), ...
    'Main','oceanlon','oceanlat');
load(fullfile(functionDir,'Areaall.mat'),'Areaall');

EddyRegion={'Coastal Peru','North Equatorial Pacific','Australia Coast', ...
    'Kuroshi Extension','Alaska Gyre'};
EddyLine(:,:,1)=[280 -7;280 -35;290 -35;290 -7;280 -7];
EddyLine(:,:,2)=[195 14;195 5;285 5;285 14;195 14];
EddyLine(:,:,3)=[154 -15;154 -43;173 -43;173 -15;154 -15];
EddyLine(:,:,4)=[144 45;144 33;188 33;188 45;144 45];
EddyLine(:,:,5)=[180 60;180 47;230 47;230 60;180 60];

BASIN=nan(size(oceanlon));
for n=1:size(EddyLine,3)
    inside=oceanlon>=EddyLine(1,1,n) & oceanlon<=EddyLine(3,1,n) ...
        & oceanlat>=EddyLine(2,2,n) & oceanlat<=EddyLine(1,2,n);
    BASIN(inside)=n;
end

BASIN_ed=BASIN;
for n=1:size(EddyLine,3)
    BASIN_ed(BASIN==n & Main~=4)=nan;
end

outputFile=fullfile(resultDir,'BASINed_maxR.mat');
save(outputFile,'BASIN','EddyRegion','EddyLine','BASIN_ed','-append');

m=matfile(fullfile(resultDir,'OSresults_monthly.mat'));
OS=cat(4,m.OS,m.OSlavoX,m.OSlavoY,m.OSw,m.OSeddy);
nRegion=numel(EddyRegion);
R=cell(1,nRegion);PV=cell(1,nRegion);
r_m=nan(nRegion,4);pv_m=nan(nRegion,4);
POterms_NET_ed=nan(size(OS,4),size(OS,3),nRegion);
POterms_NET=nan(size(OS,4),size(OS,3),nRegion);

for g=1:nRegion
    [L1,L2]=find(BASIN_ed==g);
    r=nan(numel(L1),4);pv=nan(numel(L1),4);
    for l=1:numel(L1)
        x=squeeze(OS(L1(l),L2(l),:,1));
        for i=1:4
            y=squeeze(OS(L1(l),L2(l),:,i+1));
            [a,b]=corrcoef(x',y');
            r(l,i)=a(1,2);
            pv(l,i)=b(1,2);
        end
    end
    R{g}=r;PV{g}=pv;
    r_m(g,:)=mean(r,1,'omitnan');
    pv_m(g,:)=mean(pv,1,'omitnan');

    for i=1:size(OS,4)
        for t=1:size(OS,3)
            os=OS(:,:,t,i);
            A=os(BASIN_ed==g);A1=Areaall(BASIN_ed==g);
            POterms_NET_ed(i,t,g)=(sum(A(A>0).*A1(A>0),'omitnan') ...
                +sum(A(A<0).*A1(A<0),'omitnan'))/10^6;
            A=os(BASIN==g);A1=Areaall(BASIN==g);
            POterms_NET(i,t,g)=(sum(A(A>0).*A1(A>0),'omitnan') ...
                +sum(A(A<0).*A1(A<0),'omitnan'))/10^6;
        end
    end
end

OSbar_ed=squeeze(mean(reshape(POterms_NET_ed,[5 12 14 nRegion]),3,'omitnan'));
OSbar_sigm_ed=squeeze(std(reshape(POterms_NET_ed,[5 12 14 nRegion]),0,3,'omitnan'));
OSbar=squeeze(mean(reshape(POterms_NET,[5 12 14 nRegion]),3,'omitnan'));
OSbar_sigm=squeeze(std(reshape(POterms_NET,[5 12 14 nRegion]),0,3,'omitnan'));

save(outputFile,'OSbar_ed','OSbar_sigm_ed','OSbar','OSbar_sigm', ...
    'R','PV','r_m','pv_m','-append');
