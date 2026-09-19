%% Generate FIGS13_Ot_Section
clear
scriptDir=fileparts(mfilename('fullpath'));
projectRoot=fileparts(scriptDir);
functionDir=fullfile(projectRoot,'function');
inputDir=fullfile(projectRoot,'input_data');
resultDir=fullfile(projectRoot,'result_data');
figureDir=fullfile(projectRoot,'figs');
if ~isfolder(figureDir), mkdir(figureDir); end
addpath(genpath(functionDir));
set(groot,'defaultFigureUnits','pixels');
set(groot,'defaultFigurePosition',[100 100 1400 900]);

[DEP,LAT]=meshgrid([5:10:805],[-63.5:1:60.5]);

load(fullfile(projectRoot,'result_data','Ot_section_monclim.mat'),'Oxy','Lat','Depth');
Lat=Lat(1,26:151)';
Depth=Depth(26:151,1:41);
Oxy=Oxy(26:151,1:41,:);

for m=1:12
OXY(:,:,m)=griddata(repmat(Lat,[1 41]),Depth,double(Oxy(:,:,m)),LAT,DEP,'linear');
end
clear Oxy Depth Lat

load(fullfile(projectRoot,'result_data','Ot_section_monclim.mat'),'PRHO','oceanlat','depth');

oceanlat=oceanlat(1,52:303)';
depth=depth(1:27);
PRHO=PRHO(52:303,1:27,:);

for m=1:12
pRHO(:,:,m)=griddata(repmat(oceanlat,[1 27]),repmat(depth',[252 1]),PRHO(:,:,m),LAT,DEP,'linear');
end
clear oceanlat depth PRHO

pRHO=pRHO-1000;

[~,~,OXY_t]=gradient(cat(3,OXY(:,:,12),OXY,OXY(:,:,1))); OXY_t=OXY_t(:,:,2:13);


Oxy3=squeeze(mean(reshape(cat(3,OXY_t(:,:,12),OXY_t(:,:,1:11)),[125 81 3 4]),3,'omitnan'));

Prho=squeeze(mean(reshape(cat(3,pRHO(:,:,12),pRHO(:,:,1:11)),[125 81 3 4]),3,'omitnan'));

clear  OXY OXY_t pRHO
load(fullfile(projectRoot,'input_data','PhscData_monthly.mat'),'oceanlat','Hml','Hmax');

Hml=squeeze(mean(reshape(cat(3,Hml(:,:,168),Hml(:,:,1:167)),[720 360 3 4 14]),[3 5],'omitnan'));

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');

Hml(repmat(BASIN,[1 1 4])~=2)=nan;
Hmax(BASIN~=2)=nan;

Hml=squeeze(mean(Hml,1,'omitnan'));
Hmax=squeeze(mean(Hmax,1,'omitnan'));

oceanlat=oceanlat(1,53:302)';
Hml=Hml(53:302,:);
Hmax=Hmax(53:302)';

oceanlat1=[-63.5:0.25:60.5]';
Hmax1=interp1(oceanlat,Hmax,oceanlat1,'linear');

m=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
Sox=m.OS;
Sox=cat(3,Sox(:,:,168),Sox(:,:,1:167));
Sox=squeeze(mean(reshape(Sox,[720 360 3 4 14]),[3 5],'omitnan'));
Sox(repmat(BASIN,[1 1 4])~=2)=nan;
Sox=squeeze(mean(Sox,1,'omitnan'));
Sox=Sox(53:302,:);

SoxA=Sox-mean(Sox,2,'omitnan');
SoxA(SoxA<0)=nan;

SoxA_ID=nan(size(Hmax1,1),4);
for s=1:4
N=oceanlat(~isnan(SoxA(:,s)));
for n=1:size(N,1)
a=find(oceanlat1==N(n));
SoxA_ID(a-1:a+1,s)=Hmax1(a-1:a+1);
end
end
[all_themes, all_colors] = GetColors();

Season={'DJF','MAM','JJA','SON'};

TL={'∂[O_2]/∂t'};

h=figure;
for n=1:1
    for s=1:4
c(n)=subplot(2,2,s);
hold on; box on
pcolor(LAT,DEP,Oxy3(:,:,s,n));shading flat;
[cc1,cc2]=contour(LAT,DEP,Prho(:,:,s),[21:0.5:28],'color',[.6 .6 .6],'linewidth',0.3);
       clabel(cc1,cc2,'fontsize',8,'color',[.6 .6 .6]);

 clim([-8 8]);colormap(c(n),flipud(othercolor('RdBu11',14)));%c=colorbar;

plot(oceanlat,Hml(:,s),"Color",'k','LineStyle','--','LineWidth',0.8);

plot(oceanlat1,Hmax1,'color','b','LineStyle','-','LineWidth',0.8);
plot(oceanlat1,SoxA_ID(:,s),'color','r','LineStyle','-','LineWidth',0.8);

title([TL{n} ' - ' Season{s}]);
ylabel('Depth (m)');

set(gca,'YDir','reverse','xtick',[-60 -30 0 30 60],'XTickLabel',{'60°S','30°S','0','30°N','60°N'}, ...
    'fontsize',12,'TickDir','both');

    end
end

Num={'(a)','(b)','(c)','(d)'};
for n=1:1
    for s=1:4
c(n)=subplot(2,2,s);
hold on; box on
text(45,650,Num{n,s},'FontSize',13,'FontWeight','bold');
    end
end
h1=figure;
 clim([-8 8]);colormap(flipud(othercolor('RdBu11',14)));
 c=colorbar;c.Label.String='μmol kg^-^1 mon^-^1';
 set(gca,'fontsize',14);

fn=[fullfile(projectRoot,'figs','FIGS13_Colorbar.png')];
exportgraphics(h1,fn, 'Resolution', 450);


 fn=[fullfile(projectRoot,'figs','FIGS13_Ot_Section.pdf')];
exportgraphics(h,fn, 'Resolution', 450,'ContentType','vector');

fn=[fullfile(projectRoot,'figs','FIGS14_Ot_Section.png')];
exportgraphics(h,fn, 'Resolution', 450);

