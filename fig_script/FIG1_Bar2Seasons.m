%% Generate FIG1_Bar2Seasons
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
m=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
OS=m.OS;
OS=cat(3,OS(:,:,168),OS(:,:,1:167));
OS=squeeze(mean(reshape(OS,[720 360 3 4 14]),[3 5],'omitnan'));

oceanlon=m.oceanlon;oceanlat=m.oceanlat;

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN1','groupnum');
 groupnum=flipud(groupnum);

load(fullfile(functionDir,'Areaall.mat'),'Areaall');

lats=-89:2:89;

    for s=1:4
os=OS(:,:,s);
    os(isnan(BASIN1))=nan;
Area=Areaall;

 oslat=reshape(os,[720*4 90 1]);
 area=reshape(Area,[720*4 90 1]);

 for j=1:size(lats,2)
 line1=oslat(:,j); line2=area(:,j);
 line2(isnan(line1))=nan;
  subarea=sum(line2(line1>0),'omitnan')./sum(line2,'omitnan');
 obarea=sum(line2(line1<0),'omitnan')./sum(line2,'omitnan');
 OSArea_lat(j,s)=subarea;
 subsum=sum(line2.*line1,'omitnan')./10^6;
 OSSum_lat(j,s)=subsum; % Tmol mon^-1
  subs=sum(line2(line1>0).*line1(line1>0),'omitnan')./sum(line2(line1>0),'omitnan');
  OSV_S_lat(j,s)=subs;  % mol m^-2 mon^-1
  subo=sum(line2(line1<0).*line1(line1<0),'omitnan')./sum(line2(line1<0),'omitnan');
  OSV_O_lat(j,s)=subo;  % mol m^-2 mon^-1

 end
 end


for g=1:5
    for s=1:4
os=OS(:,:,s);

Area=Areaall;

    os(BASIN1~=groupnum{g,1})=nan;
Area(BASIN1~=groupnum{g,1})=nan;
os(isnan(Area))=nan;
Area(isnan(os))=nan;

  subarea=sum(Area(os>0),'all','omitnan')./sum(Area,'all','omitnan');
 obarea=sum(Area(os<0),'all','omitnan')./sum(Area,'all','omitnan');
 OSArea(g,s)=subarea;
 subsum=sum(Area(os>0).*os(os>0),'all','omitnan')./10^6;
 OSSum_S(g,s)=subsum; % Tmol mon^-1
 obsum=sum(Area(os<0).*os(os<0),'all','omitnan')./10^6;
 OSSum_O(g,s)=obsum; % Tmol mon^-1
OSSum(g,s)=sum(Area.*os,'all','omitnan')./10^6;

  subs=sum(Area(os>0).*os(os>0),'all','omitnan')./sum(Area(os>0),'all','omitnan');
  OSV_S(g,s)=subs;  % mol m^-2 mon^-1
  subo=sum(Area(os<0).*os(os<0),'all','omitnan')./sum(Area(os<0),'all','omitnan');
 OSV_O(g,s)=subo;
 OSV(g,s)=sum(Area.*os,'all','omitnan')./sum(Area,'all','omitnan');

 end
 end

[all_themes, all_colors] = GetColors();
cls(1,:)=all_colors(8,:);cls(2,:)=all_colors(11,:);

sz1=12;

h=figure;
subplot(3,1,2);
hold on;box on;grid on;
 yyaxis left
 p=plot(lats,OSArea_lat(:,[1 3])*100,'LineStyle','-','LineWidth',1);
set(gca,'xlim',[-62 62],'ylim',[0 100],'fontsize',sz1,...
    'xtick',[-60 -47.5 -35 -22.5 -10 0 10 25 40 50 60], ...
    'XTickLabel',{'60°S',' ','35°S',' ','10°S','0','10°N',' ','40°N',' ','60°N'},'YColor','k', ...
    'ytick',[0 25 50 75 100],'fontsize',sz1);
ylabel({['Subduction area'];['(%)']},'fontsize',sz1);
 yyaxis right
 bh=bar([-47.5 -22.5 0 25 50],OSArea(:,[1 3])*100);
for c=1:2
   bh(c).FaceColor=cls(c,:);bh(c).EdgeColor=cls(c,:); bh(c).FaceAlpha=0.85;
   p(c).Color=[cls(c,:) 1];

end
set(gca,'xlim',[-62 62],'ylim',[0 100],'YColor','k','fontsize',sz1, ...
    'ytick',[0 25 50 75 100],'YTickLabel',[]);

subplot(3,1,3);
hold on;box on;grid on;
yyaxis left
plot([-65 65],[0 0],'k-');
 p=plot(lats,OSV_S_lat(:,[1 3]),'LineStyle','-','Marker','none','LineWidth',1);
 p1=plot(lats,OSV_O_lat(:,[1 3]),'LineStyle','--','Marker','none','LineWidth',1);

for c=1:2
     p1(c).Color=[cls(c,:) 1];
   p(c).Color=[cls(c,:) 1];
end
ax1=gca;
set(ax1,'xlim',[-62 62],'ylim',[-8 8],'fontsize',sz1,...
   'xtick',[-60 -47.5 -35 -22.5 -10 0 10 25 40 50 60], ...
    'XTickLabel',{'60°S',' ','35°S',' ','10°S','0','10°N',' ','40°N',' ','60°N'},'YColor','k');

ylabel({['S^o^x and O^o^x rate'];['(mol m^-^2 mon^-^1)']},'fontsize',sz1);

yyaxis right
  bh=bar([-47.5 -22.5 0 25 50],OSV_S(:,[1 3]));
  b1=bar([-47.5 -22.5 0 25 50],OSV_O(:,[1 3]));

for c=1:2
    bh(c).FaceColor=cls(c,:);bh(c).EdgeColor=cls(c,:); bh(c).FaceAlpha=0.85;
     b1(c).FaceColor='none';b1(c).EdgeColor=cls(c,:);
end
set(gca,'xlim',[-62 62],'ylim',[-8 8],'YTickLabel',[],'fontsize',sz1,'YColor','k');

subplot(3,1,1);
hold on;box on;grid on;
 yyaxis left
plot([-65 65],[0 0],'k-');
 p=plot(lats,OSSum_lat(:,[1 3]),'LineStyle','-','LineWidth',1);
 set(gca,'xlim',[-62 62],'ylim',[-28 28],'fontsize',sz1,...
   'xtick',[-60 -47.5 -35 -22.5 -10 0 10 25 40 50 60], ...
    'XTickLabel',{'60°S',' ','35°S',' ','10°S','0','10°N',' ','40°N',' ','60°N'},'YColor','k');

ylabel({['Zonally intergrated S^o^x'];['(Tmol mon^-^1)']},'fontsize',sz1);

 yyaxis right
  bh=bar([-47.5 -22.5 0 25 50],OSSum(:,[1 3]));
for c=1:2
    bh(c).FaceColor=cls(c,:);bh(c).EdgeColor=cls(c,:); bh(c).FaceAlpha=0.85;
   p(c).Color=[cls(c,:) 1];
end
ylabel({['Sub-regional S^o^x'];['(Tmol mon^-^1)']},'fontsize',sz1);
set(gca,'xlim',[-62 62],'ylim',[-45 45],'fontsize',sz1,'YColor','k');

legend(bh,{'DJF','JJA'},'fontsize',12,'Orientation','horizontal');
h;
for i=1:3
subplot(3,1,i);
hold on;
plot([-35 -35],[-50 120],'k--');
plot([-10 -10],[-50 120],'k--');
plot([10 10],[-50 120],'k--');
plot([40 40],[-50 120],'k--');

end
xlabel({['Latitude']},'fontsize',sz1);
subplot(3,1,1);
hold on;
for g=1:5
    text(-50,0,groupnum{g,5},"FontSize",12,"FontName",'Times New Roman');
end

legend(bh,{'DJF','JJA'},'fontsize',12,'Orientation','horizontal');

Num={' f)',' g)',' h)'};
for n=1:3
subplot(3,1,n);
hold on; box on
text(-60,10,Num{n},'FontSize',13,'FontWeight','bold');

end

fn=[fullfile(projectRoot,'figs','FIG1_Bar2Seasons')];
 exportgraphics(h, [fn '.png'], 'Resolution', 450);

