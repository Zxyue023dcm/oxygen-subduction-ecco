%% Generate FIG3_hmplot_ek
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
load(fullfile(projectRoot,'input_data','PhscData_monthly.mat'),'oceanlon','oceanlat');

load(fullfile(projectRoot,'function','PO_subbasins.mat'),'BASIN');

load(fullfile(projectRoot,'result_data','Wek.mat'),'WEKm');
WEKm=repmat(-WEKm,[1 2]);
Wek_a=WEKm-mean(WEKm,2,'omitnan');

load(fullfile(projectRoot,'function','Areaall.mat'),'Areaall');
Areaall(BASIN~=2)=nan;
Areaall=repmat(Areaall,[1 1 12]);

m0=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
Sox=mean(reshape(m0.OS,[720 360 12 14]),4,'omitnan');
Basin=repmat(BASIN,[1 1 12]);
Sox(Basin~=2)=nan;
Areaall(isnan(Sox))=nan;Sox(isnan(Areaall))=nan;

lats=-89:2:89;
sox=reshape(Sox,[720*4 90 12]);
Areaall=reshape(Areaall,[720*4 90 12]);
sox=squeeze(sum(sox.*Areaall,1,'omitnan')./10^6);
sox(sox==0)=nan;
sox=repmat(sox,[1 2]);

sox_anom=sox-mean(sox,2,'omitnan');

load(fullfile(projectRoot,'input_data','PhscData_monthly.mat'),'wvelb2');

wb=mean(reshape(wvelb2,[720 360 12 14]),4,'omitnan');
Basin=repmat(BASIN,[1 1 12]);
wb(Basin~=2)=nan;
wb(isnan(Sox))=nan;
clear wvelb2
latstr='[-89:2:89]';
eval(['lats=' latstr ';']);
wb=reshape(wb,[720*4 90 12]);
wb=squeeze(mean(wb,1,'omitnan'));
wb(wb==0)=nan;
wb=repmat(wb,[1 2]);

wb_anom=wb-mean(wb,2,'omitnan');

[X,Y]=meshgrid([1:24],lats);

ft=12;

h=figure;
c1=subplot(1,4,1);
hold on;box on;
 pcolor(X,Y,sox); shading interp
 contour(X,Y,wb,[0 0],'color',[.7 .7 .7],'Linewidth',1,'ShowText','off');
contour(X,Y,wb,[1:1:3].*10^(-6),'LineStyle','-','color',[.7 .7 .7]);
contour(X,Y,wb,[3:2:10].*10^(-6),'LineStyle','-','color',[.7 .7 .7]);
 contour(X,Y,wb,-[1:1:3].*10^(-6),'LineStyle','--','color',[.7 .7 .7]);
ylabel('Latitude');
set(gca,'fontsize',ft,'ylim',[-65 60],'xlim',[1 24],'TickDir','both',...
   'YTick',[-60:5:60],'YTickLabel',{'60°S','','','45°S','','','30°S','','','15°S','','','0','','','15°N','','','30°N','','','45°N','','','60°N'},...
    'XTick',[2:3:24],'XTickLabel',{'Feb','May','Aug','Nov','Feb','May','Aug','Nov'});

 caxis([-15 15]);colormap(c1,flipud(othercolor('RdYlBu10',18)));
cc1=colorbar('location','southoutside','Orientation', 'horizontal',...
    'Ticks',[-15:5:15],'Limits',[-15 15]);
cc1.Label.String={'F^o^x';'(Tmol mon^-^1)'};
grid on;box on;

cc1.FontSize=ft;

c2=subplot(1,4,2);
hold on;box on;
 pcolor(X,Y,sox_anom); shading interp
 contour(X,Y,wb_anom,[0 0],'color',[.7 .7 .7],'Linewidth',1,'ShowText','off');
 contour(X,Y,wb_anom,[0.5:0.5:3].*10^(-6),'LineStyle','-','color',[.7 .7 .7]);
 contour(X,Y,wb_anom,-[0.5:0.5:4].*10^(-6),'LineStyle','--','color',[.7 .7 .7]);

set(gca,'fontsize',ft,'ylim',[-65 60],'xlim',[1 24],'TickDir','both',...
   'YTick',[-60:5:60],'YTickLabel',{'60°S','','','45°S','','','30°S','','','15°S','','','0','','','15°N','','','30°N','','','45°N','','','60°N'},...
    'XTick',[2:3:24],'XTickLabel',{'Feb','May','Aug','Nov','Feb','May','Aug','Nov'});

 caxis([-7 7]);colormap(c2,polarmap(28));
cc2=colorbar('location','southoutside','Orientation', 'horizontal',...
    'Ticks',[-6:2:6],'Limits',[-7 7]);
cc2.Label.String={'F^o^x anomalies';'(Tmol mon^-^1)'};
cc2.FontSize=ft;

Wscale=10^6;
c3=subplot(1,4,3);
hold on;box on;
 pcolor(X,Y,WEKm*Wscale); shading interp
 contour(X,Y,wb*Wscale,[0 0],'color',[.7 .7 .7],'Linewidth',1,'ShowText','off');
contour(X,Y,wb*Wscale,[1:1:3].*10^(-6)*Wscale,'LineStyle','-','color',[.7 .7 .7]);
contour(X,Y,wb*Wscale,[3:2:10].*10^(-6)*Wscale,'LineStyle','-','color',[.7 .7 .7]);
 contour(X,Y,wb*Wscale,-[1:1:3].*10^(-6)*Wscale,'LineStyle','--','color',[.7 .7 .7]);

set(gca,'fontsize',ft,'ylim',[-65 60],'xlim',[1 24],'TickDir','both',...
   'YTick',[-60:5:60],'YTickLabel',{'60°S','','','45°S','','','30°S','','','15°S','','','0','','','15°N','','','30°N','','','45°N','','','60°N'},...
    'XTick',[2:3:24],'XTickLabel',{'Feb','May','Aug','Nov','Feb','May','Aug','Nov'});
 caxis([-3 3]);colormap(c3,flipud(othercolor('BrBG5',14)));

cc3=colorbar('location','southoutside','Orientation', 'horizontal',...
    'Ticks',[-3:1.5:3],'Limits',[-3 3]);
cc3.Label.String={'-ω_E_k';'(× 10^-^6 m s^-^1)'};cc3.FontSize=ft;


c4=subplot(1,4,4);
hold on;
 pcolor(X,Y,Wek_a*Wscale); shading interp
 contour(X,Y,wb_anom*Wscale,[0 0],'color',[.7 .7 .7],'Linewidth',1,'ShowText','off');
 contour(X,Y,wb_anom*Wscale,[0.5:0.5:4].*10^(-6)*Wscale,'LineStyle','-','color',[.7 .7 .7]);
 contour(X,Y,wb_anom*Wscale,-[0.5:0.5:4].*10^(-6)*Wscale,'LineStyle','--','color',[.7 .7 .7]);

set(gca,'fontsize',ft,'ylim',[-65 60],'xlim',[1 24],'TickDir','both',...
   'YTick',[-60:5:60],'YTickLabel',{'60°S','','','45°S','','','30°S','','','15°S','','','0','','','15°N','','','30°N','','','45°N','','','60°N'},...
    'XTick',[2:3:24],'XTickLabel',{'Feb','May','Aug','Nov','Feb','May','Aug','Nov'});

 caxis([-3 3]);colormap(c4,polarmap(28));
cc4=colorbar('location','southoutside','Orientation', 'horizontal',...
    'Ticks',[-3:1:3],'Limits',[-3 3]);
cc4.Label.String={['-ω_E_k anomalies'];['(× 10^-^6 m s^-^1)']};
cc4.Label.FontSize = ft;
grid on;box on;

h;
for c=1:4
subplot(1,4,c);hold on;
line([5 5],[-11 -33],'color','r','linestyle','-','linewidth',1);
line([5 11],[-11 -11],'color','r','linestyle','-','linewidth',1);
line([5 11],[-33 -33],'color','r','linestyle','-','linewidth',1);
line([11 11],[-11 -33],'color','r','linestyle','-','linewidth',1);
line([11 11],[11 35],'color','r','linestyle','-','linewidth',1);
line([11 16],[11 11],'color','r','linestyle','-','linewidth',1);
line([11 16],[35 35],'color','r','linestyle','-','linewidth',1);
line([16 16],[11 35],'color','r','linestyle','-','linewidth',1);

end

for c=1:4
subplot(1,4,c);hold on;
line([5 5],[35 60],'color','b','linestyle','-','linewidth',1);
line([5 11],[35 35],'color','b','linestyle','-','linewidth',1);
line([5 11],[60 60],'color','b','linestyle','-','linewidth',1);
line([11 11],[35 60],'color','b','linestyle','-','linewidth',1);
line([11 11],[-33 -49],'color','b','linestyle','-','linewidth',1);
line([11 16],[-33 -33],'color','b','linestyle','-','linewidth',1);
line([11 16],[-49 -49],'color','b','linestyle','-','linewidth',1);
line([16 16],[-33 -49],'color','b','linestyle','-','linewidth',1);

end


 fn=[fullfile(projectRoot,'figs','FIG3_hmplot_ek')];
exportgraphics(h,[fn '.pdf'], 'Resolution', 450,'ContentType','vector');

 exportgraphics(h, [fn '.png'], 'Resolution', 450);

