%% Generate FIGS10_DJF_JJA_map
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
load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');

load(fullfile(functionDir,'Areaall.mat'));

m=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
oceanlon=m.oceanlon;oceanlat=m.oceanlat;
Sox=m.OS;
Sox=cat(3,Sox(:,:,168),Sox(:,:,1:167));
Sox=squeeze(mean(reshape(Sox,[720 360 3 4 14]),[3 5],'omitnan'));
Sox=Sox(:,:,[1 3]);
Sox(repmat(BASIN,[1 1 2])~=2)=nan;

Mon={'DJF','JJA'};
f= figure;
ax = zeros(2,4);

for k = 1:8
    ax(k) = subplot(2,4,k);
    plot(ax(k), rand(10,1));
    title(ax(k), sprintf('子图%d', k));
    pos2(k,:) = get(ax(k), 'Position');
end

close(f);

sz=12;

colors = colororder;
firstColor = colors(1,:);

h1=figure;
for i=1:2
c= subplot(2,2,i*2);
      m_proj('miller','lat',[-63 62],'long',[122 292]);
        hold on
       m_pcolor(oceanlon,oceanlat,Sox(:,:,i)); shading interp
       clim([-7 7]);colormap(c,flipud(othercolor('RdYlBu10',14)));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
                m_grid('linestyle','none','tickdir','both','fontsize',sz,'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135:45:290],'ytick',[-60 -30 0 30 60]);

 title([Mon{i}],'fontsize',11);
end

lats=[-89:2:89];

for i=1:2
 oslat=reshape(Sox(:,:,i),[720*4 90 1]);
 area=reshape(Areaall,[720*4 90 1]);

 for j=1:size(lats,2)
 line1=oslat(:,j); line2=area(:,j);
 line2(isnan(line1))=nan;
  OSArea_lat(j,i)=sum(line2(line1>0),'omitnan')./sum(line2,'omitnan');
Area_lat(j,i)=sum(line2,'omitnan');
 OSSum_lat(j,i)=sum(line2.*line1,'omitnan')./10^6;  % Tmol mon^-1
 Sub_lat(j,i)=sum(line2(line1>0).*line1(line1>0),'omitnan')./10^6;
SubA_lat(j,i)=sum(line2(line1>0),'omitnan');

 end
 OSArea_lat(OSArea_lat==0)=nan;
OSSum_lat(OSSum_lat==0)=nan;
 Area_lat(Area_lat==0)=nan;
Sub_lat(Sub_lat==0)=nan;
SubA_lat(SubA_lat==0)=nan;
end

h1;
for i=1:2
ax1 = axes('Position', pos2(4*i-2,:) ,...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right',...
           'XColor', 'b', ...
           'YColor', 'b');
hold on
line([50 50],[-65 62], 'Parent', ax1, 'Color','k','LineStyle','--');
line(OSArea_lat(:,i)*100, lats', 'Parent', ax1, 'Color', 'b', 'LineWidth', 1);

xlabel(ax1, 'Subduction Area (%)', 'Color', 'b');
ylabel(ax1, ' ', 'Color', 'b');
set(ax1, 'XLim', [0 100], 'YLim', [-63 62]);

set(ax1,'FontSize',11,'YTick',[-60 -30  0  30 60],'YTickLabel',{' '});

ax2 = axes('Position', get(ax1, 'Position'), ...
           'XAxisLocation', 'bottom', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'XColor', 'k', ...
           'YColor', 'k');

line(OSSum_lat(:,i), lats', 'Parent', ax2, 'Color', 'k', 'LineWidth', 1);

xlabel(ax2, 'S^o^x (Tmol mon^-^1)', 'Color', 'k');
set(ax2, 'XLim', [-25 25], 'YLim', [-63 62]);

set(ax2,'FontSize',11,'YTick',[-60 -30  0  30 60],'YTickLabel',{'60°S','30°S','0','30°N','60°N'});

end


h=figure;
b1=subplot(1,2,1);
  colormap(b1,flipud(othercolor('RdYlBu10',14))); clim([-7 7]);
c=colorbar('location','eastoutside');
c.Label.String='mol m^-^2 mon^-^1';c.FontSize=13;

fn=[fullfile(projectRoot,'figs','FIGS10_DJF_JJA_map.png')];
 exportgraphics(h1, fn, 'Resolution', 450); % PNG

fn=[fullfile(projectRoot,'figs','FIGS10_Colorbar.png')];
 exportgraphics(h, fn, 'Resolution', 450); % PNG

