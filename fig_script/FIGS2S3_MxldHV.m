%% Generate FIGS1S2_MxldHV
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

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');
load(fullfile(projectRoot,'input_data','PhscData_monthly.mat'),'Hmax','Hml','Hst');

Hmax(BASIN~=2)=nan;
Hml=mean(Hml,[3 4],'omitnan');Hst=mean(Hst,[3 4],'omitnan');
Hml(BASIN~=2)=nan;Hst(BASIN~=2)=nan;

load(fullfile(projectRoot,'input_data','PhscData_monthly.mat'),'evel1','nvel1');
evel1=mean(evel1,[3 4],'omitnan');nvel1=mean(nvel1,[3 4],'omitnan');
evel1(BASIN~=2)=nan;nvel1(BASIN~=2)=nan;
load(fullfile(projectRoot,'input_data','PhscData_monthly.mat'),'evel2','nvel2');
evel2=mean(evel2,[3 4],'omitnan');nvel2=mean(nvel2,[3 4],'omitnan');
evel2(BASIN~=2)=nan;nvel2(BASIN~=2)=nan;
ft=12;
p1=-1.5;p2=1.35;p3=0;

h=figure;
c1=subplot(2,3,1);

hold on;
 m_proj('miller','lat',[-65 65],'long',[122 294]);
       m_pcolor(oceanlon,oceanlat,Hmax);shading interp
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
         colormap(c1,jet(20)); % Paired6 YlGn9
clim([0 600]);
 cc1=colorbar;    cc1.Label.String='Depth (m)'; cc1.Label.FontSize=ft-0.5;
 set(gca,'fontsize',ft-0.5);
T= title(['a) H_m_a_x'],'fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

c11=subplot(2,3,2);

hold on;
 m_proj('miller','lat',[-65 65],'long',[122 294]);
       m_pcolor(oceanlon,oceanlat,Hml);shading interp
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
         colormap(c11,jet(20)); % Paired6 YlGn9
clim([0 250]);
 cc1=colorbar;    cc1.Label.String='Thickness (m)'; cc1.Label.FontSize=ft-0.5;
 set(gca,'fontsize',ft-0.5);
T= title(['b) H_m_l'],'fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

c12=subplot(2,3,3);

hold on;
 m_proj('miller','lat',[-65 65],'long',[122 294]);
       m_pcolor(oceanlon,oceanlat,Hst);shading interp
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
         colormap(c12,jet(20)); % Paired6 YlGn9
clim([0 450]);
 cc1=colorbar;    cc1.Label.String='Thickness (m)'; cc1.Label.FontSize=ft-0.5;
 set(gca,'fontsize',ft-0.5);
T= title(['c) H_s_t'],'fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

data=Hmax(~isnan(Hmax));% n=63977;
        subplot(2,3,4); hold on;box on;grid on;
hist(data,100);xlabel('Depth (m)');ylabel('Grid number (× 10^3)');
set(gca,'fontsize',ft-2,'YTick',0:1000:5000,'YTickLabel',0:1:5);
 title('d) Frequency distribution of H_m_a_x','fontsize',ft,'FontWeight','bold');

 data=Hml(~isnan(Hml));% n=63977;
        subplot(2,3,5); hold on;box on;grid on;
hist(data,100);xlabel('Thickness (m)');ylabel('Grid number (× 10^3)');
set(gca,'fontsize',ft-2,'YTick',0:1000:5000,'YTickLabel',0:1:5);
 title('e) Frequency distribution of H_m_l','fontsize',ft,'FontWeight','bold');

 data=Hst(~isnan(Hst));% n=63977;
        subplot(2,3,6); hold on;box on;grid on;
hist(data,100);xlabel('Thickness (m)');ylabel('Grid number (× 10^3)');
set(gca,'fontsize',ft-2,'YTick',0:1000:5000,'YTickLabel',0:1:5);
 title('f) Frequency distribution of H_s_t','fontsize',ft,'FontWeight','bold');

             fn=[fullfile(projectRoot,'figs','FIGS2_Mxld.png')];%h=figure(3);
 exportgraphics(h, fn, 'Resolution', 450); % PNG

 [lon_l,lat_l] = distcoordinate(oceanlon(:,1),oceanlat(1,:)');
 ans=diff(lon_l,1,1);lon_l=[lon_l(1,:)-ans(1,:);lon_l];lat_l=[lat_l(1,:);lat_l];
 lon_l=[nan(size(lon_l,1),1) lon_l]; lat_l=[nan(size(lat_l,1),1) lat_l];

     hml=[Hml(end,:);Hml];hml=[nan(size(hml,1),1) hml];
[~,gHmlx,gHmly]= divdffs(hml, lon_l,hml,lat_l,'backward'); % unit: m/m

for a=1:360
gHmlx(2*a-1:2*a,:)=repmat(mean(gHmlx(2*a-1:2*a,:),1,'omitnan'),[2,1]);
end

for b=1:180
gHmly(:,2*b-1:2*b)=repmat(mean(gHmly(:,2*b-1:2*b),2,'omitnan'),[1,2]);
end

Hslpx=gHmlx*1000;
Hslpy=gHmly*1000;

 Hslpx(isnan(Hmax))=nan;
  Hslpy(isnan(Hmax))=nan;

  clear gHmlx gHmlx1 gHmly gHmly1 lon_l lat_l  hml

 [lon_l,lat_l] = distcoordinate(oceanlon(:,1),oceanlat(1,:)');
 ans=diff(lon_l,1,1);lon_l=[lon_l(1,:)-ans(1,:);lon_l];lat_l=[lat_l(1,:);lat_l];
 lon_l=[nan(size(lon_l,1),1) lon_l]; lat_l=[nan(size(lat_l,1),1) lat_l];

     hml=[Hst(end,:);Hst];hml=[nan(size(hml,1),1) hml];
[~,gHmlx,gHmly]= divdffs(hml, lon_l,hml,lat_l,'backward'); % unit: m/m

for a=1:360
gHmlx(2*a-1:2*a,:)=repmat(mean(gHmlx(2*a-1:2*a,:),1,'omitnan'),[2,1]);
end

for b=1:180
gHmly(:,2*b-1:2*b)=repmat(mean(gHmly(:,2*b-1:2*b),2,'omitnan'),[1,2]);
end

Hslpx1=gHmlx*1000;
Hslpy1=gHmly*1000;

 Hslpx1(isnan(Hmax))=nan;
  Hslpy1(isnan(Hmax))=nan;

  clear gHmlx gHmlx1 gHmly gHmly1 lon_l lat_l  hml


h1=figure;
 c3=subplot(2,2,1);
  m_proj('miller','lat',[-65 65],'long',[122 294]);
        hold on
       m_pcolor(oceanlon,oceanlat,Hslpx); clim([-0.1 0.1]); % slop x
      colormap(c3,flipud(othercolor('RdBu5',25)));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
               c=colorbar('Ticks',-0.1:0.05:0.1);
       c.Label.String='m km^-^1'; c.Label.FontSize=ft-0.5;
        set(gca,'fontsize',ft-0.5);
T= title('a) Zonal gradient (H_m_l)','fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

 c5=subplot(2,2,3);
  m_proj('miller','lat',[-65 65],'long',[122 294]);
        hold on;
       m_pcolor(oceanlon,oceanlat,evel1); clim([-0.3 0.3]); % slop x
      colormap(c5,flipud(othercolor('RdBu5',25)));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
               c=colorbar('Ticks',-0.3:0.1:0.3);
       c.Label.String='m s^-^1';  c.Label.FontSize=ft-0.5;
        set(gca,'fontsize',ft-0.5);
T= title('b) Zonal velocity (H_m_l)','fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

 c3=subplot(2,2,2);
  m_proj('miller','lat',[-65 65],'long',[122 294]);
        hold on
       m_pcolor(oceanlon,oceanlat,Hslpx1); clim([-0.15 0.15]); % slop x
      colormap(c3,flipud(othercolor('RdBu5',25)));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
               c=colorbar('Ticks',-0.15:0.05:0.15);
       c.Label.String='m km^-^1'; c.Label.FontSize=ft-0.5;
        set(gca,'fontsize',ft-0.5);
T= title('e) Zonal gradient (H_s_t)','fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

 c5=subplot(2,2,4);
  m_proj('miller','lat',[-65 65],'long',[122 294]);
        hold on;
       m_pcolor(oceanlon,oceanlat,evel2); clim([-0.3 0.3]); % slop x
      colormap(c5,flipud(othercolor('RdBu5',25)));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
               c=colorbar('Ticks',-0.3:0.1:0.3);
       c.Label.String='m s^-^1';  c.Label.FontSize=ft-0.5;
        set(gca,'fontsize',ft-0.5);
T= title('f) Zonal velocity (H_s_t)','fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

h2=figure;
       c4=subplot(2,2,1);
 m_proj('miller','lat',[-65 65],'long',[122 294]);
        hold on
  m_pcolor(oceanlon,oceanlat,Hslpy); clim([-0.2 0.2]); % slop y
      colormap(c4,flipud(othercolor('RdBu5',25)));%RdYlBu10
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
         c=colorbar('Ticks',-0.2:0.1:0.2);
       c.Label.String='m km^-^1';  c.Label.FontSize=ft-0.5;
        set(gca,'fontsize',ft-0.5);
T= title('c) Meridional gradient (H_m_l)','fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

       c6=subplot(2,2,3);
 m_proj('miller','lat',[-65 65],'long',[122 294]);
        hold on  ;
  m_pcolor(oceanlon,oceanlat,nvel1); clim([-0.1 0.1]); % slop y
      colormap(c6,flipud(othercolor('RdBu5',25)));%RdYlBu10
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
         c=colorbar('Ticks',-0.1:0.05:0.1);
       c.Label.String='m s^-^1';  c.Label.FontSize=ft-0.5;
        set(gca,'fontsize',ft-0.5);
T= title('d) Meridional velocity (H_m_l)','fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

       c4=subplot(2,2,2);
 m_proj('miller','lat',[-65 65],'long',[122 294]);
        hold on
  m_pcolor(oceanlon,oceanlat,Hslpy1); clim([-0.3 0.3]); % slop y
      colormap(c4,flipud(othercolor('RdBu5',25)));%RdYlBu10
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
         c=colorbar('Ticks',-0.3:0.1:0.3);
       c.Label.String='m km^-^1';  c.Label.FontSize=ft-0.5;
        set(gca,'fontsize',ft-0.5);
T= title('g) Meridional gradient (H_s_t)','fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);

       c6=subplot(2,2,4);
 m_proj('miller','lat',[-65 65],'long',[122 294]);
        hold on  ;
  m_pcolor(oceanlon,oceanlat,nvel2); clim([-0.1 0.1]); % slop y
      colormap(c6,flipud(othercolor('RdBu5',25)));%RdYlBu10
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 35 60],'fontsize',ft-1);
         c=colorbar('Ticks',-0.1:0.05:0.1);
       c.Label.String='m s^-^1';  c.Label.FontSize=ft-0.5;
        set(gca,'fontsize',ft-0.5);
T= title('h) Meridional velocity (H_s_t)','fontsize',ft,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [p1,p2,p3]);


        fn=[fullfile(projectRoot,'figs','FIGS3_MxldHV_a.png')];%h=figure(3);
 exportgraphics(h1, fn, 'Resolution', 450); % PNG

         fn=[fullfile(projectRoot,'figs','FIGS3_MxldHV_b.png')];%h=figure(3);
 exportgraphics(h2, fn, 'Resolution', 450); % PNG

