%% Generate FIG2_decoms_clim_seas_map
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
oceanlon=m.oceanlon;oceanlat=m.oceanlat;

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');

load(fullfile(projectRoot,'result_data','seas_corr_maxR.mat'));

OSc=cat(4,m.OS,m.OSlavoX,m.OSlavoY,m.OSw,m.OSeddy);
OSc=squeeze(mean(OSc,3,'omitnan'));
OSc(repmat(BASIN,[1 1 5])~=2)=nan;

load(fullfile(functionDir,'Areaall.mat'),'Areaall');
Areaall(BASIN~=2)=nan;
A=sum(Areaall,'all','omitnan');
P_main=nan(4,1);
for i=1:4
    P_main(i)=sum(Areaall(Main==i),'all','omitnan')/A*100;
end


tname{2,1}='a) Zonal lateral induction';tname{3,1}='b) Meridional lateral induction';
tname{4,1}='c) Vertical velocity';tname{5,1}='d) Eddy-induced subduction';

h=figure;
  for i=1:4
c(i)= subplot(2,2,i);
      m_proj('miller','lat',[-63 62],'long',[122 292]);
        hold on
       m_pcolor(oceanlon,oceanlat,OSc(:,:,i+1)); shading interp
     clim([-7 7]);colormap(c(i),flipud(othercolor('RdYlBu10',28)));

       DOTS=[lon(maxR==i) lat(maxR==i)];
       s=   m_scatter(DOTS(:,1),DOTS(:,2),2,'b','filled');
     s.MarkerFaceAlpha = 0.7;
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'fontsize',10,'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -40 -20 0 20 40 60],'fontsize',13);
T= title([tname{i+1} ' (' num2str(P_main(i),'%.1f') '%)'],'fontsize',14,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [-1.5, 1.25, 0]);
  end


h1=figure;
b1=subplot(1,2,1);
  colormap(b1,flipud(othercolor('RdYlBu10',28))); clim([-7 7]);
c=colorbar;
c.Label.String='mol O_2 m^-^2 mon^-^1';c.FontSize=12;
ylabel('O^o^x              S^o^x','FontSize',12);

fn=[fullfile(projectRoot,'figs','FIG2_climmaxR_maps.png')];%h=figure(3);
 exportgraphics(h, fn, 'Resolution', 450);

 fn=[fullfile(projectRoot,'figs','FIG2_colorbar.png')];
 exportgraphics(h1, fn, 'Resolution', 450);

