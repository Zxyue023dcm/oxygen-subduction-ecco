%% Generate FIG1_SoxSpatial
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

load(fullfile(functionDir,'Areaall.mat'),'Areaall');

m0=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
OSp=m0.OS;
oceanlon=m0.oceanlon;
oceanlat=m0.oceanlat;

OSp(repmat(BASIN,[1 1 168])~=2)=nan;

OSc=mean(OSp,3,'omitnan');

[~, all_colors] = GetColors();
cls(1,:)=all_colors(8,:);cls(2,:)=all_colors(11,:);

mj=matfile(fullfile(projectRoot,'result_data','Straj_particle.mat'));
LON=mj.oceanlon;LAT=mj.oceanlat;
OTCP=mj.OTCP;

ft=12;

h=figure;

 c2= subplot(3,1,1);% c2= figure;
  m_proj('miller','lat',[-62 62],'long',[122 294]);
        hold on
       m_pcolor(oceanlon,oceanlat,OSc); shading interp
       caxis([-7 7]);colormap(c2,flipud(othercolor('RdYlBu10',14)));
  m_contour(LON,LAT,OTCP,[50 50],'color',all_colors(20,:),'linewidth',0.75);
     c=colorbar('location','westoutside');
     c.Label.String={['O^o^x                                   S^o^x'];['(mol m^-^2 mon^-^1)']};c.FontSize=ft;

     m_coast('patch',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 40 60],'fontsize',ft-1.5);

 m_line([122 292],[-35 -35],'linewidth',1,'linestyle','--','color','k');
 m_line([122 292],[-10 -10],'linewidth',1,'linestyle','--','color','k');
 m_line([122 292],[10 10],'linewidth',1,'linestyle','--','color','k');
 m_line([122 292],[40 40],'linewidth',1,'linestyle','--','color','k');

OSstd=mean(reshape(OSp,[720 360 12 14]),4,'omitnan');

 stdp=std(OSstd,0,3,'omitnan');

h;
c4=subplot(3,1,2);
  m_proj('miller','lat',[-62 62],'long',[122 294]);
        hold on
       m_pcolor(oceanlon,oceanlat,stdp);
       colormap(c4,parula(12));
  caxis([0 3]);
   c=colorbar('location','westoutside');
c.Label.String = { '\sigma_{season\_clim}', '(mol m^{-2} month^{-1})' };c.FontSize=ft;
        m_coast('patch',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 40 60],'fontsize',ft-1.5); % title(['(d)'],'fontsize',12,'FontName','Time New Roman','FontWeight','bold');

OSia=OSp; % mon-1

OSsc=repmat(mean(reshape(OSp,[720 360 12 14]),4,'omitnan'),[1 1 14]);

for i=1:720
    for j=1:360
         x = squeeze(OSia(i,j,:));
y =squeeze(OSsc(i,j,:));
[a, b] = corrcoef(x', y');

R(i,j)=a(1,2);pValue(i,j)=b(1,2);
    end
end

R(BASIN~=2)=nan;
R(pValue>0.05)=nan;
graydot=[oceanlon(pValue>0.05 & BASIN==2) oceanlat(pValue>0.05& BASIN==2)];

h;c6=subplot(3,1,3);%c6=figure;
  m_proj('miller','lat',[-62 62],'long',[122 294]);
        hold on
        m_scatter(graydot(:,1),graydot(:,2),5,[.2 .2 .2],'filled');
       m_pcolor(oceanlon,oceanlat,R.^2);
    c=colorbar('location','westoutside');
        c.Label.String='r^2_s_e_a_s_o_n';c.FontSize=ft;
        caxis([0 1]);
 colormap(c6,othercolor('PuBuGn7'));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 40 60],'fontsize',ft-1.5);

h;
Num={' c)',' d)',' e)'};
for n=1:3
subplot(3,1,n);
hold on; box on
text(-1.5,1.1,Num{n},'FontSize',ft+1,'FontWeight','bold');

end

h;
subplot(3,1,1);
hold on;
text(1.5,1.1,'SSP','FontSize',ft,'fontName','Times New Roman');
text(1.5,0.9,'NSTP','FontSize',ft,'fontName','Times New Roman');
text(1.5,0.7,'EP','FontSize',ft,'fontName','Times New Roman');
text(1.5,0.5,'SSTP','FontSize',ft,'fontName','Times New Roman');
text(1.5,0.3,'SAP','FontSize',ft,'fontName','Times New Roman');


 fn=[fullfile(projectRoot,'figs','FIG1_SoxSpatial')];
 exportgraphics(h, [fn '.png'], 'Resolution', 450);

