%% Generate FIGS11_4TermsMaxR2
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

load(fullfile(projectRoot,'result_data','seas_corr_maxR.mat'),'R','pValue','oceanlon','oceanlat','Main');
load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');
R(repmat(BASIN,[1 1 4])~=2)=nan;
    R(pValue>0.05)=nan;


tname{1,1}='a) Zonal La.';tname{2,1}='b) Merid. La.';
tname{3,1}='c) Vert. Vel.';tname{4,1}='d) Edd.';

h=figure;
for cpn=1:4

c(cpn)=subplot(3,2,cpn);
hold on;
      m_proj('miller','lat',[-63 62],'long',[122 292]);
       m_pcolor(oceanlon,oceanlat,R(:,:,cpn));shading interp
clim([-1 1]);% colormap(c(cpn),othercolor('P6B6',21));
 set(gca,'fontsize',12);colormap(c(cpn),othercolor('PuOr9',21));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'fontsize',11,'xtick',[135 180:45:300],'ytick',[-60 -30 0 30 60]);
T= title([tname{cpn,1}],'fontsize',12,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [-1.5, 1.25, 0]);
end


Main(BASIN~=2)=nan;

data=[oceanlon(~isnan(Main)) oceanlat(~isnan(Main)) Main(~isnan(Main))];

 h;
 c(5)=subplot(3,2,5);
        m_proj('miller','lat',[-65 65],'long',[122 294]);hold on;
m_scatter(data(:,1),data(:,2),4,data(:,3),'filled');
colormap(c(5),osclormap1new());clim([0.5 4.5]);
 set(gca,'fontsize',12);
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'fontsize',11,'xtick',[135 180:45:300],'ytick',[-60 -30 0 30 60]);
T= title('e) Max Correlation Coefficient','fontsize',12,'FontWeight','bold');
set(T, 'HorizontalAlignment', 'left');
set(T, 'Position', [-1.5, 1.35, 0]);

for i=1:4
A=find(data(:,3)==i);
s(i)=m_scatter(data(A(1),1),data(A(1),2),5,data(A(1),3),'filled');
end
 legend(s,tname{1:4},'box','off');

h1=figure;

b2=subplot(1,2,2);
  colormap(b2,othercolor('PuOr9',21)); clim([-1 1]);
c=colorbar;
c.Label.String={'Correlation coefficient (r)'; 'with total F^o^x'};c.FontSize=12;


 fn=[fullfile(projectRoot,'figs','FIGS12_4TermsMaxR2.png')];
 exportgraphics(h, fn, 'Resolution', 450); % PNG

  fn=[fullfile(projectRoot,'figs','FIGS12_Colorbar.png')];
 exportgraphics(h1, fn, 'Resolution', 450); % PNG

