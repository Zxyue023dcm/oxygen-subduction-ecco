%% Generate FIGS5_eccofig2
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
m=matfile(fullfile(projectRoot,'result_data','OSresults_monthly_ecco.mat'));
oceanlon=m.oceanlon;oceanlat=m.oceanlat;

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');
load(fullfile(functionDir,'Areaall.mat'));

OSc=cat(4,m.OSp,m.OSlavoX,m.OSlavoY,m.OSw,m.OSeddy);
OSc=squeeze(mean(OSc,3,'omitnan'));
OSc(repmat(BASIN,[1 1 5])~=2)=nan;

OS(:,:,:,1)=m.OSp;
OS(:,:,:,2)=m.OSlavoX;
OS(:,:,:,3)=m.OSlavoY;
OS(:,:,:,4)=m.OSw;
OS(:,:,:,5)=m.OSeddy;

for cpn=1:size(OS,4)
for i=1:720
    for j=1:360
OS(i,j,:,cpn)=detrend(squeeze(OS(i,j,:,cpn)));
    end
end
end

MainPV=zeros(720,360);
for i=1:720
    for j=1:360
        aver=[];
for cpn=2:5
x=squeeze(OS(i,j,:,1));
y=squeeze(OS(i,j,:,cpn));
[a, b] = corrcoef(x', y');
R(i,j,cpn-1)=a(1,2);pValue(i,j,cpn-1)=b(1,2);
aver(cpn-1,1)=median(y);
end

nn=find(R(i,j,:)==max(R(i,j,:)));
if size(nn,1)==1
Main(i,j)=nn;
elseif size(nn,1)>1
Main(i,j)=find(abs(aver(nn))==max(abs(aver(nn))));
else
Main(i,j)=nan;
end

if sum(pValue(i,j,:)>0.05)==4
 MainPV(i,j)=1;
end

    end
end

Main(MainPV==1)=nan;
Main(BASIN~=2)=nan;

[lat,lon]=meshgrid([-89:2:89],[1:2:359]);
maxR=griddata(oceanlon,oceanlat,Main,lon,lat,'nearest');

Areaall(BASIN~=2)=nan;
A=sum(Areaall,'all','omitnan');
 for i=1:4
P_main(i,1)=sum(Areaall(Main==i),'all','omitnan')/A *100;
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


fn=[fullfile(projectRoot,'figs','FIGS6_eccofig2.png')];%h=figure(3);
 exportgraphics(h, fn, 'Resolution', 450);

 fn=[fullfile(projectRoot,'figs','FIGS6_colorbar.png')];
 exportgraphics(h1, fn, 'Resolution', 450);

