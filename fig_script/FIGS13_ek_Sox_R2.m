%% Generate FIGS12_ek_Sox_R2
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
load(fullfile(projectRoot,'result_data','Wek.mat'),'Wek');
Wek(repmat(BASIN,[1 1 168])~=2)=nan;
m0=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
OSp=m0.OS;oceanlon=m0.oceanlon;oceanlat=m0.oceanlat;
OSp(repmat(BASIN,[1 1 168])~=2)=nan;

for i=1:720
for j=1:360

if ~isnan(OSp(i,j,1))

x=squeeze(Wek(i,j,:));
y1=squeeze(OSp(i,j,:));

[beta,~, ~, ~, stats ] = regress(y1, [ones(size(x)), x]);
R_squared(i,j)= stats(1);
p_value(i,j)= stats(3);
else
R_squared(i,j) = nan;
p_value(i,j)= nan;
end

end
end

for b=1:180
R_squared1(:,2*b-1:2*b)=repmat(mean(R_squared(:,2*b-1:2*b),2,'omitnan'),[1,2]);
end
for b=1:179
R_squared1(:,2*b:2*b+1)=repmat(mean(R_squared1(:,2*b:2*b+1),2,'omitnan'),[1,2]);
end

R_squared1(isnan(R_squared))=nan;

pv=p_value;
pv(pv>=0.05)=nan;
pv(pv<0.05)=1;
[lat,lon]=meshgrid([-89:4:89],[1:4:359]);
pv=griddata(oceanlon,oceanlat,pv,lon,lat,'nearest');


h2=figure;

    data=[oceanlon(p_value>0.05) oceanlat(p_value>0.05)];
hold on;
      m_proj('miller','lat',[-63 62],'long',[122 292]);
        hold on
       m_pcolor(oceanlon,oceanlat,R_squared1*100); shading interp
             m_contour(oceanlon,oceanlat,R_squared1*100,[34 34],'w','linewidth',1);
       caxis([0 100]);colormap(othercolor('RdPu9'));

              DOTS=[lon(pv==1) lat(pv==1)];
       s=   m_scatter(DOTS(:,1),DOTS(:,2),10,'k','Marker','X');

        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'fontsize',8,'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -40 -15 0 15 40 60],'fontsize',12);

c=colorbar;
c.Label.String='R^2 (%)';c.FontSize=12;

load(fullfile(functionDir,'Areaall.mat'));

Areaall(isnan(R_squared))=nan;

Area_ALL=sum(Areaall,"all",'omitnan');

Area_r50=sum(Areaall(R_squared>0.5),"all",'omitnan')/Area_ALL;
Area_r40=sum(Areaall(R_squared>0.4),"all",'omitnan')/Area_ALL;
Area_r34=sum(Areaall(R_squared>0.34),"all",'omitnan')/Area_ALL;


fn=[fullfile(projectRoot,'figs','FIGS13_ek_Sox_R2.png')];
exportgraphics(h2,fn, 'Resolution', 450);

