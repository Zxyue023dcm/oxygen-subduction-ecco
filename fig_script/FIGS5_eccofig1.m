%% Generate FIGS4_eccofig1
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

m0=matfile(fullfile(projectRoot,'result_data','OSresults_monthly_ecco.mat'));
OSp=m0.OSp; % mon-1
oceanlon=m0.oceanlon;
oceanlat=m0.oceanlat;

OSp(repmat(BASIN,[1 1 168])~=2)=nan;

OSc=mean(OSp,3,'omitnan');

ft=12;

h=figure;

 c2= subplot(3,1,1);% c2= figure;
  m_proj('miller','lat',[-62 62],'long',[122 294]);
        hold on
       m_pcolor(oceanlon,oceanlat,OSc); shading interp
       caxis([-7 7]);colormap(c2,flipud(othercolor('RdYlBu10',14)));
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
c.Label.String = { '\sigma_{season\_clim}', '(mol m^{-2} mon^{-1})' };c.FontSize=ft;
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
Num={'(a)','(b)','(c)'};
for n=1:3
subplot(3,1,n);
hold on; box on
text(-1.5,1.1,Num{n},'FontSize',12,'FontWeight','bold');

end


 fn=[fullfile(projectRoot,'figs','FIGS5_a')];
 exportgraphics(h, [fn '.png'], 'Resolution', 450);


OS=m0.OSp; % month^-1

OS=cat(3,OS(:,:,168),OS(:,:,1:167));
OS=squeeze(mean(reshape(OS,[720 360 3 4 14]),[3 5],'omitnan'));

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN1','groupnum');
 groupnum=flipud(groupnum);

load(fullfile(functionDir,'Areaall.mat'));

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

h1=figure;
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
h1;
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

Num={' (d)',' (e)',' (f)'};
for n=1:3
subplot(3,1,n);
hold on; box on
text(-60,10,Num{n},'FontSize',13,'FontWeight','bold');

end

fn=[fullfile(projectRoot,'figs','FIGS5_b')];
 exportgraphics(h1, [fn '.png'], 'Resolution', 450);

