%% Generate FIGS6_Straj
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
m0=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
lon=m0.oceanlon;lat=m0.oceanlat;
load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN','BASIN1');

load(fullfile(functionDir,'Areaall.mat'),'Areaall');

OSp=m0.OS;
OSm=mean(OSp,3,'omitnan');
OSp=reshape(OSp,[720 360 12 14]);

groupnum={'SSTP (40°S-10°S)','NSTP (10°N-40°N)'};

LAT_min=[-40 10];LAT_max=[-10 40];
lon_min=122;lon_max=292;

for b=1:size(groupnum,2)

lat_min=LAT_min(b);lat_max=LAT_max(b);

for yrs=1:13
    for mm=1:12
       os=OSp(:,:,mm,yrs);
osi=os(lon>=lon_min & lon<=lon_max & lat>=lat_min & lat<=lat_max);
arei=Areaall(lon>=lon_min & lon<=lon_max & lat>=lat_min & lat<=lat_max);
arei(isnan(osi))=nan;
sub(mm,yrs)=sum(osi(osi>0).*arei(osi>0),'omitnan')/10^6; % Tmol/annual
    end
end
sub(sub==0)=nan;

SUB_miu(b,:)=mean(sub,2,'omitnan');
SUB_sigm(b,:)=std(sub,[],2,'omitnan');
SUB(b,:)=reshape(sub,[1 mm*yrs]);

end

load(fullfile(projectRoot,'result_data','Straj_particle.mat'));

a1=find(lon(:,1) >= oceanlon(1,1) & lon(:,1) <= oceanlon(end,1) );
a2=find(lat(1,:) >= oceanlat(1,1) & lat(1,:) <= oceanlat(1,end) );

OSm=OSm(a1,a2);

OSm(OSm<=0)=nan;
OSm(OSm>0)=1;
[LAT,LON]=meshgrid([-59:2:61],[123:2:291]);
osm=griddata(oceanlon,oceanlat,OSm,LON,LAT,'nearest');

Data=SUB;

for g=1:size(groupnum,2)
[a, b] = corrcoef(Nmonth(g,:)', Data(g,:)', 'rows', 'complete');

R(g,1)=a(1,2);pValue(g,1)=b(1,2);
r2(g,1)=a(1,2).^2;
end

m1=matfile(fullfile(projectRoot,'result_data','Straj_particle.mat'));

sub_month=m1.sub_month;

for g=1:size(groupnum,2)
[a, b] = corrcoef(sub_month(g,:)', Data(g,:)', 'rows', 'complete');

R1(g,1)=a(1,2);pValue1(g,1)=b(1,2);
r21(g,1)=a(1,2).^2;
end

colors =[255 255 255; 206 235 254 ;235 235 235; 244 240 64; 254 0 0 ; 151 0 130];
pos=[0,0.025,0.05,0.1,0.35,0.5]*2;
my_cmap1=custom_colormap(colors,'positions',[0,0.025,0.05,0.1,0.35,0.5]*2);
load(fullfile(functionDir,'colorData.mat'),'RdYlBu10');
colors =flipud(RdYlBu10(1:6,:));
pos=linspace(0,1,6);
my_cmap2=custom_colormap(colors,'positions',pos);

ft=12;
hf=figure;H=tight_subplot(3,2, [0.03 0.04], [0.06 0.04], [0.08 0.08]);
 axes(H(1));
 axis off
 m_proj('miller','lat',[min(oceanlat,[],'all')-0.25 max(oceanlat,[],'all')],'long',[min(oceanlon,[],'all') max(oceanlon,[],'all')]);
         hold on
       m_pcolor(oceanlon,oceanlat,OTCP); shading interp;
        colormap(H(1),my_cmap1);clim([0 140]);
        m_contour(oceanlon,oceanlat,OTCP,[50 50],'k','linewidth',1);
                    DOTS=[LON(osm==1) LAT(osm==1)];
      m_scatter(DOTS(:,1),DOTS(:,2),10,'w','Marker','.');
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',90:45:300,'ytick',[-60 -35 -10  10 35 60],'fontsize',ft);

cc1=colorbar;set(gca,'fontsize',13);title(cc1, 'days','fontsize',ft);

for b=1:size(groupnum,2)

hf;
 axes(H(2+b));
    hold on;box on;grid on;
    yyaxis left
    if b==1
    set(gca, 'XTick', 1:12, 'XTickLabel', 1:12, 'YTick', 0:10:60, 'YTickLabel',0:10:60, 'YColor', 'k','fontsize',ft);
    ylabel('Frequency (%)');
    elseif b==2
set(gca, 'XTick', 1:12, 'XTickLabel', 1:12, 'YTick', 0:10:60, 'YColor', 'k','fontsize',ft);
    end
set(gca, 'XTickMode', 'manual', 'YTickMode', 'manual');
 hb=bar(1:12,F1M(:,b), 'Facecolor',[.8 .8 .8],'EdgeColor','k','BarWidth',1);
   eb=  errorbar(1:12,F1M(:,b), F1S(:,b), ...
            'Color', 'k', 'LineStyle', 'none', ...
            'LineWidth', 1, 'CapSize', 4);  eb.CapSize = 6;
[values, indices] = maxk(F1M(:,b), 4);
offset =3;
for m = indices'
    y_top = F1M(m,b) + F1S(m,b);
    text(m, y_top + offset, sprintf('%.1f%%', F1M(m,b)), ...
         'HorizontalAlignment', 'center', 'FontSize', ft, 'FontWeight', 'normal');
end
ylim([0 60]);
if pValue(b)<0.01
text(1,55,['r = ' num2str(R(b),'%.2f') ' , p < 0.01'],"Color",'b', 'FontSize', ft);
else
text(1,55,['r = ' num2str(R(b),'%.2f') ' , p = ' num2str(pValue(b),'%.2f')],"Color",'b', 'FontSize', ft);
end

   yyaxis right
   e1 =  plot(1:12,SUB_miu(b,:), 'Marker','none',"Color",'b','LineStyle','-','LineWidth',1);
set(gca,'fontsize',ft,'xtick',1:1:12,'YColor','b');
if b==2
ylabel('S^o^x (Tmol yr^-^1)');end
xlim([0.5 12.5]);

end

Sr=mean(m1.Sr,3,'omitnan')*12;
Sr=smoothdata(Sr,1,'movmean',4,'omitnan');
 hf; %s2=subplot(3,2,2);
  Sr(:,111:130)=nan;
 axes(H(2));
m_proj('miller','lat',[-60 62],'long',[122 292]);         hold on
       m_pcolor(oceanlon,oceanlat,Sr); shading interp;
        colormap(H(2),my_cmap2);  clim([0 80]);
    m_contour(oceanlon,oceanlat,OTCP,[50 50],'k','linewidth',1);
m_scatter(DOTS(:,1),DOTS(:,2),10,'w','Marker','.');
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',90:45:300,'ytick',[-60 -35 -10  10 35 60],'fontsize',ft);

cc2=colorbar;set(gca,'fontsize',ft);title(cc2, 'm/yr','fontsize',ft);

for b=1:size(groupnum,2)

hf;
 axes(H(4+b));
    hold on;box on;grid on;
    yyaxis left
    if b==1
    set(gca, 'XTick', 1:12, 'XTickLabel', 1:12, 'YTick', 1:0.5:2.5, 'YTickLabel',1:0.5:2.5, 'YColor', 'k','fontsize',ft);
  ylabel('S_t_r_a_j (m mon^-^1)');
    elseif b==2
set(gca, 'XTick', 1:12, 'XTickLabel', 1:12, 'YTick', 1:0.5:2.5, 'YColor', 'k','fontsize',ft);
    end

set(gca, 'XTickMode', 'manual', 'YTickMode', 'manual');
 hb= plot(1:12,m1.Sm(b,:), 'Marker','none',"Color",'k','LineStyle','-','LineWidth',1);

   set(gca,'fontsize',ft,'xtick',1:1:12,'YColor','k');
xlabel('Month');

ym=find(m1.Sm(b,:)==max(m1.Sm(b,:)));
ytop=m1.Sm(b,ym);
if pValue1(b)<0.01
text(ym,ytop+.1,['r = ' num2str(R1(b),'%.2f') ' , p < 0.01'],"Color",'b', 'FontSize', ft);
else
text(ym,ytop+.1,['r = ' num2str(R1(b),'%.2f') ' , p = ' num2str(pValue1(b),'%.2f')],"Color",'b', 'FontSize', ft);
end

   yyaxis right
   e1 =  plot(1:12,SUB_miu(b,:)/12, 'Marker','none',"Color",'b','LineStyle','-','LineWidth',1);
set(gca,'fontsize',ft,'YColor','b');
if b==2
ylabel('S^o^x (Tmol mon^-^1)');end
xlim([0.5 12.5]);

end


fn2=[fullfile(projectRoot,'figs','FIGS7_Straj')];
 exportgraphics(hf, [fn2 '.png'], 'Resolution', 450);

