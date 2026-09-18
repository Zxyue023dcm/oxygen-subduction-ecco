%% Generate FIG4_eddy_regions
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
m1=matfile(fullfile(functionDir,'PO_subbasins.mat'));

load(fullfile(functionDir,'Areaall.mat'));
Areaall(m1.BASIN~=2)=nan;


m=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
oceanlon=m.oceanlon;oceanlat=m.oceanlat;
OS=m.OSeddy;
OS=mean(reshape(OS,[720 360 12 14]),[3 4],'omitnan');
OS(m1.BASIN~=2)=nan;

OS_sub=OS;OS_sub(OS<-0.05)=nan;
OS_ob=OS;OS_ob(OS>0.05)=nan;

  mf=matfile(fullfile(projectRoot,'result_data','BASINed_maxR.mat'));
R=[4 2 1];
EddyRegion=mf.EddyRegion;BASIN=mf.BASIN_ed;
for r=1:length(R)
EddyLine(:,:,r)=mf.EddyLine(:,:,R(r));
Num{r}=EddyRegion{R(r)};
end

load(fullfile(projectRoot,'result_data','seas_corr_maxR.mat'),'Main','maxR','lon','lat');

Main_sub=Main;Main_sub(OS<0)=nan;
Main_ob=Main;Main_ob(OS>0)=nan;

Main_sub=griddata(oceanlon,oceanlat,Main_sub,lon,lat,'nearest');
Main_ob=griddata(oceanlon,oceanlat,Main_ob,lon,lat,'nearest');

colorlist=slanCM('pride');
    mymap_sub = custom_colormap(colorlist(129:end,:), 'apply', 'off');
       mymap_ob = custom_colormap(colorlist(1:128,:), 'apply', 'off');
co=2;
ci=0.5;
ft=12;

H=figure;
c1=subplot(3,2,1);
      m_proj('miller','lat',[-60 62],'long',[122 292]);         hold on
       m_pcolor(oceanlon,oceanlat,OS_sub); shading interp;
        clim([0 co]); % colormap(c1,slanCM('pride'));
        colormap(c1,mymap_sub);
              DOTS=[lon(Main_sub==4) lat(Main_sub==4)];
          m_scatter(DOTS(:,1),DOTS(:,2),2,'b','filled');
        m_coast('patch',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
        m_grid('linestyle','none','tickdir','both','fontsize',ft,'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135:45:290],'ytick',[-60 -30 0 30 60]);
title(['a) Eddy-induced S^o^x'],'fontsize',ft,'FontWeight','bold');
 cl=colorbar;cl.Label.String='mol O_2 m^-^2 mon^-^1';cl.FontSize=ft;
 cl.Ticks=[0:ci:co];

a=0.5;
EddyLine1([2 3],2,:)=EddyLine([2 3],2,:)-a;
EddyLine1([1 4 5],2,:)=EddyLine([1 4 5],2,:)+a;
EddyLine1([1 2 5],1,:)=EddyLine([1 2 5],1,:)-a;
EddyLine1([3 4],1,:)=EddyLine([3 4],1,:)+a;

for n=1:3
hold on;
for i=1:4
m_line(EddyLine1(i:i+1,1,n),EddyLine1(i:i+1,2,n),'color','k','linewidth',1);
end
if n==3
    m_text(EddyLine1(4,1,n)-2,EddyLine1(1,2,n)-3,Num{n}, ...
        'color','k','fontsize',10,'HorizontalAlignment','right');
else
    m_text(EddyLine1(1,1,n)-8,EddyLine1(1,2,n)-3,Num{n}, ...
        'color','k','fontsize',10);
end
end

c2=subplot(3,2,2);
      m_proj('miller','lat',[-60 62],'long',[122 292]);         hold on
       m_pcolor(oceanlon,oceanlat,OS_ob); shading interp;
          clim([-co 0]); % colormap(c2,slanCM('pride'));
           colormap(c2,mymap_ob);
          DOTS=[lon(Main_ob==4) lat(Main_ob==4)];
          m_scatter(DOTS(:,1),DOTS(:,2),2,'b','filled');
        m_coast('patch',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
        m_grid('linestyle','none','tickdir','both','fontsize',ft,'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135:45:290],'ytick',[-60 -30 0 30 60]);
title(['b) Eddy-induced O^o^x'],'fontsize',ft,'FontWeight','bold');
 cl=colorbar; cl.Label.String='mol O_2 m^-^2 mon^-^1';cl.FontSize=ft;
 cl.Ticks=[-co:ci:0];


    mf=matfile(fullfile(projectRoot,'result_data','BASINed_maxR.mat'));

EddyRegion=mf.EddyRegion;BASIN=mf.BASIN_ed;
R=[4 2 1];
for r=1:length(R)
EddyLine(:,:,r)=mf.EddyLine(:,:,R(r));
Num{r}=EddyRegion{R(r)};
OSbar(:,:,r)=mf.OSbar(:,:,R(r));
OSbar_ed(:,:,r)=mf.OSbar_ed(:,:,R(r));
end

RR=mf.R;PV=mf.PV;
for r=1:length(R)
rr=RR{R(r)};pv=PV{R(r)};
for i=1:4
rr1=rr(pv(:,i)<0.05,i);
 corr_data{r,i}=rr1.^2;
r_m(r,i)=mean(rr1.^2,'omitnan');
end
end
ft=12;

[all_themes, all_colors] = GetColors();
 cls(1,:)=[.4 .4 .4]; cls(2,:)=all_colors(14,:); cls(3,:)=all_colors(13,:);
cls(4,:)=all_colors(17,:); cls(5,:)=all_colors(16,:);

n_sites = length(Num);

term_names = {'Zonal La.','Merid. La.','Vert. Vel.','Edd.'};
n_terms = length(term_names);

offset = linspace(-0.3, 0.3, n_terms);
colors = [all_colors(14,:);
          all_colors(13,:);
          all_colors(17,:);
          all_colors(16,:)];
box_alpha = 0.5;

H;
subplot(3,1,2);
hold on; grid on;

box_handles = gobjects(n_sites, n_terms);
for i_site = 1:n_sites
    for i_term = 1:n_terms
        x_pos = i_site + offset(i_term);
        y_data = corr_data{i_site, i_term};
        h(i_term) =  boxplot_custom(y_data,'linecolor','k','fillcolor',colors(i_term,:),...
              'x_position',x_pos, 'Width', 0.18,...
              'outliercolor',colors(i_term,:),'add_mean',1,...
              'mean_marker','o','mean_face_color',colors(i_term,:),...
             'mean_edge_color','k','mean_size',10);

    end
end

xlim([0.5, n_sites+0.5]);
xticks(1:n_sites);
xticklabels(Num);
ylabel({'Squared';'correlation coefficients (r^2)'; 'with total F^o^x'}, 'FontSize', ft);

legend_handles = gobjects(1, n_terms);
for i = 1:n_terms
    legend_handles(i) = patch(NaN, NaN, colors(i,:), ...
                               'EdgeColor', 'k', ...
                               'DisplayName', term_names{i});
end

set(gca, 'FontSize', ft, 'Box', 'on', 'TickDir', 'out');
grid on;
hold off;

mon=1:12;

H;

for g=1:size(EddyLine,3)

   subplot(3,3,g+6);
    hold on;box on;
hb=bar(mon,OSbar_ed(2:5,:,g));
for c=1:4
   hb(c).FaceColor=cls(c+1,:);hb(c).EdgeColor='none';
end
h1=plot(mon,OSbar_ed(1,:,g),"Color",cls(1,:),'LineStyle','-','Marker','*','LineWidth',1);
set(gca,'xlim',[0.5 12.5],'XTick',1:1:12,'TickDir','both','fontsize',ft);%, ...

    if g==1
        ylabel({'F^o^x','(Tmol O_2 mon^-^1)'},'fontsize',ft);
    end
    xlabel('Month','fontsize',ft)
end
leg1 = legend([h1 hb],{'Total','Zonal La.','Merid. La.','Vert. Vel.','Edd.'},'Orientation','horizon', 'FontSize', ft-1);
leg1.AutoUpdate = 'off';

fn=[fullfile(projectRoot,'figs','FIG4_EddyRegions.png')];
exportgraphics(H,fn, 'Resolution', 450);

