%% Generate FIGS3_o2_uncer
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


OXY={'oxy_in_mml','oxy_in_mml10'}; clb_o=[0 300];
UNC={'unc_in_mml','unc_in_mml10'}; clb_u=[0 25];
Tit={'above the H_m_a_x','10m below the H_m_a_x'};

load(fullfile(projectRoot,'result_data','UncO2_Hmax.mat'));

m0=matfile(fullfile(functionDir,'PO_subbasins.mat'));
BASIN=griddata(m0.LON,m0.LAT,m0.BASIN,double(LON),double(LAT),'nearest');

ft=11;
hf=figure;H=tight_subplot(2,3, [0.04 0.04], [0.06 0.04], [0.08 0.08]);

for i=1:2
 axes(H(i));
 axis off
   m_proj('miller','lat',[-60 62],'long',[122 292]);
         hold on ; eval(['data=mean(' OXY{i} '(:,:,1:168),3,"omitnan");']);
data(BASIN~=2)=nan;
       m_pcolor(LON,LAT,data); shading interp;
         clim(clb_o); colormap(H(i),slanCM('haline'));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'fontsize',8,'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',90:45:300,'ytick',[-60 -35 -10  10 35 60],'fontsize',ft);
   title(Tit{i},'fontsize',ft+1);

 axes(H(i+3));
 axis off
   m_proj('miller','lat',[-60 62],'long',[122 292]);
         hold on ; eval(['data=mean(' UNC{i} '(:,:,1:168),3,"omitnan");']);
data(BASIN~=2)=nan;
       m_pcolor(LON,LAT,data); shading interp;
      clim(clb_u); colormap(H(i+3),slanCM('haline'));
        m_coast('patch',[.7 .7 .7],'edgecolor',[.3 .3 .3]);
        m_grid('linestyle','none','tickdir',[],'fontsize',8,'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',90:45:300,'ytick',[-60 -35 -10  10 35 60],'fontsize',ft);

end

for i=1:2
axes(H(3*i))
axis off
if i==1
colormap(H(3*i),slanCM('haline'));
clim(clb_o);
elseif i==2
colormap(H(3*i),slanCM('haline'));
clim(clb_u);
end
cbar(i) = colorbar;
cbar(i).Label.FontSize = ft;
cbar(i).FontSize = 10;
cbar(i).YAxisLocation = 'right';
cbar(i).Label.HorizontalAlignment = 'right';
if i==1
cbar(i).Label.String = {'GOBAI-O_2 (μmol kg^-^1)'};
cbar(i).Label.Position(1) = 2.5;
cbar(i).Label.Position(2) = 200;
elseif i==2
   cbar(i).Label.String = {'Uncertainty (μmol kg^-^1)'};
cbar(i).Label.Position(1) = 2.5;
cbar(i).Label.Position(2) = 17;
end
pos = get(H(3*i), 'Position');

cb_width  = 0.1 * pos(3);  % 5% of subplot width
cb_height = 0.9 * pos(4);   % 60% of subplot height
cb_left   = pos(1) + 0.01 * (pos(3) - cb_width);   % Center horizontally
cb_bottom = pos(2) + 0.06 * pos(4);                % Slightly lower
set(cbar(i), 'Position', [cb_left, cb_bottom, cb_width, cb_height]);

end


fn2=[fullfile(projectRoot,'figs','FIGS3_o2_uncer')];
 exportgraphics(hf, [fn2 '.png'], 'Resolution', 450);

