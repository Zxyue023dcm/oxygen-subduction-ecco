%% Generate FIGS7_TimeSeries
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
OSp=m0.OS;
oceanlon=m0.oceanlon;oceanlat=m0.oceanlat;

Station=[391,211;506,184;413,92];
SNum={'a)','b)','c)'};
for i=1:size(Station,1)
    lon=round(oceanlon(Station(i,1),Station(i,2)));
     lat=round(oceanlat(Station(i,1),Station(i,2)));
     if lat>=0 && lon<=180
LonLat{i,1}=[num2str(lat) '°N, ' num2str(lon) '°E'];
     elseif lat>=0 && lon>180
         lon=360-lon;
LonLat{i,1}=[num2str(lat) '°N, ' num2str(lon) '°W'];
     elseif lat<0 && lon<=180
         lat=abs(lat);
LonLat{i,1}=[num2str(lat) '°S, ' num2str(lon) '°E'];
     elseif lat<0 && lon>180
         lat=abs(lat);
         lon=360-lon;
LonLat{i,1}=[num2str(lat) '°S, ' num2str(lon) '°W'];
     end

end

[all_themes, all_colors] = GetColors();
cls=all_colors(5,:);
ft=11;

h=figure;
for i=1:size(Station,1)
subplot(size(Station,1),1,i);hold on;box on;
TSeries=squeeze(OSp(Station(i,1),Station(i,2),:));
p(1)=plot([2004:1/12:2018-1/12],TSeries,'LineStyle','-','color',cls,'LineWidth',1);
RSeries=repmat(mean(reshape(TSeries,[12 14]),2),[14 1]);
p(2)=plot([2004:1/12:2018-1/12],RSeries,'LineStyle','--','color',cls,'LineWidth',1);
ylabel({['F^o^x'];[ '(mol m^-^2 mon^-^1)']});
xlabel('Year');
title([SNum{i} ' ' LonLat{i} '' ],'FontSize',ft,'FontWeight','bold');
set(gca,'fontsize',ft);
end
legend(p,{'Original series','Seasonal cycles'},'fontsize',ft,'box','off');

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');

OSp(repmat(BASIN,[1 1 168])~=2)=nan;

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

   for  i=1:size(Station,1)
       r(i)=R(Station(i,1),Station(i,2));
 r2(i)=R(Station(i,1),Station(i,2)).^2;
 pv(i)=pValue(Station(i,1),Station(i,2));
   end

   h;
   for i=1:size(Station,1)
subplot(size(Station,1),1,i);hold on;box on;
text(2006,0,['r^2_s_e_a_s_o_n = ' num2str(r2(i),'%.2f') ' , p < 0.01'],'fontsize',ft);

   end
fn=[fullfile(projectRoot,'figs','FIG8_TimeSeries.png')];%h=figure(3);
 exportgraphics(h, fn, 'Resolution', 450);

