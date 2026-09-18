%% Generate FIGS9_SeaonalLines
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


load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN','BASIN1','groupnum');

load(fullfile(functionDir,'Areaall.mat'));

m0=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
OSp=m0.OS;oceanlon=m0.oceanlon;oceanlat=m0.oceanlat;
OSp=reshape(OSp,[720 360 12 14]);

for g=1:size(groupnum,1)
    bnum=groupnum{g,1};
    for y=1:14
    for s=1:12
        os=OSp(:,:,s,y);
    data=[];
for b=1:size(bnum,1)
data=[data;oceanlon(BASIN1==bnum(b)) oceanlat(BASIN1==bnum(b)) os(BASIN1==bnum(b)) Areaall(BASIN1==bnum(b)) BASIN1(BASIN1==bnum(b))];
end

 line1=data(:,3); line2=data(:,4);
 subarea=sum(line2(line1>0),'omitnan')./sum(line2,'omitnan');
obarea=sum(line2(line1<0),'omitnan')./sum(line2,'omitnan');
regions_AreaPPT(s,g,y)= subarea;
 subsum=sum(line2(line1>0).*line1(line1>0),'omitnan')./10^6;
obsum=sum(line2(line1<0).*line1(line1<0),'omitnan')./10^6;
regions_NET(s,g,y)= (subsum+obsum); % T mol /mon

 regions_net(s,g,y)=sum(line2.*line1,'omitnan')./sum(line2(~isnan(line2.*line1)),'omitnan'); % mol m-2 mon-1

    end
    end
end

regions_AreaPPT_miu=mean(regions_AreaPPT,3,'omitnan');regions_AreaPPT_sigm=std(regions_AreaPPT,0,3,'omitnan');
 regions_NET_miu=mean(regions_NET,3,'omitnan');regions_NET_sigm=std(regions_NET,0,3,'omitnan');
 regions_net_miu=mean(regions_net,3,'omitnan');regions_net_sigm=std(regions_net,0,3,'omitnan');

regions_AreaPPT1=permute(regions_AreaPPT,[2 1 3]);
regions_NET1=permute(regions_NET,[2 1 3]);

regions_AreaPPT1=reshape(regions_AreaPPT1,[5 12*14]);
regions_NET1=reshape(regions_NET1,[5 12*14]);

for g=1:5
[a, b] = corrcoef(regions_AreaPPT1(g,:)', regions_NET1(g,:)');

R(g,1)=a(1,2);pValue(g,1)=b(1,2);
r2(g,1)=a(1,2).^2;
end

ft=11;

subnum=[2,4,5,3,1];
NO={'a)','b)','c)','d)','e)'};

h=figure;
for g=1:5
subplot(3,2,subnum(g));
hold on;box on;
yyaxis left
e1=errorbar([1:12],regions_NET_miu(:,g),regions_NET_sigm(:,g),'Marker','none',...
    "MarkerEdgeColor",'k',"MarkerFaceColor",'k','LineStyle','-','LineWidth',1);
ylabel({' Tmol mon^-^1'});
e1.Color = 'k';
e1.CapSize = 2.5;
axis padded
set(gca,'YColor','k','xlim',[0.5 12.5],'fontsize',ft);
if g==1 || g==2 ||g==4 || g==5
    plot([0.5 12.5],[0 0],'k--');
end

yyaxis right
e2=errorbar([1:12],regions_AreaPPT_miu(:,g)*100,regions_AreaPPT_sigm(:,g)*100,'Marker','none',...
    "MarkerEdgeColor",'b',"MarkerFaceColor",'b','LineStyle','-','LineWidth',1);
e2.Color = 'b';
e2.CapSize = 2.5;
axis padded
set(gca,'YColor','b','xlim',[0.5 12.5],'fontsize',ft);
ylabel('%');xlabel('Month');

xticks([1:1:12]);
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'});

title([NO{subnum(g)} ' ' groupnum{g,5}],'fontsize',ft,'FontWeight','bold');

text(6,55,['r = ' num2str(R(g),'%.2f')],"Color",'b','fontsize',ft);
end

legend([e1,e2],{'Integrated S^o^x','Subduction Area'},'Box','off','fontsize',ft);


fn=[fullfile(projectRoot,'figs','FIGS9_SeasonalLines.png')];
 exportgraphics(h, fn, 'Resolution', 450);

