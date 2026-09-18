%% Generate FIGS8_ENSO
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
load(fullfile(projectRoot,'input_data','nino34.mat'));

  ifenso=zeros(172,1);
for i=3:length(nino34)-2
data=[nino34(i-2) nino34(i-1) nino34(i) nino34(i+1) nino34(i+2)] ;

if sum(data>0.5)==5
    ifenso(i-2:i+2)=1;
elseif sum(data<-0.5)==5
 ifenso(i-2:i+2)=-1;
else
     ifenso(i)=0;
end
end

ifenso=ifenso(3:170,:);
nino34=nino34(:,3:170);

load(fullfile(functionDir,'PO_subbasins.mat'),'BASIN');

load(fullfile(functionDir,'Areaall.mat'));

m0=matfile(fullfile(projectRoot,'result_data','OSresults_monthly.mat'));
OSp=m0.OS;
oceanlon=m0.oceanlon;
oceanlat=m0.oceanlat;

for i=1:168
os=OSp(:,:,i);area=Areaall;
os(BASIN~=2)=nan;  os(oceanlat>2.5 | oceanlat<-2.5 | oceanlon<170 | oceanlon>280)=nan; % 2.5-2.5 170E-80W
area(isnan(os))=nan;
os(isnan(area))=nan;

fullseries(i,1)=sum(os.*area,'all','omitnan')/10^6;

end

repeatseries=repmat(mean(reshape(fullseries,[12 14]),2),[14 1]);
deSox=fullseries-repeatseries;

cutoff_period = 24;
fs_monthly = 12;
cutoff_frequency = 1 / (cutoff_period / 12);
Wn = cutoff_frequency / (fs_monthly / 2);
fprintf('Filter parameters:\n');
fprintf('Wn = %.6f\n', Wn);

order = 4;
[b, a] = butter(order, Wn, 'low');

x = deSox;

y = filtfilt(b, a, x);

[a, b] = corrcoef(y, nino34);
R2=a(1,2);
pValue=b(1,2);

ft=12;
year=[2004:1/12:2018-1/12];
h=figure;

subplot(2,1,1);hold on;box on;
yyaxis left
plot([2004 2018],[0 0],'k--');

for i=1:length(ifenso)
v2 = [year(i)-1/24 -40; year(i)-1/24 40; year(i)+1/24 40; year(i)+1/24 -40];
f2 = [1 2 3 4];

if ifenso(i)==1
patch('Faces',f2,'Vertices',v2,'FaceColor','r','FaceAlpha',.3,'EdgeColor','none');
elseif ifenso(i)==-1
 patch('Faces',f2,'Vertices',v2,'FaceColor','b','FaceAlpha',.3,'EdgeColor','none');
end
end
p1(1)=plot([2004:1/12:2017+11/12],x,'color',[.4 .4 .4],'LineStyle','-','LineWidth',1);
p1(2)=plot([2004:1/12:2017+11/12],y,'k-','LineWidth',1.5);
ylabel({['Deseasonalized F^o^x'];['(Tmol mon^-^1)']}); xlabel('Year');
set(gca,'YColor','k','ylim',[-40 40],'xlim',[2004 2018],'fontsize',ft);

yyaxis right
p1(2)=plot(2004:1/12:2018-1/12,nino34,'r-','LineWidth',1);
ylabel(['Nino3.4 Index']);
set(gca,'YColor','r','ylim',[-3 3],'xlim',[2004 2018],'fontsize',ft);

text(2005,0,['r = ' num2str(R2,'%.2f') ' , p < 0.01'],'fontsize',ft);

cutoffs = [12, 18, 24, 30, 36, 48];
n_cut = length(cutoffs);
rs = zeros(n_cut, 1);
fs_monthly = 12;
for k = 1:n_cut
    cutoff_period = cutoffs(k);
    cutoff_frequency = 1 / (cutoff_period / 12);
    Wn = cutoff_frequency / (fs_monthly / 2);
    [b, a] = butter(4, Wn, 'low');
    y = filtfilt(b, a, deSox);
    R = corrcoef(y, nino34);
    rs(k) = R(1,2);
end

h;
subplot(2,2,3);grid on;
hold on; box on;

plot(cutoffs, rs, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('Low-pass filter cutoff period (months)');
ylabel('Correlation coefficient (r)');
line([24 24], [0.4, 1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
text(24, min([0.4, 1])+0.05, 'Original cutoff (24 months)', 'HorizontalAlignment', 'center', 'Color', 'r','fontsize',ft);
set(gca,'xlim',[min(cutoffs)-2, max(cutoffs)+2],'ylim',[0.4, 1],'fontsize',ft);

cutoff_period = 24;
fs_monthly = 12;
cutoff_frequency = 1 / (cutoff_period / 12);
Wn = cutoff_frequency / (fs_monthly / 2);

orders = [2, 4, 6];
n_ord = length(orders);
rs = zeros(n_ord, 1);

for k = 1:n_ord
    [b, a] = butter(orders(k), Wn, 'low');
    y = filtfilt(b, a, deSox);
    R = corrcoef(y, nino34);
    rs(k) = R(1,2);
end

h;
subplot(2,2,4);
hold on; box on;grid on;

bar(orders, rs, 'FaceColor', [0.6 0.6 0.8]);
xlabel('Butterworth filter order');
ylabel('Correlation coefficient (r)');
set(gca,'ylim',[0.5 1],'fontsize',ft);
for k = 1:n_ord
    text(orders(k), rs(k)+0.02, sprintf('%.3f', rs(k)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold','fontsize',ft-1);
end
line([0.5 6.5], [rs(2) rs(2)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
legend('Correlation', sprintf('4th-order reference (r=%.3f)', rs(2)), 'Location', 'best','fontsize',ft-1,'box','off');


fn=[fullfile(projectRoot,'figs','FIGS8_ENSO.png')];
exportgraphics(h,fn, 'Resolution', 450);

