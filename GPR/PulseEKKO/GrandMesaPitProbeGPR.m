%% Grand Mesa Pit, Probe, and QPGPR Data Anaylsis
clear;% close all; clc;
workDir = pwd;
dataDir = 'D:\git-repository\GrandMesaMarch2019\';
QPdataDir = 'D:\GrandMesaGPR\';
addpath '.\functions';
addpath '.\colormaps';
set(0,'DefaultAxesFontName','Serif')
set(0,'DefaultTextFontName','Serif')
yetBlack = load('yetBlack.txt');
%% Load Data
% Load All Depth Probe Data
probeFilename = 'GM_March2019_depths_sitevar_11Jul2019.txt';
probeData = readtable([dataDir,probeFilename]);
% Convert Lon,Lat to UTM
[Easting,Northing] = deg2utm(probeData.Latitude,probeData.Longitude);
probeData.Easting = Easting; probeData.Northing = Northing;
probeData = probeData(:,[[1:4],[18:19],[5:17]]);
% Load Summaray Pit Data
pitFilename = 'PitSummary.csv';
pitData = readtable([dataDir,pitFilename]);
% Load All GPR QP Data
qpFilename = 'PulseEKKO_QP_25March2019.csv';
qpData25 = readtable([QPdataDir,qpFilename]);
qpFilename = 'PulseEKKO_QP_26March2019.csv';
qpData26 = readtable([QPdataDir,qpFilename]);
qpFilename = 'PulseEKKO_QP_27March2019.csv';
qpData27 = readtable([QPdataDir,qpFilename]);
qpFilename = 'PulseEKKO_QP_28March2019.csv';
qpData28 = readtable([QPdataDir,qpFilename]);
% Concatentate Daily Files
qpData = [qpData25;qpData26;qpData27;qpData28];
clear('qpData25','qpData26','qpData27','qpData28')

%% Search for Probe Depths
qpEasting = qpData.Easting;
qpNorthing = qpData.Northing;
qpDepth = qpData.Depth;
gprDepth = zeros(size(probeData,1),1);distDepth = gprDepth;
probeDepth = probeData.Depthcm;prbDepth = gprDepth;
for kk = 1:size(probeData,1)
    tmp = sqrt((Easting(kk)-qpEasting).^2+(Northing(kk)-qpNorthing).^2);
    tmpPrb = sqrt((Easting(kk)-Easting).^2+(Northing(kk)-Northing).^2);
    qpIx = find(tmp<3);
%     prbIx = find(tmpPrb<10);
    [~,prbIx] = sort(tmpPrb);
    prbIx = prbIx(1:10);
%     prbIx = prbIx(1:5);
    prbDepth(kk) = mean(probeDepth(prbIx));
% prbDepth(kk)  = probeDepth(kk);
    w = 1./tmp(qpIx).^2;
    gprDepth(kk) = sum(w.*qpDepth(qpIx))./sum(w); 
    distDepth(kk) =  sum(w.*tmp(qpIx))./sum(w); 
end
nanIx = isnan(gprDepth);
gprDepth(nanIx) = [];
prbDepth(nanIx) = [];
distDepth(nanIx) = [];
zedIx = find(prbDepth == 0);
prbDepth(zedIx) = [];
gprDepth(zedIx) = [];
distDepth(zedIx) = [];
% Correlation
[R,P] = corrcoef(gprDepth,prbDepth);
% RMSE
RMSEdepth = sqrt(mean((gprDepth-prbDepth).^2));
MAEdepth = mean(abs(gprDepth-prbDepth));
% binDepth(length(binDepth)) = [];
% prbDepth(length(prbDepth)) = [];
% Scatter
b = colormap(bone);
b = b(1:63,:);
figure();
subplot(1,2,1)
plot([100,300],[100,300],'k','linewidth',2);hold on;
scatter(prbDepth,gprDepth,15,distDepth,'filled');colormap((b));caxis([0,3])
c = colorbar; c.Location = 'northoutside';c.Label.String = 'Distance from Probe (m)';
c.FontSize = 12; c.Label.FontSize = 12;
xlim([100,300]);ylim([100,300]);
axis square
xlabel('Probe Depth (cm)')
ylabel('GPR Depth (cm)')
annotation('textbox',[.31,.225,.1,.1],'FitBoxToText','on','linestyle','none','string',['       R^2 = ',num2str(round(R(1,2).^2,2))],'fontweight','bold')
annotation('textbox',[.31,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['RMSE = ',num2str(round(RMSEdepth,1))],'fontweight','bold')
annotation('textbox',[.18,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['MAE = ',num2str(round(MAEdepth,1))],'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','ytick',[100,200,300])


%% SWE
pitEasting = pitData.Easting;
pitNorthing = pitData.Northing;
gprSWE = zeros(size(pitData,1),1);
qpSWE = qpData.SWE; distSWE = gprSWE;

for kk = 1:size(pitData,1)
        tmp = sqrt((pitEasting(kk)-qpEasting).^2+(pitNorthing(kk)-qpNorthing).^2);
            qpIx = find(tmp<10);
    w = 1./tmp(qpIx).^2;
    distSWE(kk) =  sum(w.*tmp(qpIx))./sum(w); 
    gprSWE(kk) = sum(w.*qpSWE(qpIx))./sum(w);  
end

nanIx = isnan(gprSWE);
gprSWE(nanIx) = [];
distSWE(nanIx) = [];
pitSWE = pitData.SWE_mm_;
pitSWE(nanIx) = [];
% Correlation
[R,P] = corrcoef(gprSWE,pitSWE);
% RMSE
RMSESWE = sqrt(mean((gprSWE-pitSWE).^2));
MAESWE = mean(abs(gprSWE-pitSWE));

% Scatter
% figure();
subplot(1,2,2)
plot([400,1000],[400,1000],'k','linewidth',2);hold on;
scatter(pitSWE,gprSWE,35,distSWE,'filled');colormap((b));
xlim([400,1000]);ylim([400,1000]);caxis([0,10]);%caxis([min(distSWE),max(distSWE)]);
c = colorbar; c.Location = 'northoutside';c.Label.String = 'Distance from Pit (m)';
c.FontSize = 12; c.Label.FontSize = 12;c.Ticks = [0:5:10];
axis square
xlabel('Pit SWE (mm)')
ylabel('GPR SWE (mm)')
set(gca,'fontsize',12,'fontweight','bold','xtick',[400:200:1000])
annotation('textbox',[.75,.225,.1,.1],'FitBoxToText','on','linestyle','none','string',['       R^2 = ',num2str(round(R(1,2).^2,2))],'fontweight','bold')
annotation('textbox',[.75,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['RMSE = ',num2str(round(RMSESWE,1))],'fontweight','bold')
annotation('textbox',[.62,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['MAE = ',num2str(round(MAESWE,1))],'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold')


% daspect([1,1,1])
% uistack(h,'bottom')

figure();
scatter(qpEasting,qpNorthing,15,qpSWE,'filled');colormap(b)
xlabel('Easting (km)')
ylabel('Northing (km)')
c = colorbar; c.Location = 'northoutside';c.Label.String = 'SWE (mm)';
c.FontSize = 12; c.Label.FontSize = 14;caxis([100,1000])
xlim([740.000,749.500])
ylim([4321.000,4327.500])
hold on
I = imread([QPdataDir,'pit_locations_clusters_crop.png']); 
Eaxis = linspace(740.000,749.500,size(I,2));
Naxis = linspace(4321.000,4327.500,size(I,1));
% h = imagesc(Eaxis,Naxis,I); 
image(Eaxis,Naxis,flipud(I))
% uistack(h,'bottom')
scatter(qpEasting./1000,qpNorthing./1000,2,qpSWE,'filled');colormap(b)
set(gca,'fontsize',14,'fontweight','bold')
daspect([1,1,1])
% ax = gca;
% ax.XRuler.Exponent = 0;
% ax.YRuler.Exponent = 0;

figure();
scatter(qpEasting,qpNorthing,15,qpDepth,'filled');colormap(b)
xlabel('Easting (km)')
ylabel('Northing (km)')
c = colorbar; c.Location = 'northoutside';c.Label.String = 'Depth (cm)';
c.FontSize = 12; c.Label.FontSize = 14;caxis([50,300])
xlim([740.000,749.500])
ylim([4321.000,4327.500])
hold on
I = imread([QPdataDir,'pit_locations_clusters_crop.png']); 
Eaxis = linspace(740.000,749.500,size(I,2));
Naxis = linspace(4321.000,4327.500,size(I,1));
% h = imagesc(Eaxis,Naxis,I); 
image(Eaxis,Naxis,flipud(I))
% uistack(h,'bottom')
scatter(qpEasting./1000,qpNorthing./1000,2,qpDepth,'filled');colormap(b)
set(gca,'fontsize',14,'fontweight','bold')
daspect([1,1,1])