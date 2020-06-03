%% Multi-channel Quad-Pole GPR (MxQP)
clear; close all; clc;
workDir = pwd;
dataDir = 'D:\git-repository\GrandMesaMarch2019\';
QPdataDir = 'D:\GrandMesaGPR\';
addpath '.\functions';
addpath '.\colormaps';
set(0,'DefaultAxesFontName','Serif')
set(0,'DefaultTextFontName','Serif')
yetBlack = load('yetBlack.txt');
%% Load All GPR QP Data
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
% Load Version1 data
qpFilename = 'PulseEKKO_QP_25March2019_v1.csv';
qpData25 = readtable([QPdataDir,qpFilename]);
qpFilename = 'PulseEKKO_QP_26March2019_v1.csv';
qpData26 = readtable([QPdataDir,qpFilename]);
qpFilename = 'PulseEKKO_QP_27March2019_v1.csv';
qpData27 = readtable([QPdataDir,qpFilename]);
qpFilename = 'PulseEKKO_QP_28March2019_v1.csv';
qpData28 = readtable([QPdataDir,qpFilename]);
% Concatentate Daily Files
qpDataV1 = [qpData25;qpData26;qpData27;qpData28];
clear('qpData25','qpData26','qpData27','qpData28')
%% Correlate Depth and SWE
qpEasting = qpData.Easting;
qpNorthing = qpData.Northing;
qpDepth = qpData.Depth;
qpSWE = qpData.SWE;
qpTWT = qpData.TWT;

qpEastingV1 = qpDataV1.Easting;
qpNorthingV1 = qpDataV1.Northing;
qpDepthV1 = qpDataV1.Depth;
qpSWEV1 = qpDataV1.SWE;
qpTWTV1 = qpDataV1.TWT;

% find nearest locations to V1
radius = 3;
gprSWE = zeros(size(qpSWEV1,1),1);
gprDepth = gprSWE; distDepth = gprSWE; gprTWT = gprSWE;
for kk = 1:size(qpSWEV1,1)
    tmp = sqrt((qpEastingV1(kk)-qpEasting).^2+(qpNorthingV1(kk)-qpNorthing).^2);
    qpIx = find(tmp<radius);
% Weighted Average
    w = 1./tmp(qpIx).^2;
    gprSWE(kk) = sum(w.*qpSWE(qpIx))./sum(w);
    gprDepth(kk) = sum(w.*qpDepth(qpIx))./sum(w); 
    distDepth(kk) =  sum(w.*tmp(qpIx))./sum(w); 
    gprTWT(kk) = sum(w.*qpTWT(qpIx))./sum(w); 
end
nanIx = isnan(gprDepth);
gprDepth(nanIx) = [];
gprSWE(nanIx) = [];
qpSWEV1(nanIx) = [];
qpDepthV1(nanIx) = [];
distDepth(nanIx) = [];
distSWE = distDepth;
gprTWT(nanIx) = [];
qpTWTV1(nanIx) = [];
% Correlation
[Rdepth,Pdepth] = corrcoef(gprDepth(:),qpDepthV1);
% RMSE
RMSEdepth = sqrt(mean((gprDepth(:)-qpDepthV1).^2));
MAEdepth = mean(abs(gprDepth(:)-qpDepthV1));

% Correlation
[Rswe,Pswe] = corrcoef(gprSWE(:),qpSWEV1);
% RMSE
RMSESWE = sqrt(mean((gprSWE(:)-qpSWEV1).^2));
MAESWE = mean(abs(gprSWE(:)-qpSWEV1));

% Correlation
[Rtwt,Ptwt] = corrcoef(gprTWT(:),qpTWTV1);
% RMSE
RMSETWT = sqrt(mean((gprTWT(:)-qpTWTV1).^2));
MAETWT = mean(abs(gprTWT(:)-qpTWTV1));
%% Scatter
b = colormap(bone);
b = b(1:63,:);
figure();
subplot(1,2,1)
plot([100,300],[100,300],'k','linewidth',2);hold on;
scatter(qpDepthV1,gprDepth,15,distDepth,'filled');colormap((b));caxis([0,radius])
c = colorbar; c.Location = 'northoutside';c.Label.String = 'Distance from Probe (m)';
c.FontSize = 12; c.Label.FontSize = 12;
xlim([100,300]);ylim([100,300]);
axis square
xlabel('Probe Depth (cm)')
ylabel('GPR Depth (cm)')
annotation('textbox',[.31,.225,.1,.1],'FitBoxToText','on','linestyle','none','string',['       R^2 = ',num2str(round(Rdepth(1,2).^2,2))],'fontweight','bold')
annotation('textbox',[.31,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['RMSE = ',num2str(round(RMSEdepth,1))],'fontweight','bold')
annotation('textbox',[.18,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['MAE = ',num2str(round(MAEdepth,1))],'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','ytick',[100,200,300])

subplot(1,2,2)
plot([400,1000],[400,1000],'k','linewidth',2);hold on;
scatter(qpSWEV1,gprSWE,35,distSWE,'filled');colormap((b));
xlim([400,1000]);ylim([400,1000]);caxis([0,radius]);%caxis([min(distSWE),max(distSWE)]);
c = colorbar; c.Location = 'northoutside';c.Label.String = 'Distance from Pit (m)';
c.FontSize = 12; c.Label.FontSize = 12;c.Ticks = [0:5:10];
axis square
xlabel('Pit SWE (mm)')
ylabel('GPR SWE (mm)')
set(gca,'fontsize',12,'fontweight','bold','xtick',[400:200:1000])
annotation('textbox',[.75,.225,.1,.1],'FitBoxToText','on','linestyle','none','string',['       R^2 = ',num2str(round(Rswe(1,2).^2,2))],'fontweight','bold')
annotation('textbox',[.75,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['RMSE = ',num2str(round(RMSESWE,1))],'fontweight','bold')
annotation('textbox',[.62,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['MAE = ',num2str(round(MAESWE,1))],'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold')

figure();
plot([0,30],[0,30],'k','linewidth',2);hold on;
scatter(qpTWTV1,gprTWT,15,distDepth,'filled');colormap((b));caxis([0,radius])
c = colorbar; c.Location = 'northoutside';c.Label.String = 'Distance from Probe (m)';
c.FontSize = 12; c.Label.FontSize = 12;
xlim([0,30]);ylim([0,30]);
axis square
xlabel('Probe TWT (ns)')
ylabel('GPR TWT (ns)')
annotation('textbox',[.31,.225,.1,.1],'FitBoxToText','on','linestyle','none','string',['       R^2 = ',num2str(round(Rtwt(1,2).^2,2))],'fontweight','bold')
annotation('textbox',[.31,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['RMSE = ',num2str(round(RMSETWT,1))],'fontweight','bold')
annotation('textbox',[.18,.18,.1,.1],'FitBoxToText','on','linestyle','none','string',['MAE = ',num2str(round(MAETWT,1))],'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','ytick',[10,20,30])
