%% Process Raw .GPS
% Tate Meehan
clear;close all;clc;
%% Read .GPS data file
dataDir = 'C:\Users\snowfield\Desktop\GrandMesaGPR\PulseEKKO_26March2019\PulseEKKO\';
contents = dir(dataDir);
gpsfilenames = dir(fullfile(dataDir,'*.GPS'));
% [filename,pathname] = uigetfile([dataDir,'*.GPS'] );
% Open the Selected File
% fid=fopen([pathname,filename]);
fid=fopen([gpsfilenames.folder,gpsfilenames.name]);

% Read the whole lines of .GPS file into a cell array
rawGPS=textscan(fid,'%s','Delimiter','\n');
% Get the Indicies of the important .GPS Rows
hdix = 1:3:length(rawGPS{1});
GPGGAix = 2:3:length(rawGPS{1});
GPZDAix = 3:3:length(rawGPS{1});
for ii = 1:length(hdix)
% Get the Date and Time of the Trace
% GPZDA string
tmp = rawGPS{1}{GPZDAix(ii)};
time(ii,:) = str2num(tmp(8:13));
date(ii,:) = datestr(datenum(tmp([18,19,21,22,24,25,26,27]),'ddmmyyyy'),'mm/dd/yyyy');

% Get the Trace Number and Position
tmp = strsplit(rawGPS{1}{hdix(ii)});
trcno(ii,:) = str2num(tmp{2}(2:end));
posx(ii,:) = str2num(tmp{5});

% Get the Longitude DDDMM.MMMMM
tmp = rawGPS{1}{GPGGAix(ii)};
lon(ii,:) = str2num(tmp(32:34));
lonm(ii,:) = str2num(tmp(35:43));
% Convert to Decimal Degees
lonm(ii,:) = lonm(ii)./60;
lon(ii,:) = lon(ii)+lonm(ii);
% Check Sign
if (tmp(45)) == 'W'
    lon(ii,:) = -lon(ii,:); 
end

% Get the Latitude DDMM.MMMMM
lat(ii,:) = str2num(tmp(18:19));
latm(ii,:) = str2num(tmp(20:28));
% Convert to Decimal Degees
latm(ii,:) = latm(ii)./60;
lat(ii,:) = lat(ii)+latm(ii);
% Check Sign
if (tmp(30)) == 'S'
    lat(ii,:) = -lat(ii,:); 
end
end
% Convert to UTM
[x,y,utmzone] = deg2utm(lat,lon);
% Moving Least Squares Interpolation



