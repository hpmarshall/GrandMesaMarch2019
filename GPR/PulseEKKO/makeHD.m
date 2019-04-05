%% Replicate .HD
% Tate Meehan - Oct. 2018
clear; close all; clc
addpath '/sonichome/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
isWrite = 1;

%% From Data Directory Read and Sort SubDirectories
workingDirectory = pwd;
directory = '/SNOWDATA/NSF_GREENTRACS/GreenTrACS2017/PulseEKKO/500MHz';
folders = dir(directory);
folders(1:2) = []; % Remove Hidden Directories
datearray = cat(1,folders.date);
daymonthyear = (datearray(:,1:11));
num = datenum(daymonthyear,'dd-mmm-yyyy');
[~,sortIx] = sort(num);
folders = folders(sortIx);

% Loop Over each day of Acquisition
for ff = 16%:length(folders)%1:length(folders)
    dataDir = [directory,'/',folders(ff).name];
    % Get .HD Files
    hdfilenames = dir(fullfile(dataDir,'*.HD'));
    hdfiles = struct2cell(hdfilenames);
    hdfiles(2:end,:) = [];
    hdfiles = cell2mat(hdfiles');
    % This IS Robust for 2016 w/GPS files
    % Get .DT1 Files
    dt1filenames = dir(fullfile(dataDir,'*.DT1'));
    dt1files = struct2cell(dt1filenames);
    dt1files(2:end,:) = [];
    dt1files = cell2mat(dt1files');
    dt1bytes = [dt1filenames.bytes];
    testbytes = dt1bytes == 0;
    if any(testbytes)
        dt1files(find(testbytes)) = [];
    end
    
    if size(dt1files,1) == size(hdfiles,1)
        fprintf([folders(ff).name, ' is OK','\n'])
    else
        hdnum = str2num(hdfiles(:,5:6));
        nhd = hdnum+1;
        dt1num = str2num(dt1files(:,5:6));
        missingHD = dt1num(find(~ismember(dt1num,hdnum)));        
        files = dt1files(:,1:6);
        nFiles = length(dt1num);
%     end
    trhd = cell(nFiles,1);
    hdr1 = cell(nFiles,1);
    for ii = 1:length(hdnum)%hdnum+1
%         tic
        %------------------------------------------------------------------
        % Multiplexed Channel Record
        filename = files(nhd(ii),:);
        
        % Read Data
        filepath = fullfile(dataDir,filename);
        [~,hdr1{nhd(ii)},trhd{nhd(ii)},~,~,~,~,~] = readSensorsSoftwareData( filepath );
         ntrcs(ii) = [size(trhd{nhd(ii)},2)];
    end
    
    % Least-Squares Approximation for Number of Traces in Missing .HD files
    % Not Robust for less than 2 existing files..
       G = [ones(length(hdnum),1),dt1bytes(nhd)'];
       d = ntrcs';
       m = G\d;
       for ii = 1:length(missingHD)
           foundTraces(ii) = m(2).*dt1bytes(missingHD(ii)+1) + m(1);
       end
       
    end
    % This is a Work in Progress, Fixed bug For now with testbytes. 
    % Must Complete .HD File Writing and foundTrace Storage for the case of 
    % multiple missing files and multiple directories.
end