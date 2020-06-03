%% sensorsSoftwareZeroLevel.m
% This code unpacks the binary S&S data file and writes the raw traces and
% meta data to ncdf

% Tate Meehan - Oct. 2018 - Updated April 2019 and March 2020
clear; %close all; clc
% addpath '/sonichome/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
addpath './functions'
isWrite = 1;
isDumpGPS = 1;
isGeode = 0;
gong = load('gong');
% % Create Error Log
%  f = fopen( 'ErrorLog.txt', 'w' );  
%  fclose(f);
% % Create Report Log
%  f = fopen('ReportLog.txt','w');
%% From Data Directory Read and Sort SubDirectories
workingDirectory = pwd;
% Enter the Appropriate Data Directory
% directory = '/SNOWDATA/GrandMesa2019/GPR/';
% directory = 'D:\GrandMesaGPR';
% directory = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\LDP\122019';
% directory = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\GPR';
% ReProcess 2019
directory = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2019\GrandMesaGPR';
% workingDirectory = directory;
% Create Error Log
 f = fopen(fullfile(directory, 'ErrorLog.txt'), 'w' );  
 fclose(f);
% Create Report Log
 f = fopen(fullfile(directory,'ReportLog.txt'),'w');
folders = dir(directory);
lf = length(folders);
% folders(1:2) = []; % Remove Hidden Directories
% folders([1:3,5:end]) = []; % Remove Extra Folders
% folders([1:2,lf-1:lf]) = []; % Remove Extra Folders
% folders(3:5) = []; % Remove CMP directories
% folders(5:end) = []; % 06022020 run
folders([[1:11],[16:lf]]) = [];

datearray = cat(1,folders.date);
daymonthyear = (datearray(:,1:11));
num = datenum(daymonthyear,'dd-mmm-yyyy');
[~,sortIx] = sort(num);
folders = folders(sortIx);

%% Loop Over each day of Acquisition
for ff = 2:length(folders)
%     dataDir = [directory,'/',folders(ff).name,'/PulseEKKO/'];
    dataDir = [directory,'\',folders(ff).name,'\PulseEKKO\'];
%         if strcmp(dataDir,'D:\GrandMesaGPR\PulseEKKO_28March2019\PulseEKKO\')
%     if strcmp(dataDir,'/SNOWDATA/GrandMesa2019/GPR/PulseEKKO_28March2019/PulseEKKO/')
        % Move into one additional directory for this day.
%         dataDir = '/SNOWDATA/GrandMesa2019/GPR/PulseEKKO_28March2019/PulseEKKO/D04/';
%         dataDir = 'D:\GrandMesaGPR\PulseEKKO_28March2019\PulseEKKO\D04\';
        if strcmp(dataDir,'D:\CRREL_SnowCompaction\CRREL\SnowEx2019\GrandMesaGPR\PulseEKKO_28March2019\PulseEKKO\')
            dataDir = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2019\GrandMesaGPR\PulseEKKO_28March2019\PulseEKKO\D04\';
        end
    % This IS Robust for 2016 w/GPS files
    % Get .DT1 Files
    dt1filenames = dir(fullfile(dataDir,'*.DT1'));
    dt1files = struct2cell(dt1filenames);
    dt1files(2:end,:) = [];
    dt1files = cell2mat(dt1files');
    dt1bytes = [dt1filenames.bytes];
    % Remove Files Without Data from Query
    testbytes = dt1bytes == 0;
    if any(testbytes)
        rmDT1ix = find(testbytes);
        for ii = 1:length(any(testbytes))
%             msg = [dataDir,'/',dt1files(rmDT1ix,:),' has 0 bytes.','\n'];
            msg = [dataDir,'\',dt1files(rmDT1ix,:),' has 0 bytes.','\n'];
            warning(msg)
            WarnUser(msg,'ErrorLog.txt',workingDirectory)
        end
        dt1files(rmDT1ix,:) = [];
    end
    % Get .HD Files
    hdfilenames = dir(fullfile(dataDir,'*.HD'));
    hdfiles = struct2cell(hdfilenames);
    hdfiles(2:end,:) = [];
    hdfiles = cell2mat(hdfiles');
    
    % Get .GPS Files ...
    gpsfilenames = dir(fullfile(dataDir,'*.GPS'));
    
    % Construct Files
    files = dt1files(:,1:6);
    nFiles = size(files,1);
    
        if size(dt1files,1) ~= size(hdfiles,1)
            hdnum = str2num(hdfiles(:,5:6));
            nhd = hdnum+1;
            dt1num = str2num(dt1files(:,5:6));
            missingHD = dt1num(find(~ismember(dt1num,hdnum)));
            for ii = 1:length(missingHD)
%                 msg = ['Missing ', dataDir,'/',files(missingHD(ii),:),'.HD','\n'];
                msg = ['Missing ', dataDir,'\',files(missingHD(ii),:),'.HD','\n'];
                warning(msg);
                WarnUser(msg,'ErrorLog.txt',workingDirectory)
            end
        else
        end
       % Load Post Processed GPS
       sound(gong.y,gong.Fs)
       isPostProcessed = questdlg('Load a Post-Processed GPS File?','Post-Processed GPS','Yes','No','Yes');
       if strcmp(isPostProcessed,'Yes')
           disp('Loading Post-Processed GPS')
           tic
           if strcmp(dataDir,'D:\CRREL_SnowCompaction\CRREL\SnowEx2019\GrandMesaGPR\PulseEKKO_28March2019\PulseEKKO\D04\')
               [fname,pname] = uigetfile([dataDir(1:end-14),'GeoXH\*.csv']);
           else
               [fname,pname] = uigetfile([dataDir(1:end-10),'GeoXH\*.csv']);
           end
           %           opts = detectImportOptions([pname,fname])
           ppGPS = readtable([pname,fname]);
           clear('pname','fname')
           tmp = ppGPS.Var7;
         for jj = 1:length(tmp)
            ppTime(jj,:) = (tmp{jj}(10:end));
            % Correct time to seconds of day
            [~,~,~,H,M,S] = datevec(ppTime(jj,:));
            ppSec(jj,:) = H*3600+M*60+S;
            ppDate(jj,:) = datestr(tmp{jj}(1:8),'mm/dd/yyyy');
         end
         clear tmp
         disp('Post-Procesed GPS Loaded')
         toc
         disp(' ')
         
       end
      
    
    for ii = 1 : nFiles
        tic
        %------------------------------------------------------------------
        % Multiplexed Channel Record
        filename = files(ii,:);
        
        % Read Data
        filepath = fullfile(dataDir,filename);
        [Rad,hdr1,trhd,~,f0,~,~,~] = readSensorsSoftwareData( filepath );
        % Append f0 to Trace Header
        trhd(9,:) = ones(1,size(trhd,2)).*f0;
        nChan = max(unique(trhd(23,:)));
        % Nominal Frequency GHz
        f0GHz = f0/1000;
        % No. Traces in Multiplexed Data
        [~, multiplexNtrcs] = size(Rad);
        
        % Read .GPS file
        fid=fopen([dataDir,gpsfilenames(ii).name]);
        % Read the whole lines of .GPS file into a cell array
        rawGPS=textscan(fid,'%s','Delimiter','\n');
        fclose(fid);
        % Get the Indicies of the important .GPS Rows
        hdix = 1:3:length(rawGPS{1});
        GPGGAix = 2:3:length(rawGPS{1});
        GPZDAix = 3:3:length(rawGPS{1});
%         sec = zeros(length(hdix),1); time = sec; date = sec; posx = sec;
%         lat = sec; latm = sec; lon = sec; lonm = sec; z = sec;
        kk = 0;
        for jj = 1:length(hdix)
            % Break Loop if Index Exceeds
            if jj > length(rawGPS{1})/3%jj + kk > length(hdix)
                break
            end
            % Get the Date and Time of the Trace
            % GPZDA string
            tmp = rawGPS{1}{GPZDAix(jj)};
            tmpGPZDA = strsplit([rawGPS{1}{GPZDAix(jj)}],',');
            if strcmp(tmpGPZDA{1},'$GPZDA')
            else
                % Error Catch for Duplicated Trace ID
                rawGPS{1}(GPZDAix(jj)) = [];
                tmp = rawGPS{1}{GPZDAix(jj)};
            end
            time(jj,:) = str2num(tmp(8:13));
            % Correct time to seconds of day
            timestr = [tmp(8:9),':',tmp(10:11),':',tmp(12:13)];
            [~,~,~,H,M,S] = datevec(timestr);
            sec(jj,:) = H*3600+M*60+S;
            date(jj,:) = datestr(datenum(tmp([18,19,21,22,24,25,26,27]),'ddmmyyyy'),'mm/dd/yyyy');
                        
            % Get the Trace Number and Position
            tmp = strsplit(rawGPS{1}{hdix(jj)});
            trcno(jj,:) = str2num(tmp{2}(2:end));
            posx(jj,:) = str2num(tmp{5});
            
            % Include Post Processed GPS
            if strcmp(isPostProcessed,'Yes')
                ppIx = find(sec(jj,:) == ppSec(:));
                if ~isempty(ppIx)
                    lon(jj,:) = ppGPS.Var1(ppIx);
                    lat(jj,:) = ppGPS.Var2(ppIx);
                    z(jj,:) = ppGPS.Var6(ppIx);
                else
                    % Get the Longitude DDDMM.MMMMM
                    tmp = rawGPS{1}{GPGGAix(jj)};
                    try
                        lon(jj,:) = str2num(tmp(32:34));
                        lonm(jj,:) = str2num(tmp(35:43));
                        % Convert to Decimal Degees
                        lonm(jj,:) = lonm(jj)./60;
                        lon(jj,:) = lon(jj)+lonm(jj);
                        % Check Sign
                        if (tmp(45)) == 'W'
                            lon(jj,:) = -lon(jj,:);
                        end
                    catch
                        lon(jj,:) = NaN;
                    end
                    
                    % Get the Latitude DDMM.MMMMM
                    try
                        lat(jj,:) = str2num(tmp(18:19));
                        latm(jj,:) = str2num(tmp(20:28));
                        % Convert to Decimal Degees
                        latm(jj,:) = latm(jj)./60;
                        lat(jj,:) = lat(jj)+latm(jj);
                        % Check Sign
                        if (tmp(30)) == 'S'
                            lat(jj,:) = -lat(jj,:);
                        end
                    catch
                        lat(jj,:) = NaN;
                    end
                    % Post Processed Elevation Data is unavailable
%                     z(jj,:) = -999999;%NaN;
                    % Get Raw Elevation from GPGGA dummy!
                    % Get the Elevation cm precision
                    try
                        z(jj,:) = round(str2num(tmp(56:62)),3,'decimal');
                    catch
                        z(jj,:) = NaN;
                    end
                end
                   
            else
            if isGeode
                % Get the Longitude DDDMM.MMMMMMM
            tmp = rawGPS{1}{GPGGAix(jj)};
            try
            lon(jj,:) = str2num(tmp(33:35));
            lonm(jj,:) = str2num(tmp(36:45));
            % Convert to Decimal Degees
            lonm(jj,:) = lonm(jj)./60;
            lon(jj,:) = lon(jj)+lonm(jj);
            % Check Sign
            if (tmp(47)) == 'W'
                lon(jj,:) = -lon(jj,:);
            end
            catch
                lon(jj,:) = NaN;
            end
                
            else
            % Get the Longitude DDDMM.MMMMM
            tmp = rawGPS{1}{GPGGAix(jj)};
            try
            lon(jj,:) = str2num(tmp(32:34));
            lonm(jj,:) = str2num(tmp(35:43));
            % Convert to Decimal Degees
            lonm(jj,:) = lonm(jj)./60;
            lon(jj,:) = lon(jj)+lonm(jj);
            % Check Sign
            if (tmp(45)) == 'W'
                lon(jj,:) = -lon(jj,:);
            end
            catch
                lon(jj,:) = NaN;
            end
            end
            
            % Get the Latitude DDMM.MMMMM
            if isGeode
            try
            lat(jj,:) = str2num(tmp(18:19));
            latm(jj,:) = str2num(tmp(20:29));
            % Convert to Decimal Degees
            latm(jj,:) = latm(jj)./60;
            lat(jj,:) = lat(jj)+latm(jj);
            % Check Sign
            if (tmp(31)) == 'S'
                lat(jj,:) = -lat(jj,:);
            end
            catch
                lat(jj,:) = NaN;
            end
            else
            % Less Precision
            % Get the Latitude DDMM.MMMMM
            try
            lat(jj,:) = str2num(tmp(18:19));
            latm(jj,:) = str2num(tmp(20:28));
            % Convert to Decimal Degees
            latm(jj,:) = latm(jj)./60;
            lat(jj,:) = lat(jj)+latm(jj);
            % Check Sign
            if (tmp(30)) == 'S'
                lat(jj,:) = -lat(jj,:);
            end
            catch
                lat(jj,:) = NaN;
            end
            end
            if isGeode
            % Get the Elevation mm precision
            try
            z(jj,:) = round(str2num(tmp(58:66)),3,'decimal');
            catch
                z(jj,:) = NaN;
            end
            else % Get Trimble Raw Elevation
                try % Change these indicies
                    z(jj,:) = round(str2num(tmp(56:62)),3,'decimal');
                catch
                    z(jj,:) = NaN;
                end
            end
            end
        end
        % Remove Non Unique Values
        [time,unIx] = unique(time);
        sec = sec(unIx);
        lon = lon(unIx);
        lat = lat(unIx);
        trcno = trcno(unIx);
%         if strcmp(isPostProcessed,'Yes')
        z = z(unIx);
%         end
        % Remove NaN Values
        nanIx = find(isnan(lon));
        lon(nanIx) = [];
        lat(nanIx) = [];
        trcno(nanIx) = [];
        time(nanIx) = [];
        sec(nanIx) = [];
%         if strcmp(isPostProcessed,'Yes')
            z(nanIx) = [];
%         end
        
        % Least Squares Estimate of Traces per Second
        G = [ones(length(time),1),trcno]; d = sec;
        m = G\d; % 1./m(2) = Trace Rate;
        
        % Reconfigure GPS Trace Number
        dsec = sec - sec(1);
        trcno = round(dsec./m(2))+1;
        
        % Ensure that we haven't Created Extraneous Trace Numbers
        % Due to Free Run Acquisition this Traces will be Removed Later
            % This does not Introduce any Permanent or Correlated Errors
        trcendIx = find(trcno > multiplexNtrcs);
        if ~isempty(trcendIx)
            trcno(trcendIx(end):-1:trcendIx(1)-1) = ...
                multiplexNtrcs:-1:(multiplexNtrcs-(length(trcendIx)));
            % Algorithm may fail (Recurssion Catch)
            if trcno(trcendIx(1)- 1)  ==  trcno(trcendIx(1)-2)
                dtrc = trcno(trcendIx(1)- 2) - trcno(trcendIx(1)- 3) ;
                if dtrc > 1
                    trcno(trcendIx(1)- 2) = trcno(trcendIx(1)- 2) - 1;
                else
                    disp('Algortihm is not robust')
                    keyboard
                end
            end
            
        end
        
        % Linear Time Sample Interpolation
        s = sec; % Copy
        sec = interp1(trcno,sec,1:multiplexNtrcs);
        secNanIx = find(isnan(sec));
        secGoodIx = find(~isnan(sec));
        dsec = mean(diff(sec(secGoodIx)));
        % Recursively Repair NaN Values caused by Extrapolation
        for jj = 1:length(secNanIx)
            sec(secNanIx(jj)) = sec(secNanIx(jj) - 1) + dsec;
        end
        decsec = sec./(24*3600);
        
        % Convert to Time of Day
        tmp = datestr(decsec,'HH:MM:SS.FFF');
        tmp(:,[3,6]) = [];
        time = str2num(tmp);
   
        % Convert to UTM
        [E,N,utmzone] = deg2utm(lat,lon);
        % Correction for Antenna Position
        
        % PCHIP Interpolation (OutPerforms MLS)
        X = pchip(s,E,sec);
        Y = pchip(s,N,sec);
%         if strcmp(isPostProcessed,'Yes')
            % Include Elevations
            antennaHeight = 1.5; %[m]
        Z = pchip(s,z,sec)-antennaHeight;
%         end
        
        % Dead Reckoning occurs in DeMux
        
        
        % Append GPS Data to Trace Header
        trhd(4,:) = time;
        % Antenna Midpoint
        trhd(13,:) = X;
        trhd(14,:) = Y;
%         if strcmp(isPostProcessed,'Yes')
            trhd(15,:) = Z;
%         end
        
        % Install Transmitter and Receiver Sequencing and Geometry
        % 500 MHz Offsets
        if f0 == 500
            error('Undefined Offset Array!')
            % 1 GHz Offsets
        elseif f0 == 1000
            % QuadPol Array
            txGeo = [-0.2,-.2,-.6, -.6 ]; % Tx Sequence & Absolute Position
            rxGeo = [0, -.4,0,-.4]; % Rx Sequence & Absolute Position
%             % SkiPole Array
%             txGeo = [0 0]; % Tx Sequence & Absolute Position
%             rxGeo = [0.3, -0.3]; % Rx Sequence & Absolute Position
            % Grand Mesa Sled Array
%             txGeo = [0 0]; % Tx Sequence & Absolute Position
%             rxGeo = [0.25, -0.25]; % Rx Sequence & Absolute Position
        else
            error('Undefined Offset Array!')
        end
        
        offsetArray = abs(rxGeo - txGeo);     % Offset Array for CMP Gather
        % Append the Offset Array for NetCDF Export
        channelArray = trhd(23,:);
        offsetAppend = offsetArray(channelArray);
        trhd(22,:) = offsetAppend;
        
        %% Recycled Code Block from Preprocess.m
            X = X(:); Y = Y(:); Z = Z(:);
            GPS(:,[1:3]) = [X,Y,Z];
            % Filter Smoothing Distance
            R = 51;
            % Interpolated Sites
            dX = [0;diff(X)];
            dX = medfilt1(dX,R);
            dY = [0;diff(Y)];
            dY = medfilt1(dY,R);
            dZ = diff(Z);
            dZ = [dZ(1);dZ];
            dZ = medfilt1(dZ,R);
            % Raw Sites
            dx = [0;diff(E)];
            dy = [0;diff(N)];
            dz = diff(z);
            dz = [dz(1);dz];
            
            % Compute Distance
            if quantile(Z,0.5 < 0)
                % If elevation is unavailable
                dS = sqrt(dX.^2+dY.^2);
                ds = sqrt(dx.^2+dy.^2);
            else
                dS = sqrt(dX.^2+dY.^2+dZ.^2);
                ds = sqrt(dx.^2+dy.^2+dz.^2);
            end
            if ii>2
                try
                    endDist;
                catch
                    endDist = 0;
                end
                Distance = endDist+cumsum(dS);
                distance = endDist+cumsum(ds);
            else
                Distance = cumsum(dS); %[m]
                distance = cumsum(ds);
            end
            % Smooth Delta Distance
            dS = medfilt1(dS,R);
            endDist = Distance(end);
            % DistanceKm = Distance./1000; % [km]
            trhd(2,:) = Distance; % Configure Distance [m]
            % Compute Velocity
            tmpTime = trhd(4,:);
            tmpTime = num2str(tmpTime');
            Time = zeros(length(tmpTime),1);
            for kk = 1:length(tmpTime)
            spaceIx = find(isspace(tmpTime(kk,:)));
            spaceTime = tmpTime(kk,:);
            spaceTime(spaceIx) = [];
            Time(kk) = str2num(strcat(spaceTime([1,2]))).*3600+...
                str2num(strcat(spaceTime([3,4]))).*60+...
                str2num(strcat(spaceTime(5:end)));
            end
            dT = [mean(diff(Time));diff(Time)];
            dt = [mean(diff(s));diff(s)];
            Velocity = dS./dT; % [m/s]
            velocity = ds./dt;
%             Velocity = medfilt1(Velocity,R);
            Velocity = movmean(Velocity,R);
            velocityThreshold = 0.25;
            dumIx = find(velocity<velocityThreshold);

            
            % Heading
            dH = [0,0;diff([X,Y])];
            dh = [0,0;diff([E,N])];
            Heading = mod(atan2(dH(:,1), dH(:,2))*180/pi, 360);
            heading = mod(atan2(dh(:,1), dh(:,2))*180/pi, 360);
            dheading = [0;diff(heading)];
           %  heading degree threshold
           headThresh = 300;
            wrapIx = find(abs(dheading)>headThresh);
            rmvIx = find(ismember(wrapIx,dumIx));
            % Remover Static Positions
            wrapIx(rmvIx) = [];
            wrapIx = wrapIx - 1;
            for kk = 1:length(wrapIx)
            [~,interpIx(kk)] = min(abs(Time-s(wrapIx(kk))));
            Heading(interpIx(kk)+1:end) = Heading(interpIx(kk)+1:end)-360;
            end
            
%             Heading = medfilt1(Heading,R);
%             Heading = movmean(Heading,R);
%             Tailing = mod(Heading+180,360);
            
            % Slope
            Slope = atand(dZ./dS);
%             Slope = medfilt1(Slope,R);
            Slope = movmean(Slope,R);

            % Append GPS to Trace Header
            trhd(17:19,:) = [Slope';Velocity';Heading'];
            
%             GPS = [GPS,Distance,Slope,Velocity,Heading];
            clear('GPSix','GPSixEdges','dS','dX','dY','dZ''dT','ds','dx','dy','dz''dt',...
                'Distance','Slope','Velocity','Heading','Tailing','Time','tmpTime','GPS',...
                'interpIx','heading','dheading','wrapIx','distance');
 
        % Nominal Frequency GHz
%         D.f0GHz = D.f0/1000;
        % No. Traces in Multiplexed Data
%         [~, multiplexNtrcs] = size(trhd);
        
%         % Load Transmitter and Receiver Sequencing and Geometry        
%         % 1 GHz Offsets
%         % Subrtact 5m Distance from GPS to Rx1;
%         txGeo = [-.2, -.2, -.6, -.6] - 5; % Tx Sequence & Absolute Position
%         rxGeo = [0, -.4, 0, -.4] - 5; % Rx Sequence & Absolute Position
%         
%         % Input Polaraztion Configureation
%         % QuadPol
%         D.Polarization = {'HV','HH','VV','VH'};
%         % SkiPol
%         D.Polarization = {'HH','HV'};
% 
%         
%         D.offsetArray = D.trhd{ii}(22,chan);      % Offset Array for CMP Gather
        midPointArray = mean([txGeo;rxGeo],1);  % Midpoint Locations Relative to Tx 1
        trhd(21,:) = midPointArray([trhd(23,:)]);
        trhd(16,:) = midPointArray([trhd(23,:)])+trhd(2,:); % Append Midpoint Locations
        nChan = length(offsetArray);           % Refresh nChan before kill
        polOffset = min(offsetArray);          % Determine Near Offset
        pol = 1;
        
        % Determine Data Acquisition  Method
        if all(trhd(5,:) == 0)
            isFreeRun = 1;
            isOdometer = 0;
        else
            isOdometer = 1;
            isFreeRun = 0;        
        end
        
        % Remove Static Traces if is Free Run Acquisition
        if isFreeRun
            % Process Near Offset Channel
            disp(' ')
            fprintf('Begin Static Trace Removal \n')
%             tic
%             [polRad, polL2] = processTraceRemoval... % Processes Near Offset Data
%                 (D.Rad{ii}(:,pol:nChan:end), D.f0, D.dt );
            % Remove Duplicate Static Traces
%             threshold = quantile(D.trhd{ii}(18,:),.75);
%             dupIx = find(D.trhd{ii}(18,:) < 1);
%             dupIx = removeStaticTrace( polRad, polL2, D.trhd{ii}, multiplexNtrcs, pol, nChan );
%             [GPSix,edgeIx] = removeStaticGPS(D.trhd{1}(13,:),D.trhd{1}(14,:),500);
%             clear('nearRad');
            
%             figure();
%             plot(D.trhd{ii}(16,:),1:multiplexNtrcs,'.k');hold on
%             plot(D.trhd{ii}(16,dupIx),dupIx,'.r');
%             figure();plot(D.trhd{ii}(16,:),D.trhd{ii}(18,:),'.k')
            % Remove Static Trace Headers from Multiplexed Record
            dupIx = removeStaticPositions(trhd,velocityThreshold);
%             dupIx = [];
            % Manually Pluck a few
            if (strcmp(dataDir,'D:\CRREL_SnowCompaction\CRREL\SnowEx2019\GrandMesaGPR\PulseEKKO_26March2019\PulseEKKO\')...
                    && strcmp(filename,'LINE01'))
                allIx = [1:multiplexNtrcs]';
                pluckIx = allIx(~ismember(allIx,dupIx));
                dupIx = [dupIx;pluckIx(1:165)];
            end
            if (strcmp(dataDir,'D:\CRREL_SnowCompaction\CRREL\SnowEx2019\GrandMesaGPR\PulseEKKO_27March2019\PulseEKKO\')...
                    && strcmp(filename,'LINE04'))
                allIx = [1:multiplexNtrcs]';
                pluckIx = allIx(~ismember(allIx,dupIx));
                dupIx = [dupIx;pluckIx(4751:end)];
                %                 trhd(:,4751:end) = [];
            end
            trhd(:,dupIx) = []; 
            % Configure Trace Indicies
            trhd(1,:) = 1:length(trhd);
            % re-Incorporate GPS here
            X = trhd(13,:)';Y = trhd(14,:)';Z = trhd(15,:)';
            
            dX = [0;diff(X)];
            dX = medfilt1(dX,R);
            dY = [0;diff(Y)];
            dY = medfilt1(dY,R);
            dZ = diff(Z);
            dZ = [dZ(1);dZ];
            dZ = medfilt1(dZ,R);
            % Compute Distance
            if quantile(Z,0.5 < 0)
                % If elevation is unavailable
                dS = sqrt(dX.^2+dY.^2);
            else
                dS = sqrt(dX.^2+dY.^2+dZ.^2);
            end
            dS = medfilt1(dS,R);
            if ii>1
                Distance = endDist+cumsum(dS);
            else
                Distance = cumsum(dS); %[m]
            end
            endDist = Distance(end);
            
            trhd(2,:) = Distance; % Configure Distance
            Rad(:,dupIx) = [];      % Remove Static Traces from Multiplexed Data
            xArray = trhd(2,:);         % Define Configured Distance xArray
            
            fprintf('Static Trace Removal Done \n')
%             toc
            display(' ')
            
        end
        
        % Interpolate GPS
        % GPS = [X,Y,Z,Distance,Slope,Speed,Heading,Tailing];
%         nGPS = size(GPS,1);
%         nTrcs = size(D.Rad{ii},2);
%         xq = linspace(1,nGPS,nTrcs);
%         % Piecewise Cubic Hermite Interpolating Polynomials
%         Si = pchip(1:nGPS,GPS(:,4),xq);
%         Xi = pchip(GPS(:,4),GPS(:,1),Si);
%         Yi = pchip(GPS(:,4),GPS(:,2),Si);
%         Zi = pchip(GPS(:,4),GPS(:,3),Si);
%         Sxi= pchip(GPS(:,4),GPS(:,5),Si);
%         Vi = pchip(GPS(:,4),GPS(:,6),Si);
%         Hi = pchip(GPS(:,4),GPS(:,7),Si);
%         Ti = pchip(GPS(:,4),GPS(:,8),Si);
%         % Append GPS to Trace Header
%         D.trhd{ii}(13:20,:) = [Xi;Yi;Zi;Si;Sxi;Vi;Hi;Ti];
        clear('GPS');%,'nGPS','nTrcs','Si','Xi','Yi','Zi','Sxi','Vi','Hi','Ti','xq');
        
        % Remove Skip Traces if Wheel Odometer Acquisition
        if isOdometer
            disp('Somethings Wrong,Because Data is Free Run')
        dupIx = find(~(diff(trhd(23,:)))); % Find Skipped Traces
        
        for jj = dupIx
            trhd(1,jj:end) = trhd(1,jj:end) - 1; % ReConfigure Trace Indicies
        end
        
        trhd(:,dupIx) = []; % Remove Skipped Traces from Trace Header
        Rad(:,dupIx) = []; % Remove Skipped Traces from Multiplexed Data
        xArray = trhd(2,:); % Define ReConfigured Distance as xArray
        end
                
%         clear('dupIx');
        
        % Allocation Here
%         if ii == 1
%             Radar = cell(nChan,MD.nFiles); traceIx = cell(nChan,MD.nFiles);
% %             Array = cell(nChan,MD.nFiles);
%         end

            % GPS DeadReckoning Within Demux
                % Calculate Truer Antenna Positions using Array Geometry
                % delta is the correction from Tx1 antenna center
                delta = [.35,5.95,0];
                tmptrhd = trhd;
                for kk = 1:nChan
                    [trhd] = deadReckon(trhd,kk,delta);
                end
                % Smooth Positions post-Reckoning
                R = 101;
                for kk = 13:15
                    trhd(kk,:) = movmean(trhd(kk,:),R);
                end
        
        multiplexNtrcs = size(trhd,2);% Length of Multiplex
        traceMod = mod(multiplexNtrcs,nChan);% Count Unbinned Traces
        trhd(:,(multiplexNtrcs - traceMod)+1 :multiplexNtrcs ) = []; % Remove Unbinned Trace Headers
        Rad(:,(multiplexNtrcs - traceMod)+1 :multiplexNtrcs ) = []; % Remove Unbinned Shots
        
        % Average Positions for Bin Centers (The GPS Location and Distance)
        for jj = 1:nChan:size(trhd,2)
            % Store Bin Centers [m]
            trhd(10:12,jj:jj+nChan-1) = mean(trhd(13:15,jj:jj+nChan-1),2)*ones(1,nChan);
            % Overwrite Distance with Bin Center Position [m]
            tmp = mean(trhd(16,jj:jj+nChan-1),2)*ones(1,nChan);
            trhd(2,jj:jj+nChan-1) = tmp;
            % Overwrite Tailing with Average Bin Center Heading
            trhd(20,jj:jj+nChan-1) = mean(trhd(19,jj:jj+nChan-1),2)*ones(1,nChan);
        end
        % Smooth Heading
%             D.trhd{ii}(20,:) = nonParametricSmooth(1:length(D.trhd{ii}),D.trhd{ii}(20,:),1:length(D.trhd{ii}),251.*nChan);
        % Zero Starting Distance
        tmp = trhd(2,1);
        trhd(2,:) = trhd(2,:) - tmp;
        trhd(16,:) = trhd(16,:)-tmp;
        
        clear('tmp')
        
        % Remove Every Nth Trace for Data Reduction
%         rmNtrc = 2;
        
%         parfor (jj =  1:nChan, nWorkers)
% %         for jj = 1:nChan
% 
%             % Extract Full-fold Traces & Sort Antenna Positions
%             gatherLength = size(Radar{jj,ii},2); % Length of Each Channel
%             
%             % Flag Un-binned Traces
%             if jj <= traceMod
%                 xTrc = 1;
%             else
%                 xTrc = 0;
%             end
%                 
%                 % Remove un-Binned Common Offset Traces
%                 Radar{jj,ii} = Radar{jj,ii}(:,1:gatherLength - xTrc);
%                 tmpDistance{jj,ii} = D.Distance{jj,ii}(:,1:gatherLength - xTrc);
%                 
%                 % Remove un-Binned Trace Indicies
%                 traceIx{jj,ii} = traceIx{jj,ii}(:,1:gatherLength - xTrc);
%                 
%                 % Reduce Data Volume
%                 if isReduceData
%                     Radar{jj,ii} = Radar{jj,ii}(:,1:rmNtrc:end);
%                     traceIx{jj,ii} = traceIx{jj,ii}(1:rmNtrc:end);
%                 end
                %%
        % Reorganize Trace Header
        % GPS time
        tmp4 = trhd(4,:);
        % frequency MHz
        tmp9 = trhd(9,:);
        % Move GPS time to Geolocation block
        trhd(4,:) = tmp9;
        % Move Frequency
        trhd(9,:) = tmp4;
        % Eliminate Irrelavant array
        trhd(24,:) = zeros(1,length(trhd(24,:)));
        % channel number
        tmp23 = trhd(23,:);
        % midpoint geometery
        tmp21 = trhd(21,:);
        % offset
        tmp22 = trhd(22,:);
        % geolocation
        tmpgeo = trhd(9:20,:);
        % shot gather bin center distance
        tmp2 = trhd(2,:);
        % Antenna Midpoint Distance
        tmp16 = trhd(16,:);
        % Move Channel Number
        trhd(2,:) = tmp23;
        % nsamps
        tmp3 = trhd(3,:);
        % Move Frequency
        trhd(3,:) = tmp9;
        % dt
        tmp7 = trhd(7,:);
        % Move dt
        trhd(4,:) = tmp7./1000;
        % Move nsamps
        trhd(5,:) = tmp3;
        % Move midpoint Geometery
        trhd(6,:) = tmp21;
        % Move Offset array
        trhd(7,:) = tmp22;
        % Move Bin Center Geolocation
        trhd(8:11,:) = trhd(9:12,:);
        % Move Bin Center Distance
        trhd(12,:) = tmp2;
        % Move Midpoint Geoloctaion
        trhd(21:23,:) = trhd(13:15,:); 
        % Move Antenna Midpoint Distance
        trhd(24,:) = tmp16;
        % Move Bin Center Heading
%         trhd(13,:) = trhd(20,:);
        % Move Antenna Center Heading
        trhd(25,:) = trhd(19,:);
        % Move Velocity
        trhd(13,:) = trhd(18,:);
        % Move Slope
        trhd(15,:) = trhd(17,:);
        % Move Heading
        trhd(14,:) = trhd(20,:);
        % Move Midpoint Geolocations
        trhd(16:20,:) = trhd(21:25,:);
        % Remove Extra Rows
        trhd(21:25,:) = []; 
        
        % Attach trace header to data for NetCDF Export
        DATA = [trhd;Rad];
        
        clear('time','s','sec','date','trcno','posx','lat','latm','lon','lonm','z')
        
        %         GPR = struct('metaData',{hdr1},'traceHeader',trhd,'traces',Rad,...
        %             'offsets',offsetArray,'frequency',f0,'timeSample',dt,'spaceSample',dx);
%         toc
        %% Write to NetCDF
        if isWrite
            % Create the output Directory if it does not exist
%             if exist([workingDirectory,'/',folders(ff).name],'dir') ~= 7
%                 mkdir(workingDirectory,folders(ff).name)
%             end
            % Change to outdir for File Write
%             cd([workingDirectory,'/',folders(ff).name])
               cd(dataDir)
            
            % Name the NetCDF.nc file
            lfile = length(filename);
            fileDate = ['03',folders(ff).name(11:12),folders(ff).name(18:21)];
            ncdfName = ['BSU_pE_GPR_',fileDate,'_',filename(lfile-1:lfile),'.nc'];
%             ncdfName = ['BSU_pE_GPR_',folders(ff).name,'_',filename(lfile-1:lfile),'.nc'];
%             ncdfName = ['pulseEKKO-1GHz-GPR-HH-HV-multiplexed-radargram-',folders(ff).name,'-',filename(lfile-1:lfile),'.nc'];
%             ncdfName = [folders(ff).name,'-',filename,'.nc'];
%             wName = 'pulseEKKO-122019';
%             ncdfName = [wName,'-',filename,'.nc'];

            
            fprintf(['Writing ',ncdfName, '\n'])
            
            numrow = size(DATA,1);
            
            numcol = size(DATA,2);
            
            % Create the NetCDF-4 file with Dimensions of DATA
            netcdf.setDefaultFormat('FORMAT_NETCDF4') ; 
            
            ncid = netcdf.create(ncdfName,'CLOBBER');
            
            dimidrow = netcdf.defDim(ncid,'rows',numrow);
            
            dimidcol = netcdf.defDim(ncid,'length',numcol);
            
            % Define and put Variable DATA
            varid = netcdf.defVar(ncid,'DATA','NC_DOUBLE',[dimidrow dimidcol]);
            
            % Increase NetCDF Compression
            netcdf.defVarDeflate(ncid,varid,true,true,5);
            
            netcdf.endDef(ncid);
            
            netcdf.putVar(ncid,varid,DATA);
                        
            % Close the .nc file after writing
            netcdf.close(ncid);
            % Test for Proper Writing
            %     ncid2 = netcdf.open(ncdfName,'NC_NOWRITE');
            %
            %     data_copy = netcdf.getVar(ncid2,0);
            %
            %     if isequal(DATA,data_copy)
            %
            %         disp('Data match');
            %
            %     else
            %
            %         disp('Data mis-match');
            %     end
            
            % Return to the working directory for next iteration
            cd(workingDirectory)
            fprintf(['Done writing ',folders(ff).name,'-',files(ii,:),'\n'])
            reportMessage = ['Wrote ',folders(ff).name,'-',filename,'.nc','\n'];
            Report2User(reportMessage,'ReportLog.txt',workingDirectory)
        end
        % Extract GPS Positions
        if isDumpGPS
        if ii == 1 && ff == 1
        csvGPS = [trhd(9,1:nChan:end)',trhd(10,1:nChan:end)',trhd(11,1:nChan:end)'];
        end
        try csvGPS;
            csvGPS = [csvGPS;[trhd(9,1:nChan:end)',trhd(10,1:nChan:end)',trhd(11,1:nChan:end)']];
        catch
            csvGPS = [trhd(9,1:nChan:end)',trhd(10,1:nChan:end)',trhd(11,1:nChan:end)'];
        end
        end
        % Clear Memory
        clear('DATA','trhd','Rad')
        toc
        display(' ')
    end
    if strcmp(isPostProcessed, 'Yes')
        clear('ppGPS','ppTime','ppDate','ppSec')
    end
end
%%
if isDumpGPS
% Write GPS to .csv
cd(directory)
writematrix(csvGPS,'BSU_pE_GPR_UTM.csv')
writematrix(csvGPS(1:100:end,:),'BSU_pE_GPR_UTM_decimated.csv')
writematrix([utm2ll(csvGPS(:,1),csvGPS(:,2),12),csvGPS(:,3)],'BSU_pE_GPR_LL.csv');
writematrix([utm2ll(csvGPS(1:100:end,1),csvGPS(1:100:end,2),12),csvGPS(1:100:end,3)],'BSU_pE_GPR_LL_decimated.csv');
cd(workingDirectory)
end
