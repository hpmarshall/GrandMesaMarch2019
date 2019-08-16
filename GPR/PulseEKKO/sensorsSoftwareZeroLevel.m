%% sensorsSoftwareZeroLevel.m
% This code unpacks the binary S&S data file and writes the raw traces and
% meta data to ncdf

% Tate Meehan - Oct. 2018 - Updated April 2019
clear; close all; clc
% addpath '/sonichome/tatemeehan/GreenTracs2017/GPR_Processing/MultiOffset/TM'
addpath './functions'
isWrite = 1;

% Create Error Log
 f = fopen( 'ErrorLog.txt', 'w' );  
 fclose(f);
% Create Report Log
 f = fopen('ReportLog.txt','w');
%% From Data Directory Read and Sort SubDirectories
workingDirectory = pwd;
% Enter the Appropriate Data Directory
% directory = '/SNOWDATA/GrandMesa2019/GPR/';
directory = 'D:\GrandMesaGPR';
folders = dir(directory);
% folders(1:2) = []; % Remove Hidden Directories
folders([1:4,6:end]) = []; % Remove Extra Folders
datearray = cat(1,folders.date);
daymonthyear = (datearray(:,1:11));
num = datenum(daymonthyear,'dd-mmm-yyyy');
[~,sortIx] = sort(num);
folders = folders(sortIx);

% Loop Over each day of Acquisition
for ff = 1:length(folders)
%     dataDir = [directory,'/',folders(ff).name,'/PulseEKKO/'];
    dataDir = [directory,'\',folders(ff).name,'\'];%,'\PulseEKKO\'];
        if strcmp(dataDir,'D:\GrandMesaGPR\PulseEKKO_28March2019\PulseEKKO\')
%     if strcmp(dataDir,'/SNOWDATA/GrandMesa2019/GPR/PulseEKKO_28March2019/PulseEKKO/')
        % Move into one additional directory for this day.
%         dataDir = '/SNOWDATA/GrandMesa2019/GPR/PulseEKKO_28March2019/PulseEKKO/D04/';
        dataDir = 'D:\GrandMesaGPR\PulseEKKO_28March2019\PulseEKKO\D04\';
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
       isPostProcessed = questdlg('Load a Post-Processed GPS File?','Post-Processed GPS','Yes','No','Yes');
       if strcmp(isPostProcessed,'Yes')
           disp('Loading Post-Processed GPS')
           tic
           if strcmp(dataDir,'D:\GrandMesaGPR\PulseEKKO_28March2019\PulseEKKO\D04\')
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
        [Rad,hdr1,trhd,dt,f0,~,dx,~] = readSensorsSoftwareData( filepath );
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
%         fclose(fid);
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
                tmp = rawGPS{1}{GPZDAix(jj)};;
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
                    z(jj,:) = -999999;%NaN;
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
        end
        % Remove Non Unique Values
        [time,unIx] = unique(time);
        sec = sec(unIx);
        lon = lon(unIx);
        lat = lat(unIx);
        trcno = trcno(unIx);
        if strcmp(isPostProcessed,'Yes')
        z = z(unIx);
        end
        % Remove NaN Values
        nanIx = find(isnan(lon));
        lon(nanIx) = [];
        lat(nanIx) = [];
        trcno(nanIx) = [];
        time(nanIx) = [];
        sec(nanIx) = [];
        if strcmp(isPostProcessed,'Yes')
            z(nanIx) = [];
        end
        
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
        
        % PCHIP Interpolation (OutPerforms MLS)
        X = pchip(s,E,sec);
        Y = pchip(s,N,sec);
        if strcmp(isPostProcessed,'Yes')
            % Include Elevations
        Z = pchip(s,z,sec);
        end
        
        % Append GPS Data to Trace Header
        trhd(4,:) = time;
        trhd(13,:) = X;
        trhd(14,:) = Y;
        if strcmp(isPostProcessed,'Yes')
            trhd(15,:) = Z;
        end
        
        % Install Transmitter and Receiver Sequencing and Geometry
        % 500 MHz Offsets
        if f0 == 500
            display('Undefined Offset Array!')
            % 1 GHz Offsets
        elseif f0 == 1000
            txGeo = [-0.2,-.2,-.6, -.6 ]; % Tx Sequence & Absolute Position
            rxGeo = [0, -.4,0,-.4]; % Rx Sequence & Absolute Position
        else
            display('Undefined Offset Array!')
        end
        
        offsetArray = abs(rxGeo - txGeo);     % Offset Array for CMP Gather
        % Append the Offset Array for NetCDF Export
        channelArray = trhd(23,:);
        offsetAppend = offsetArray(channelArray);
        trhd(22,:) = offsetAppend;
        
        % Attach trace header to data for NetCDF Export
        DATA = [trhd;Rad];
        
        clear('time','s','sec','date','trcno','posx','lat','latm','lon','lonm','z')
        
        %         GPR = struct('metaData',{hdr1},'traceHeader',trhd,'traces',Rad,...
        %             'offsets',offsetArray,'frequency',f0,'timeSample',dt,'spaceSample',dx);
        
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
            ncdfName = [folders(ff).name,'-',filename,'.nc'];
            
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
        toc
        display(' ')
    end
end