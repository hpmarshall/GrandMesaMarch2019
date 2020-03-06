%% Read Sensors and Software Data
% This Script Reads the Binary Multi-channel GPR Data and Performs the
% pre-processing and signal processing routines. 
% GPS Locations are incorporated in NSIDC Sea Ice Polarstereographic North
%
% The Array Geometry is Installed
% The Acquisition Method is Decided
% The Near Channel (7) is Killed
% The Time Window can be trimmed
% Data is de-multiplex (grouped into channels) and coarse trace shifts
% are applied as the first cut at correcting for digital time-sampling 
% errors.
% Inside the function wrapper processCommonOffset.m de-WOW filtering and
% trace stacking are applied. There are Various filter parameters decided
% in this function.

    D.Rad = cell(1,MD.nFiles);
    D.trhd = cell(1,MD.nFiles);
%     GPS = cell(1,MD.nFiles);
%     Year = cell(1,MD.nFiles);
    TimeAxis = cell(1,MD.nFiles);
    tmpDistance = cell(nChan,MD.nFiles);

    
    for ii = 1 : MD.nFiles
        tic
        %------------------------------------------------------------------
        % Multiplexed Channel Record
        filename = MD.fileNames(MD.lineNo(ii)+1).name;
        filepath = fullfile(MD.dataDir,filename);
        % Read netCDF data file
        disp(' ')
        fprintf('Reading .nc File \n')
        tic
        ncRad = ncread(filepath,'DATA');
        D.trhd{ii} = ncRad(1:25,:);
        D.Rad{ii} = ncRad(26:end,:);
        clear('ncRad');
        D.f0 = (D.trhd{ii}(9,1)); % [MHz]
        D.dt = D.trhd{ii}(7,1)/1000; % [ns]
%         D.dx = diff(D.trhd{ii}(2,1:2)); % [m]
        % Need to Automate Offset Array from .nc File
%         offsetArray = unique(D.trhd{ii}(22,1:100));

        % Configure GPS
            disp(' ')
            fprintf('Configuring GPS \n')
            tic
            % GPS = [X,Y,Z,Distance,Slope,Speed,Heading,Tailing];
            X = D.trhd{ii}(13,:);
            Y = D.trhd{ii}(14,:);
            Z = D.trhd{ii}(15,:);
            X = X(:); Y = Y(:); Z = Z(:);
            GPS(:,[1:3]) = [X,Y,Z];
            
            dX = [0;diff(X)];
            dY = [0;diff(Y)];
            dZ = diff(Z);
            dZ = [dZ(1);dZ];
            
            % Compute Distance
            if quantile(Z,0.5 < 0)
                % If elevation is unavailable
                dS = sqrt(dX.^2+dY.^2);
            else
                dS = sqrt(dX.^2+dY.^2+dZ.^2);
            end
            if ii>1
                Distance = endDist+cumsum(dS);
            else
                Distance = cumsum(dS); %[m]
            end
            endDist = Distance(end);
            % DistanceKm = Distance./1000; % [km]
            D.trhd{ii}(2,:) = Distance; % Configure Distance [m]
            % Filter Smoothing Distance
            R = 11;
            % Compute Velocity
            tmpTime = D.trhd{ii}(4,:);
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
            Velocity = dS./dT; % [m/s]
            Velocity = medfilt1(Velocity,R);
            
            % Heading
            dH = [0,0;diff([X,Y])];
            Heading = mod(atan2(dH(:,1), dH(:,2))*180/pi, 360);
            Heading = medfilt1(Heading,R);
%             Tailing = mod(Heading+180,360);
            
            % Slope
            Slope = atand(dZ./dS);
            Slope = medfilt1(Slope,R);

            % Append GPS to Trace Header
            D.trhd{ii}(17:19,:) = [Slope';Velocity';Heading'];
            
%             GPS = [GPS,Distance,Slope,Velocity,Heading];
            clear('GPSix','GPSixEdges','dS','dX','dY','dZ''dT',...
                'Distance','Slope','Velocity','Heading','Tailing','Time','tmpTime','R','GPS');
 
        % Nominal Frequency GHz
        D.f0GHz = D.f0/1000;
        % No. Traces in Multiplexed Data
        [~, multiplexNtrcs] = size(D.Rad{ii});
        
        % Load Transmitter and Receiver Sequencing and Geometry        
        % 1 GHz Offsets
        % Subrtact 5m Distance from GPS to Rx1;
        txGeo = [-.2, -.2, -.6, -.6]; % Tx Sequence & Absolute Position
        rxGeo = [0, -.4, 0, -.4]; % Rx Sequence & Absolute Position
        
        % Input Polaraztion Configureation
        D.Polarization = {'HV','HH','VV','VH'};
        
        D.offsetArray = D.trhd{ii}(22,chan);      % Offset Array for CMP Gather
        midPointArray = mean([txGeo;rxGeo],1);  % Midpoint Locations Relative to Tx 1
        D.trhd{ii}(21,:) = midPointArray([D.trhd{ii}(23,:)]);
        D.trhd{ii}(16,:) = midPointArray([D.trhd{ii}(23,:)])+D.trhd{ii}(2,:); % Append Midpoint Locations
        nChan = length(D.offsetArray);           % Refresh nChan before kill
        polOffset = min(D.offsetArray);          % Determine Near Offset
        pol = 1;
        
        % Determine Data Acquisition  Method
        if all(D.trhd{ii}(5,:) == 0)
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
            tic
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
            vThreshold = .25;
            dupIx = removeStaticPositions(D.trhd{ii},vThreshold);
            D.trhd{ii}(:,dupIx) = [];         
            % Configure Trace Indicies
            D.trhd{ii}(1,:) = 1:length(D.trhd{ii});
            % re-Incorporate GPS here
            X = D.trhd{ii}(13,:)';Y = D.trhd{ii}(14,:)';Z = D.trhd{ii}(15,:)';
            
            dX = [0;diff(X)];
            dY = [0;diff(Y)];
            dZ = diff(Z);
            dZ = [dZ(1);dZ];
            
            % Compute Distance
            if quantile(Z,0.5 < 0)
                % If elevation is unavailable
                dS = sqrt(dX.^2+dY.^2);
            else
                dS = sqrt(dX.^2+dY.^2+dZ.^2);
            end
            if ii>1
                Distance = endDist+cumsum(dS);
            else
                Distance = cumsum(dS); %[m]
            end
            endDist = Distance(end);
            
            D.trhd{ii}(2,:) = Distance; % Configure Distance
            D.Rad{ii}(:,dupIx) = [];      % Remove Static Traces from Multiplexed Data
            xArray = D.trhd{ii}(2,:);         % Define Configured Distance xArray
            
            fprintf('Static Trace Removal Done \n')
            toc
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
        dupIx = find(~(diff(D.trhd{ii}(23,:)))); % Find Skipped Traces
        
        for jj = dupIx
            D.trhd{ii}(1,jj:end) = D.trhd{ii}(1,jj:end) - 1; % ReConfigure Trace Indicies
        end
        
        D.trhd{ii}(:,dupIx) = []; % Remove Skipped Traces from Trace Header
        D.Rad{ii}(:,dupIx) = []; % Remove Skipped Traces from Multiplexed Data
        xArray = D.trhd{ii}(2,:); % Define ReConfigured Distance as xArray
        end
                
%         clear('dupIx');

        % Pad Data with Instrument Zero
        padding = 0;
        instrumentPad = zeros(padding,size(D.Rad{ii},2));
        if padding ~= 0
            for jj = 1:size(D.Rad{ii},2)
                instrumentZero = D.Rad{ii}(1,jj);
                instrumentPad(:,jj) = ones(padding,1).*instrumentZero;
            end
            D.Rad{ii} = [instrumentPad;D.Rad{ii}];
        end  
        
        % Trim Time Window
        if isTrimTWT
            reSample = 575 + padding;   % Number of Wanted Samples
            D.Rad{ii} = D.Rad{ii}(1:reSample,:);
        end
        
        % Allocation Here
        if ii == 1
            Radar = cell(nChan,MD.nFiles); traceIx = cell(nChan,MD.nFiles);
%             Array = cell(nChan,MD.nFiles);
        end
        
        for jj = chan
            % DeMux Sequential Data
            [Radar{jj,ii},D.trhd{ii},traceIx{jj,ii},D.Distance{jj,ii}] = DeMux(D.Rad{ii},D.trhd{ii},chan(jj));
            % GPS DeadReckoning After Demux
            delta = [0.35,5.0,0]; % Reckoning Perturbations
            % Calculate Truer Antenna Positions using Array Geometry
            [D.trhd{ii}] = deadReckon( D.trhd{ii}, jj, delta );
        end
        % Smooth Positions post-Reckoning
        R = 51;
        for kk = 13:15
            D.trhd{ii}(kk,:) = movmean(D.trhd{ii}(kk,:),R);
        end
        
        multiplexNtrcs = size(D.trhd{ii},2);% Length of Multiplex
        traceMod = mod(multiplexNtrcs,nChan);% Count Unbinned Traces
        D.trhd{ii}(:,(multiplexNtrcs - traceMod)+1 :multiplexNtrcs ) = []; % Remove Unbinned Trace Headers
        
        % Average Positions for Bin Centers (The GPS Location and Distance)
        for jj = 1:nChan:size(D.trhd{ii},2)
            % Store Bin Centers [m]
            D.trhd{ii}(10:12,jj:jj+nChan-1) = mean(D.trhd{ii}(13:15,jj:jj+nChan-1),2)*ones(1,nChan);
            % Overwrite Distance with Bin Center Position [m]
            tmp = mean(D.trhd{ii}(16,jj:jj+nChan-1),2)*ones(1,nChan);
            D.trhd{ii}(2,jj:jj+nChan-1) = tmp;
            % Overwrite Tailing with Average Bin Center Heading
            D.trhd{ii}(20,jj:jj+nChan-1) = mean(D.trhd{ii}(19,jj:jj+nChan-1),2)*ones(1,nChan);
        end
        % Smooth Heading
%             D.trhd{ii}(20,:) = nonParametricSmooth(1:length(D.trhd{ii}),D.trhd{ii}(20,:),1:length(D.trhd{ii}),251.*nChan);
        % Zero Starting Distance
        tmp = D.trhd{1}(2,1);
        D.trhd{ii}(2,:) = D.trhd{ii}(2,:) - tmp;
        
        clear('tmp')
        
        % Remove Every Nth Trace for Data Reduction
        rmNtrc = 2;
        dt = D.dt;
        parfor (jj =  1:nChan, nWorkers)
%         for jj = 1:nChan

            % Extract Full-fold Traces & Sort Antenna Positions
            gatherLength = size(Radar{jj,ii},2); % Length of Each Channel
            
            % Flag Un-binned Traces
            if jj <= traceMod
                xTrc = 1;
            else
                xTrc = 0;
            end
                
                % Remove un-Binned Common Offset Traces
                Radar{jj,ii} = Radar{jj,ii}(:,1:gatherLength - xTrc);
                tmpDistance{jj,ii} = D.Distance{jj,ii}(:,1:gatherLength - xTrc);
                
                % Remove un-Binned Trace Indicies
                traceIx{jj,ii} = traceIx{jj,ii}(:,1:gatherLength - xTrc);
                
                % Reduce Data Volume
                if isReduceData
                    Radar{jj,ii} = Radar{jj,ii}(:,1:rmNtrc:end);
                    traceIx{jj,ii} = traceIx{jj,ii}(1:rmNtrc:end);
                end
                
                
                % Process Common Offset Channels
%                 disp(' ')           
%                 fprintf(['Begin Signal Processing in Common-Offset Domain ',...
%                     filename, ' CHAN0', num2str(jj),'\n'])
            % Store Travel-Time Axis
            TWT = [0:dt:(dt.*(size(Radar{jj,ii},1)-1))]';

            [Radar{jj,ii}, TimeAxis{jj,ii}] = processCommonOffset(Radar{jj,ii}, ...
                D.f0, D.dt, TWT, D.offsetArray(jj) );

        end
        if isReduceData
            tmpIx = sort(cat(2,traceIx{:,ii}));
            D.trhd{ii} = D.trhd{ii}(:,tmpIx);
        end
        
        % Trim all time windows to same length
        for jj = chan
            if jj == 1
            minIx = length(TimeAxis{jj,ii});
            minChan = jj;
            elseif minIx > length(TimeAxis{jj,ii})
                minIx = length(TimeAxis{jj,ii});
                minChan = jj;
            end
        end
        clear('TimeAxis')
        % Store Travel-Time Axis
            TimeAxis = [0:dt:(dt.*(minIx-1))]';
        % Trim Channels to Consistent Travel Time Axis
        for jj = chan
            Radar{jj,ii} = Radar{jj,ii}(1:minIx,:);            
        end
        D.Radar = Radar;
        D.DistanceAxis{ii} = mean(cat(1,tmpDistance{:,ii}));
        D.DistanceAxis{ii} = D.DistanceAxis{ii} - D.DistanceAxis{ii}(1);
        kk = 0;
        D.X{ii} = zeros(length(1:nChan:size(D.trhd{ii},2)),1);
        D.Y{ii} = D.X{ii}; D.Z{ii} = D.X{ii};
        for jj = 1:nChan:size(D.trhd{ii},2)
            kk = kk + 1;
            D.X{ii}(kk) = D.trhd{ii}(10,jj);
            D.Y{ii}(kk) = D.trhd{ii}(11,jj);
            D.Z{ii}(kk) = D.trhd{ii}(12,jj);            
        end
        D.TimeAxis{ii} = TimeAxis;
        clear('TimeAxis')
        
        fprintf('Signal Processing Done \n')
        toc
        display(' ')
    end
    MD.nChan = nChan; MD.chan = chan;
    D = rmfield(D,'Distance');
    clear('Rad','tmpIx','Distance','endDist','chan','nChan','TWT','dt','TimeAxis','minIx','minChan',...
        'dH','dS','dT','dupIx','dX','dY','dZ','gatherLength','instrumentPad',...
        'padding','midPointArray','multiplexNtrcs','pol','polOffset','Radar','rmNtrc',...
        'rxGeo','spaceIx','spaceTime','traceIx','traceMod','txGeo','X','xArray','xTrc','Y','Z','vThreshold');