function [ Rad, l2Rad ] = processTraceRemoval( Rad, f0, dt )
% processTraceRemoval is a data processing Subroutine explicitly created
% for the removal of static traces.
%   
%   Written by Tate Meehan, Boise State University, GreenTrACS 2016
%
% Process Data
% Filter Parameters are Established within Function
isDisplay = 0;       % Control Text Display to Screen
isNormalize = 0;     % Control Flag to Normalize Data
isDeWoW = 0;         % Control Flag to De-Wow Data
isDeTrend = 1;       % Control Flag to De-Trend Data
isMedFilt = 1;       % Control Flag for Median Subtraction Filter
isSVDfilt = 0;       % Control Flag for SVD Component Subtraction Filter
isDecon = 0;         % Control Flag for Deconvolution of Data
isBandPass = 1;      % Control Flag to Band-Pass Filter Data
isTimeZero = 0;      % Control Flag for Time-Zero Correction
isExpGain = 1;       % Control Flag for Ramped Gain of Data
isTraceBalance = 0;  % Control Flag for Trace Balanced Gain of Data
isACGGain = 0;       % Control Fla for AGC Gain of Data
isCleanNaN = 0;      % Control Flag for NaN Clean-up of Data
isStak = 1;          % Control Flag to Stack Data
isSumEnergy = 1;

    %----------------------------------------------------------------------      
    % Convert Units
    f0Hz = f0 * 1e6;        % [Hz]
    dtSec = dt * 1e-9;      % [s]
    [nsamp, ntrcs] = size(Rad);

    %----------------------------------------------------------------------
    % Normalize Data
    if isNormalize
        if isDisplay
        display( 'Begin Normalize')
        tic
        end
        
        Rad = normalize_dylan( Rad );
        
        if isDisplay
        display( 'Normalize Done')
        toc
        display(' ')
        end
    end    
    
    %----------------------------------------------------------------------      
    % De-WoW
    if isDeWoW
        if isDisplay
        display( 'Begin De-WOW Filter')
        tic
        end
        
        Rad = T8WoW(Rad,dt,f0);
        
        if isDisplay
        display( 'De-WOW Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------    
    % Remove Trends (De-WoW)
    if isDeTrend
        if isDisplay
        display( 'Begin De-Trend')
        tic
        end
        
        % Parameters
        isMeanFilt = 1;
        isDeStripe = 0;
        
        % Function
        if isMeanFilt
            Rad = detrend( Rad, 'constant' ); % Subtract Mean of Trace
        end
        if isDeStripe
            Rad = detrend( Rad, 'linear' ); % Subtract Linear Trend of Trace
        end
        
        if isDisplay
        display( 'De-Trend Done')
        toc
        display(' ')
        end
    end
        
    %----------------------------------------------------------------------
    % Median Subtraction Filter
    if isMedFilt
        if isDisplay
        display( 'Begin Median Subtraction Filter')
        tic
        end
        
        % Parameters
        MedFiltR = 2.*ceil(1/(f0Hz*dtSec))+1;% Rank of Median Subtraction Filter
        
        % Function
        Rad = medfilt1( Rad, MedFiltR, [], 2 );
        
        if isDisplay
        display( 'Median Subtraction Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % SVD First Principle Component Subtraction Filter
    if isSVDfilt
        if isDisplay
        display( 'Begin Singular Value Decomposition Filter')
        tic
        end
        
        Rad = SVDSfilter( Rad, .7 );
        
        if isDisplay
        display( 'Singular Value Decomposition Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Predictive Decon (Package Unavailable)
    if isDecon 
        % Paramters
        L = 107;                % [sample] prediction distance
        mu = 0.01;              % prewhitening
        NF = round(L*1.75);     % length of the inverse filter
        Rad2ude = Rad .* 0;     % allocate
        
        % Function
        for jj = 1 : size(Rad, 2)
            [f,o] = predictive( Rad(:,jj), NF, L, mu );
            Rad2ude(:,jj) = o;
        end
        Rad = Rad2ude;
    end

    %----------------------------------------------------------------------    
    % Band-Pass Filter
    if isBandPass
        if isDisplay
        display( 'Begin Band-Pass Filter')
        tic
        end
        
        % Parameters
        % Build Two Octave Filter About Nominal Frequency
        fMin = f0Hz/2; % [Hz]   
        fMax = f0Hz*2; % [Hz]
        
        % Function
        Rad = filter_dylan( Rad, dtSec, fMin, fMax, 2 );
        
        if isDisplay
        display( 'Band-Pass Filter Done')
        toc
        display(' ')
        end
    end
            
    %----------------------------------------------------------------------
    % Time-Zero Correction
    if isTimeZero
        if isDisplay
        display( 'Begin Time Zero Correction')
        tic
        end
        
        % Parameters
        % Establish Trigger Amplitude for Time Zero Correction
        if isCommon
            timeR = 1000;
        end
        if isMultiPlex
            timeR = 5;
        end
        
        % Function
        Rad = timeZero(Rad, timeR);
        
        if isDisplay
        display( 'Time Zero Correction Done')
        toc
        display(' ')
        end
    end

    %----------------------------------------------------------------------
    % Exponential Time Dependant Gain
    if isExpGain
        if isDisplay
        display( 'Begin Power Gain')
        tic
        end
        
        % Parameters
        tpow = 2.25;%1.5; % Filter Order, 1 is Linear
        
        % Function
        Rad = gain_bradford(Rad, tpow );
        
        if isDisplay
        display( 'Power Gain Done')
        toc
        display(' ')
        end
    end

    %----------------------------------------------------------------------
    % Trace Balanced Gain
    if isTraceBalance
        if isDisplay
        tic
        display( 'Begin Trace Balance')
        end
        
        % Parameters
        pow = .05; % Power of Gain
        
        % Function
        Rad = gainT8(Rad,pow,f0,dt);
        
        if isDisplay
        display( 'Trace Balance Done')
        toc
        display(' ')
        end
    end
    
    % ---------------------------------------------------------------------
    % Automatic Gain Control
    if isACGGain
        if isDisplay
            tic
            display( 'Begin AGC')
        end
        
        % Parameters
%         R = f0/4.*dtSec; % Quarter Wavelength Convolution Filter
        R = ceil(nsamp/2);
        type = 1; % Trace Normalize: 0 = None, 1 = RMSnorm, 2 = amplitude

        Rad = AGCgain(Rad,R,type);
        
        if isDisplay
            display( 'AGC Done')
            toc
            display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Clean Non-Existant Values
    if isCleanNaN
        if isDisplay
        tic
        display( 'Begin Clean NAN')
        end
        
        % Function
        [Rad,NaNno,datum] = cleanNaN(Rad);
        
        if isDisplay
        display( 'Clean NaN Done')
        toc
        cleanRecord = sprintf('Of %10.f Datum %10.f NaN Values Cleansed'...
            ,datum,NaNno);
        disp(cleanRecord)
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Stacking Filter
    if isStak
        if isDisplay
        display( 'Begin Stack')
        tic
        end
        
        % Parameters
        StakFiltR = 10; % Filter Rank
        
        % Function
        Rad = meanStakR(Rad,StakFiltR);
        
        if isDisplay
        display( 'Stacking Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Sum Energy
    if isSumEnergy
        l2Rad = zeros(ntrcs,1);
%         energy = Rad.^2;
        for kk = 1:ntrcs
%             sumEnergy(kk) = sum(energy(:,kk));
            l2Rad(kk) = norm(Rad(:,kk),2);
        end
    end
end