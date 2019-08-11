function [ meanStak ] = meanStakR( data, R )
% meanStakR avereages adjacent traces determined by R; the stacking window rank.
%   Window Size is length 2R+1 and is centered about column(:,ii)
%   Inputs:     data            - Data Matrix
%               R               - Rank of Stacking Window
%
%   Outputs:    meanStak       - Stacked Data Matrix
%
% Written by Tate Meehan, Boise State University, GreenTrACS 2017

[nr, nc] = size(data);

% Appropiate Window Size
if R > nc/2+1
    R = floor(nc/2);
end
% Allocate Stacked Data Matrix
meanStak = zeros(nr,nc);

% Stack Traces
for ii = 1: nc
    if ii <= R
        % Leading Stak Edge Taper
        meanStak(:,ii) = (mean(data(:,1:ii+R),2));
    elseif ii > R && ii <= (nc-R)
        % Stack Traces
        meanStak(:,ii) = mean(data(:,ii-R:ii+R),2);
    else
        % Trailing Stak edge Taper
        meanStak(:,ii) = mean(data(:,ii-R:nc),2);
        
    end
end

end

    