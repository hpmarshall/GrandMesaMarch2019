function [removeTraces] = removeStaticPositions(trhd,v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

multiplexNtrcs = length(trhd);
nChan = length(unique(trhd(23,:)));

binIx = find(trhd(23,:) == 1);
rmIx = find(trhd(18,binIx) < v);
staticPlexTraces = binIx(rmIx);
% Infill Multiplex Trace Indicies For Removal
% removeBin = zeros(length(staticPlexTraces),nChan);
removeBin = [staticPlexTraces(:),staticPlexTraces(:)+1,...
    staticPlexTraces(:)+2,staticPlexTraces(:)+3];


% for ii = 1:nChan
%     if ii == nChan
%         removeBin(:,ii) = staticPlexTraces;
%     else
%     removeBin(:,ii) = staticPlexTraces - (nChan - ii);
%     end
% end
    % Indicies of Static Traces
    removeTraces = removeBin(removeBin(:) <= multiplexNtrcs);

end

