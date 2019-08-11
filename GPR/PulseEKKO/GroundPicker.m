%% Automatic Ground Picker
deconIx = 10;
for ii = 1:MD.nFiles
    groundPick = zeros(size(D.Coherence{ii},2),1);
    for jj = 1:size(D.Coherence{ii},2)
        % Adaptive Treshold
        trc = D.Coherence{ii}(:,jj);
        threshold = quantile(trc,.99);
        pickIx = find(trc > threshold);
%         filtPickIx = find(pickIx > quantile(pickIx,.05) & pickIx < quantile(pickIx,0.95));
        % Cluster Weighting
%         if std(pickIx) > 1/D.dt
        groundPick(jj) = mean(pickIx)-deconIx;
    end
    D.groundIx{ii} = groundPick;
    D.Time2Ground{ii} = groundPick.*D.dt;
end
% Noise Region Picker
% for ii = 1:MD.nFiles
%     topPick = zeros(size(D.Coherence{ii},2));
%     for jj = 1:size(D.Coherence{ii},2)
%         % Adaptive Treshold
%         trc = D.Coherence{ii}(:,jj);
%         threshold = quantile(trc,.99);
%         pickIx = find(trc > threshold);
% %         filtPickIx = find(pickIx > quantile(pickIx,.05) & pickIx < quantile(pickIx,0.95));
%         % Cluster Weighting
% %         if std(pickIx) > 1/D.dt
%         groundPick(jj) = mean(pickIx);
%     end
%     D.
%     D.Time2Ground = groundPick.*D.dt;
% end

    clear('groundPick','threshold','pickIx','trc','topPick','deconIx')
