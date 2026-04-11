function updateAcrossOffsetSummaries()

d = dir(fullfile('Summaries','kernelAverage_probe*.mat'));
if isempty(d); return; end

probeDirs = [];
scales = [];

for i = 1:numel(d)
    S = load(fullfile(d(i).folder, d(i).name));

    probeDirs(end+1) = S.probeDirDeg; %#ok<AGROW>

    % example metric
    if isfield(S,'compStats') && isfield(S.compStats,'scale')
        scales(end+1) = S.compStats.scale(1); %#ok<AGROW>
    else
        scales(end+1) = NaN;
    end
end

% sort
[probeDirs, idx] = sort(probeDirs);
scales = scales(idx);

% --- plot ---
figure(200); clf;
plot(probeDirs, scales, 'o-','LineWidth',1.5);
xlabel('Probe offset (deg)');
ylabel('Scale');
title('Scale vs Probe Offset');
box off;

% --- save ---
save(fullfile('Summaries','readoutAcrossOffsets.mat'), ...
    'probeDirs','scales');

end