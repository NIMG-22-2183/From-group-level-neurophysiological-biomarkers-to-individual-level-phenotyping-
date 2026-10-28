clc
close all
clear all

% Local folders: edit these paths to match your data.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(repoRoot, 'data');
resultsRoot = fullfile(repoRoot, 'results');


% =========================================================================
% Grand average +- STD log PSD overlays for each channel and group
% + pairwise frequency-bin Welch t-tests with significance shown on plots
%
% Groups:
%   Group 1 = PSP 
%   Group 2 = bvFTD
%   Group 3 = CON
%
% For each subject:
%   - load SPM MEEG object D from MAT file
%   - compute PSD for each epoch and channel using pwelch
%   - convert PSD to log10 scale
%   - average across epochs within subject
%
% For each group:
%   - compute grand average across subjects
%   - compute STD across subjects
%
% Statistics:
%   - pairwise Welch two-sample t-tests at each frequency bin
%   - comparisons:
%       1) PSP vs bvFTD
%       2) PSP vs CON
%       3) bvFTD   vs CON
%   - significant bins are drawn on the corresponding channel plot
%
% Output:
%   - one overlaid PNG figure per channel: group mean +- STD log PSD
%   - significance markers for pairwise group comparisons
%   - MAT files with per-group PSD data
%   - ONE combined Excel workbook summarising stats across ALL channels
% =========================================================================


%% ==================== USER SETTINGS ====================
excelFolder = dataRoot;
excelFile = 'participants.xlsx';
excelPath   = fullfile(excelFolder, excelFile);

% Folder containing the .mat files
matFolder = fullfile(dataRoot, 'sources');

% PSD settings
fs             = 500;       % sampling rate (Hz)
useLogPower    = true;      % plot/store log10(PSD)
epsPow         = 1e-30;     % avoids log10(0)
psdFreqRange   = [1 64];    % display AND statistics frequency range
makePSDPlots   = true;

% Statistics settings
alphaSig       = 0.05;      % significance threshold
useFDR         = false;     % true = apply BH-FDR across frequencies per comparison/channel
minNPerGroup   = 2;         % minimum valid subjects per bin per group to run t-test

% Fixed channel labels in order
fixedChanLabels = {'RAI','ACC','LAI','MPFC','PCC','LAG','RAG'};

%% ==================== STEP 1: LOAD EXCEL ====================
if ~isfile(excelPath)
    error('Excel file not found: %s', excelPath);
end

data = readtable(excelPath);

% Clean IDs (remove "sub-" prefix)
IDs    = regexprep(data{:,1}, '^sub-', '');
Groups = data{:,2};

% Define group labels
group1_labels = {'PSP'};
group2_labels = {'bvFTD'};
group3_labels = {'CON'};

%% ==================== STEP 2: LIST MAT FILES ====================
matFiles = dir(fullfile(matFolder, 'XX_sub-Sub0*_DSALL.mat'));
matNames = {matFiles.name}';

% Extract ID part (XX_sub-Sub0XXX) from MAT filenames
matIDs = regexp(matNames, '^(XX_sub-Sub0\d+)_DSALL\.mat$', 'tokens', 'once');
validMask = ~cellfun(@isempty, matIDs);
matNames = matNames(validMask);
matIDs   = matIDs(validMask);
matIDs   = cellfun(@(c) c{1}, matIDs, 'UniformOutput', false);

%% ==================== STEP 3: MATCH MAT FILES TO GROUPS ====================
Group1_Files = {};
Group2_Files = {};
Group3_Files = {};

for i = 1:numel(matIDs)
    thisIDF = matIDs{i};
    thisID  = regexp(thisIDF, 'Sub\d+', 'match', 'once');
    idx     = strcmp(IDs, thisID);

    if any(idx)
        label = Groups{find(idx,1,'first')};

        if ismember(label, group1_labels)
            Group1_Files{end+1} = fullfile(matFolder, matNames{i}); 
        elseif ismember(label, group2_labels)
            Group2_Files{end+1} = fullfile(matFolder, matNames{i}); 
        elseif ismember(label, group3_labels)
            Group3_Files{end+1} = fullfile(matFolder, matNames{i}); 
        end
    else
        warning('No group found for MAT file: %s', matNames{i});
    end
end

%% ==================== STEP 4: DISPLAY RESULTS ====================
disp('Group 1 MAT files (PSP only):'); disp(Group1_Files');
disp('Group 2 MAT files (bvFTD):');    disp(Group2_Files');
disp('Group 3 MAT files (CON):');      disp(Group3_Files');

%% ==================== SETUP ====================
if isstring(Group1_Files), Group1_Files = cellstr(Group1_Files); end
if isstring(Group2_Files), Group2_Files = cellstr(Group2_Files); end
if isstring(Group3_Files), Group3_Files = cellstr(Group3_Files); end

g1 = Group1_Files(:);   % PSP only
g2 = Group2_Files(:);   % bvFTD
g3 = Group3_Files(:);   % CON

groups     = {g1, g2, g3};
groupNames = {'PSP','bvFTD','CON'};
nGroups    = numel(groups);

% Output folder
saveFolder = fullfile(resultsRoot, 'psd_overlays');
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

% Containers for group-level PSD results
groupPSD        = cell(nGroups,1);   % each cell: [nChans x nFreq x nSubjects]
groupF          = cell(nGroups,1);   % frequency vector
groupChanLabels = cell(nGroups,1);   % channel labels
groupSubjIDs    = cell(nGroups,1);   % subject IDs
groupFilesUsed  = cell(nGroups,1);   % files actually used

%% ==================== LOAD & AGGREGATE PER GROUP ====================
for g = 1:nGroups
    fileList = groups{g};
    N = numel(fileList);

    fprintf('\n============================================================\n');
    fprintf('Processing group: %s\n', groupNames{g});
    fprintf('Number of candidate files: %d\n', N);
    fprintf('============================================================\n');

    PSD_subjects = [];      % [nChans x nFreq x K]
    f_ref        = [];
    subjIDs      = strings(1, N);
    filesUsed    = cell(1, N);
    chanLabels   = {};
    k            = 0;

    for fIdx = 1:N
        matFull = fileList{fIdx};
        fprintf('\nProcessing: %s\n', matFull);

        % ---- Load MAT file ----
        try
            S = load(matFull, 'D');
        catch ME
            warning('Could not load %s. Reason: %s', matFull, ME.message);
            continue;
        end

        if ~isfield(S, 'D')
            warning('Variable D not found in %s. Skipping.', matFull);
            continue;
        end

        D = S.D;
        Y = D.data;

        % ---- Repair stale .dat path if needed ----
        try
            needsRepair = ~exist(Y.fname, 'file');
        catch
            needsRepair = false;
        end

        if needsRepair
            [mPath,mName,~] = fileparts(matFull);
            datGuess = fullfile(mPath, [mName '.dat']);

            if exist(datGuess,'file')
                try
                    Y = file_array(datGuess, Y.dim, Y.dtype, Y.offset, Y.scl_slope, Y.scl_inter);
                    D.data = Y;
                    fprintf('Repaired .dat path to: %s\n', datGuess);
                catch ME
                    warning('Failed to rebuild file_array for %s: %s', matFull, ME.message);
                    continue;
                end
            else
                warning('Missing .dat for %s. Expected: %s. Skipping subject.', matFull, datGuess);
                continue;
            end
        end

        % ---- Safe size handling ----
        sz = size(Y);    % expected [nChans nSamples nEpochs]
        if numel(sz) < 2
            warning('Unexpected data dims for %s (size=%s). Skipping.', matFull, mat2str(sz));
            continue;
        end

        nChans   = sz(1);
        nSamples = sz(2);
        nEpochs  = 1;
        if numel(sz) >= 3
            nEpochs = max(1, sz(3));
        end

        if nSamples < 8
            warning('Too few samples in %s. Skipping.', matFull);
            continue;
        end

        % ---- Frequency vector from representative segment ----
        try
            tmpSeg = detrend(double(Y(1,:,1)));
            [~, fWelch] = pwelch(tmpSeg, hamming(nSamples), 0, nSamples, fs);
        catch ME
            warning('Could not compute representative PSD for %s: %s', matFull, ME.message);
            continue;
        end
        nFreq = numel(fWelch);

        % ---- Channel labels ----
        if isempty(chanLabels)
            if nChans == numel(fixedChanLabels)
                chanLabels = fixedChanLabels(:);
            else
                chanLabels = arrayfun(@(x) sprintf('Ch%d', x), 1:nChans, 'UniformOutput', false)';
            end
        end

        % ---- Subject PSD storage ----
        subjPSD = nan(nChans, nFreq);

        % ---- Compute PSD channel by channel ----
        for ch = 1:nChans
            psd_epochs = nan(nFreq, nEpochs);

            for e = 1:nEpochs
                try
                    seg = detrend(double(Y(ch,:,e)));
                    [Pxx, ~] = pwelch(seg, hamming(nSamples), 0, nSamples, fs);
                    psd_epochs(:,e) = Pxx(:);
                catch
                    % leave epoch as NaN
                end
            end

            valid_psd = all(isfinite(psd_epochs),1) & all(psd_epochs >= 0, 1);
            psd_epochs_valid = psd_epochs(:, valid_psd);

            if isempty(psd_epochs_valid)
                continue;
            end

            if useLogPower
                subjPSD(ch,:) = mean(log10(psd_epochs_valid + epsPow), 2, 'omitnan');
            else
                subjPSD(ch,:) = mean(psd_epochs_valid, 2, 'omitnan');
            end
        end

        % ---- Accept subject ----
        if all(all(isnan(subjPSD)))
            warning('All PSD values NaN for %s. Skipping subject.', matFull);
            continue;
        end

        k = k + 1;

        if k == 1
            PSD_subjects = nan(nChans, nFreq, 0);
            f_ref = fWelch(:);
        else
            if numel(fWelch) ~= numel(f_ref) || any(abs(fWelch(:) - f_ref(:)) > 1e-10)
                warning('Frequency vector mismatch in %s. Skipping subject.', matFull);
                k = k - 1;
                continue;
            end
            if nChans ~= size(PSD_subjects,1)
                warning('Channel count mismatch in %s. Expected %d, got %d. Skipping subject.', ...
                    matFull, size(PSD_subjects,1), nChans);
                k = k - 1;
                continue;
            end
        end

        PSD_subjects(:,:,k) = subjPSD;

        tok = regexp(matFull, '(Sub0\d+)', 'tokens', 'once');
        if ~isempty(tok)
            subjIDs(k) = string(tok{1});
        else
            subjIDs(k) = sprintf("sub%03d", k);
        end
        filesUsed{k} = matFull;
    end

    if k == 0
        warning('No valid subjects found in group %s.', groupNames{g});
        groupPSD{g}        = [];
        groupF{g}          = [];
        groupChanLabels{g} = {};
        groupSubjIDs{g}    = strings(0,1);
        groupFilesUsed{g}  = {};
        continue;
    end

    PSD_subjects = PSD_subjects(:,:,1:k);
    subjIDs      = subjIDs(1:k);
    filesUsed    = filesUsed(1:k);

    groupPSD{g}        = PSD_subjects;
    groupF{g}          = f_ref;
    groupChanLabels{g} = chanLabels;
    groupSubjIDs{g}    = subjIDs;
    groupFilesUsed{g}  = filesUsed;

    % Save raw per-group PSD results
    results = struct();
    results.groupName    = groupNames{g};
    results.filesUsed    = filesUsed;
    results.subjIDs      = subjIDs;
    results.PSD_subjects = PSD_subjects;   % [nChans x nFreq x nSubjects]
    results.f            = f_ref;
    results.chanLabels   = chanLabels;
    results.fs           = fs;
    results.useLogPower  = useLogPower;
    results.psdFreqRange = psdFreqRange;

    outMat = fullfile(saveFolder, sprintf('Group_%s_LogPSD_results.mat', regexprep(groupNames{g}, '\W', '')));
    save(outMat, 'results');
    fprintf('Saved group PSD results: %s\n', outMat);
end

%% ==================== PLOT OVERLAID GRAND AVERAGE +- STD + STATS ====================
allStatsTables = {};
allSummaryTables = {};

if makePSDPlots
    refg = find(~cellfun(@isempty, groupPSD), 1, 'first');

    if isempty(refg)
        warning('No valid group PSD data available for plotting.');
    else
        fPlot = groupF{refg};
        nChansRef = size(groupPSD{refg}, 1);

        % Restrict BOTH plotting and stats to 1-64 Hz
        fMask = (fPlot >= psdFreqRange(1)) & (fPlot <= psdFreqRange(2));
        fShow = fPlot(fMask);

        chanLabelsRef = groupChanLabels{refg};
        if isempty(chanLabelsRef)
            if nChansRef == numel(fixedChanLabels)
                chanLabelsRef = fixedChanLabels;
            else
                chanLabelsRef = arrayfun(@(x) sprintf('Ch%d', x), 1:nChansRef, 'UniformOutput', false);
            end
        end

        % fixed colors for groups
        groupColors = [
            0.0000    0.4470    0.7410   % PSP
            0.8500    0.3250    0.0980   % bvFTD
            0.4660    0.6740    0.1880   % CON
            ];

        % pairwise comparisons
        pairIdx = [
            1 2
            1 3
            2 3
            ];

        pairNames = {
            'PSP vs bvFTD'
            'PSP vs CON'
            'bvFTD vs CON'
            };

        pairColors = [
            0.4940    0.1840    0.5560
            0.6350    0.0780    0.1840
            0.3010    0.7450    0.9330
            ];

        for ch = 1:nChansRef
            fig = figure('Color','w', ...
                'Name', sprintf('Overlay_GrandAvgSTD_LogPSD_Channel_%d', ch), ...
                'Units','pixels', 'Position', [100 100 1100 760]);
            hold on;

            hLine = gobjects(0);
            legTxt = {};

            allMu = [];
            allSd = [];

            % -------- Plot grand mean +- STD for groups --------
            for g = 1:nGroups
                if isempty(groupPSD{g})
                    continue;
                end

                PSDg = groupPSD{g};   % [nChans x nFreq x nSubjects]
                fg   = groupF{g};

                if size(PSDg,1) < ch
                    continue;
                end

                if numel(fg) ~= numel(fPlot) || any(abs(fg(:) - fPlot(:)) > 1e-10)
                    warning('Skipping group %s for channel %d due to frequency mismatch.', groupNames{g}, ch);
                    continue;
                end

                X = squeeze(PSDg(ch,:,:))';   % [nSubjects x nFreq]
                if isempty(X)
                    continue;
                end

                mu = mean(X, 1, 'omitnan');
                sd = std(X, 0, 1, 'omitnan');

                mu = mu(fMask);
                sd = sd(fMask);

                allMu = [allMu, mu(:)']; %#ok<AGROW>
                allSd = [allSd, sd(:)']; %#ok<AGROW>

                c = groupColors(min(g,size(groupColors,1)),:);

                fill([fShow; flipud(fShow)], ...
                     [mu(:)-sd(:); flipud(mu(:)+sd(:))], ...
                     c, 'FaceAlpha', 0.18, 'EdgeColor', 'none');

                h = plot(fShow, mu, 'LineWidth', 2.5, 'Color', c);
                hLine(end+1) = h; %#ok<SAGROW>
                legTxt{end+1} = sprintf('%s (n=%d)', groupNames{g}, size(X,1)); %#ok<SAGROW>
            end

            grid off;
            box on;
            xlabel('Frequency (Hz)');

            if useLogPower
                ylabel('Grand average log_{10}(PSD)');
            else
                ylabel('Grand average PSD');
            end

            if ch <= numel(chanLabelsRef)
                chLabel = chanLabelsRef{ch};
            else
                chLabel = sprintf('Ch%d', ch);
            end

            xlim([psdFreqRange(1) psdFreqRange(2)]);

            % -------- Statistics --------
            statsOut = struct();
            statsOut.channel      = ch;
            statsOut.channelLabel = chLabel;
            statsOut.f            = fShow;
            statsOut.alpha        = alphaSig;
            statsOut.useFDR       = useFDR;
            statsOut.comparisons  = struct([]);

            % determine y-limits based on curves
            yDataMin = min(allMu - allSd, [], 'omitnan');
            yDataMax = max(allMu + allSd, [], 'omitnan');
            if isempty(yDataMin) || ~isfinite(yDataMin), yDataMin = -1; end
            if isempty(yDataMax) || ~isfinite(yDataMax), yDataMax =  1; end
            if yDataMax <= yDataMin
                yDataMax = yDataMin + 1;
            end
            yRange = yDataMax - yDataMin;

            % reserve some headroom for significance rows
            topPadFrac = 0.18;
            ylim([yDataMin - 0.05*yRange, yDataMax + topPadFrac*yRange]);

            ySigBase = yDataMax + 0.05*yRange;
            ySigStep = 0.045*yRange;

            sigHandles = gobjects(0);
            sigLegTxt  = {};

            for p = 1:size(pairIdx,1)
                gA = pairIdx(p,1);
                gB = pairIdx(p,2);

                if isempty(groupPSD{gA}) || isempty(groupPSD{gB})
                    continue;
                end

                XAfull = squeeze(groupPSD{gA}(ch,:,:))';   % [nSubA x nFreq]
                XBfull = squeeze(groupPSD{gB}(ch,:,:))';   % [nSubB x nFreq]

                if isempty(XAfull) || isempty(XBfull)
                    continue;
                end

                XA = XAfull(:, fMask);
                XB = XBfull(:, fMask);

                nFreqShow = numel(fShow);
                pvals  = nan(1, nFreqShow);
                tstats = nan(1, nFreqShow);
                dfvals = nan(1, nFreqShow);
                nAvals = zeros(1, nFreqShow);
                nBvals = zeros(1, nFreqShow);

                for kf = 1:nFreqShow
                    xa = XA(:,kf);
                    xb = XB(:,kf);

                    xa = xa(isfinite(xa));
                    xb = xb(isfinite(xb));

                    nAvals(kf) = numel(xa);
                    nBvals(kf) = numel(xb);

                    if numel(xa) < minNPerGroup || numel(xb) < minNPerGroup
                        continue;
                    end

                    try
                        [~, ptmp, ~, statsTmp] = ttest2(xa, xb, ...
                            'Vartype', 'unequal', ...
                            'Tail', 'both');

                        pvals(kf) = ptmp;

                        if isfield(statsTmp, 'tstat')
                            tstats(kf) = statsTmp.tstat;
                        end
                        if isfield(statsTmp, 'df')
                            dfvals(kf) = statsTmp.df;
                        end
                    catch
                        % leave as NaN
                    end
                end

                if useFDR
                    sigMask = bh_fdr_mask(pvals, alphaSig);
                else
                    sigMask = pvals < alphaSig;
                end
                sigMask(~isfinite(pvals)) = false;

                ySig = ySigBase + (p-1)*ySigStep;
                cSig = pairColors(p,:);

                % show significant frequency bins as dots
                if any(sigMask)
                    hs = plot(fShow(sigMask), ySig*ones(1,sum(sigMask)), '.', ...
                        'Color', cSig, 'MarkerSize', 12);
                    sigHandles(end+1) = hs; 
                    if useFDR
                        sigLegTxt{end+1} = sprintf('%s (q<%.3f)', pairNames{p}, alphaSig); 
                    else
                        sigLegTxt{end+1} = sprintf('%s (p<%.3f)', pairNames{p}, alphaSig); 
                    end
                else
                    hs = plot(nan, nan, '.', 'Color', cSig, 'MarkerSize', 12);
                    sigHandles(end+1) = hs; 
                    if useFDR
                        sigLegTxt{end+1} = sprintf('%s (no q<%.3f bins)', pairNames{p}, alphaSig); 
                    else
                        sigLegTxt{end+1} = sprintf('%s (no p<%.3f bins)', pairNames{p}, alphaSig); 
                    end
                end

                statsOut.comparisons(p).name        = pairNames{p};
                statsOut.comparisons(p).groupA      = groupNames{gA};
                statsOut.comparisons(p).groupB      = groupNames{gB};
                statsOut.comparisons(p).pvals       = pvals;
                statsOut.comparisons(p).tstats      = tstats;
                statsOut.comparisons(p).df          = dfvals;
                statsOut.comparisons(p).nA          = nAvals;
                statsOut.comparisons(p).nB          = nBvals;
                statsOut.comparisons(p).sigMask     = sigMask;
                statsOut.comparisons(p).f           = fShow;
                statsOut.comparisons(p).alpha       = alphaSig;
                statsOut.comparisons(p).useFDR      = useFDR;
                statsOut.comparisons(p).nSigBins    = sum(sigMask);
            end

            title(sprintf('%s — Grand average \\pm STD log PSD + pairwise t-tests', chLabel), ...
                'Interpreter', 'tex');

            % Combined legend
            allHandles = [hLine(:); sigHandles(:)];
            allLegTxt  = [legTxt(:); sigLegTxt(:)];
            if ~isempty(allHandles)
                legend(allHandles, allLegTxt, 'Location', 'best');
            end

            safeLbl = regexprep(chLabel, '[^a-zA-Z0-9]', '');
            if isempty(safeLbl)
                safeLbl = sprintf('Ch%d', ch);
            end

            outPng  = fullfile(saveFolder, sprintf('Overlay_GrandAvgSTD_LogPSD_%s.png', safeLbl));
            outMat  = fullfile(saveFolder, sprintf('Overlay_GrandAvgSTD_LogPSD_%s_data.mat', safeLbl));

            savepng(fig, outPng);
            save(outMat, 'statsOut');

            Tstats = stats_struct_to_table(statsOut);
            Tsum   = stats_summary_table(statsOut);

            allStatsTables{end+1}   = Tstats; 
            allSummaryTables{end+1} = Tsum;   

            fprintf('Saved overlay plot: %s\n', outPng);
            close(fig);
        end
    end
end

%% ==================== SAVE ONE COMBINED EXCEL FOR ALL CHANNELS ====================
allStatsXlsx = fullfile(saveFolder, 'AllChannels_LogPSD_Stats_Summary.xlsx');

try
    if isempty(allStatsTables)
        T_all = table();
    else
        T_all = vertcat(allStatsTables{:});
    end

    if isempty(allSummaryTables)
        T_sum = table();
    else
        T_sum = vertcat(allSummaryTables{:});
    end

    writetable(T_all, allStatsXlsx, 'Sheet', 'AllStats_LongFormat');
    writetable(T_sum, allStatsXlsx, 'Sheet', 'Summary_ByChannel_Comparison');

    fprintf('Saved combined Excel summary: %s\n', allStatsXlsx);
catch ME
    warning('Could not write combined Excel summary: %s', ME.message);
end

fprintf('\nAll done.\n');
fprintf('Grand-average +- STD log PSD PNG figures and combined stats saved in:\n%s\n', saveFolder);

%% ==================== HELPER FUNCTION ====================
function savepng(fig, filename)
    [fpath,~,~] = fileparts(filename);
    if ~exist(fpath,'dir')
        mkdir(fpath);
    end
    try
        exportgraphics(fig, filename, 'Resolution', 300);
    catch
        print(fig, '-dpng', '-r300', filename);
    end
end

%% ==================== HELPER FUNCTION ====================
function sigMask = bh_fdr_mask(pvals, q)
    % Benjamini-Hochberg FDR mask
    pvals = pvals(:);
    sigMask = false(size(pvals));

    valid = isfinite(pvals);
    pv = pvals(valid);

    if isempty(pv)
        sigMask = reshape(sigMask, size(pvals));
        return;
    end

    [ps, ~] = sort(pv);
    m = numel(ps);
    thresh = (1:m)'/m * q;
    below = ps <= thresh;

    if any(below)
        k = find(below, 1, 'last');
        cutoff = ps(k);
        sigMask(valid) = pv <= cutoff;
    end

    sigMask = sigMask(:).';
end

%% ==================== HELPER FUNCTION ====================
function T = stats_struct_to_table(statsOut)
    % Convert statsOut for one channel into long-form table
    rows = [];

    for i = 1:numel(statsOut.comparisons)
        C = statsOut.comparisons(i);
        nF = numel(C.f);

        Ti = table();
        Ti.Channel         = repmat(statsOut.channel, nF, 1);
        Ti.ChannelLabel    = repmat(string(statsOut.channelLabel), nF, 1);
        Ti.Comparison      = repmat(string(C.name), nF, 1);
        Ti.GroupA          = repmat(string(C.groupA), nF, 1);
        Ti.GroupB          = repmat(string(C.groupB), nF, 1);
        Ti.FrequencyHz     = C.f(:);
        Ti.nGroupA         = C.nA(:);
        Ti.nGroupB         = C.nB(:);
        Ti.tStatistic      = C.tstats(:);
        Ti.DF              = C.df(:);
        Ti.pValue          = C.pvals(:);
        Ti.Significant     = C.sigMask(:);
        Ti.Alpha           = repmat(C.alpha, nF, 1);
        Ti.UseFDR          = repmat(C.useFDR, nF, 1);

        if isempty(rows)
            rows = Ti;
        else
            rows = [rows; Ti]; 
        end
    end

    if isempty(rows)
        T = table();
    else
        T = rows;
    end
end

%% ==================== HELPER FUNCTION ====================
function Tsum = stats_summary_table(statsOut)
    % Summary table: one row per comparison for one channel
    rows = [];

    for i = 1:numel(statsOut.comparisons)
        C = statsOut.comparisons(i);

        sigIdx = find(C.sigMask);
        if isempty(sigIdx)
            firstSig = NaN;
            lastSig  = NaN;
            sigFreqs = "";
        else
            firstSig = C.f(sigIdx(1));
            lastSig  = C.f(sigIdx(end));
            sigFreqs = strjoin(string(round(C.f(sigIdx), 6)), ', ');
        end

        Ti = table();
        Ti.Channel            = statsOut.channel;
        Ti.ChannelLabel       = string(statsOut.channelLabel);
        Ti.Comparison         = string(C.name);
        Ti.GroupA             = string(C.groupA);
        Ti.GroupB             = string(C.groupB);
        Ti.Alpha              = C.alpha;
        Ti.UseFDR             = C.useFDR;
        Ti.NumTestedBins      = numel(C.f);
        Ti.NumSignificantBins = sum(C.sigMask);
        Ti.FirstSigFreqHz     = firstSig;
        Ti.LastSigFreqHz      = lastSig;
        Ti.SignificantFreqHz  = sigFreqs;

        if isempty(rows)
            rows = Ti;
        else
            rows = [rows; Ti]; 
        end
    end

    if isempty(rows)
        Tsum = table();
    else
        Tsum = rows;
    end
end
