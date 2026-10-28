function CLASSIFY_FEATURE_TABLE_FROM_STATS_REPORT_LOO()
% CLASSIFY_FEATURE_TABLE_FROM_STATS_REPORT_LOO
%
% Classification from precomputed subject-feature table and feature-statistics
% report, using leave-one-out (LOO) validation.
% Cohorts: CON, bvFTD and PSP. PC is not merged with PSP.
%
% INPUT FILES
%   1) Features_AllSubjects_Optimized.xlsx
%      Sheet: 'Features'
%
%   2) Feature_Statistics_Report_AllComparisons.xlsx
%      Sheets used explicitly:
%         'CON_vs_bvFTD'
%         'CON_vs_PSP'
%         'PSP_vs_bvFTD'
%         'CON_vs_Patient'   % Patient = bvFTD + PSP only
%         'TopHits_AllComparisons'
%
% OUTPUT
%   results/psd_probe/
%
% FEATURE SELECTION RULES
%   P_TOP5 / P_TOP10         -> smallest p
%   D_TOP5 / D_TOP10         -> largest |Cohen d|
%   COMBINED_TOP5 / TOP10    -> z(-log10 p) + z(|d|)
%   STABLE_TOP5 / TOP10      -> combined + recurrence bonus
%
% -------------------------------------------------------------------------

clc; close all; clearvars;

% Local folders: edit these paths to match your data.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(repoRoot, 'data');
resultsRoot = fullfile(repoRoot, 'results');


%% ====================== USER SETTINGS ===========================

featureFileCandidates = {fullfile(resultsRoot, 'features', 'Features_AllSubjects_Optimized.xlsx')};

statsFileCandidates = {fullfile(resultsRoot, 'feature_statistics', 'Feature_Statistics_Report_AllComparisons.xlsx')};

% Four contrasts from the feature-statistics workbook.

comparisonSheets = cell(4,1);
comparisonSheets{1} = 'CON_vs_bvFTD';
comparisonSheets{2} = 'CON_vs_PSP';
comparisonSheets{3} = 'PSP_vs_bvFTD';
comparisonSheets{4} = 'CON_vs_Patient';  % Patient = bvFTD + PSP only

topHitsSheet = 'TopHits_AllComparisons';

rngSeed   = 1;
nTreesRF  = 300;
nTreesAda = 200;

panelSizes = [5 10];
minPanelSize = 2;

outRoot = fullfile(resultsRoot, 'psd_probe');
if ~exist(outRoot, 'dir'), mkdir(outRoot); end

%% ====================== RESOLVE FILE PATHS ======================

featureFile = resolve_existing_file(featureFileCandidates, 'Features_AllSubjects_Optimized.xlsx');
statsFile   = resolve_existing_file(statsFileCandidates, 'Feature_Statistics_Report_AllComparisons.xlsx');

fprintf('Using feature file:\n%s\n', featureFile);
fprintf('Using stats file:\n%s\n', statsFile);

%% ====================== LOAD FEATURE TABLE ======================

Tfeat = readtable(featureFile, 'Sheet', 'Features');

requiredFeatureVars = {'SubjectID','Group'};
for i = 1:numel(requiredFeatureVars)
    if ~ismember(requiredFeatureVars{i}, Tfeat.Properties.VariableNames)
        error('Feature table is missing required column: %s', requiredFeatureVars{i});
    end
end

if isstring(Tfeat.SubjectID), Tfeat.SubjectID = cellstr(Tfeat.SubjectID); end
if ischar(Tfeat.SubjectID),   Tfeat.SubjectID = cellstr(Tfeat.SubjectID); end

if isstring(Tfeat.Group), Tfeat.Group = cellstr(Tfeat.Group); end
if ischar(Tfeat.Group),   Tfeat.Group = cellstr(Tfeat.Group); end

Tfeat.SubjectID = regexprep(Tfeat.SubjectID, '^sub-', '');
Tfeat.SubjectID = regexprep(Tfeat.SubjectID, '.*?(Sub0*\d+).*', '$1');
Tfeat.Group     = cellfun(@char, Tfeat.Group, 'UniformOutput', false);

% Keep PSP separate from PC.
% Exact PSP labels become 'PSP'; exact PC labels become 'PC' and are not used
% in any PSP comparison. Combined PSP/PC labels are marked as excluded,
% because they cannot be separated into PSP-only subjects here.
Tfeat.Group     = standardize_group_labels_psp_only(Tfeat.Group);

fprintf('Group counts in feature table:\n');
disp(groupcounts(categorical(Tfeat.Group)));

allVarNames = Tfeat.Properties.VariableNames;
featureVars = setdiff(allVarNames, {'SubjectID','Group'}, 'stable');

if isempty(featureVars)
    error('No feature columns found in Features sheet.');
end

Xall      = table2array(Tfeat(:, featureVars));
allLabels = Tfeat.Group(:);
allIDs    = Tfeat.SubjectID(:);

%% ====================== LOAD TOPHITS IF PRESENT =================

TopHits = [];
try
    TopHits = readtable(statsFile, 'Sheet', topHitsSheet);
catch
    warning('Could not read TopHits_AllComparisons sheet. Continuing without recurrence bonus.');
    TopHits = [];
end

%% ====================== BUILD GLOBAL RECURRENCE BONUS ===========

recurrenceBonus = build_recurrence_bonus(TopHits, featureVars);

%% ====================== RUN EACH COMPARISON =====================

summaryMetricsAll = table();

for s = 1:numel(comparisonSheets)

    sheetName = comparisonSheets{s};

    try
        Tstats = readtable(statsFile, 'Sheet', sheetName);
    catch ME
        warning('Skipping %s: %s', sheetName, ME.message);
        continue;
    end

    Tstats = clean_stats_table_feature_report(Tstats);

    if isempty(Tstats)
        warning('Skipping %s because stats table is empty or malformed.', sheetName);
        continue;
    end

    g1 = unique(Tstats.GroupA, 'stable');
    g2 = unique(Tstats.GroupB, 'stable');

    if numel(g1) ~= 1 || numel(g2) ~= 1
        warning('Skipping %s because GroupA / GroupB are not uniquely defined.', sheetName);
        continue;
    end

    group1FromStats = char(g1{1});
    group2FromStats = char(g2{1});

    % Resolve actual row-selection groups. This is especially important for
    % CON_vs_Patient, where the stats sheet may say GroupB='Patient', but
    % the feature table contains the real labels bvFTD and PSP.
    [group1List, group2List, group1, group2] = resolve_comparison_groups_psp_only( ...
        sheetName, group1FromStats, group2FromStats);

    fprintf('\n============================================================\n');
    fprintf('Processing comparison: %s\n', sheetName);
    fprintf('Group1: %s\n', group1);
    fprintf('Group2: %s\n', group2);
    fprintf('============================================================\n');

    idx1 = ismember(allLabels, group1List);
    idx2 = ismember(allLabels, group2List);

    if sum(idx1) < 2 || sum(idx2) < 2
        warning('%s skipped: need at least 2 subjects in each group.', sheetName);
        continue;
    end

    [tf, idxFeat] = ismember(Tstats.Feature, featureVars);
    Tstats  = Tstats(tf,:);
    idxFeat = idxFeat(tf);

    if isempty(Tstats)
        warning('%s skipped: no overlapping features found in feature table.', sheetName);
        continue;
    end

    p  = Tstats.p_Welch;
    d  = Tstats.Cohen_d;
    lp = -log10(max(p, realmin));
    ad = abs(d);

    score_p        = -p;
    score_d        = ad;
    score_combined = safe_zscore(lp) + safe_zscore(ad);
    score_stable   = score_combined + 0.50 * recurrenceBonus(idxFeat);

    panelDefs = {};

    for kk = 1:numel(panelSizes)
        K = panelSizes(kk);

        panelDefs(end+1,:) = {sprintf('P_TOP%d', K),        rank_select(idxFeat, score_p,        K, 'descend')}; %#ok<AGROW>
        panelDefs(end+1,:) = {sprintf('D_TOP%d', K),        rank_select(idxFeat, score_d,        K, 'descend')}; %#ok<AGROW>
        panelDefs(end+1,:) = {sprintf('COMBINED_TOP%d', K), rank_select(idxFeat, score_combined, K, 'descend')}; %#ok<AGROW>
        panelDefs(end+1,:) = {sprintf('STABLE_TOP%d', K),   rank_select(idxFeat, score_stable,   K, 'descend')}; %#ok<AGROW>
    end

    compDir = fullfile(outRoot, matlab.lang.makeValidName(sheetName));
    if ~exist(compDir, 'dir'), mkdir(compDir); end

    panelLog = table( ...
        cell(0,1), ...
        cell(0,1), ...
        zeros(0,1), ...
        cell(0,1), ...
        'VariableNames', {'Comparison','Panel','NumFeatures','FeatureList'});

    for pp = 1:size(panelDefs,1)
        panelName = panelDefs{pp,1};
        featIdx = unique(panelDefs{pp,2}, 'stable');
        featIdx = featIdx(~isnan(featIdx));

        if numel(featIdx) < minPanelSize
            warning('%s / %s skipped: too few selected features.', sheetName, panelName);
            continue;
        end

        useRows = idx1 | idx2;
        Xcmp    = Xall(useRows, featIdx);
        yText   = allLabels(useRows);
        subjIDs = allIDs(useRows);

        Y = double(ismember(yText, group2List));

        keepRows = ~any(isnan(Xcmp),2);
        Xcmp    = Xcmp(keepRows,:);
        Y       = Y(keepRows);
        yText   = yText(keepRows);
        subjIDs = subjIDs(keepRows);

        nPos = sum(Y==1);
        nNeg = sum(Y==0);

        if nPos < 2 || nNeg < 2
            warning('%s / %s skipped after NaN removal: too few subjects.', sheetName, panelName);
            continue;
        end

        selectedNames = featureVars(featIdx);

        %% Save selected features table
        featTable = table;
        featTable.SubjectID = subjIDs(:);
        featTable.Group     = yText(:);

        for j = 1:numel(selectedNames)
            safeName = matlab.lang.makeValidName(selectedNames{j});
            featTable.(safeName) = Xcmp(:,j);
        end

        featXlsx = fullfile(compDir, sprintf('SelectedFeatures_%s.xlsx', panelName));
        writetable(featTable, featXlsx, 'Sheet', 'SelectedFeatures');

        featInfo = table((1:numel(selectedNames))', featIdx(:), selectedNames(:), ...
            'VariableNames', {'FeatureOrder','FeatureColumnIndex','FeatureName'});
        writetable(featInfo, featXlsx, 'Sheet', 'FeatureInfo');

        featSummary = table;
        featSummary.Comparison  = {sheetName};
        featSummary.Panel       = {panelName};
        featSummary.Group1      = {group1};
        featSummary.Group2      = {group2};
        featSummary.NumSubjects = height(featTable);
        featSummary.NumGroup1   = sum(Y==0);
        featSummary.NumGroup2   = sum(Y==1);
        featSummary.NumFeatures = numel(selectedNames);
        writetable(featSummary, featXlsx, 'Sheet', 'Summary');

        newRow = table({sheetName}, {panelName}, numel(selectedNames), {strjoin(selectedNames, ', ')}, ...
            'VariableNames', {'Comparison','Panel','NumFeatures','FeatureList'});
        panelLog = [panelLog; newRow];

        %% LOO classification
        [metricsTable, predTable, rocStruct] = run_models_loo_feature_table( ...
            Xcmp, Y, yText, subjIDs, panelName, sheetName, group1, group2, rngSeed, nTreesRF, nTreesAda);

        classXlsx = fullfile(compDir, sprintf('Classification_%s.xlsx', panelName));
        writetable(metricsTable, classXlsx, 'Sheet', 'Metrics');
        writetable(predTable, classXlsx, 'Sheet', 'Predictions');
        writetable(featInfo, classXlsx, 'Sheet', 'SelectedFeatures');

        modelNames = fieldnames(rocStruct);
        for mm = 1:numel(modelNames)
            modelName = modelNames{mm};
            rocTab = table(rocStruct.(modelName).fpr(:), rocStruct.(modelName).tpr(:), ...
                'VariableNames', {'FPR','TPR'});
            writetable(rocTab, classXlsx, 'Sheet', safe_sheet_name(modelName));
        end

        figDir = fullfile(compDir, 'figures', panelName);
        if ~exist(figDir, 'dir'), mkdir(figDir); end
        save_roc_and_confusion_figures_feature_table(rocStruct, predTable, figDir, panelName, sheetName);

        tmp = metricsTable;
        tmp.Comparison = repmat({sheetName}, height(tmp), 1);
        tmp.Panel      = repmat({panelName}, height(tmp), 1);
        summaryMetricsAll = [summaryMetricsAll; tmp];
    end

    if ~isempty(panelLog)
        writetable(panelLog, fullfile(compDir, 'PanelLog.xlsx'), 'Sheet', 'Panels');
    end
end

%% ====================== GLOBAL SUMMARY ==========================

if ~isempty(summaryMetricsAll)
    summaryXlsx = fullfile(outRoot, 'Summary_AllComparisons.xlsx');
    writetable(summaryMetricsAll, summaryXlsx, 'Sheet', 'AllMetrics');

    [~, ord] = sortrows(summaryMetricsAll(:, {'Comparison','AUC'}), {'Comparison','AUC'}, {'ascend','descend'});
    Tsorted = summaryMetricsAll(ord,:);

    bestRows = false(height(Tsorted),1);
    comps = unique(Tsorted.Comparison, 'stable');
    for i = 1:numel(comps)
        idx = strcmp(Tsorted.Comparison, comps{i});
        ii = find(idx, 1, 'first');
        if ~isempty(ii), bestRows(ii) = true; end
    end

    Tbest = Tsorted(bestRows,:);
    writetable(Tbest, summaryXlsx, 'Sheet', 'BestByComparison');
end

fprintf('\nAll done.\nOutputs in:\n%s\n', outRoot);

end

%% =====================================================================
%% ============================ HELPERS =================================
%% =====================================================================

function f = resolve_existing_file(candidates, labelName)
f = '';
for i = 1:numel(candidates)
    if isfile(candidates{i})
        f = candidates{i};
        return;
    end
end
error('Could not find required file: %s', labelName);
end

function gOut = standardize_group_labels_psp_only(gIn)
% Standardise feature-table group labels for PSP-only analysis.
% PC is intentionally kept separate from PSP and is not included in PSP contrasts.

if isstring(gIn)
    gIn = cellstr(gIn);
elseif ischar(gIn)
    gIn = cellstr(gIn);
end

gOut = cell(size(gIn));

for i = 1:numel(gIn)
    raw = strtrim(char(gIn{i}));
    s = lower(raw);
    s = strrep(s, '+', '');
    s = strrep(s, '/', '');
    s = strrep(s, '-', '');
    s = strrep(s, '_', '');
    s = regexprep(s, '\s+', '');

    if any(strcmp(s, {'con','control','controls','healthycontrol','healthycontrols'}))
        gOut{i} = 'CON';
    elseif strcmp(s, 'bvftd')
        gOut{i} = 'bvFTD';
    elseif strcmp(s, 'psp')
        gOut{i} = 'PSP';
    elseif strcmp(s, 'pc')
        gOut{i} = 'PC';
    elseif any(strcmp(s, {'psppc','pcpsp','psppccombined','psppcgroup'}))
        gOut{i} = 'PSPPC_EXCLUDED';
    else
        gOut{i} = raw;
    end
end

gOut = gOut(:);
end

function [group1List, group2List, group1Label, group2Label] = resolve_comparison_groups_psp_only(sheetName, group1FromStats, group2FromStats)
% Convert comparison names / stats labels into actual feature-table labels.
% This keeps PSP separate from PC. Pooled Patient = bvFTD + PSP only.

switch lower(strrep(sheetName, '/', '_'))
    case 'con_vs_bvftd'
        group1List = {'CON'};
        group2List = {'bvFTD'};
        group1Label = 'CON';
        group2Label = 'bvFTD';

    case 'con_vs_psp'
        group1List = {'CON'};
        group2List = {'PSP'};
        group1Label = 'CON';
        group2Label = 'PSP';

    case 'psp_vs_bvftd'
        group1List = {'PSP'};
        group2List = {'bvFTD'};
        group1Label = 'PSP';
        group2Label = 'bvFTD';

    case {'con_vs_patient','con_vs_pooledpatient','con_vs_pat'}
        group1List = {'CON'};
        group2List = {'bvFTD','PSP'};
        group1Label = 'CON';
        group2Label = 'Patient_bvFTD_PSP';

    otherwise
        group1List = {standardize_single_group_label_psp_only(group1FromStats)};
        group2List = {standardize_single_group_label_psp_only(group2FromStats)};
        group1Label = strjoin(group1List, '+');
        group2Label = strjoin(group2List, '+');
end
end

function g = standardize_single_group_label_psp_only(g)
raw = strtrim(char(g));
s = lower(raw);
s = strrep(s, '+', '');
s = strrep(s, '/', '');
s = strrep(s, '-', '');
s = strrep(s, '_', '');
s = regexprep(s, '\s+', '');

if any(strcmp(s, {'con','control','controls','healthycontrol','healthycontrols'}))
    g = 'CON';
elseif strcmp(s, 'bvftd')
    g = 'bvFTD';
elseif strcmp(s, 'psp')
    g = 'PSP';
elseif strcmp(s, 'pc')
    g = 'PC';
elseif any(strcmp(s, {'patient','patients','pooledpatient','noncon'}))
    g = 'Patient_bvFTD_PSP';
elseif any(strcmp(s, {'psppc','pcpsp','psppccombined','psppcgroup'}))
    g = 'PSPPC_EXCLUDED';
else
    g = raw;
end
end

function T = clean_stats_table_feature_report(T)

T = standardise_variable_names(T);

requiredVars = {'Feature','GroupA','GroupB','p_Welch','Cohen_d'};
for i = 1:numel(requiredVars)
    if ~ismember(requiredVars{i}, T.Properties.VariableNames)
        warning('Stats table missing variable: %s', requiredVars{i});
        T = table();
        return;
    end
end

if isstring(T.Feature), T.Feature = cellstr(T.Feature); end
if ischar(T.Feature),   T.Feature = cellstr(T.Feature); end

if isstring(T.GroupA), T.GroupA = cellstr(T.GroupA); end
if ischar(T.GroupA),   T.GroupA = cellstr(T.GroupA); end

if isstring(T.GroupB), T.GroupB = cellstr(T.GroupB); end
if ischar(T.GroupB),   T.GroupB = cellstr(T.GroupB); end

bad = ismissing(T.Feature) | isnan(T.p_Welch) | isnan(T.Cohen_d);
T = T(~bad,:);

T = T(T.p_Welch > 0,:);
end

function T = standardise_variable_names(T)
vn = T.Properties.VariableNames;

mapFrom = {
    {'Feature','FeatureName','Parameter','Variable'}
    {'GroupA','Group1','ComparisonA'}
    {'GroupB','Group2','ComparisonB'}
    {'p_Welch','p','pValue','p_value','p_uncorrected','Welch_p'}
    {'Cohen_d','Cohens_d','EffectSize','d','cohen_d'}
    };

mapTo = {'Feature','GroupA','GroupB','p_Welch','Cohen_d'};

for k = 1:numel(mapTo)
    target = mapTo{k};
    candidates = mapFrom{k};
    idxFound = find(ismember(lower(vn), lower(candidates)), 1, 'first');
    if ~isempty(idxFound)
        vn{idxFound} = target;
    end
end

T.Properties.VariableNames = matlab.lang.makeUniqueStrings(vn);
end

function recurrenceBonus = build_recurrence_bonus(TopHits, featureVars)

recurrenceBonus = zeros(numel(featureVars),1);

if isempty(TopHits)
    return;
end

TopHits = standardise_tophits_names(TopHits);

if ~ismember('Feature', TopHits.Properties.VariableNames), return; end
if ~ismember('Comparison', TopHits.Properties.VariableNames), return; end
if ~ismember('RankScore', TopHits.Properties.VariableNames), return; end

if isstring(TopHits.Feature), TopHits.Feature = cellstr(TopHits.Feature); end
if ischar(TopHits.Feature),   TopHits.Feature = cellstr(TopHits.Feature); end

if isstring(TopHits.Comparison), TopHits.Comparison = cellstr(TopHits.Comparison); end
if ischar(TopHits.Comparison),   TopHits.Comparison = cellstr(TopHits.Comparison); end

comparisons = unique(TopHits.Comparison, 'stable');
recurrence = zeros(numel(featureVars),1);

for c = 1:numel(comparisons)
    idxC = strcmp(TopHits.Comparison, comparisons{c});
    Tc = TopHits(idxC,:);

    [tf, idxF] = ismember(Tc.Feature, featureVars);
    Tc = Tc(tf,:);
    idxF = idxF(tf);

    if isempty(Tc), continue; end

    [~, ord] = sort(Tc.RankScore, 'descend', 'MissingPlacement', 'last');
    ord = ord(1:min(20, numel(ord)));
    recurrence(idxF(ord)) = recurrence(idxF(ord)) + 1;
end

if max(recurrence) > 0
    recurrenceBonus = recurrence / max(recurrence);
end
end

function T = standardise_tophits_names(T)
vn = T.Properties.VariableNames;

if any(strcmpi(vn,'FeatureName')) && ~any(strcmp(vn,'Feature'))
    vn{find(strcmpi(vn,'FeatureName'),1)} = 'Feature';
end
if any(strcmpi(vn,'Comp')) && ~any(strcmp(vn,'Comparison'))
    vn{find(strcmpi(vn,'Comp'),1)} = 'Comparison';
end
if any(strcmpi(vn,'Score')) && ~any(strcmp(vn,'RankScore'))
    vn{find(strcmpi(vn,'Score'),1)} = 'RankScore';
end

T.Properties.VariableNames = matlab.lang.makeUniqueStrings(vn);
end

function z = safe_zscore(x)
x = x(:);
mu = mean(x, 'omitnan');
sd = std(x, 0, 'omitnan');
if isnan(sd) || sd == 0
    z = zeros(size(x));
else
    z = (x - mu) ./ sd;
end
end

function featIdx = rank_select(idxRef, score, K, direction)

score = score(:);
idxRef = idxRef(:);

valid = ~isnan(score) & ~isnan(idxRef);
score = score(valid);
idxRef = idxRef(valid);

if isempty(score)
    featIdx = [];
    return;
end

switch lower(direction)
    case 'descend'
        [~, ord] = sort(score, 'descend', 'MissingPlacement', 'last');
    otherwise
        [~, ord] = sort(score, 'ascend', 'MissingPlacement', 'last');
end

ord = ord(1:min(K, numel(ord)));
featIdx = idxRef(ord);
end

function [metricsTable, predTable, rocStruct] = run_models_loo_feature_table( ...
    X, Y, yText, subjIDs, panelName, compName, group1, group2, rngSeed, nTreesRF, nTreesAda)

rng(rngSeed);

cvp = cvpartition(Y, 'LeaveOut');
nFolds = cvp.NumTestSets;

modelNames = {'LinearSVM','RBFSVM','RandomForest','AdaBoost'};

metricsTable = table();
predTable = table();
rocStruct = struct();

for m = 1:numel(modelNames)
    modelName = modelNames{m};

    yTrue  = nan(numel(Y),1);
    yPred  = nan(numel(Y),1);
    yScore = nan(numel(Y),1);

    for fold = 1:nFolds
        trIdx = training(cvp, fold);
        teIdx = test(cvp, fold);

        Xtr = X(trIdx,:);
        Xte = X(teIdx,:);
        Ytr = Y(trIdx);
        Yte = Y(teIdx);

        if numel(unique(Ytr)) < 2
            yTrue(teIdx)  = Yte;
            yPred(teIdx)  = mode(Ytr);
            yScore(teIdx) = double(mode(Ytr));
            continue;
        end

        mu = mean(Xtr, 1);
        sd = std(Xtr, 0, 1);
        sd(sd==0) = 1;

        Xtrz = (Xtr - mu) ./ sd;
        Xtez = (Xte - mu) ./ sd;

        switch modelName

            case 'LinearSVM'
                mdl = fitcsvm(Xtrz, Ytr, ...
                    'KernelFunction', 'linear', ...
                    'Standardize', false, ...
                    'ClassNames', [0 1], ...
                    'BoxConstraint', 1);
                mdl = fitPosterior(mdl);
                [pred, score] = predict(mdl, Xtez);
                posScore = get_positive_score(score);

            case 'RBFSVM'
                mdl = fitcsvm(Xtrz, Ytr, ...
                    'KernelFunction', 'rbf', ...
                    'KernelScale', 'auto', ...
                    'Standardize', false, ...
                    'ClassNames', [0 1], ...
                    'BoxConstraint', 1);
                mdl = fitPosterior(mdl);
                [pred, score] = predict(mdl, Xtez);
                posScore = get_positive_score(score);

            case 'RandomForest'
                mdl = TreeBagger(nTreesRF, Xtr, Ytr, ...
                    'Method', 'classification', ...
                    'OOBPrediction', 'off', ...
                    'MinLeafSize', 1, ...
                    'NumPredictorsToSample', 'all');

                [predCell, score] = predict(mdl, Xte);
                pred = str2double(predCell);

                classNamesRF = mdl.ClassNames;
                if iscell(classNamesRF)
                    posCol = find(strcmp(classNamesRF, '1'), 1);
                    if isempty(posCol), posCol = size(score,2); end
                else
                    posCol = find(classNamesRF == 1, 1);
                    if isempty(posCol), posCol = size(score,2); end
                end
                posScore = score(:, posCol);

            case 'AdaBoost'
                t = templateTree('MaxNumSplits', 20);
                mdl = fitcensemble(Xtr, Ytr, ...
                    'Method', 'AdaBoostM1', ...
                    'Learners', t, ...
                    'NumLearningCycles', nTreesAda);

                [pred, score] = predict(mdl, Xte);
                posScore = get_positive_score(score);
        end

        yTrue(teIdx)  = Yte;
        yPred(teIdx)  = pred;
        yScore(teIdx) = posScore;
    end

    valid = ~isnan(yTrue) & ~isnan(yPred) & ~isnan(yScore);
    yTrue_v  = yTrue(valid);
    yPred_v  = yPred(valid);
    yScore_v = yScore(valid);

    M = compute_binary_metrics(yTrue_v, yPred_v, yScore_v);
    M.Model = {modelName};
    M.NumEvaluated = numel(yTrue_v);

    metricsTable = [metricsTable; struct2table(M)];

    Tpred = table;
    Tpred.SubjectID          = subjIDs(:);
    Tpred.GroupText          = yText(:);
    Tpred.TrueLabel          = yTrue;
    Tpred.PredictedLabel     = yPred;
    Tpred.Score_Positive     = yScore;
    Tpred.CorrectPrediction  = (yTrue == yPred);
    Tpred.Model              = repmat({modelName}, numel(Y), 1);
    Tpred.Comparison         = repmat({compName}, numel(Y), 1);
    Tpred.Panel              = repmat({panelName}, numel(Y), 1);
    Tpred.Group1Label        = repmat({group1}, numel(Y), 1);
    Tpred.Group2Label        = repmat({group2}, numel(Y), 1);

    predTable = [predTable; Tpred];

    [fpr, tpr, ~, auc] = perfcurve(yTrue_v, yScore_v, 1);
    rocStruct.(modelName).fpr = fpr;
    rocStruct.(modelName).tpr = tpr;
    rocStruct.(modelName).auc = auc;
end
end

function posScore = get_positive_score(score)
if size(score,2) >= 2
    posScore = score(:,2);
else
    posScore = score(:,1);
end
end

function M = compute_binary_metrics(yTrue, yPred, yScore)

yTrue = yTrue(:);
yPred = yPred(:);
yScore = yScore(:);

TP = sum((yTrue==1) & (yPred==1));
TN = sum((yTrue==0) & (yPred==0));
FP = sum((yTrue==0) & (yPred==1));
FN = sum((yTrue==1) & (yPred==0));

acc = (TP + TN) / numel(yTrue);

if (TP+FN) > 0
    sens = TP / (TP + FN);
else
    sens = NaN;
end

if (TN+FP) > 0
    spec = TN / (TN + FP);
else
    spec = NaN;
end

if (TP+FP) > 0
    prec = TP / (TP + FP);
else
    prec = NaN;
end

if (prec + sens) > 0
    f1 = 2 * prec * sens / (prec + sens);
else
    f1 = NaN;
end

if (TP+FP) > 0 && (TP+FN) > 0 && (TN+FP) > 0 && (TN+FN) > 0
    mcc = ((TP*TN) - (FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
else
    mcc = NaN;
end

try
    [~,~,~,auc] = perfcurve(yTrue, yScore, 1);
catch
    auc = NaN;
end

M = struct();
M.Accuracy    = acc;
M.Sensitivity = sens;
M.Specificity = spec;
M.Precision   = prec;
M.F1score     = f1;
M.MCC         = mcc;
M.AUC         = auc;
M.TP          = TP;
M.TN          = TN;
M.FP          = FP;
M.FN          = FN;
end

function save_roc_and_confusion_figures_feature_table(rocStruct, predTable, figDir, panelName, compName)

modelNames = fieldnames(rocStruct);

for i = 1:numel(modelNames)
    modelName = modelNames{i};

    h = figure('Color','w','Position',[100 100 700 550]);
    plot(rocStruct.(modelName).fpr, rocStruct.(modelName).tpr, 'LineWidth', 2); hold on;
    plot([0 1], [0 1], '--k', 'LineWidth', 1.2);
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title(sprintf('%s | %s | AUC = %.3f', compName, modelName, rocStruct.(modelName).auc), ...
        'Interpreter','none');
    grid off; axis square;
    saveas(h, fullfile(figDir, sprintf('ROC_%s.png', modelName)));
    close(h);

    idx = strcmp(predTable.Model, modelName);
    yTrue = predTable.TrueLabel(idx);
    yPred = predTable.PredictedLabel(idx);

    valid = ~isnan(yTrue) & ~isnan(yPred);
    yTrue = yTrue(valid);
    yPred = yPred(valid);

    h2 = figure('Color','w','Position',[100 100 600 500]);
    cm = confusionchart(yTrue, yPred, ...
        'RowSummary','row-normalized', ...
        'ColumnSummary','column-normalized');
    cm.Title = sprintf('%s | %s | %s', compName, panelName, modelName);
    saveas(h2, fullfile(figDir, sprintf('Confusion_%s.png', modelName)));
    close(h2);
end

h3 = figure('Color','w','Position',[100 100 750 600]); hold on;
for i = 1:numel(modelNames)
    modelName = modelNames{i};
    plot(rocStruct.(modelName).fpr, rocStruct.(modelName).tpr, ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('%s (AUC=%.3f)', modelName, rocStruct.(modelName).auc));
end
plot([0 1],[0 1],'--k','HandleVisibility','off');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(sprintf('%s | %s | ROC comparison', compName, panelName), 'Interpreter', 'none');
legend('Location','southeast');
grid off; axis square;
saveas(h3, fullfile(figDir, 'ROC_AllModels.png'));
close(h3);

end

function s = safe_sheet_name(s)
s = regexprep(s, '[:\\/?*\[\]]', '_');
if numel(s) > 31
    s = s(1:31);
end
end
