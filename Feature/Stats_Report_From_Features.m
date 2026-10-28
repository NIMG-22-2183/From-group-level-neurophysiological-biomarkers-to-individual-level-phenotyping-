function Stats_Report_From_Features()
% =========================================================================
% Stats_Report_From_Features
%
% Reads the saved optimized feature table and performs comprehensive
% between-group statistical testing for:
%   1) CON vs bvFTD
%   2) CON vs PSP
%   3) PSP vs bvFTD
%   4) CON vs Patient (bvFTD + PSP)
%
% Outputs:
%   - One Excel workbook with all comparison reports (one sheet per contrast)
%   - One summary sheet with top hits across all contrasts
%   - One subfolder per comparison containing plots:
%       * Heatmap of Hedges g
%       * Heatmap of -log10(q)
%       * Heatmap of signed score = sign(g)*(-log10(q))
%       * Heatmap of significant cells
%       * Volcano plot
%       * Top-N bar plot
%
% Statistical tests / metrics:
%   - Welch's unequal-variance two-sample t-test (two-sided)
%   - Cohen's d (pooled SD)
%   - Hedges' g (bias-corrected)
%   - Benjamini-Hochberg FDR q-values
%
% Extra outputs:
%   - Means, SDs, Ns for each group
%   - Mean difference
%   - Welch t-statistic and Welch df
%   - RankScore = abs(Hedges_g) * max(0,-log10(q))
%
% =========================================================================

clc; close all; clearvars;

% Local folders: edit these paths to match your data.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(repoRoot, 'data');
resultsRoot = fullfile(repoRoot, 'results');


%% ============================ USER SETTINGS ============================
featureFolder = fullfile(resultsRoot, 'features');
featureExcel  = fullfile(featureFolder, 'Features_AllSubjects_Optimized.xlsx');
featureSheet  = 'Features';

outBase = fullfile(resultsRoot, 'feature_statistics');
if ~exist(outBase, 'dir'), mkdir(outBase); end

excelOut = fullfile(outBase, 'Feature_Statistics_Report_AllComparisons.xlsx');
matOut   = fullfile(outBase, 'Feature_Statistics_Report_AllComparisons.mat');

topNBar        = 30;
qSigThreshold  = 0.05;
minNonNaNEach  = 3;
saveFigFiles   = true;
makeHeatmaps   = true;
makeVolcano    = true;
makeTopBar     = true;

%% ============================== LOAD DATA ==============================
assert(isfile(featureExcel), 'Feature Excel not found:\n%s', featureExcel);

T = readtable(featureExcel, 'Sheet', featureSheet, 'PreserveVariableNames', true);
assert(ismember('SubjectID', T.Properties.VariableNames), 'SubjectID column not found.');
assert(ismember('Group',     T.Properties.VariableNames), 'Group column not found.');

T.Group = standardize_group_labels(T.Group);

allVars = T.Properties.VariableNames;
metaVars = {'SubjectID','Group'};
featVars = setdiff(allVars, metaVars, 'stable');

isNum = false(size(featVars));
for i = 1:numel(featVars)
    isNum(i) = isnumeric(T.(featVars{i})) || islogical(T.(featVars{i}));
end
featVars = featVars(isNum);
assert(~isempty(featVars), 'No numeric feature columns found.');

fprintf('Loaded %d subjects and %d numeric features.\n', height(T), numel(featVars));
disp(groupcounts(categorical(T.Group)));

%% =========================== DEFINE CONTRASTS ==========================
contrasts = struct([]);

contrasts(1).name   = 'CON_vs_bvFTD';
contrasts(1).groupA = 'CON';
contrasts(1).groupB = 'bvFTD';
contrasts(1).idxA   = strcmpi(T.Group, 'CON');
contrasts(1).idxB   = strcmpi(T.Group, 'bvFTD');

contrasts(2).name   = 'CON_vs_PSP';
contrasts(2).groupA = 'CON';
contrasts(2).groupB = 'PSP';
contrasts(2).idxA   = strcmpi(T.Group, 'CON');
contrasts(2).idxB   = strcmpi(T.Group, 'PSP');

contrasts(3).name   = 'PSP_vs_bvFTD';
contrasts(3).groupA = 'PSP';
contrasts(3).groupB = 'bvFTD';
contrasts(3).idxA   = strcmpi(T.Group, 'PSP');
contrasts(3).idxB   = strcmpi(T.Group, 'bvFTD');

contrasts(4).name   = 'CON_vs_Patient';
contrasts(4).groupA = 'CON';
contrasts(4).groupB = 'Patient';
contrasts(4).idxA   = strcmpi(T.Group, 'CON');
contrasts(4).idxB   = strcmpi(T.Group, 'PSP') | strcmpi(T.Group, 'bvFTD');

%% ========================== RUN COMPARISONS ============================
allResults = cell(numel(contrasts),1);
allTop     = table();

for ci = 1:numel(contrasts)
    C = contrasts(ci);

    fprintf('\n=============================================================\n');
    fprintf('Running comparison: %s\n', C.name);
    fprintf('Group A = %s (n=%d)\n', C.groupA, sum(C.idxA));
    fprintf('Group B = %s (n=%d)\n', C.groupB, sum(C.idxB));
    fprintf('=============================================================\n');

    compFolder = fullfile(outBase, C.name);
    if ~exist(compFolder, 'dir'), mkdir(compFolder); end

    Res = run_one_comparison(T, featVars, C.groupA, C.groupB, C.idxA, C.idxB, minNonNaNEach);

    Res = sortrows(Res, {'q_BH','Abs_Hedges_g','Abs_Cohen_d'}, {'ascend','descend','descend'});
    Res.Rank = (1:height(Res))';

    if makeHeatmaps
        make_comparison_heatmaps(Res, C.name, compFolder, qSigThreshold, saveFigFiles);
    end
    if makeVolcano
        make_volcano_plot(Res, C.name, compFolder, qSigThreshold, saveFigFiles);
    end
    if makeTopBar
        make_top_features_barplot(Res, C.name, compFolder, topNBar, saveFigFiles);
    end

    allResults{ci} = Res;

    tmp = Res;
    tmp.Comparison = repmat(string(C.name), height(tmp), 1);
    allTop = [allTop; tmp]; %#ok<AGROW>
end

%% ============================== SAVE EXCEL =============================
README = build_readme_table(T, featVars, contrasts, qSigThreshold, minNonNaNEach);
safe_write_sheet(excelOut, README, 'README');

for ci = 1:numel(contrasts)
    safe_write_sheet(excelOut, allResults{ci}, contrasts(ci).name);
end

allTop = sortrows(allTop, {'q_BH','Abs_Hedges_g','Abs_Cohen_d'}, {'ascend','descend','descend'});
safe_write_sheet(excelOut, allTop, 'TopHits_AllComparisons');

save(matOut, 'T', 'featVars', 'contrasts', 'allResults', 'allTop', ...
    'qSigThreshold', 'minNonNaNEach', 'featureExcel', 'excelOut');

fprintf('\nDone.\n');
fprintf('Excel report:\n  %s\n', excelOut);
fprintf('MAT summary:\n  %s\n', matOut);
fprintf('Plots saved under:\n  %s\n', outBase);

end % main


%% =======================================================================
function Res = run_one_comparison(T, featVars, groupAName, groupBName, idxA, idxB, minNonNaNEach)

nFeat = numel(featVars);

Feature          = strings(nFeat,1);
FeatureFamily    = strings(nFeat,1);
Channel          = strings(nFeat,1);
GroupA           = strings(nFeat,1);
GroupB           = strings(nFeat,1);
nA               = nan(nFeat,1);
nB               = nan(nFeat,1);
meanA            = nan(nFeat,1);
meanB            = nan(nFeat,1);
sdA              = nan(nFeat,1);
sdB              = nan(nFeat,1);
meanDiff_AminusB = nan(nFeat,1);
tWelch           = nan(nFeat,1);
dfWelch          = nan(nFeat,1);
pWelch           = nan(nFeat,1);
cohenD           = nan(nFeat,1);
hedgesG          = nan(nFeat,1);
absD             = nan(nFeat,1);
absG             = nan(nFeat,1);
neglog10p        = nan(nFeat,1);
rankScore        = nan(nFeat,1);

for j = 1:nFeat
    vname = featVars{j};
    Feature(j) = string(vname);
    [fam, ch] = split_feature_and_channel(vname);
    FeatureFamily(j) = string(fam);
    Channel(j) = string(ch);
    GroupA(j) = string(groupAName);
    GroupB(j) = string(groupBName);

    xA = T{idxA, vname};
    xB = T{idxB, vname};

    xA = double(xA(:));
    xB = double(xB(:));

    xA = xA(isfinite(xA));
    xB = xB(isfinite(xB));

    nA(j) = numel(xA);
    nB(j) = numel(xB);

    if numel(xA) < minNonNaNEach || numel(xB) < minNonNaNEach
        continue;
    end

    meanA(j) = mean(xA);
    meanB(j) = mean(xB);
    sdA(j)   = std(xA, 0);
    sdB(j)   = std(xB, 0);
    meanDiff_AminusB(j) = meanA(j) - meanB(j);

    [tWelch(j), dfWelch(j), pWelch(j)] = welch_ttest_manual(xA, xB);

    cohenD(j) = cohens_d_pooled(xA, xB);
    hedgesG(j) = hedges_g_from_d(cohenD(j), numel(xA), numel(xB));

    absD(j) = abs(cohenD(j));
    absG(j) = abs(hedgesG(j));

    if isfinite(pWelch(j)) && pWelch(j) > 0
        neglog10p(j) = -log10(pWelch(j));
    end
end

qBH = bh_fdr(pWelch);

for j = 1:nFeat
    if isfinite(qBH(j)) && qBH(j) > 0 && isfinite(absG(j))
        rankScore(j) = absG(j) * max(0, -log10(qBH(j)));
    else
        rankScore(j) = NaN;
    end
end

Res = table( ...
    Feature, FeatureFamily, Channel, GroupA, GroupB, ...
    nA, nB, meanA, meanB, sdA, sdB, meanDiff_AminusB, ...
    tWelch, dfWelch, pWelch, qBH, ...
    cohenD, absD, hedgesG, absG, neglog10p, rankScore, ...
    'VariableNames', { ...
    'Feature','FeatureFamily','Channel','GroupA','GroupB', ...
    'N_A','N_B','Mean_A','Mean_B','SD_A','SD_B','MeanDiff_AminusB', ...
    'Welch_t','Welch_df','p_Welch','q_BH', ...
    'Cohen_d','Abs_Cohen_d','Hedges_g','Abs_Hedges_g','NegLog10_p','RankScore'});

end


%% =======================================================================
function g = standardize_group_labels(g)

if iscell(g)
    g = string(g);
elseif iscategorical(g)
    g = string(g);
elseif ischar(g)
    g = string(cellstr(g));
elseif ~isstring(g)
    g = string(g);
end

g = strtrim(g);

for i = 1:numel(g)
    gi = lower(strrep(strrep(g(i), ' ', ''), '_', ''));

    if any(strcmpi(gi, ["con","control","controls"]))
        g(i) = "CON";
    elseif any(strcmpi(gi, ["bvftd","bv-ftd","bv_ftd","b v f t d"]))
        g(i) = "bvFTD";
    elseif any(strcmpi(gi, ["psp","progressivesupranuclearpalsy"]))
        g(i) = "PSP";
    else
        g(i) = string(g(i));
    end
end

end


%% =======================================================================
function [fam, ch] = split_feature_and_channel(vname)

parts = split(string(vname), "__");
if numel(parts) >= 2
    fam = parts(1);
    ch  = join(parts(2:end), "__");
else
    fam = string(vname);
    ch  = "";
end

fam = char(fam);
ch  = char(ch);

end


%% =======================================================================
function [tstat, df, p] = welch_ttest_manual(x1, x2)

x1 = x1(isfinite(x1));
x2 = x2(isfinite(x2));

n1 = numel(x1);
n2 = numel(x2);

if n1 < 2 || n2 < 2
    tstat = NaN; df = NaN; p = NaN;
    return;
end

m1 = mean(x1);
m2 = mean(x2);
s1 = var(x1, 0);
s2 = var(x2, 0);

se2 = s1/n1 + s2/n2;
if ~isfinite(se2) || se2 <= 0
    tstat = NaN; df = NaN; p = NaN;
    return;
end

tstat = (m1 - m2) / sqrt(se2);

num = se2^2;
den = (s1^2 / (n1^2 * (n1-1))) + (s2^2 / (n2^2 * (n2-1)));
if ~isfinite(den) || den <= 0
    df = NaN;
    p = NaN;
    return;
end

df = num / den;
p  = 2 * tcdf(-abs(tstat), df);

end


%% =======================================================================
function d = cohens_d_pooled(x1, x2)

x1 = x1(isfinite(x1));
x2 = x2(isfinite(x2));

n1 = numel(x1);
n2 = numel(x2);

if n1 < 2 || n2 < 2
    d = NaN;
    return;
end

m1 = mean(x1);
m2 = mean(x2);
s1 = std(x1, 0);
s2 = std(x2, 0);

sp2 = ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2);
if ~isfinite(sp2) || sp2 <= 0
    d = NaN;
    return;
end

sp = sqrt(sp2);
d  = (m1 - m2) / sp;

end


%% =======================================================================
function g = hedges_g_from_d(d, n1, n2)

if ~isfinite(d) || (n1 + n2 - 2) <= 0
    g = NaN;
    return;
end

df = n1 + n2 - 2;
J  = 1 - (3 / (4*df - 1));
g  = J * d;

end


%% =======================================================================
function q = bh_fdr(p)

p = p(:);
q = nan(size(p));

good = isfinite(p);
pg = p(good);

if isempty(pg)
    return;
end

[ps, ord] = sort(pg);
m = numel(ps);

qs = (m ./ (1:m)') .* ps;
qs = flipud(cummin(flipud(qs)));
qs(qs > 1) = 1;

tmp = nan(m,1);
tmp(ord) = qs;
q(good) = tmp;

end


%% =======================================================================
function make_comparison_heatmaps(Res, compName, outDir, qThr, saveFigFiles)

if isempty(Res)
    return;
end

use = Res.Channel ~= "";
R = Res(use,:);

if isempty(R)
    return;
end

families = unique(R.FeatureFamily, 'stable');
channels = unique(R.Channel, 'stable');

% ---- force natural channel order like Ch1, Ch2, ..., Ch7
channels = natural_channel_sort(channels);

Gmat = nan(numel(families), numel(channels));
Qmat = nan(numel(families), numel(channels));
Smat = nan(numel(families), numel(channels));
SigMat = nan(numel(families), numel(channels));

for i = 1:height(R)
    r = find(strcmp(families, R.FeatureFamily(i)), 1, 'first');
    c = find(strcmp(channels, R.Channel(i)), 1, 'first');
    if isempty(r) || isempty(c), continue; end

    Gmat(r,c) = R.Hedges_g(i);

    if isfinite(R.q_BH(i)) && R.q_BH(i) > 0
        qv = R.q_BH(i);
        qlog = -log10(qv);
        Qmat(r,c) = qlog;
        Smat(r,c) = sign(R.Hedges_g(i)) * qlog;
        if qv <= qThr
            SigMat(r,c) = 1;
        else
            SigMat(r,c) = 0;
        end
    end
end

% ---- Hedges g heatmap
f1 = figure('Color','w','Position',[100 100 1200 700],'Visible','off');
imagesc(Gmat);
axis tight;
set(gca,'XTick',1:numel(channels),'XTickLabel',channels, ...
        'YTick',1:numel(families),'YTickLabel',families, ...
        'TickLabelInterpreter','none','FontSize',10);
xtickangle(45);
colorbar;
title([compName ': Hedges g (GroupA - GroupB)'], 'Interpreter','none', 'FontWeight','bold');
xlabel('Channel');
ylabel('Feature family');
colormap(f1, parula);
saveas(f1, fullfile(outDir, 'Heatmap_HedgesG.png'));
if saveFigFiles, savefig(f1, fullfile(outDir, 'Heatmap_HedgesG.fig')); end
close(f1);

% ---- -log10(q) heatmap
f2 = figure('Color','w','Position',[100 100 1200 700],'Visible','off');
imagesc(Qmat);
axis tight;
set(gca,'XTick',1:numel(channels),'XTickLabel',channels, ...
        'YTick',1:numel(families),'YTickLabel',families, ...
        'TickLabelInterpreter','none','FontSize',10);
xtickangle(45);
colorbar;
title([compName ': -log10(q_BH)'], 'Interpreter','none', 'FontWeight','bold');
xlabel('Channel');
ylabel('Feature family');
colormap(f2, hot);
saveas(f2, fullfile(outDir, 'Heatmap_NegLog10Q.png'));
if saveFigFiles, savefig(f2, fullfile(outDir, 'Heatmap_NegLog10Q.fig')); end
close(f2);

% ---- signed score heatmap
f3 = figure('Color','w','Position',[100 100 1200 700],'Visible','off');
imagesc(Smat);
axis tight;
set(gca,'XTick',1:numel(channels),'XTickLabel',channels, ...
        'YTick',1:numel(families),'YTickLabel',families, ...
        'TickLabelInterpreter','none','FontSize',10);
xtickangle(45);
colorbar;
title([compName ': Signed score = sign(g) x -log10(q)'], 'Interpreter','none', 'FontWeight','bold');
xlabel('Channel');
ylabel('Feature family');
colormap(f3, turbo);
saveas(f3, fullfile(outDir, 'Heatmap_SignedScore.png'));
if saveFigFiles, savefig(f3, fullfile(outDir, 'Heatmap_SignedScore.fig')); end
close(f3);

% ---- significant cells heatmap
f4 = figure('Color','w','Position',[100 100 1200 700],'Visible','off');
imagesc(SigMat);
axis tight;
set(gca,'XTick',1:numel(channels),'XTickLabel',channels, ...
        'YTick',1:numel(families),'YTickLabel',families, ...
        'TickLabelInterpreter','none','FontSize',10);
xtickangle(45);
colorbar;
title(sprintf('%s: Significant cells (q <= %.3f)', compName, qThr), ...
    'Interpreter','none', 'FontWeight','bold');
xlabel('Channel');
ylabel('Feature family');
colormap(f4, gray);
caxis([0 1]);
saveas(f4, fullfile(outDir, 'Heatmap_SignificantCells.png'));
if saveFigFiles, savefig(f4, fullfile(outDir, 'Heatmap_SignificantCells.fig')); end
close(f4);

end


%% =======================================================================
function make_volcano_plot(Res, compName, outDir, qThr, saveFigFiles)

R = Res;
good = isfinite(R.Hedges_g) & isfinite(R.p_Welch) & R.p_Welch > 0;
R = R(good,:);

if isempty(R), return; end

x = R.Hedges_g;
y = -log10(R.p_Welch);
sig = isfinite(R.q_BH) & R.q_BH <= qThr;

f = figure('Color','w','Position',[100 100 1000 720],'Visible','off');
hold on;

scatter(x(~sig), y(~sig), 28, 'filled', 'MarkerFaceAlpha', 0.45);
scatter(x(sig),  y(sig),  40, 'filled', 'MarkerFaceAlpha', 0.90);

yl = yline(-log10(0.05), '--', 'p = 0.05');
xl1 = xline(-0.5, ':', 'g = -0.5');
xl2 = xline(0.5,  ':', 'g = 0.5');

yl.LabelHorizontalAlignment = 'right';
yl.LabelVerticalAlignment   = 'bottom';
xl1.LabelHorizontalAlignment = 'left';
xl2.LabelHorizontalAlignment = 'right';

yl.Color  = [0.3 0.3 0.3];
xl1.Color = [0.4 0.4 0.4];
xl2.Color = [0.4 0.4 0.4];

title([compName ': Volcano plot'], 'Interpreter','none', 'FontWeight','bold');
xlabel('Hedges g (GroupA - GroupB)');
ylabel('-log10(p)');
grid on;
box on;

% annotate top 10
R2 = sortrows(R, {'q_BH','Abs_Hedges_g'}, {'ascend','descend'});
nLab = min(10, height(R2));
for i = 1:nLab
    text(R2.Hedges_g(i), -log10(R2.p_Welch(i)), ['  ' char(R2.Feature(i))], ...
        'FontSize', 8, 'Interpreter', 'none', 'Clipping','on');
end

saveas(f, fullfile(outDir, 'Volcano_HedgesG_vs_p.png'));
if saveFigFiles, savefig(f, fullfile(outDir, 'Volcano_HedgesG_vs_p.fig')); end
close(f);

end


%% =======================================================================
function make_top_features_barplot(Res, compName, outDir, topN, saveFigFiles)

R = Res;
good = isfinite(R.RankScore);
R = R(good,:);

if isempty(R), return; end

R = sortrows(R, {'q_BH','Abs_Hedges_g','Abs_Cohen_d'}, {'ascend','descend','descend'});
nShow = min(topN, height(R));
R = R(1:nShow,:);

labels = cellstr(R.Feature);
vals   = R.Hedges_g;

f = figure('Color','w','Position',[100 100 1400 950],'Visible','off');
barh(vals, 'FaceColor', [0.12 0.47 0.71], 'EdgeColor', 'k', 'LineWidth', 0.6);
set(gca, 'YDir','reverse', 'YTick',1:nShow, 'YTickLabel',labels, ...
    'TickLabelInterpreter','none', 'FontSize',9);
xlabel('Hedges g (GroupA - GroupB)');
ylabel('Feature');
title(sprintf('%s: Top %d features', compName, nShow), 'Interpreter','none', 'FontWeight','bold');
grid on;
box on;
xline(0, 'k-', 'LineWidth', 0.75);

saveas(f, fullfile(outDir, sprintf('Top%d_Features_BarPlot.png', nShow)));
if saveFigFiles, savefig(f, fullfile(outDir, sprintf('Top%d_Features_BarPlot.fig', nShow))); end
close(f);

end


%% =======================================================================
function README = build_readme_table(T, featVars, contrasts, qSigThreshold, minNonNaNEach)

txt = strings(0,1);

txt(end+1,1) = "Feature statistics report";
txt(end+1,1) = " ";
txt(end+1,1) = "Input summary:";
txt(end+1,1) = "  Number of subjects = " + string(height(T));
txt(end+1,1) = "  Number of numeric features = " + string(numel(featVars));
txt(end+1,1) = " ";
txt(end+1,1) = "Groups present:";

ug = unique(string(T.Group), 'stable');
for i = 1:numel(ug)
    txt(end+1,1) = "  " + ug(i) + " : n = " + string(sum(strcmpi(T.Group, ug(i))));
end

txt(end+1,1) = " ";
txt(end+1,1) = "Contrasts analysed:";
for i = 1:numel(contrasts)
    txt(end+1,1) = "  " + string(contrasts(i).name);
end

txt(end+1,1) = " ";
txt(end+1,1) = "Statistics reported:";
txt(end+1,1) = "  - Welch unequal-variance two-sample t-test (two-sided)";
txt(end+1,1) = "  - Cohen d (pooled SD)";
txt(end+1,1) = "  - Hedges g (bias-corrected)";
txt(end+1,1) = "  - Benjamini-Hochberg FDR q-values";

txt(end+1,1) = " ";
txt(end+1,1) = "Feature direction:";
txt(end+1,1) = "  MeanDiff_AminusB, Cohen_d, and Hedges_g are GroupA - GroupB";

txt(end+1,1) = " ";
txt(end+1,1) = "Sorting within each sheet:";
txt(end+1,1) = "  q_BH ascending, then Abs_Hedges_g descending, then Abs_Cohen_d descending";

txt(end+1,1) = " ";
txt(end+1,1) = "Additional ranking score:";
txt(end+1,1) = "  RankScore = abs(Hedges_g) * max(0,-log10(q_BH))";

txt(end+1,1) = " ";
txt(end+1,1) = "Thresholds:";
txt(end+1,1) = "  q significance threshold = " + string(qSigThreshold);
txt(end+1,1) = "  minimum non-NaN per group per feature = " + string(minNonNaNEach);

README = table((1:numel(txt))', txt(:), 'VariableNames', {'Line','Description'});

end


%% =======================================================================
function safe_write_sheet(excelFile, Tbl, sheetName)

sheetName = char(sheetName);
if numel(sheetName) > 31
    sheetName = sheetName(1:31);
end

try
    writetable(Tbl, excelFile, 'Sheet', sheetName, 'WriteMode', 'overwritesheet');
catch
    writetable(Tbl, excelFile, 'Sheet', sheetName);
end

end


%% =======================================================================
function chOut = natural_channel_sort(chIn)
% Sort channel labels like Ch1, Ch2, ..., Ch10 naturally.

if isempty(chIn)
    chOut = chIn;
    return;
end

chIn = string(chIn(:));
numVals = nan(numel(chIn),1);

for i = 1:numel(chIn)
    tok = regexp(chIn(i), '(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        numVals(i) = str2double(tok{1});
    end
end

hasNum = isfinite(numVals);
if all(hasNum)
    [~, ord] = sort(numVals);
else
    [~, ord] = sort(lower(chIn));
end

chOut = chIn(ord);

end
