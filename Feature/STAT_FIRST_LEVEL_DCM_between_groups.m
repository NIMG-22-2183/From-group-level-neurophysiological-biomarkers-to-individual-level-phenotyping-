function STAT_FIRST_LEVEL_DCM_between_groups()
% STAT_FIRST_LEVEL_DCM_between_groups
%
% Between-group stats on estimated DCM parameters in DCM.Ep.
% For each subject DCM: Ep = spm_vec(DCM.Ep); names = spm_fieldindices(DCM.Ep, 1:np)
% Then compute Welch t-tests + effect sizes + multiple-comparison correction.
%
% Outputs: one Excel per comparison in results/dcm_statistics/
%   - Patient group = bvFTD + PSP.

clc; close all; clearvars;

% Local folders: edit these paths to match your data.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(repoRoot, 'data');
resultsRoot = fullfile(repoRoot, 'results');


%% ====================== USER SETTINGS ===========================
excelFolder = dataRoot;
excelFile = 'participants.xlsx';
excelPath   = fullfile(excelFolder, excelFile);

matFolder = fullfile(resultsRoot, 'dcm');
matPattern  = '*_DCM_DSALL.mat';

% Group definitions
groupA_labels = {'CON'};        % controls
groupB_labels = {'bvFTD'};      % bvFTD
groupC_labels = {'PSP'};        % PSP 
% Comparisons to run (you can add/remove lines)
comparisons = { ...
    'CON_vs_bvFTD',      groupA_labels, groupB_labels; ...
    'CON_vs_PSP',        groupA_labels, groupC_labels; ...
    'bvFTD_vs_PSP',      groupB_labels, groupC_labels; ...
    'CON_vs_PATIENTS',   groupA_labels, [groupB_labels, groupC_labels]; ...
    };

% Output folder in the CURRENT working directory
outDir = fullfile(resultsRoot, 'dcm_statistics');
if ~exist(outDir), mkdir(outDir); end

%% ====================== STEP 1: LOAD EXCEL ======================
if ~isfile(excelPath)
    error('Excel file not found: %s', excelPath);
end
data = readtable(excelPath);

IDs    = data{:,1};
Groups = data{:,2};

% Ensure cellstr
if isstring(IDs),    IDs    = cellstr(IDs);    end
if isnumeric(IDs),   IDs    = arrayfun(@(n) sprintf('Sub%04d', n), IDs, 'UniformOutput', false); end
if isstring(Groups), Groups = cellstr(Groups); end

% Canonicalise IDs to 'Sub####'
IDs = regexprep(IDs, '^sub-', '');
IDs = regexprep(IDs, '.*?(Sub0*\d+).*', '$1');

% Normalise group labels to char
Groups = cellfun(@(x) char(x), Groups, 'UniformOutput', false);

%% ====================== STEP 2: LIST DCM FILES ==================
allFiles = dir(fullfile(matFolder, matPattern));
matNames = {allFiles.name}';


tok  = regexp(matNames, '^(Sub0\d+)_DCM_DSALL\.mat$', 'tokens', 'once');
keep = ~cellfun('isempty', tok);
matNames = matNames(keep);
matIDs   = cellfun(@(c) c{1}, tok(keep), 'UniformOutput', false);

if isempty(matNames)
    error('No files matching pattern found in %s', matFolder);
end

%% ====================== STEP 3: MAP SUBJECT -> LABEL -> FILE =====
subj = struct('id', {}, 'label', {}, 'file', {});
k = 0;
for i = 1:numel(matIDs)
    thisID = matIDs{i};
    idx = find(strcmpi(IDs, thisID), 1, 'first');
    if isempty(idx)
        warning('No group found for ID %s (file %s). Skipping.', thisID, matNames{i});
        continue;
    end

    k = k + 1;
    subj(k).id    = thisID;
    subj(k).label = Groups{idx};
    subj(k).file  = fullfile(matFolder, matNames{i});
end

if isempty(subj)
    error('No subjects could be matched between Excel and MAT files.');
end

%% ====================== STEP 4: EXTRACT Ep + NAMES FOR ALL ==============
refNames = [];
np_ref   = [];
refOK    = false;

% First pass: determine reference parameterisation
for i = 1:numel(subj)
    if ~isfile(subj(i).file), continue; end
    S = load(subj(i).file);
    if ~isfield(S, 'DCM') || ~isfield(S.DCM, 'Ep')
        warning('No DCM.Ep in %s. Skipping.', subj(i).file);
        continue;
    end
    DCM = S.DCM;

    Ep_vec = spm_vec(DCM.Ep);
    np     = numel(Ep_vec);

    try
        nm = spm_fieldindices(DCM.Ep, 1:np);
    catch
        %  at least store numeric index names
        nm = arrayfun(@(j) sprintf('p_%04d', j), 1:np, 'UniformOutput', false);
        warning('spm_fieldindices failed for %s. Using fallback names p_####.', subj(i).id);
    end

    % Ensure names is a cellstr
    if isstring(nm), nm = cellstr(nm); end
    if ~iscell(nm),  nm = cellstr(string(nm)); end

    refNames = nm(:);
    np_ref   = np;
    refOK    = true;
    break;
end

if ~refOK
    error('Could not load any valid DCM.Ep to establish reference names.');
end

% Second pass: build matrix X (subjects x params) with NaNs for mismatches
Nsub = numel(subj);
X    = nan(Nsub, np_ref);

for i = 1:Nsub
    if ~isfile(subj(i).file)
        warning('Missing file: %s', subj(i).file);
        continue;
    end
    S = load(subj(i).file);
    if ~isfield(S, 'DCM') || ~isfield(S.DCM, 'Ep')
        warning('No DCM.Ep in %s', subj(i).file);
        continue;
    end
    DCM = S.DCM;

    Ep_vec = spm_vec(DCM.Ep);
    np     = numel(Ep_vec);

    if np ~= np_ref
        % If parameter vector length differs, we keep overlapping indices only
        warning('Param length mismatch for %s: np=%d, ref=%d. Trunc/pad with NaN.', subj(i).id, np, np_ref);
        m = min(np, np_ref);
        X(i,1:m) = Ep_vec(1:m);
    else
        X(i,:) = Ep_vec(:)';
    end
end

% Remove subjects with all-NaN
allNaN = all(isnan(X), 2);
if any(allNaN)
    warning('Removing %d subjects with no parameters extracted.', sum(allNaN));
    subj = subj(~allNaN);
    X    = X(~allNaN, :);
    Nsub = numel(subj);
end

%% ====================== STEP 5: RUN COMPARISONS ==================
allLabels = {subj.label};

for c = 1:size(comparisons,1)
    compName   = comparisons{c,1};
    labels_1   = comparisons{c,2};
    labels_2   = comparisons{c,3};

    idx1 = ismember(allLabels, labels_1);
    idx2 = ismember(allLabels, labels_2);

    X1 = X(idx1,:);
    X2 = X(idx2,:);

    n1 = sum(idx1); n2 = sum(idx2);
    fprintf('\n%s: n1=%d vs n2=%d\n', compName, n1, n2);

    if n1 < 2 || n2 < 2
        warning('%s skipped: need at least 2 subjects per group.', compName);
        continue;
    end

    % Preallocate
    p  = nan(np_ref,1);
    t  = nan(np_ref,1);
    df = nan(np_ref,1);

    mu1 = nan(np_ref,1); sd1 = nan(np_ref,1);
    mu2 = nan(np_ref,1); sd2 = nan(np_ref,1);

    d_cohen  = nan(np_ref,1);
    g_hedges = nan(np_ref,1);

    % Loop each parameter
    for j = 1:np_ref
        a = X1(:,j); b = X2(:,j);
        a = a(~isnan(a)); b = b(~isnan(b));

        if numel(a) < 2 || numel(b) < 2
            continue;
        end

        mu1(j) = mean(a); sd1(j) = std(a, 0);
        mu2(j) = mean(b); sd2(j) = std(b, 0);

        % Welch t-test
        [t(j), df(j), p(j)] = welch_ttest(a, b);

        % Effect sizes
        [d_cohen(j), g_hedges(j)] = effect_sizes(a, b);
    end

    % Multiple comparison corrections
    p_fdr  = fdr_bh(p);
    p_bonf = min(p * np_ref, 1);

    % Assemble output table
    T = table;
    T.ParamIndex     = (1:np_ref)';
    T.ParamName      = refNames;
    T.N1             = repmat(n1, np_ref, 1);
    T.N2             = repmat(n2, np_ref, 1);
    T.Mean1          = mu1;
    T.SD1            = sd1;
    T.Mean2          = mu2;
    T.SD2            = sd2;
    T.t_Welch        = t;
    T.df_Welch       = df;
    T.p_uncorrected  = p;
    T.p_FDR_BH       = p_fdr;
    T.p_Bonferroni   = p_bonf;
    T.Cohens_d       = d_cohen;
    T.Hedges_g       = g_hedges;

    % Sort by uncorrected p (smallest first)
    [~, ord] = sort(T.p_uncorrected, 'ascend', 'MissingPlacement','last');
    T = T(ord,:);

    % Write excel
    outXlsx = fullfile(outDir, [compName '.xlsx']);
    writetable(T, outXlsx, 'Sheet', 'Stats');

    % Add a short summary sheet
    summary = table;
    summary.Comparison = {compName};
    summary.Group1_Labels = {strjoin(labels_1, ',')};
    summary.Group2_Labels = {strjoin(labels_2, ',')};
    summary.N1 = n1;
    summary.N2 = n2;
    summary.NumParams = np_ref;
    summary.NumValidTests = sum(~isnan(p));
    summary.NumSig_FDR_0p05 = sum(p_fdr < 0.05, 'omitnan');
    summary.NumSig_Bonf_0p05 = sum(p_bonf < 0.05, 'omitnan');
    writetable(summary, outXlsx, 'Sheet', 'Summary');

    fprintf('Saved: %s\n', outXlsx);
end

fprintf('\nAll done. Outputs in: %s\n', outDir);

end

%% ====================== HELPERS ======================

function [t, df, p] = welch_ttest(a, b)
% Welch's two-sample t-test (unequal variances), two-sided
na = numel(a); nb = numel(b);
ma = mean(a);  mb = mean(b);
va = var(a, 0); vb = var(b, 0);

se2 = va/na + vb/nb;
t   = (ma - mb) / sqrt(se2);

df  = (se2^2) / ((va^2)/(na^2*(na-1)) + (vb^2)/(nb^2*(nb-1)));

% Two-sided p-value
p = 2 * tcdf(-abs(t), df);
end

function [d, g] = effect_sizes(a, b)
% Cohen's d (pooled SD) and Hedges g (bias-corrected)
na = numel(a); nb = numel(b);
sa = std(a,0); sb = std(b,0);
sp = sqrt(((na-1)*sa^2 + (nb-1)*sb^2) / (na+nb-2));

if sp == 0 || isnan(sp)
    d = NaN; g = NaN; return;
end

d = (mean(a) - mean(b)) / sp;

% Hedges correction
J = 1 - (3 / (4*(na+nb) - 9));
g = J * d;
end

function p_adj = fdr_bh(p)
% Benjamini-Hochberg FDR adjustment (q-values)
p_adj = nan(size(p));
valid = ~isnan(p);
pv = p(valid);

[m,~] = size(pv); %#ok<ASGLU>
[ps, idx] = sort(pv(:), 'ascend');
rank = (1:numel(ps))';

q = ps .* numel(ps) ./ rank;
% enforce monotonicity
q = flipud(cummin(flipud(q)));
q(q>1) = 1;

tmp = nan(size(ps));
tmp(idx) = q;

p_adj(valid) = tmp;
end
