function Extraction_Features_optimised()
% Extraction_Features_optimised
% This extracts features for CON, PSP, and bvFTD only.
%   results/features/Features_AllSubjects_Optimized.xlsx
%   Sheet: Features
% with Group labels: CON, PSP, bvFTD

clc; close all; clearvars;

% Local folders: edit these paths to match your data.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(repoRoot, 'data');
resultsRoot = fullfile(repoRoot, 'results');


%% ====================== USER SETTINGS ===========================
% Excel with ID (col1) and Group (col2). IDs look like Sub0xxx or sub-Sub0xxx.
excelFolder = dataRoot;
excelFile = 'participants.xlsx';
excelPath   = fullfile(excelFolder, excelFile);

% Folder containing the SPM MAT files (with D.data as file_array)
matFolder = fullfile(dataRoot, 'sources');

% Output folder and filenames are unchanged for downstream compatibility
outFolder   = fullfile(resultsRoot, 'features');

if ~exist(outFolder,'dir'), mkdir(outFolder); end
excelOut    = fullfile(outFolder, 'Features_AllSubjects_Optimized.xlsx');
matOut      = fullfile(outFolder, 'Features_AllSubjects_Optimized.mat');

% ========= Analysis params (saved to .mat for reproducibility)
fs          = 500;                           % (Hz) sampling rate

% Canonical bands
bands.Delta = [1 3];                % (Hz)
bands.Theta = [3 7];                % (Hz)
bands.Alpha = [8 13];               % (Hz)
bands.Beta  = [13 30];              % (Hz)
bands.Gamma = [30 45];              % (Hz)

% PSD settings
psd.minHz        = 1.0;             % analyze PSD from 1 Hz
psd.maxHz        = 45;              % upper limit
psd.nfft_mode    = 'fullseg';       % 'fullseg' uses nSamples; or set numeric NFFT
psd.window       = 'hamming';       % window type for pwelch()
psd.overlap      = 0;               % samples overlapped

% Use same PSD domain for relative powers
denBand = [psd.minHz psd.maxHz];

% Aperiodic (1/f) fit on log10 P vs log10 f
oneoverf.fitBand      = [2 40];        % Hz
oneoverf.excludeBands = [8 13; 13 30]; % exclude alpha/beta when fitting 1/f
oneoverf.robust       = false;         % if true, uses robustfit (Stats TB)
oneoverf.eps          = 1e-30;         % epsilon for logs

% Peak finding in alpha/beta on raw PSD
peak.minProminence   = 0.0;         % power units
peak.minWidthHz      = 0.0;         % Hz optional minimum width

% CFC params (alpha phase -> beta amplitude)
cfc.nBins         = 18;             % Tort MI histogram bins
cfc.minLenCycles  = 3;              % epochs need min alpha cycles
cfc.eps           = 1e-12;
cfc.enable        = true;           % set false if you want to skip PAC

% Random seed
rng(1);

%% ================== LOAD EXCEL & MATCH FILES ====================
assert(isfile(excelPath), 'Excel file not found: %s', excelPath);
T      = readtable(excelPath, 'PreserveVariableNames', true);
IDs    = T{:,1};
Groups = T{:,2};

IDs    = stringify_cell_column(IDs);
Groups = stringify_cell_column(Groups);

% Canonicalise IDs to 'Sub####'
IDs = regexprep(IDs, '^sub-', '');
IDs = regexprep(IDs, '.*?(Sub0*\d+).*', '$1');

m = dir(fullfile(matFolder, 'XX_sub-Sub0*_DSALL.mat'));
matNames = {m.name}';
matIDs   = regexp(matNames, '^(XX_sub-Sub0\d+)_DSALL\.mat$', 'tokens', 'once');
validMat = ~cellfun(@isempty, matIDs);
matNames = matNames(validMat);
matIDs   = matIDs(validMat);
matIDs   = cellfun(@(c) c{1}, matIDs, 'UniformOutput', false);

G_CON   = {};
G_PSP   = {};
G_bvFTD = {};

skipped_PC = 0;
skipped_PSPPC_merged = 0;
skipped_other = 0;

for i = 1:numel(matIDs)
    thisID  = regexp(matIDs{i}, 'Sub\d+', 'match', 'once');
    idx     = find(strcmpi(IDs, thisID), 1, 'first');

    if ~isempty(idx)
        labelRaw = Groups{idx};
        label = canonical_group_label(labelRaw);

        switch label
            case 'CON'
                G_CON{end+1} = fullfile(matFolder, matNames{i}); 
            case 'PSP'
                G_PSP{end+1} = fullfile(matFolder, matNames{i}); 
            case 'bvFTD'
                G_bvFTD{end+1} = fullfile(matFolder, matNames{i}); 
            case 'PC'
                skipped_PC = skipped_PC + 1;
            case 'PSPPC_EXCLUDED'
                skipped_PSPPC_merged = skipped_PSPPC_merged + 1;
                warning('Subject %s has merged label "%s". Excluded because PSP vs PC cannot be separated from a merged label.', thisID, labelRaw);
            otherwise
                skipped_other = skipped_other + 1;
        end
    else
        warning('No group found for MAT file: %s', matNames{i});
    end
end

G_CON   = G_CON(:);
G_PSP   = G_PSP(:);
G_bvFTD = G_bvFTD(:);

fprintf('Matched for feature extraction: CON=%d, PSP=%d, bvFTD=%d\n', numel(G_CON), numel(G_PSP), numel(G_bvFTD));
fprintf('Excluded: PC=%d, merged PSP/PC labels=%d, other/unmatched labels=%d\n', skipped_PC, skipped_PSPPC_merged, skipped_other);

if numel(G_PSP) == 0
    warning(['No PSP files were matched. Check the exact labels in column 2 of %s. ', ...
             'They must identify PSP separately from PC; a merged PSP/PC label cannot be split automatically.'], excelPath);
end

%% =============== FEATURE EXTRACTION (per subject) ===============
fprintf('\nExtracting optimized features...\n');
[FeatTbl, chanLabels, featureList] = extract_all_features_only_optimized( ...
    struct('CON',{G_CON}, 'PSP',{G_PSP}, 'bvFTD',{G_bvFTD}), ...
    fs, bands, denBand, psd, oneoverf, peak, cfc);

assert(~isempty(FeatTbl) && height(FeatTbl)>0, 'No features extracted.');

% Save outputs with unchanged downstream-compatible names
writetable(FeatTbl, excelOut, 'Sheet','Features', 'WriteMode','overwritesheet');
save(matOut, 'FeatTbl', 'chanLabels', 'featureList', 'fs', 'bands', 'denBand', 'psd', 'oneoverf', 'peak', 'cfc');

fprintf('\nSaved feature table to:\n  %s\nSaved .mat bundle to:\n  %s\n', excelOut, matOut);
fprintf('\nGroups saved in feature table:\n');
disp(groupcounts(categorical(FeatTbl.Group)));

fprintf('\nKept feature families per channel:\n');
disp(featureList(:));

end % === main ===


%% ========================== HELPERS =============================
function [FeatTbl, chanLabels, featureList] = extract_all_features_only_optimized(groupFiles, fs, bands, denBand, psd, oneoverf, peak, cfc)

    featureList = { ...
        'RelDelta', ...
        'RelTheta', ...
        'RelAlpha', ...
        'RelBeta', ...
        'RelGamma', ...
        'Ratio_Delta_Alpha', ...
        'Ratio_Theta_Alpha', ...
        'OneOverF_slope', ...
        'Alpha_peak_freq', ...
        'Beta_peak_freq', ...
        'Alpha_residual_max', ...
        'Beta_residual_max', ...
        'SpecEntropy', ...
        'MedianFreq', ...
        'PAC_MI_alphaPhase_betaAmp'};

    FeatTbl = table();
    chanLabels = {};

    % Downstream-compatible group names. No PC/PSP combined group is created.
    groups = {'CON','PSP','bvFTD'};

    for gi = 1:numel(groups)
        gname = groups{gi};

        switch gname
            case 'CON'
                files = groupFiles.CON;
            case 'PSP'
                files = groupFiles.PSP;
            case 'bvFTD'
                files = groupFiles.bvFTD;
        end

        for f = 1:numel(files)
            matFull = files{f};

            S = safe_load_D(matFull);
            if isempty(S), continue; end
            if ~isfield(S,'D')
                warning('No D in %s (skip)', matFull);
                continue;
            end

            D = S.D;
            Y = D.data;
            Y = maybe_repair_file_array(Y, matFull);

            if isempty(Y)
                warning('Missing data for %s (skip)', matFull);
                continue;
            end

            sz = size(Y); % [nChans nSamples nEpochs]
            if numel(sz) < 2
                warning('Bad dims for %s: %s', matFull, mat2str(sz));
                continue;
            end

            nChans   = sz(1);
            nSamples = sz(2);
            if numel(sz) >= 3
                nEpochs = sz(3);
            else
                nEpochs = 1;
            end

            if isempty(chanLabels)
                chanLabels = get_channel_labels_safe(D, nChans);
            else
                nChans = min(nChans, numel(chanLabels));
            end

            tok = regexp(matFull, '(Sub0\d+)', 'tokens', 'once');
            if ~isempty(tok)
                sid = tok{1};
            else
                sid = sprintf('sub_%04d', height(FeatTbl)+1);
            end

            % Critical downstream-compatible labels: CON, PSP, bvFTD
            row = table(string(sid), string(gname), 'VariableNames', {'SubjectID','Group'});

            for ch = 1:nChans
                relDe = nan(nEpochs,1);
                relTh = nan(nEpochs,1);
                relAl = nan(nEpochs,1);
                relBe = nan(nEpochs,1);
                relGa = nan(nEpochs,1);

                ratio_Delta_Alpha = nan(nEpochs,1);
                ratio_Theta_Alpha = nan(nEpochs,1);

                aperiodic_slope = nan(nEpochs,1);

                alpha_peak_f = nan(nEpochs,1);
                beta_peak_f  = nan(nEpochs,1);

                alpha_res_max = nan(nEpochs,1);
                beta_res_max  = nan(nEpochs,1);

                spec_entropy = nan(nEpochs,1);
                spec_median_f = nan(nEpochs,1);

                pac_AB = nan(nEpochs,1);

                for e = 1:nEpochs
                    seg = detrend(double(Y(ch,:,e)));

                    [Pxx, fv] = compute_psd(seg, fs, psd, nSamples);

                    keep = isfinite(fv) & isfinite(Pxx) & (fv >= psd.minHz) & (fv <= psd.maxHz);
                    f1 = fv(keep);
                    P1 = Pxx(keep);

                    if isempty(f1) || isempty(P1) || numel(f1) < 2
                        continue;
                    end

                    % ---- safer PSD integration instead of bandpower(...,'psd')
                    de = integrate_psd_bandpower(P1, f1, bands.Delta);
                    th = integrate_psd_bandpower(P1, f1, bands.Theta);
                    al = integrate_psd_bandpower(P1, f1, bands.Alpha);
                    be = integrate_psd_bandpower(P1, f1, bands.Beta);
                    ga = integrate_psd_bandpower(P1, f1, bands.Gamma);
                    dTotal = integrate_psd_bandpower(P1, f1, denBand);

                    if isfinite(dTotal) && dTotal > 0
                        relDe(e) = de / dTotal;
                        relTh(e) = th / dTotal;
                        relAl(e) = al / dTotal;
                        relBe(e) = be / dTotal;
                        relGa(e) = ga / dTotal;
                    end

                    ratio_Delta_Alpha(e) = safe_div(de, al);
                    ratio_Theta_Alpha(e) = safe_div(th, al);

                    [b0, b1, maskFit] = fit_aperiodic_loglog(f1, P1, oneoverf); 
                    aperiodic_slope(e) = b1;

                    [~, fA] = peak_in_band(f1, P1, bands.Alpha, peak);
                    [~, fB] = peak_in_band(f1, P1, bands.Beta,  peak);
                    alpha_peak_f(e) = fA;
                    beta_peak_f(e)  = fB;

                    if any(maskFit) && isfinite(b0) && isfinite(b1)
                        logP = log10(P1 + oneoverf.eps);
                        logF = log10(f1 + eps);
                        fitLine = b0 + b1 .* logF;
                        res = logP - fitLine;

                        [alpha_res_max(e), ~, ~] = band_residual_stats(f1, res, bands.Alpha);
                        [beta_res_max(e),  ~, ~] = band_residual_stats(f1, res, bands.Beta);
                    end

                    [spec_entropy(e), medTmp, ~] = spectral_summaries(f1, P1);
                    spec_median_f(e) = medTmp;

                    if cfc.enable
                        try
                            pac_AB(e) = compute_pac_alpha_beta(seg, fs, bands, cfc);
                        catch
                            pac_AB(e) = NaN;
                        end
                    end
                end

                tag = matlab.lang.makeValidName(chanLabels{ch});

                row.(['RelDelta__' tag]) = mean(relDe, 'omitnan');
                row.(['RelTheta__' tag]) = mean(relTh, 'omitnan');
                row.(['RelAlpha__' tag]) = mean(relAl, 'omitnan');
                row.(['RelBeta__'  tag]) = mean(relBe, 'omitnan');
                row.(['RelGamma__' tag]) = mean(relGa, 'omitnan');

                row.(['Ratio_Delta_Alpha__' tag]) = mean(ratio_Delta_Alpha, 'omitnan');
                row.(['Ratio_Theta_Alpha__' tag]) = mean(ratio_Theta_Alpha, 'omitnan');

                row.(['OneOverF_slope__' tag]) = mean(aperiodic_slope, 'omitnan');

                row.(['Alpha_peak_freq__' tag]) = mean(alpha_peak_f, 'omitnan');
                row.(['Beta_peak_freq__'  tag]) = mean(beta_peak_f,  'omitnan');

                row.(['Alpha_residual_max__' tag]) = mean(alpha_res_max, 'omitnan');
                row.(['Beta_residual_max__'  tag]) = mean(beta_res_max,  'omitnan');

                row.(['SpecEntropy__' tag]) = mean(spec_entropy, 'omitnan');
                row.(['MedianFreq__'  tag]) = mean(spec_median_f, 'omitnan');

                if cfc.enable
                    row.(['PAC_MI_alphaPhase_betaAmp__' tag]) = mean(pac_AB, 'omitnan');
                end
            end

            FeatTbl = [FeatTbl; row]; 
        end
    end
end


%% ---------------- group-label utilities ----------------
function c = stringify_cell_column(x)
    if iscell(x)
        c = cell(size(x));
        for i = 1:numel(x)
            if ischar(x{i})
                c{i} = x{i};
            elseif isstring(x{i})
                c{i} = char(x{i});
            elseif isnumeric(x{i})
                c{i} = num2str(x{i});
            elseif iscategorical(x{i})
                c{i} = char(string(x{i}));
            else
                c{i} = char(string(x{i}));
            end
        end
    elseif isstring(x)
        c = cellstr(x);
    elseif ischar(x)
        c = cellstr(x);
    elseif isnumeric(x)
        c = arrayfun(@(n) char(string(n)), x, 'UniformOutput', false);
    elseif iscategorical(x)
        c = cellstr(string(x));
    else
        c = cellstr(string(x));
    end
    c = c(:);
end

function label = canonical_group_label(rawLabel)
    s = lower(strtrim(char(string(rawLabel))));
    s = strrep(s, '+', '');
    s = strrep(s, '-', '');
    s = strrep(s, '_', '');
    sNoSlash = strrep(s, '/', '');
    sCompact = regexprep(sNoSlash, '\s+', '');

    if any(strcmp(sCompact, {'con','control','controls','healthycontrol','healthycontrols'}))
        label = 'CON';
    elseif any(strcmp(sCompact, {'bvftd','bvfrontotemporaldementia'}))
        label = 'bvFTD';
    elseif any(strcmp(sCompact, {'psp','progressivesupranuclearpalsy'}))
        label = 'PSP';
    elseif any(strcmp(sCompact, {'pc'}))
        label = 'PC';
    elseif any(strcmp(sCompact, {'pcpsp','psppc','psppccombined','pcpspcombined'})) || contains(s, 'psp/pc') || contains(s, 'pc/psp')
        label = 'PSPPC_EXCLUDED';
    else
        label = char(strtrim(string(rawLabel)));
    end
end

%% ---------------- low-level utilities ----------------
function S = safe_load_D(matFull)
    S = struct();
    try
        S = load(matFull, 'D');
    catch
        warning('Could not load %s', matFull);
        S = [];
    end
end

function Y = maybe_repair_file_array(Y, matFull)
    try
        hasFname = isprop(Y,'fname');
    catch
        try
            hasFname = isfield(struct(Y),'fname');
        catch
            hasFname = false;
        end
    end

    if hasFname
        try
            needRepair = ~isempty(Y.fname) && ~exist(Y.fname,'file');
        catch
            needRepair = false;
        end
    else
        needRepair = false;
    end

    if needRepair
        [mPath,mName,~] = fileparts(matFull);
        datGuess = fullfile(mPath, [mName '.dat']);
        if exist(datGuess,'file')
            try
                Y = file_array(datGuess, Y.dim, Y.dtype, Y.offset, Y.scl_slope, Y.scl_inter);
                fprintf('Repaired .dat for %s\n', mName);
            catch
                Y = [];
            end
        else
            Y = [];
        end
    end
end

function lbls = get_channel_labels_safe(D, nChans)
    try
        lbls = cell(1,nChans);
        for ch = 1:nChans
            try
                lbl = D.chanlabels(ch);
            catch
                lbl = chanlabels(D, ch);
            end
            if iscell(lbl), lbl = lbl{1}; end
            if isempty(lbl), lbl = sprintf('Ch%d', ch); end
            lbls{ch} = char(lbl);
        end
    catch
        lbls = arrayfun(@(c) sprintf('Ch%d',c), 1:nChans, 'UniformOutput', false);
    end
end

function [Pxx,f] = compute_psd(seg, fs, psd, nSamples)
    switch psd.nfft_mode
        case 'fullseg'
            nfft = nSamples;
        otherwise
            if isnumeric(psd.nfft_mode)
                nfft = psd.nfft_mode;
            else
                nfft = nSamples;
            end
    end

    switch lower(psd.window)
        case 'hamming'
            win = hamming(nSamples);
        case 'hann'
            win = hann(nSamples);
        otherwise
            win = hamming(nSamples);
    end

    [Pxx,f] = pwelch(seg, win, psd.overlap, nfft, fs);
end

function pwr = integrate_psd_bandpower(P, f, band)


    pwr = NaN;

    if isempty(P) || isempty(f) || numel(P) ~= numel(f) || numel(f) < 2
        return;
    end

    P = P(:);
    f = f(:);

    good = isfinite(P) & isfinite(f);
    P = P(good);
    f = f(good);

    if numel(f) < 2
        return;
    end

    [f, ord] = sort(f);
    P = P(ord);

    lo = max(band(1), min(f));
    hi = min(band(2), max(f));

    if ~(isfinite(lo) && isfinite(hi)) || hi <= lo
        return;
    end

    mask = (f >= lo) & (f <= hi);

    if nnz(mask) >= 2
        pwr = trapz(f(mask), P(mask));
        return;
    end

    fi = f;
    Pi = P;

    try
        gridf = unique([lo; fi(fi > lo & fi < hi); hi]);
        if numel(gridf) < 2
            return;
        end
        gridP = interp1(fi, Pi, gridf, 'linear', 'extrap');
        pwr = trapz(gridf, gridP);
    catch
        pwr = NaN;
    end
end

function val = safe_div(a, b)
    if ~isfinite(a) || ~isfinite(b) || b == 0
        val = NaN;
    else
        val = a / b;
    end
end

function [b0, b1, maskFit] = fit_aperiodic_loglog(f, P, oneoverf)
    mask = (f >= oneoverf.fitBand(1)) & (f <= oneoverf.fitBand(2));

    for k = 1:size(oneoverf.excludeBands,1)
        band = oneoverf.excludeBands(k,:);
        mask = mask & ~(f >= band(1) & f <= band(2));
    end

    maskFit = mask & isfinite(P) & P > 0;
    if ~any(maskFit)
        b0 = NaN;
        b1 = NaN;
        return;
    end

    xf = log10(f(maskFit));
    yf = log10(P(maskFit) + oneoverf.eps);

    if oneoverf.robust && exist('robustfit','file') == 2
        B = robustfit(xf, yf);
        b0 = B(1);
        b1 = B(2);
    else
        p  = polyfit(xf, yf, 1);
        b1 = p(1);
        b0 = p(2);
    end
end

function [pkPower, pkFreq] = peak_in_band(f, P, band, peak)
    inB = (f >= band(1) & f <= band(2));
    if ~any(inB)
        pkPower = NaN;
        pkFreq  = NaN;
        return;
    end

    Pf = P(inB);
    ff = f(inB);

    try
        if peak.minWidthHz > 0
            minWidthSamples = max(1, round(peak.minWidthHz / max(eps, median(diff(ff)))));
            [pks,locs] = findpeaks(Pf, ff, ...
                'MinPeakProminence', peak.minProminence, ...
                'MinPeakWidth', minWidthSamples);
        else
            [pks,locs] = findpeaks(Pf, ff, ...
                'MinPeakProminence', peak.minProminence);
        end

        if isempty(pks)
            [pkPower, idx] = max(Pf);
            pkFreq = ff(idx);
        else
            [pkPower, mx] = max(pks);
            pkFreq = locs(mx);
        end
    catch
        [pkPower, idx] = max(Pf);
        pkFreq = ff(idx);
    end
end

function [r_max, f_at_max, r_mean] = band_residual_stats(f, res, band)
    m = (f >= band(1) & f <= band(2));
    if ~any(m)
        r_max = NaN;
        f_at_max = NaN;
        r_mean = NaN;
        return;
    end

    rr = res(m);
    ff = f(m);

    [r_max, idx] = max(rr);
    f_at_max = ff(idx);
    r_mean   = mean(rr, 'omitnan');
end

function [H, f_med, f_edge] = spectral_summaries(f, P)
    P = P(:);
    f = f(:);
    P(P < 0) = 0;

    if ~any(isfinite(P))
        H = NaN;
        f_med = NaN;
        f_edge = NaN;
        return;
    end

    den = sum(P, 'omitnan');
    if ~isfinite(den) || den <= 0
        H = NaN;
        f_med = NaN;
        f_edge = NaN;
        return;
    end

    Pn = P / den;

    H = -nansum(Pn .* log2(Pn + eps));
    H = H / log2(numel(Pn));

    c = cumsum(P);
    tot = c(end);
    if ~isfinite(tot) || tot <= 0
        f_med = NaN;
        f_edge = NaN;
    else
        f_med  = interp1(c, f, 0.5 * tot, 'linear', 'extrap');
        f_edge = interp1(c, f, 0.9 * tot, 'linear', 'extrap');
    end
end

%% ---------------- PAC helper ----------------
function pacAB = compute_pac_alpha_beta(x, fs, bands, p)
    minSamples = max(round(p.minLenCycles * fs / max(bands.Alpha(1), eps)), 20);
    if numel(x) < minSamples
        error('Not enough samples for CFC analysis');
    end

    xa = bandpass_zp(x, fs, bands.Alpha);
    xb = bandpass_zp(x, fs, bands.Beta);

    ha = hilbert(xa);
    hb = hilbert(xb);

    phi_a = angle(ha);
    amp_b = abs(hb);

    pacAB = tort_mi(phi_a, amp_b, p.nBins, p.eps);
end

function y = bandpass_zp(x, fs, band)
    Wn = band(:) / (fs/2);
    Wn = max(min(Wn, 0.999), 0.001);
    [b,a] = butter(4, Wn, 'bandpass');
    y = filtfilt(b, a, double(x));
end

function mi = tort_mi(phase, amp, nBins, epsv)
    if nargin < 4
        epsv = 1e-12;
    end

    phase = phase(:);
    amp   = amp(:);

    good = isfinite(phase) & isfinite(amp);
    phase = phase(good);
    amp   = amp(good);

    if isempty(phase) || isempty(amp)
        mi = NaN;
        return;
    end

    edges = linspace(-pi, pi, nBins+1);
    [~,bin] = histc(phase, edges); 

    P = zeros(nBins,1);
    for k = 1:nBins
        ak = amp(bin == k);
        if isempty(ak)
            P(k) = 0;
        else
            P(k) = mean(ak, 'omitnan');
            if ~isfinite(P(k))
                P(k) = 0;
            end
        end
    end

    P = P / (sum(P) + epsv);
    U = (1/nBins) * ones(nBins,1);

    mi = nansum(P .* log((P + epsv) ./ (U + epsv))) / log(nBins);
    if ~isfinite(mi)
        mi = NaN;
    end
end
