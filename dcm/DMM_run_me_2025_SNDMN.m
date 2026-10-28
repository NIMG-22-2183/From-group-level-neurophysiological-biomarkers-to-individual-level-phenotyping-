clc
clear all;

% Local folders: edit these paths to match your data.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(repoRoot, 'data');
resultsRoot = fullfile(repoRoot, 'results');
sourceDir = fullfile(dataRoot, 'sources');
dcmDir = fullfile(resultsRoot, 'dcm');
if ~exist(dcmDir,'dir'), mkdir(dcmDir); end
timeWindowMs = [];  % Set [start end] in milliseconds for the source epoch.
assert(numel(timeWindowMs)==2 && all(isfinite(timeWindowMs)) && timeWindowMs(2)>timeWindowMs(1), ...
    'Set timeWindowMs at the top of this script before fitting DCM.');

close all;


% Initialise SPM.
%==========================================================================
spm('defaults','EEG');
%==========================================================================

BASE_BIDS_DIR = fullfile(dataRoot, 'bids');
OUTDIR                   = fullfile(BASE_BIDS_DIR, 'derivatives', 'meg_derivatives');
participants = readtable(fullfile(dataRoot, 'participants.xlsx'));
participant_id = cellstr(string(participants{:,1}));
participant_id = regexprep(participant_id, '^(?!sub-)', 'sub-');
group = cellstr(string(participants{:,2}));
do_subjs = find(ismember(lower(strtrim(string(group))), ["con","control","controls","bvftd","psp"]))';
n_subjs                  = length(participant_id);
modelling_folder = BASE_BIDS_DIR;
M_folder = fullfile(BASE_BIDS_DIR, 'derivatives', 'meg_derivatives');
ana = BASE_BIDS_DIR;

subjects_name = do_subjs ; 
Dx3 = {} ;


for i = 1: length(subjects_name)
    Subname              = participant_id{subjects_name(i)};
    sub_dir              = sourceDir;
    %     file_path      = spm_select('FPList', sub_dir,[Subname,'_DMN.mat']); % for DMN only
    file_path            = spm_select('FPList', sub_dir,append('XX_',[Subname,'_DSALL.mat'])); 
    file_path            = strtrim(file_path);
    sub                  =  Subname;
    sub = regexprep(char(sub), '^sub-', '');
    Dx3{i}               =  DCM_CSD_DSALL(sub,file_path,timeWindowMs,dcmDir);
end

function DCM = DCM_CSD_DSALL(sub,meg_data,timeWindowMs,dcmDir)

% Specify the seven-source DCM.
%==========================================================================

DCM.options.analysis = 'CSD';       % DCM for cross-spectral density 
DCM.options.spatial  = 'LFP';       % virtual electrode data as data 
DCM.options.model    = 'CMM_NMDA';  % generative model is  spm_fx_cmc
DCM.options.trials   = 1;        
DCM.options.Tdcm = timeWindowMs;

DCM.options.Fdcm     = [1   64];    % frequency range  
DCM.M.Hz             = (1 : 64)';
DCM.options.D        = 1;            % downsampling

DCM.Sname            = {'RAI' ,  'ACC' ,   'LAI' ,   'MPFC'  ,  'PCC'  , 'LAG',   'RAG'}; 
% Source data
%==========================================================================
 
DCM.xY.Dfile         = spm_eeg_load(meg_data); 

% if model one LFP channel 
%==========================================================================

DCM.xY.modality      = 'LFP';
DCM.xY.Ic            =  1:7;
DCM.xY.name          = { 'RAI' ,  'ACC' ,   'LAI' ,   'MPFC'  ,  'PCC'  , 'LAG',   'RAG'}; 

% A, B ,C parameters spesficiation in the generative model  spm_fx_cmc
% for one region source  A =[0]. For one condtion B = [0] 
%==========================================================================
Ns                  = 3+4; 
LAI                  = 1 ; 
ACC                 = 2; 
RAI                  = 3; 
MPFC                = 4;
PCC                 = 5;
LAG                 = 6;
RAG                 = 7;
%-------------------------------------------------------------------------
DCM.A{1}                = ones(Ns,Ns);
DCM.A{1}                = DCM.A{1} - diag(diag(DCM.A{1}));
DCM.A{2}                = ones(Ns,Ns);
DCM.A{2}                = DCM.A{2} - diag(diag(DCM.A{2}));
DCM.A{3}                = zeros(Ns,Ns);
DCM.G                   = ones(Ns,1);
DCM.B{1}                = 0.*DCM.A{1};
DCM.C                   = sparse(length(DCM.A{1}), 0);
DCM.options.D           = 1;
DCM.M.dipfit.Nm         = 1;
DCM.M.dipfit.model      = DCM.options.model;
DCM.M.dipfit.type       = DCM.options.spatial;
DCM.M.dipfit.Nc         = Ns;
DCM.M.dipfit.Ns         = Ns;
DCM.M.U                 = eye(Ns,Ns);
DCM.name = fullfile(dcmDir, [char(sub) '_DCM_DSALL']);
DCM.xU.X                = [1];
DCM.xU.name             = {'F1'};

[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC); % 2024
[pE,pC]  = spm_ssr_priors(pE,pC);
DCM.M.pE      = pE;
DCM.M.pC      = pC;


DCM.M.Nmax           = 2000;
DCM.options.DATA     = 0;
DCM                  = spm_dcm_csd_data(DCM);
   
DCM = spm_dcm_csd(DCM);
save(fullfile(dcmDir, [char(sub) '_DCM_DSALL.mat']), 'DCM', '-v7.3');
end
