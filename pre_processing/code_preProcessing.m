% Local folders: edit these paths to match your data.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(repoRoot, 'data');
resultsRoot = fullfile(repoRoot, 'results');



% Process MaxFilter output: conversion, filtering, artefact detection and epoching. 
%==========================================================================

BASE_BIDS_DIR = fullfile(dataRoot, 'bids');
OUTDIR                               = fullfile(BASE_BIDS_DIR, 'derivatives', 'meg_derivatives');
participants = readtable(fullfile(dataRoot, 'participants.xlsx'));
participant_id = cellstr(string(participants{:,1}));
participant_id = regexprep(participant_id, '^(?!sub-)', 'sub-');
group = cellstr(string(participants{:,2}));
do_subjs = find(ismember(lower(strtrim(string(group))), ["con","control","controls","bvftd","psp"]))';
n_subjs                              = length(participant_id); 

cd (OUTDIR)
%=======================================================================

chanfile = fullfile(dataRoot, 'chan_select_MEG_EEG_STI101.mat');
chan              = load(chanfile);%ref_old_names = {'EOG061','EOG062','ECG063'};
ref_old_names     = {'EOG001','EOG002','ECG003'};
ref_new_names     = {'HEOG';'VEOG';'ECG'};
%==========================================================================
for sub = do_subjs
%==========================================================================

       sub_dir           = fullfile(OUTDIR, participant_id{sub});
       % change folder to the subject directory folder
        cd(fullfile(sub_dir, 'ses-meg1', 'meg'));
       dir               = fullfile(sub_dir, 'ses-meg1', 'meg');
 
       
       Subname           = participant_id{sub};
       files             = spm_select('FPList', dir, '.*Rest_proc-sss\.fif$');
        S                = [];
        S.outfile        = fullfile(sub_dir, 'ses-meg1', 'meg','C');
        S.dataset        = spm_select('FPList', dir, '.*Rest_proc-sss\.fif$');
        S.save           = 0;  % Save after channel names and units are set.
        S.reviewtrials   = 0;
        S.channels       = chan.label;
        S.continuous     = 1;
        S.checkboundary  = 0;
        D                = spm_eeg_convert(S);
        tc               = find(strcmp(D.chanlabels,'STI101'));
        D(tc,:)          = D(tc,:)/1e6;  % 1e6 new scaling introduced by SPM12!
        D                = units(D,tc,'V');
        
        for c = 1:length(ref_old_names)
                    ch = indchannel(D,ref_old_names{c});
                    D = chanlabels(D,ch,ref_new_names{c});
        end
        D.save;
%==============================================downsample
        S = [];
        S.D = spm_select('FPList', dir, '^C\.mat$');
        S.method = 'downsample';
        S.fsample_new = 500;
        D = spm_eeg_downsample(S); 
        D.save;
% delete(spm_select('FPList', dir, '^C.mat')); % Optional: remove an intermediate after checking the output.
% delete(spm_select('FPList', dir, '^C.dat')); % Optional: remove an intermediate after checking the output.
        
%=============================================Filter=======================
        filtwindow    = [0.1 100]; % Filter window
        low          = filtwindow(2); 
        high         = filtwindow(1);
        S            = [];
        S.D          = spm_select('FPList', dir, '^dC\.mat$');
        S.type       = 'butterworth';
        S.order      = 5;
        S.band       = 'high';
        S.freq       = high;
        D            = spm_eeg_filter(S);
% delete(spm_select('FPList', dir, '^dC.mat')); % Optional: remove an intermediate after checking the output.
% delete(spm_select('FPList', dir, '^dC.dat')); % Optional: remove an intermediate after checking the output.
        
%==========================================================================
        S= [];
        S.D         = spm_select('FPList', dir, '^fdC\.mat$');
        S.type      = 'butterworth';
        S.order     = 5;
        S.band      = 'low';
        S.freq      =  low;
        D          = spm_eeg_filter(S);
% delete(spm_select('FPList', dir, '^fdC.mat')); % Optional: remove an intermediate after checking the output.
% delete(spm_select('FPList', dir, '^fdC.dat')); % Optional: remove an intermediate after checking the output.
%==========================================================================
        S = [];
        S.D         = spm_select('FPList', dir, '^ffdC\.mat$');
        S.type      = 'butterworth';
        S.order     = 5;
        S.band      = 'stop';
        S.freq      = [48 52];
        D           = spm_eeg_filter(S);
% delete(spm_select('FPList', dir, '^ffdC.mat')); % Optional: remove an intermediate after checking the output.
% delete(spm_select('FPList', dir, '^ffdC.dat')); % Optional: remove an intermediate after checking the output.
        
        S = [];
        S.D = spm_select('FPList', dir, '^fffdC\.mat$');
        S.type = 'butterworth';
        S.order = 5;
        S.band = 'stop';
        S.freq = [98 102];
        D = spm_eeg_filter(S);
% delete(spm_select('FPList', dir, '^fffdC.mat')); % Optional: remove an intermediate after checking the output.
% delete(spm_select('FPList', dir, '^fffdC.dat')); % Optional: remove an intermediate after checking the output.

%================================OSL Artefact detection====================
        add        = spm_select('FPList', dir, '^ffffdC\.mat$');
        D          = spm_eeg_load(add); 
        D          = D.badchannels(1:D.nchannels,0);
        try
            D2=osl_detect_artefacts(D,'modalities',{'EEG','MEGMAG','MEGPLANAR'});
        catch
            D2=osl_detect_artefacts(D,'modalities',{'MEGMAG','MEGPLANAR'});
        end
        D=D2;
        D.save;
%==========================================================================
         % Epoch the continuous object marked by OSL above.
 
         S= [];
         S.D          = D;
         S.trialength = 1000;
         S.bc         = 0;
         D            = spm_eeg_epochs(S);    
%======================================== Keep the continous---------------
%data======================================================================
%          delete(spm_select('FPList', dir, '^fffdC.mat')); 
%         delete(spm_select('FPList', dir, '^fffdC.dat'));  

%================================RUN OSL Artefact on Epoched data----------
  D           = spm_eeg_load(spm_select('FPList', dir, '.*effffdC.mat$')); 
  D2          = osl_detect_artefacts(D,'modalities',{'MEGMAG','MEGPLANAR'});
  D           = D2;
  D.save
  S           =[]; 
  S.D         = spm_select('FPList', dir, '.*effffdC.mat$');
  D          = spm_eeg_remove_bad_trials(S); 
  D.save
%=========================Saved cleaned Epoched data======================
%   file_path1 = fullfile(dir, [Subname '_EP.mat']);
%              save(file_path1,'D','-v7.3');
%   file_path2 = fullfile(dir, [Subname '_EP.dat']);
%             save(file_path2,'D','-v7.3');
%============================================Delete all redundant files
% delete(spm_select('FPList', dir, '^effffdC.mat')); % Optional: remove an intermediate after checking the output.
% delete(spm_select('FPList', dir, '^effffdC.dat')); % Optional: remove an intermediate after checking the output.
%  delete(spm_select('FPList', dir, '^reffffdC.mat')); 
%  delete(spm_select('FPList', dir, '^reffffdC.dat'));  
end 
 
%==For info only =This section is optional in case that data is not=========
%=====reday yet. These are extra removal options in
%SPM12.====================
% D = spm_eeg_artefact_peak2peak(S);

% return;

% S =[];                
% S.mode                = 'mark'; 
% S.D                   = spm_eeg_load(spm_select('FPList', pwd, '.*EP.mat$'));
% S.methods(1).channels = 'EOG061'; 
% S.methods(1).fun      = 'eyeblink';
% S.append              = true;
% S.methods(1).settings.threshold = 4;
% S.methods(1).settings.excwin = 0;
% S.methods(2).channels = 'EOG062'; 
% S.methods(2).fun      = 'eyeblink';
% S.append              = true;
% S.methods(2).settings.threshold = 4;
% S.methods(2).settings.excwin = 0;
% D2 = spm_eeg_artefact(S);
% 
% S_r =[];                
% S_r.mode                = 'reject'; 
% S_r.methods(1).channels    = 'ALL'; 
% S_r.D                   = D2;
% S_r.badchanthresh       =.2;
% S_r.methods(1).fun      = 'events';
% S_r.methods(1).settings.threshold= 1;
% S_r.methods(1).settings.whatevents = 'artefact_eyeblink';
% D = spm_eeg_artefact(S_r);
% close all;
% %==========================================================================
% 
% S =[];                
% S.mode               = 'reject'; 
% S.D                  = spm_eeg_load(spm_select('FPList', pwd, 'aaEP.mat$'));
% % S.badchanthresh    =.3;
% S.methods(1).channels = 'MEG'; 
% S.methods(1).fun      = 'peak2peak';
% S.append              = true;
% S.methods(1).settings.threshold = 1e-16;
% D = spm_eeg_artefact(S);
% %==========================================================================
% S =[];                
% S.mode                = 'reject'; 
% S.append              = true;
% S.D                   = spm_eeg_load(spm_select('FPList', pwd, 'aaaEP.mat$'));
% % S.badchanthresh     =.3;
% S.methods(1).channels = 'MEG'; 
% S.methods(1).fun      = 'jump';
% S.methods(1).settings.threshold = 1;
% S.methods(2).settings.excwin = 0;
% D = spm_eeg_artefact(S);
% %==========================================================================
% S =[];                
% S.mode                = 'reject'; 
% S.append              = true;
% S.D                   = spm_eeg_load(spm_select('FPList', pwd, 'aaaaEP.mat$'));
% % S.badchanthresh     =.3;
% S.methods(1).channels = 'MEG'; 
% S.methods(1).fun      = 'flat';
% S.methods(1).settings.threshold = .001;
% S.methods(2).settings.seqlength= 1000;
% D = spm_eeg_artefact(S);
% 
% %==========================================================================
% S =[];                
% S.mode                = 'reject'; 
% S.append              = true;
% S.D                   = spm_eeg_load(spm_select('FPList', pwd, 'aaaaaEP.mat$'));
% % S.badchanthresh     =.3;
% S.methods(1).channels = 'MEG'; 
% S.methods(1).fun      = 'zscore'; %zscorediff
% S.methods(1).settings.threshold = .2;
% % S.methods(2).settings.seqlength= 00;
% D = spm_eeg_artefact(S);
% 
% %==========================================================================
% S =[];                
% S.mode                = 'reject'; 
% S.append              = true;
% S.D                   = spm_eeg_load(spm_select('FPList', pwd, 'aaaaaaEP.mat$'));
% % S.badchanthresh     =.3;
% S.methods(1).channels = 'MEG'; 
% S.methods(1).fun      = 'threshchan';
% S.methods(1).settings.threshold = .1;
% S.methods(2).settings.excwin = 0;
% D = spm_eeg_artefact(S);
% %==========================================================================
% S =[];                
% S.mode                = 'reject'; 
% S.D                   = spm_eeg_load(spm_select('FPList', pwd, 'aaaaaaaEP.mat$'));
% S.methods(1).channels = {'MEG', 'EEG'}; 
% S.methods(1).fun      = 'zscorediff'; %
% S.methods(1).settings.threshold = 2;
% D = spm_eeg_artefact(S);
% %==========================================================================

