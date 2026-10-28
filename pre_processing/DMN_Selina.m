clc
clear all, 

% Local folders: edit these paths to match your data.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(repoRoot, 'data');
resultsRoot = fullfile(repoRoot, 'results');
sourceDir = fullfile(dataRoot, 'sources');
if ~exist(sourceDir,'dir'), mkdir(sourceDir); end

close all, 
%==========================================================================
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

 %=========================================================================
 epo           = [0 1000];
 f_win         = [0.1 100];
 mesh_size     = 2;
 hanning       = 0;
 trialtype     = 'induced';
 inv_trls      = {'Undefined'};
 use_headshape = 1;
 con_win       = {epo};
 inv_typ       = 'COH';    
 val           = 1; 
%==========================================================================


error_ind =[];

for subnum = do_subjs
    try
    Subname             = participant_id{subnum};
    sub_dir             = fullfile(M_folder, participant_id{subnum});
    cd(sub_dir )
    file_path           = spm_select('FPList', fullfile(sub_dir, 'ses-meg1', 'meg'), '^reffffdC.*\.mat$');
    file_path           = strtrim(file_path);
    D                  = spm_eeg_load(file_path);
    
    anat_path = fullfile(BASE_BIDS_DIR, Subname, 'ses-meg1', 'anat');
    nii_file          = spm_select('FPList',anat_path, '.*\.nii$');
    
    inv_mods             = {'MEGPLANAR', 'MEG'};
    D.val                = val;
    D.inv{val}.date      = strvcat(date, datestr(now, 15));
    D.inv{val}.comment   = {sprintf('%s_%s', inv_typ, char(inv_mods)')};
    D.inv{val}.mesh      = [];
    D.inv{val}.datareg   = [];
    D.inv{val}.forward   = [];
    D.inv{val}.mesh.sMRI = nii_file ;
    D.inv{val}.mesh      = spm_eeg_inv_mesh(D.inv{val}.mesh.sMRI, mesh_size);
    
    D.save;
    cd(sub_dir);
    
    newmrifid           = [];
    newmrifid.fid.pnt   = D.inv{val}.mesh.fid.fid.pnt(1:3,:); % EK force selection
    newmrifid.fid.label = {'Nasion';'LPA';'RPA'};
    newmrifid.pnt       = D.inv{val}.mesh.fid.pnt;
    
    meegfid = D.fiducials;
    meegfid.pnt(find(meegfid.pnt(:,2)>0 & meegfid.pnt(:,3)<0),:) = [];
    
    D = spm_eeg_inv_datareg_ui_ek(D, val, meegfid, newmrifid,use_headshape);
    D.save;
    
    fprintf(1, 'Creating forward model\n');
    D.inv{val}.forward = struct([]);
    
    for ind = 1:length(D.inv{val}.datareg) %EK
        if strcmp(D.inv{val}.datareg(ind).modality,'MEG')
            D.inv{val}.forward(ind).voltype = 'Single Shell';
        elseif strcmp(D.inv{val}.datareg(ind).modality,'EEG')
            D.inv{val}.forward(ind).voltype = 'EEG BEM';
        end
    end
    
    fprintf(1, 'Computing leadfield\n');
    D = spm_eeg_inv_forward(D,val);
    D.save;
    
    D.inv{val}.inverse = [];
    D.inv{val}.inverse.trials = inv_trls ;
    D.inv{val}.inverse.type   = inv_typ;
    D.inv{val}.inverse.lpf    = f_win(1);
    D.inv{val}.inverse.hpf    = f_win(2);
    D.inv{val}.inverse.woi    = epo;% as much as the memory can handle
    D.inv{val}.inverse.Han    = hanning; %no hanning
    D.inv{val}.inverse.modality = inv_mods;
    D.inv{val}.inverse.dplot = 0;
    
    D = spm_eeg_invert(D,val);
    
    fileend='_rest_inverted.mat';
    
    
    save(fullfile(sourceDir, [Subname fileend]), 'D', '-v7.3');
    catch ME
        
       error_ind(end+1) = subnum ;
    end
    
end

do_subjs = setdiff(do_subjs, error_ind);

%==========================================================================



for subnum = do_subjs
Subname             = participant_id{subnum}
sub_dir = sourceDir;
cd(sub_dir )
pattern = ['^' Subname '_rest_inverted\.mat$'];
file_path           = spm_select('FPList', sub_dir, pattern);
file_path           = strtrim(file_path);
D                   = spm_eeg_load(file_path);

D.inv{val}.source.XYZ          = [38 ,22, -4 ; 0 ,24, 36  ; -32, 24 ,-2; -1 ,54, 27 ; 0 ,-55, 32  ; 49, -63 ,33 ; -46 ,-66 ,30 ];
D.inv{val}.source.label        = {'RAI', 'ACC', 'LAI','MPFC', 'PCC', 'LAG', 'RAG'};
D.inv{val}.source.rad          = 8;
D.inv{val}.source.fname = fullfile(sourceDir, ['XX_' Subname '_DSALL']);
D.inv{val}.source.type         = 'trials';
D = spm_eeg_inv_extract(D);
end

