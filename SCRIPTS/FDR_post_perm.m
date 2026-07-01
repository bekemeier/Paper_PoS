%% Statistical Analysis: Spatial Correction
% This script aggregates parcel-wise statistical results for the interaction
% models and performs correction for multiple comparisons across cortical
% parcels. For each parcel, time-resolved t-values and permutation-corrected
% p-values are loaded from the previous analysis stage, in which temporal
% multiple comparisons were controlled using a max-statistic permutation
% procedure.

% The time-corrected p-values are combined across parcels and subjected to
% Benjamini–Hochberg false discovery rate (FDR) correction (alpha = 0.05)
% to control for spatial multiple comparisons.

% Adjusted p-values are reshaped into time × parcel format and used to
% create a binary significance mask. Parcels containing at least one
% significant timepoint are identified as significant parcels, whereas
% parcels with no surviving effects are stored separately.

% The smallest t-value surviving FDR correction is retained as a threshold
% for visualization purposes.

% The script saves the corrected p-values, significance mask, significant
% and non-significant parcel indices, together with the original t- and
% p-value matrices for subsequent analysis and visualization.
clear
thisDir = mfilename("fullpath");
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1);
%addpath to the fieldtrip folder
addpath /mnt/smbdir/fieldtrip-20230402;
ft_defaults;
ft_hastoolbox('cellfunction',1)

dat_dir = fullfile(basedirectory,'MOUSDATA','Data','Derivatives','Words_V_Sent','stats_allsce');
datadir = fullfile(dat_dir,'t_stat_Interactions'); %interactions
savedir = fullfile(dat_dir); 
n_parcels = 382;

filename = fullfile(basedirectory,'MOUSDATA','Data','Preproc','time.mat');
load(filename, 'time');
t0_idx = find(time >= 0);
n_time = numel(t0_idx);

% allocate
all_t = nan(n_time, n_parcels);
all_p = nan(n_time, n_parcels);
% loop through the models
 
modelname = {'Cont'};%'Entropy', 'Index','Lexfreq','Logprob','Word';
for j=1:numel(modelname)
    model_name = modelname{j};
    disp(['Processing: ', model_name]);
     % --------------------------------------------------
    % Load the per-parcel results
    % --------------------------------------------------
    for p = 1:n_parcels
        fname = fullfile(datadir, sprintf('mscca_Interaction_parcel%03d_%s_FULLRED.mat', p, model_name));
        S = load(fname,'t_stat','p_perm');
    
        all_t(:,p) = S.t_stat;      % 109×1 → stored columnwise
        all_p(:,p) = S.p_perm;      % already time-corrected
    end
    
   %% % FDR across all parcels × timepoints
% -----------------------------
% Flatten all p-values into a vector
p_vec = all_p(:);   % [n_time * n_parcels x 1]

% Apply Benjamini-Hochberg FDR correction, assuming dependent tests
[~, ~, adj_pvals_flat] = fdr_bh(p_vec, 0.05, 'pdep', 'yes');

% Reshape back to [n_time x n_parcels]
adj_pvals = reshape(adj_pvals_flat, size(all_p));

% Create a binary mask of significant timepoints × parcels
signif_mask = adj_pvals < 0.05;

% Identify parcels with at least one significant timepoint
significant_parcels = find(any(signif_mask, 1));

% Identify parcels with no significant timepoints
null_parcels = find(all(~signif_mask, 1));
    
    % Threshold t-value (smallest t that survived)
    if ~isempty(significant_parcels)
        t_thresh = min(all_t(signif_mask));
    else
        t_thresh = NaN;
    end
    
    fprintf('Model: %s\n', model_name);
    fprintf('Significant parcels after FDR: %d\n', numel(significant_parcels));
    filename=fullfile(savedir, sprintf('Interaction_FDR_%s_FMperm.mat',model_name)); % IndPred for individual predictors, Interaction
    save(filename, 'signif_mask', 'significant_parcels', 'null_parcels', 'adj_pvals', 't_thresh', 'all_p', 'all_t');
 end
