%% Post-hoc FDR correction across parcels (PoS contrasts)
% The script aggregates post-hoc permutation test results across parcels
% showing a significant main effect of PoS and performs spatial multiple
% comparisons correction.
%
% For each selected parcel, precomputed post-hoc permutation p-values
% (verb vs noun, adjective vs noun, verb vs adjective) are loaded.
% These p-values are time-resolved (time × subjects) and derived from
% sign-flipping permutation tests on PoS regression coefficients.
%
% For each parcel and each contrast, the minimum p-value across timepoints
% is extracted to obtain a single summary statistic per parcel.
%
% Spatial multiple comparisons across parcels are controlled using the
% Benjamini–Hochberg false discovery rate (FDR) procedure applied separately
% for each PoS contrast (verb vs noun, adjective vs noun, verb vs adjective).
%
% Parcels surviving FDR correction are identified as showing significant
% PoS-related effects at at least one timepoint.
%
% For interpretability, the script additionally extracts the specific
% timepoints at which uncorrected permutation p-values fall below the
% significance threshold (α = 0.0167), separately for each surviving parcel
% and each contrast.
%
% The resulting FDR-corrected parcel sets, q-values, and time-resolved
% significance profiles are saved for subsequent visualization and
% statistical interpretation.
clear
thisDir = mfilename("fullpath")
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1)


addpath /mnt/smbdir/fieldtrip-20230402; %check!!!

 
ft_defaults;
ft_hastoolbox('cellfunction',1)
basedir= fullfile(basedirectory,'MOUSDATA','Data','Derivatives','Words_V_Sent','stats_allsce');
datadir = fullfile(basedir,'t_stat_Ind_Pred'); %individual predictors
savedir = fullfile(basedir); % 
 
% Load necessary data
filename = fullfile(basedirectory,'MOUSDATA','Data','Preproc','time.mat');
load(filename, 'time');

modelname = {'POS'};
for j=1:numel(modelname)
    model_name = modelname{j};
    disp(['Processing: ', model_name]);
alpha_uncorr = 0.0167; 
filename=fullfile(savedir, sprintf('Ind_Pred_FDR_%s_FMperm.mat',model_name)); 
load(filename,"significant_parcels")
n_parcels = length(significant_parcels);
n_time = numel(time);
p_vals_verb_vs_noun = zeros(n_time,n_parcels);
p_vals_adj_vs_noun = zeros(n_time,n_parcels);
p_vals_verb_vs_adj = zeros(n_time,n_parcels);


for k = 1:n_parcels
        % Load perm for the parcel
        parcel_num = significant_parcels(k);
        filename = sprintf('PostHoc_Ind_Pred_STparcel%03d_%s.mat', parcel_num, model_name);
        fullpath = fullfile(datadir, filename);

        if exist(fullpath, 'file')
            load(fullpath, 'p_valsVA', 'p_valsV', 'p_valsA');
            p_vals_verb_vs_noun(:, k) = p_valsV(:);
            p_vals_adj_vs_noun(:,  k) = p_valsA(:);
            p_vals_verb_vs_adj(:,  k) = p_valsVA(:);
        else
            warning('Missing file: %s', fullpath);
            p_vals_verb_vs_noun(:, k) = NaN;
            p_vals_adj_vs_noun(:,  k) = NaN;
            p_vals_verb_vs_adj(:,  k) = NaN;
        end
        
end
% --- Reduce time series to one p-value per parcel (min across time) ---
    min_p_V  = min(p_vals_verb_vs_noun, [], 1, 'omitnan');
    min_p_A  = min(p_vals_adj_vs_noun,  [], 1, 'omitnan');
    min_p_VA = min(p_vals_verb_vs_adj,  [], 1, 'omitnan');

    % --- Apply FDR correction (Benjamini-Hochberg) across parcels ---
    [h_V, q_V]   = fdr_bh(min_p_V, alpha_uncorr, 'pdep');
    [h_A, q_A]   = fdr_bh(min_p_A, alpha_uncorr, 'pdep');
    [h_VA, q_VA] = fdr_bh(min_p_VA, alpha_uncorr, 'pdep');

    % Get the surviving parcel indices (relative to significant_parcels)
    sig_parcels_V  = significant_parcels(h_V == 1);
    sig_parcels_A  = significant_parcels(h_A == 1);
    sig_parcels_VA = significant_parcels(h_VA == 1);

    % --- Also store timepoints of significance for each FDR-surviving parcel ---
    % These are timepoints with uncorrected p < 0.05 (before FDR)
    timepoints_sig_V  = cell(size(sig_parcels_V));
    timepoints_sig_A  = cell(size(sig_parcels_A));
    timepoints_sig_VA = cell(size(sig_parcels_VA));

    % V (Verb vs Noun)
    for idx = 1:numel(sig_parcels_V)
        parcel_idx = find(significant_parcels == sig_parcels_V(idx));
        timepoints_sig_V{idx} = time(p_vals_verb_vs_noun(:, parcel_idx) < alpha_uncorr);
    end

    % A (Adj vs Noun)
    for idx = 1:numel(sig_parcels_A)
        parcel_idx = find(significant_parcels == sig_parcels_A(idx));
        timepoints_sig_A{idx} = time(p_vals_adj_vs_noun(:, parcel_idx) < alpha_uncorr);
    end

    % VA (Verb vs Adj)
    for idx = 1:numel(sig_parcels_VA)
        parcel_idx = find(significant_parcels == sig_parcels_VA(idx));
        timepoints_sig_VA{idx} = time(p_vals_verb_vs_adj(:, parcel_idx) < alpha_uncorr);
    end

    % Save results
    save(fullfile(savedir, sprintf('postFDR_perm_%s_post_hoc.mat', model_name)), ...
        'sig_parcels_V', 'sig_parcels_A', 'sig_parcels_VA', ...
        'q_V', 'q_A', 'q_VA', ...
        'h_V', 'h_A', 'h_VA', ...
        'min_p_V', 'min_p_A', 'min_p_VA', ...
        'timepoints_sig_V', 'timepoints_sig_A', 'timepoints_sig_VA');
end