%% Post-hoc permutation tests (PoS contrasts)
% The script loops through each model and each significant parcel and performs
% post-hoc tests on PoS regression coefficients extracted from the encoding models.
%
% For each parcel, regression coefficients are loaded from the full model.
% PoS effects are defined by the verb and adjective coefficients, while the
% noun condition is implicitly represented as the zero baseline due to dummy coding.
%
% Three contrasts are computed at each timepoint and subject:
% (i) verb vs adjective (direct PoS category comparison),
% (ii) verb vs noun (verb coefficients tested against zero),
% (iii) adjective vs noun (adjective coefficients tested against zero).
%
% Statistical inference is performed using non-parametric permutation-based
% dependent-samples t-tests with sign-flipping across subjects (5000 permutations).
% For each contrast, a t-statistic is computed at each timepoint and a
% permutation null distribution is generated using the maximum absolute t-value
% across timepoints to control for multiple comparisons over time (FWER correction).
%
% Significance is assessed using the corresponding permutation-derived p-values,
% with Bonferroni correction applied across the three planned contrasts
% (α = 0.05/3 ≈ 0.0167). The procedure ensures control of family-wise error
% across time while testing multiple PoS contrasts within each parcel.
%
% The resulting t-statistics, p-values, and significant timepoints are
% saved for each parcel and each contrast for subsequent statistical analysis.
clear
thisDir = mfilename("fullpath")
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1)


addpath /mnt/smbdir/fieldtrip-20230402; %check!!!

 
ft_defaults;
ft_hastoolbox('cellfunction',1)
basedir = fullfile(basedirectory,'MOUSDATA','Data','Derivatives','Words_V_Sent','stats_allsce');
datadir =fullfile(basedir,'Ind_Pred');
savedir =fullfile(basedir,'t_stat_Ind_Pred');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Beginning analysis script %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_names = {'POS'};
for k =1:length(model_names)
    modelname = model_names{k};
    fprintf('\nRunning %s',modelname)
    % Vector of significant parcels
filename=fullfile(basedir, sprintf('Ind_Pred_FDR_%s_FMperm.mat',modelname)); % IndPred for individual predictors, Interaction
load(filename,"significant_parcels")

% Folder containing the .mat files
data_dir = savedir;  

    % Loop over each significant parcel
    for i = 1:length(significant_parcels)
        parcel_num = significant_parcels(i);
        fprintf('\nModel %s – parcel %d / %d\n', ...
        modelname, i, length(significant_parcels));

        % Construct the filename (with leading zeros)
        filename = sprintf('PostHoc_Ind_Pred_STparcel%03d_%s.mat', parcel_num, modelname);
        fullpath = fullfile(data_dir, filename);
        
                 
                    
                
                %load coefficients 
                filename = fullfile(datadir,sprintf('mscca_sceALL_STparcel%03d_FULL_coeffs_WL',parcel_num));
                load(filename, 'FULLcoeffs');
                POS = FULLcoeffs(:, [end-1 end], :);
                % Inputs for parcel p:
                verb_betas = squeeze(POS(:, 1, :));  % [109 × 100]
                adj_betas  = squeeze(POS(:, 2, :));  % [109 × 100]
                noun_zeros = zeros(size(verb_betas));
            
                % Run permutation test
                [t_statVA, p_valsVA, sig_timesVA] = compare_betas_fast(verb_betas, adj_betas, 5000, 0.0167); %two-tailed betw. conditions
                [t_statV, p_valsV, sig_timesV] = compare_betas_fast(verb_betas, noun_zeros,5000, 0.0167); % one-tailed (?)
                [t_statA, p_valsA, sig_timesA] = compare_betas_fast(adj_betas, noun_zeros,5000, 0.0167); % one-tailed (?)
            
                filename = fullfile(savedir, sprintf('PostHoc_Ind_Pred_STparcel%03d_%s',parcel_num,modelname));
                save(filename, 't_statVA','p_valsVA','sig_timesVA','t_statV','p_valsV','sig_timesV','t_statA','p_valsA','sig_timesA');
                clear POS verb_betas adj_betas t_statVA p_valsVA sig_timesVA t_statV p_valsV sig_timesV t_statA p_valsA sig_timesA
            
        
    end

end
