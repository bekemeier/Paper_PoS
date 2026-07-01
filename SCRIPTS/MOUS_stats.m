%% Permutation tests using t-statistic with max-stat correction across time.
% The code loads coefficients of determination of the original data
% (CODall) and the shuffled data (shufCODall) from a file for each
% parcel. Next, a t-test is performed between the original and shuffled
% data followed by t-tests for permuted data labels (5000 permutations).
% Then the p-value is computed. The p-values<alpha (here, it is 0.05) and
% the timing of the p-values<alpha are written in to variables. The script
% then saves 't_stat' and 'p_perm' values for each timepoint, and
% p-values<alpha ('sig_tp') and their timing ('sig_time') into a file.
clear
thisDir = mfilename("fullpath");
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1);
%addpath to the fieldtrip folder
addpath /mnt/smbdir/fieldtrip-20230402; %check!!!

ft_defaults;
ft_hastoolbox('cellfunction',1)

basedir = fullfile(basedirectory,'MOUSDATA','Data','Derivatives','Words_V_Sent','stats_allsce');
datadir =fullfile(basedir,'Interactions');
savedir =fullfile(basedir,'t_stat_Interactions');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Beginning analysis script %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = fullfile(basedirectory,'MOUSDATA','Data','Preproc','time.mat');
load(filename, 'time');

% define the models that will be analysed
model_names = {'Cont'}; % can add more models

for j = 1:length(model_names)
    
    modelname = model_names{j};
    disp(['Processing: ', modelname]);

    n_time    = 109;
    n_parcels = 382;
   
    for p = 1:n_parcels
        fname = fullfile(datadir, ...
            sprintf('mscca_sceALL_parcel%03d_Int_POS_%s_FULLRED.mat', p, modelname));
        load(fname, 'CODall','shufCODall');
        COD=CODall;
        meanshufCOD=shufCODall;
        clear CODall shufCODall
        % here, we remove the baseline time window from the calculations
        tp_use = find(time >= 0);
        % Run the permutation test
        COD_use        = COD(tp_use,:);
        meanshuf_use  = meanshufCOD(tp_use,:);
        clear COD meanshufCOD
     % Run the permutation test
        results = perm_test_COD_onetail(COD_use, meanshuf_use, 5000);

        t_stat = results.Tobs;
        p_perm = results.Pvals;
        max_null = results.max_null;
        sig_tp = find(p_perm < 0.05);
        sig_time = time(sig_tp);
        
      
        filename = fullfile(savedir, sprintf('mscca_Interaction_parcel%03d_%s_FULLRED.mat',p,modelname));
        save(filename, 't_stat','p_perm','sig_tp','sig_time','max_null');
        clear COD_use meanshuf_use    
    end

   
   
end
