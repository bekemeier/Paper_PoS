%% Model comparison
% In order to separate the unique variance explained by each variable of interest from that
% explained by all other variables, Huizeling et al.(2021) applied a model comparison scheme. 
% The model comparison procedure quantified the extent to which a model including a predictor of interest
% explained variance in the MEG signal, above and beyond a reduced model that did not
% include the given predictor. 

% This code computes the coefficient of determination (COD or R^2) as: 
% R^2 = 1-(errorFullModel/errorReducedModel)
% The script first loads the error matrices (dim: time x subject) for the full model and reduced
% model, computes individual CODs per timepoint and subject and saves the
% results for statistics. 
% In the next step, the script calcuates CODs for the additional random
% permutation models and saves the results for statistics.


clear
thisDir = mfilename("fullpath");
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1);

%addpath  to the fieldtrip folder
addpath /mnt/smbdir/fieldtrip-20230402; %check!!!
ft_defaults;
ft_hastoolbox('cellfunction',1)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Beginning analysis script %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for scenario = 1:6
    % Retrieve subjects for this scenario (only visual modality)
    load(fullfile(basedirectory,'MOUSDATA','Data','Preproc','subj_sce_info.mat'))
    sel = strncmp(sce, num2str(scenario), 1);
    subj = subj(sel);
    selaud = strncmp(subj, 'V', 1);  % auditory only
    subj = subj(selaud);
    clear sce

    fprintf('Scenario %d | nSub = %d\n', scenario, length(subj));

    % Directories
    dat_dir= fullfile(savedirectory,'Derivatives','Words_V_Sent',sprintf('scenario%d',scenario));
    datadir = fullfile(dat_dir,'Models');
    savedir = fullfile(dat_dir,'Stats');

    for p = 1:382  % loop over parcels
        parcel_indx = p;

        % --- Load Full Model ---
        filename = fullfile(datadir,sprintf('mscca_sce%d_parcel%03d_Interaction_POSCont_FULLRED',scenario,parcel_indx));
        load(filename,'err_full','err_red','shuf_full','shuf_red');  
        

        
        nSubj = length(subj);
        nTime = size(err_full, 1);

        % ---------- PREALLOC ----------
        COD         = zeros(109, nSubj);
        meanShufCOD = zeros(109, nSubj);

        % ---------- COMPUTE COD ----------
        for tp = 1:109
            for s = 1:nSubj

                ef  = err_full(tp,s);
                er  = err_red(tp,s);

                esf = mean(shuf_full(tp,:,s), 2);
                esr = mean(shuf_red(tp,:,s), 2);

                % numerical safety
                if er > 0 && esr > 0
                    COD(tp,s)         = 1 - ef / er;
                    meanShufCOD(tp,s)= 1 - esf / esr;
                else
                    COD(tp,s)         = NaN;
                    meanShufCOD(tp,s)= NaN;
                end
            end
        end
        % --- Save COD and shufCOD ---
        filename = fullfile(savedir, sprintf('mscca_sce%d_parcel%03d_Int_POS_Cont_FULLRED',scenario,parcel_indx));
        save(filename, 'COD','meanShufCOD','-v7.3');
        clear COD meanShufCOD err_full err_red shuf_full shuf_red
    end
end